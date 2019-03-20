#include <iostream>
#include <chrono>
#include <random>
#include <mpi.h>
#include <sys/types.h>
#include <cstring>

using namespace std::chrono;

int cmp(const void *pa, const void *pb) {
  uint32_t *a = (uint32_t *)pa;
  uint32_t *b = (uint32_t *)pb;
  return (*a < *b) ? -1 : ((*a > *b) ? 1 : 0);
}

int main(int argc, char *argv[]) {
  int r, k;
  uint32_t a = 0, i, j, c, n, m, v;
  uint32_t S1 = -1, lb, ub;
  uint32_t pivot, pivot_cand;
  uint32_t msgs = 0;
  uint32_t rnds = 1;
  uint32_t prior = 0;
  uint32_t dups = 0;
  MPI_File popl_fp;
  MPI_Offset popl_off = 0;
  MPI_Request popl_req;
  MPI_Status popl_status;
  time_point<system_clock, nanoseconds> tp_a = system_clock::now();
  const char *fname = "popl.out";
  bool popldump = false;
  int genorder = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  MPI_Comm_size(MPI_COMM_WORLD, &k);

  // Default global population.
  n = k * 256;
  // Check for user-defined global population.
  while (a < argc) {
    // Population count argument.
    if (!strcmp("-n", argv[a])) {
      a++;
      if (a < argc) sscanf(argv[a], "%u", &n);
    }
    // ith index argument.
    if (!strcmp("-i", argv[a])) {
      a++;
      if (a < argc) sscanf(argv[a], "%u", &i);
    }
    // Input filename argument.
    if (!strcmp("-f", argv[a])) {
      a++;
      if (a < argc) fname = argv[a];
    }
    // Debug argument (population dump).
    if (!strcmp("-g", argv[a])) {
      popldump = true;
    }
    if (!strcmp("--123", argv[a])) {
      genorder = 1;
    }
    if (!strcmp("--321", argv[a])) {
      genorder = -1;
    }
    if (!strcmp("-d", argv[a])) {
      a++;
      if (a < argc) sscanf(argv[a], "%u", &dups);
    }
    a++;
  }

  // Determine the local population, spreading the remainder evenly over the nodes.
  m = n / k + ((n % k) > r);
  for (j = 0; j < r; j++) {
    prior += n / k + ((n % k) > j);
  }
  popl_off = prior * sizeof(uint32_t);

  // Create a seeded PNRG based on the standard Mersenne Twister engine.
  unsigned seed = system_clock::now().time_since_epoch().count();
  std::mt19937 prng(seed);

  uint32_t *pivot_cands = new uint32_t[k];
  memset(pivot_cands, prng.min(), k * sizeof(uint32_t));
  // Populate the local set of candidates.
  uint32_t *M = new uint32_t[m];
  memset(M, prng.min(), m * sizeof(uint32_t));

  // Generate the local population.
  j = 0;
  while (j < m) {
    v = (genorder > 0) ? prior + j : ((genorder < 0) ? n - (prior + j + 1) : prng());
    M[j] = v;
    j++;
    uint32_t d = dups;
    while (d-- > 0) {
      M[j] = v;
      j++;
    }
  }

  // Pre-sort the local population for easier handling.
  qsort(M, m, sizeof(uint32_t), cmp);

  // Dump local population for post-mortem debugging.
  if (popldump) {
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &popl_fp);
    MPI_File_set_view(popl_fp, popl_off, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
    MPI_File_iwrite_at(popl_fp, 0, M, m, MPI_UNSIGNED, &popl_req);
  }

  lb = 0;
  ub = m;
  // Round-robin pivot selection counter.
  uint32_t p = 0;
  // If |S1| = i - 1 then done.
  while (S1 != i - 1) {
    pivot_cand = M[lb + ((ub - lb) >> 1)];
    // Gather pivot candidates from all nodes.
    MPI_Gather(&pivot_cand, 1, MPI_UNSIGNED, pivot_cands, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    msgs += k - 1;
    // Root chooses pivot p from candidates in round-robin fashion (poor man's rando).
    if (r == 0) {
      pivot = pivot_cands[p];
      p = (p + 1) % k;
    }
    // Broadcast pivot to all nodes.
    MPI_Bcast(&pivot, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    msgs += k - 1;
    // Set the lower bound for starting the search.
    c = lb;
    // Increment until the pivot is passed.
    while ((M[c] < pivot) && (c < ub)) { c++; }
    // Root sums all cardinalities sent and broadcasts the sum to all nodes.
    MPI_Allreduce(&c, &S1, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    msgs += (k - 1) * 2;
    // If |S1| > i - 1 then search in the range [lb,c).
    if (S1 > i - 1) {
      ub = c;
    }
    // If |S1| < i - 1 then search in the range (c,ub).
    else if (S1 < i - 1) {
      lb = c;
    }
    rnds++;
  }
  // Root announces ith member of population.
  if (r == 0) {
    std::cout << "n,k,i,p,m,r";
    for (i = 0; i < k; i++) {
      std::cout << ",t" << i;
    }
    std::cout << std::endl;
    std::cout << n << ",";
    std::cout << k << ",";
    std::cout << i << ",";
    std::cout << pivot << ",";
    std::cout << msgs << ",";
    std::cout << rnds << std::flush;
  }

  if (popldump) {
    // Wait for population dump. Close dump file.
    MPI_Wait(&popl_req, &popl_status);
    MPI_File_close(&popl_fp);
  }

  // Clean up.
  delete[] M;
  delete[] pivot_cands;
  MPI_Finalize();

  // Report execution time (wall-clock).
  time_point<system_clock, nanoseconds> tp_b = system_clock::now();
  std::cout << "," << (tp_b - tp_a).count();
  return 0;
}
