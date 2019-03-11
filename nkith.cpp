#include <iostream>
#include <chrono>
#include <random>
#include <mpi.h>
#include <sys/types.h>
#include <cstring>

using namespace std::chrono;

int main(int argc, char *argv[]) {
  int r, k;
  uint32_t a = 0, i, j, c, n, m;
  uint32_t S1 = -1, lb, ub;
  uint32_t pivot, pivot_cand;
  MPI_File popl_fp;
  MPI_Offset popl_off = 0;
  MPI_Request popl_req;
  MPI_Status popl_status;
  time_point<system_clock, nanoseconds> tp_a = system_clock::now();
  const char *fname = "popl.out";

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
    a++;
  }

  // Determine the local population, spreading the remainder evenly over the nodes.
  m = n / k + ((n % k) > r);
  for (j = 0; j < r; j++) {
    popl_off += n / k + ((n % k) > j);
  }
  popl_off *= sizeof(uint32_t);

  // Create a seeded PNRG based on the standard Mersenne Twister engine.
  unsigned seed = system_clock::now().time_since_epoch().count();
  std::mt19937 prng(seed);

  uint32_t *pivot_cands = new uint32_t[k];
  memset(pivot_cands, prng.min(), k * sizeof(uint32_t));
  // Populate the local set of candidates.
  uint32_t *M = new uint32_t[m];
  memset(M, prng.min(), m * sizeof(uint32_t));
  uint32_t local_max = prng.min();
  uint32_t global_max = prng.min();
  for (j = 0; j < m; j++) {
    uint32_t v = prng();
    std::cerr << v << std::endl;
    c = 0;
    while ((c < j) && (v >= M[c])) { c++; }
    if (c < j) {
      //memmove(M + (c + 1) * sizeof(uint32_t), M + c * sizeof(uint32_t), (j - c) * sizeof(uint32_t));
      for (int d = j; d > c; d--) {
	M[d] = M[d-1];
      }
    }
    M[c] = v;
  }

  // Dump local population for post-mortem debugging.
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &popl_fp);
  MPI_File_set_view(popl_fp, popl_off, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  MPI_File_iwrite_at(popl_fp, 0, M, m, MPI_UNSIGNED, &popl_req);

  lb = 0;
  ub = m;
  // Round-robin pivot selection counter.
  uint32_t p = 0;
  // If |S1| = i - 1 then done.
  while (S1 != i - 1) {
    pivot_cand = M[lb + ((ub - lb) >> 1)];
    // Gather pivot candidates from all nodes.
    MPI_Gather(&pivot_cand, 1, MPI_UNSIGNED, pivot_cands, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    // Root chooses pivot p from candidates in round-robin fashion (poor man's rando).
    if (r == 0) {
      pivot = pivot_cands[p];
      p = (p + 1) % k;
    }
    // Broadcast pivot to all nodes.
    MPI_Bcast(&pivot, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    // Set the lower bound for starting the search.
    c = lb;
    // Increment until the pivot is passed.
    while ((M[c] < pivot) && (c < ub)) { c++; }
    // Root sums all cardinalities sent and broadcasts the sum to all nodes.
    MPI_Allreduce(&c, &S1, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    // If |S1| > i - 1 then search in the range [lb,c).
    if (S1 > i - 1) {
      ub = c;
    }
    // If |S1| < i - 1 then search in the range (c,ub).
    else if (S1 < i - 1) {
      lb = c;
    }
  }
  // Root announces ith member of population.
  if (r == 0) {
    std::cout << i << "th = " << pivot << std::endl;
  }

  // Wait for population dump. Close dump file.
  MPI_Wait(&popl_req, &popl_status);
  MPI_File_close(&popl_fp);

  // Clean up.
  delete[] M;
  delete[] pivot_cands;
  MPI_Finalize();

  // Report execution time (wall-clock).
  time_point<system_clock, nanoseconds> tp_b = system_clock::now();
  std::cout << "t = " << (tp_b - tp_a).count() << std::endl;
  return 0;
}
