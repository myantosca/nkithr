#include <iostream>
#include <chrono>
#include <random>
#include <mpi.h>
#include <sys/types.h>
#include <cstring>

using namespace std::chrono;

int main(int argc, char *argv[]) {
  int r, k;
  uint32_t a = 0, i, j, n, m, v, lb, ub;
  uint32_t c, e;
  uint32_t s[3] = {0, 0, 0};
  uint32_t S[3] = {0, 0, 0};
  uint32_t pivot, pivot_cand;
  double pivot_weight;
  uint32_t msgs = 0;
  uint32_t rnds = 1;
  uint32_t prior = 0;
  uint32_t dups = 0;
  MPI_File popl_fp;
  MPI_Offset popl_off = 0;
  MPI_Request popl_req;
  MPI_Status popl_status;
  const char *fname = "popl.out";
  bool popldump = false;
  int genorder = 0;
  bool deterministic = false;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  MPI_Comm_size(MPI_COMM_WORLD, &k);

  time_point<system_clock, nanoseconds> tp_a = system_clock::now();

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
    if (!strcmp("-D", argv[a])) {
      deterministic = true;
    }
    a++;
  }

  // Determine the local population, spreading the remainder evenly over the nodes.
  m = n / k + ((n % k) > r);
  for (j = 0; j < r; j++) {
    prior += n / k + ((n % k) > j);
  }
  popl_off = prior * sizeof(uint32_t);

  // Put bound on dups to avoid requiring communication for duplication of random elements.
  if (genorder == 0) {
    dups = dups < n / k ? dups : n / k - 1;
  }
  // Create a seeded PNRG based on the standard Mersenne Twister engine.
  unsigned seed = system_clock::now().time_since_epoch().count();
  std::mt19937 prng(seed);

  uint32_t *pivot_cands = new uint32_t[k];
  double *pivot_weights = new double[k];
  memset(pivot_cands, prng.min(), k * sizeof(uint32_t));
  // Populate the local set of candidates.
  uint32_t *M = new uint32_t[m];
  uint32_t *Q = new uint32_t[m];
  memset(M, prng.min(), m * sizeof(uint32_t));

  // Generate the local population.
  j = 0;
  while (j < m) {
    v = (genorder > 0) ? (prior + j) / (dups + 1) : ((genorder < 0) ? (n - (prior + j + 1)) / (dups + 1) : prng());
    M[j++] = v;
    if (genorder == 0) {
      for (uint32_t d = 0; d < dups; d++) { M[j++] = v; }
    }
  }

  // Dump local population for post-mortem debugging.
  if (popldump) {
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &popl_fp);
    MPI_File_set_view(popl_fp, popl_off, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
    MPI_File_iwrite_at(popl_fp, 0, M, m, MPI_UNSIGNED, &popl_req);
    // Wait for population dump. Close dump file.
    MPI_Wait(&popl_req, &popl_status);
    MPI_File_close(&popl_fp);
  }

  lb = 0;
  ub = m;

  uint32_t p = 0;
  // Pivot selection PRNGs.
  std::knuth_b pivot_cand_gen (seed);
  std::knuth_b pivot_gen (seed);

  time_point<system_clock, nanoseconds> tp_x = system_clock::now();
  // If |S1| < i <= |S1| + |P| then done.
  while (!(S[0] < i && S[1] >= i)) {
    std::uniform_int_distribution<uint32_t> pivot_cand_select(lb, ub - 1);
    // Deterministic case: Take the middle element of the unsorted live population.
    // Random case: Select from a uniform distribution over the unsorted live population. Choose lower bound if no live population.
    // The fallback is safe since it will not change the overall set cardinalities.
    // It has either already served as a pivot or is out of range for pivoting.
    // The only cost is a wasted round.
    pivot_cand = deterministic ? M[lb + ((ub - lb) >> 1)] : (lb < ub) ? M[pivot_cand_select(pivot_cand_gen)] : M[lb];
    pivot_weight = ub - lb;
    // Gather pivot candidate weights from all nodes.
    MPI_Gather(&pivot_weight, 1, MPI_DOUBLE, pivot_weights, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    msgs += k - 1;

    // Gather pivot candidates from all nodes.
    MPI_Gather(&pivot_cand, 1, MPI_UNSIGNED, pivot_cands, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    msgs += k - 1;

    // Root chooses pivot p from candidates according to weighted discrete distribution over [0,k).
    if (r == 0) {
      if (deterministic) {
	// Round-robin selection of live pivots (poor man's uniform distribution).
	while (pivot_weights[p] == 0) { p = (p + 1) % k; }
	pivot = pivot_cands[p];
      }
      else {
	// No reason to waste messages since k should be small relative to n.
	// The sum of weights is computed at the root.
	double pivot_weight_total = 0;
	for (j = 0; j < k; j++) {
	  pivot_weight_total += pivot_weights[j];
	}
	// To make use of std::discrete_distribution, we have to copy the array to something iterable.
	// The division by live population total is done here also.
	std::vector<double> pivot_weights_vec;
	for (j = 0; j < k; j++) {
	  pivot_weights_vec.push_back(pivot_weights[j] / pivot_weight_total);
	}
	// Create the round pivot selection distribution.
	std::discrete_distribution<uint32_t> pivot_select(pivot_weights_vec.begin(), pivot_weights_vec.end());
	// Select the pivot.
	pivot = pivot_cands[pivot_select(pivot_gen)];
      }
    }
    // Broadcast pivot to all nodes.
    MPI_Bcast(&pivot, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    msgs += k - 1;

    // Initialize the set cardinality index variables.
    s[0] = lb;
    s[1] = 0;
    s[2] = ub;

    // Set the lower bound for starting the search.
    uint32_t c = lb;

    // Increment until remaining search space is exhausted.
    while (c < ub) {
      if (M[c] < pivot) {
	Q[s[0]++] = M[c];
      }
      else if (M[c] > pivot) {
	Q[--s[2]] = M[c];
      }
      else {
	s[1]++;
      }
      c++;
    }

    // s[1] = |S1| + |P|
    s[1] += s[0];

    // Fill in P.
    for (c = s[0]; c < s[1]; c++) {
      Q[c] = pivot;
    }

    // Copy back rearranged population subspace.
    if (lb < ub) memcpy(M + lb, Q + lb, (ub - lb) * sizeof(uint32_t));

    // Root sums all cardinalities sent and broadcasts the sum to all nodes.
    MPI_Allreduce(&s, &S, 3, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    msgs += (k - 1) * 2;
    // If |S1| > i - 1 then search in the range [lb,s[0]), i.e., below the pivot.
    if (S[0] > i - 1) {
      ub = s[0];
    }
    // If |S1| + |P| < i then search in the range (s[1],ub), i.e., above the pivot.
    else if (S[1] < i) {
      lb = s[1];
    }

    rnds++;
  }

  delete[] M;
  delete[] Q;
  delete[] pivot_cands;
  delete[] pivot_weights;

  time_point<system_clock, nanoseconds> tp_y = system_clock::now();
  int64_t T[2] = { (tp_y - tp_x).count(), (tp_y - tp_a).count() };
  int64_t T0[2] = { 0, 0 };
  MPI_Reduce(&T, &T0, 2, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  // Root announces ith member of population.
  if (r == 0) {
    std::cout << "n,k,i,p,m,r,t,T";
    std::cout << std::endl;
    std::cout << n << ",";
    std::cout << k << ",";
    std::cout << i << ",";
    std::cout << pivot << ",";
    std::cout << msgs << ",";
    std::cout << rnds << ",";
    std::cout << T0[0] << ",";
    std::cout << T0[1] << std::endl;
  }

  // Clean up.
  MPI_Finalize();

  return 0;
}
