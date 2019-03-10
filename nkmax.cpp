#include <iostream>
#include <chrono>
#include <random>
#include <mpi.h>
#include <sys/types.h>
#include <cstring>

using namespace std::chrono;

int main(int argc, char *argv[]) {
  int r, k;
  uint32_t a, i, n, m;
  time_point<system_clock, nanoseconds> tp_a = system_clock::now();
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
    a++;
  }

  // Determine the local population, adding 1 in cases where n is not evenly divisible by k.
  m = n / k + ((n % k) != 0);

  // Create a seeded PNRG based on the standard Mersenne Twister engine.
  unsigned seed = system_clock::now().time_since_epoch().count();
  std::mt19937 prng(seed);

  // Populate the local set of candidates.
  uint32_t *M = new uint32_t[m];
  uint32_t local_max = prng.min();
  uint32_t global_max = prng.min();
  for (i = 0; i < m; i++) {
    M[i] = prng();
    // Keep running tab of local maximum.
    local_max = M[i] > local_max ? M[i] : local_max;
  }

  // Send local max to root and determine global max.
  MPI_Reduce(&local_max, &global_max, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);

  // Root announces global max.
  if (r == 0) {
    std::cout << "max = " << global_max << std::endl;
  }

  // Clean up.
  delete[] M;
  MPI_Finalize();

  // Report execution time (wall-clock).
  time_point<system_clock, nanoseconds> tp_b = system_clock::now();
  std::cout << "t = " << (tp_b - tp_a).count() << std::endl;
  return 0;
}
