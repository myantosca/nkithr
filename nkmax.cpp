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
  uint32_t msgs = 0;
  uint32_t rnds = 1;
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
    // Input filename argument.
    if (!strcmp("-i", argv[a])) {
      a++;
      if (a < argc) fname = argv[a];
    }
    a++;
  }

  // Determine the local population, spreading the remainder evenly over the nodes.
  m = n / k + ((n % k) > r);
  for (i = 0; i < r; i++) {
    popl_off += n / k + ((n % k) > i);
  }
  popl_off *= sizeof(uint32_t);
  MPI_Datatype filetype, contig;
  MPI_Type_contiguous(m, MPI_UNSIGNED, &filetype);

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

  // Dump local population for post-mortem debugging.
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &popl_fp);
  MPI_File_set_view(popl_fp, popl_off, MPI_UNSIGNED, MPI_UNSIGNED, "native", MPI_INFO_NULL);
  MPI_File_iwrite_at(popl_fp, 0, M, m, MPI_UNSIGNED, &popl_req);
  // Send local max to root and determine global max.
  MPI_Reduce(&local_max, &global_max, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
  msgs += k - 1;
  // Root announces global max.
  if (r == 0) {
    std::cout << "n = " << n << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "u = " << global_max << std::endl;
    std::cout << "m = " << msgs << std::endl;
    std::cout << "r = " << rnds << std::endl;
  }

  // Wait for population dump. Close dump file.
  MPI_Wait(&popl_req, &popl_status);
  MPI_File_close(&popl_fp);

  // Clean up.
  delete[] M;
  MPI_Finalize();

  // Report execution time (wall-clock).
  time_point<system_clock, nanoseconds> tp_b = system_clock::now();
  std::cout << "t = " << (tp_b - tp_a).count() << std::endl;
  return 0;
}
