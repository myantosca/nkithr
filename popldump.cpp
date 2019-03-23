#include <iostream>
#include <mpi.h>
#include <sys/types.h>
#include <cstring>

int main(int argc, char *argv[]) {
  int r, k;
  uint32_t a = 0, n, i;
  const char *fname = "popl.out";
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  MPI_Comm_size(MPI_COMM_WORLD, &k);
  
  if (r == 0) {
    MPI_File popl_fp;
    MPI_Status popl_status;

    // Default global population.
    n = 0;
    // Check for user-defined global population.
    while (a < argc) {
      // Population count argument.
      if (!strcmp("-n", argv[a])) {
	a++;
	if (a < argc) sscanf(argv[a], "%u", &n);
      }
      // Input filename argument.
      if (!strcmp("-f", argv[a])) {
	a++;
	if (a < argc) fname = argv[a];
      }
      a++;
    }

    // Allocate and zero out a read buffer.
    uint32_t *N = new uint32_t[n];
    memset(N, 0x0, n * sizeof(uint32_t));

    // Read all the values from the dump file.
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &popl_fp); 
    MPI_File_read(popl_fp, N, n, MPI_UNSIGNED, &popl_status);
    MPI_File_close(&popl_fp);

    // Print out each value on a separate line.
    for (i = 0; i < n; i++) {
      std::cout << N[i] << std::endl;
    }
    // Clean up the read buffer.
    delete[] N;
  }
  // Clean up.
  MPI_Finalize();
  return 0;
}
