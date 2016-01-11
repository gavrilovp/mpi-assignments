#include "stdio.h"
#include "mpi.h"

int main(int argc, char **argv) {
  int nprocs, pid;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  printf("my id = %d\n", pid);
  if (pid == 0) {
    printf("total = %d\n", nprocs);
  }
  
  MPI_Finalize();

  return 0;
}
