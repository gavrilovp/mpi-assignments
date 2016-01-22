#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef MPI
#ifdef MPI_PCONTROL
#include <pcontrol.h>
#endif
#include <mpi.h>
#include "utils.h"
#elif OPENMP
#include <omp.h>
#endif

typedef int bool;
#define true 1
#define false 0

typedef struct {
  unsigned vertex_count;
  int** matrix;
} AdjacencyMatrix;

typedef struct {
  unsigned vertex_count;
  int* vertices;
} Clique;

int add_vertex_to_clique(int* clique, int* clique_size,
                         AdjacencyMatrix* adj_matrix, int v_i) {
  int clique_vertex_i = 0;
  for (clique_vertex_i; clique_vertex_i < *clique_size; clique_vertex_i++) {
    if (!(*adj_matrix).matrix[clique_vertex_i, v_i]) {
      return -1;
    }
  }
  clique[*clique_size] = v_i;
  (*clique_size)++;
  return 0;
}

#ifdef OPENMP
bool belongs_clique(int v_i, int* clique, int* clique_size) {
  bool found = false;
  for (int i=0; i<*clique_size; i++) {
    if (v_i == clique[i])
      if (!found) {
        found = true;
      }
  }
  return found;
}
#else
bool belongs_clique(int v_i, int* clique, int* clique_size) {
  for (int i=0; i<*clique_size; i++) {
    if (v_i == clique[i])
      return true;
  }
  return false;
}
#endif

#ifdef OPENMP
bool is_adjoinable_vertex(int* clique, int* clique_size,
                          int v_i, AdjacencyMatrix* adj_matrix) {
  bool found = true;
  for (int i=0; i<*clique_size; i++) {
    if (!(*adj_matrix).matrix[clique[i]][v_i]) {
      if (found) {
        found = false;
      }
    }
  }
  return found;
}
#else
bool is_adjoinable_vertex(int* clique, int* clique_size,
                          int v_i, AdjacencyMatrix* adj_matrix) {
  for (int i=0; i<*clique_size; i++) {
    if (!(*adj_matrix).matrix[clique[i]][v_i]) {
      return false;
    }
  }
  return true;
}
#endif

void procedure_1(int* clique, int* clique_size, AdjacencyMatrix* adj_matrix) {
  int max_rho = 0;
  int* max_clique = clique;
  int* max_clique_size = clique_size;
  
  for (int i=0; i<(*adj_matrix).vertex_count; i++) {
    if (!belongs_clique(i, clique, clique_size) &&
        is_adjoinable_vertex(clique, clique_size, i, adj_matrix)) {
      // Obtain new clique with adjoinable vertex.
      int tmp_clique_size = (*clique_size) + 1;
      int tmp_clique[tmp_clique_size];
      memcpy(tmp_clique, clique, sizeof(int) * (*clique_size));
      tmp_clique[tmp_clique_size-1] = i;

      int current_rho = 0;

      for (int j=0; j<(*adj_matrix).vertex_count; j++) {
        if (!belongs_clique(j, clique, clique_size) &&
            is_adjoinable_vertex(clique,clique_size, j, adj_matrix)) {
          current_rho++;
        }
      }
      if (current_rho > max_rho) {
        max_rho = current_rho;
        memcpy(max_clique, tmp_clique, sizeof(int) * tmp_clique_size);
        *max_clique_size = tmp_clique_size;
      }
    }
  }

  clique = max_clique;
  clique_size = max_clique_size;
}

int malloc2dint(int ***matrix, int n) {

    /* allocate the n*m contiguous items */
    int *p = (int *)malloc(n*n*sizeof(int));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*matrix) = (int **)malloc(n*sizeof(int*));
    if (!(*matrix)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<n; i++) 
       (*matrix)[i] = &(p[i*n]);

    return 0;
}

int free2dint(float ***matrix) {
    /* free the memory - the first element of the array is at the start */
    free(&((*matrix)[0][0]));

    /* free the pointers into the memory */
    free(*matrix);

    return 0;
}

int main(int argc, char** argv) {
  AdjacencyMatrix adj_matrix;

#ifdef MPI
  mpi_check(MPI_Init(&argc, &argv));

  // Set a sane error handler (return errors and don't kill our process)
  mpi_check(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN));

#ifdef MPI_PCONTROL
  // Set trace parameters
  MPI_Pcontrol(TRACELEVEL, 1, 1, 1);
  MPI_Pcontrol(TRACEFILES, "clique.trace.tmp", "clique.trace", 0);
  MPI_Pcontrol(TRACESTATISTICS, 200, 1, 1, 1, 1, 1);
  // Start trace
  MPI_Pcontrol(TRACENODE, 1024 * 1024, 1, 1);
#endif

  MPI_Datatype mpi_clique_type;
  {
    MPI_Type_vector( 16, 1, 16, MPI_INT, &mpi_clique_type); 
    MPI_Type_commit( &mpi_clique_type );
    mpi_check(MPI_Type_commit(&mpi_clique_type));
  }
  
  int pid;
  mpi_check(MPI_Comm_rank(MPI_COMM_WORLD, &pid));
  int nprocs;
  mpi_check(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));

  if (pid == 0) {
#endif

  {
    FILE *f = fopen("adjacency_matrix.txt", "r");
    if (f == NULL) {
      perror("fopen");
      exit(EXIT_FAILURE);
    }
    fscanf(f, "%u", &adj_matrix.vertex_count);
    if (fclose(f) == -1) {
      perror("close");
      exit(EXIT_FAILURE);
    }
  }

#ifdef MPI
  }
  mpi_check(MPI_Bcast(&adj_matrix.vertex_count, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD));

  int* displs = NULL;
  if (pid == 0) {
    displs = calloc(nprocs, sizeof(int));
  }

  int* recvcounts = calloc(nprocs, sizeof(int));
  mpi_check(UMPI_Scatterv_lens(adj_matrix.vertex_count, recvcounts, displs, 0, MPI_COMM_WORLD));
  int recvlen = recvcounts[pid];
  if (pid != 0) {
    free(recvcounts);
  }
  printf("Pid: %i, recvlen: %d\n", pid, recvlen);
#endif

  malloc2dint(&adj_matrix.matrix, adj_matrix.vertex_count);

#ifdef MPI
  printf("Pid %d, vertex_count %d\n", pid, adj_matrix.vertex_count);
  if (pid == 0) {
#endif

  // Second fopen because of necessary malloc for all MPI master/workers
  FILE *f = fopen("adjacency_matrix.txt", "r");
  if (f == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  fscanf(f, "%u", &adj_matrix.vertex_count);
  
  for (int i=0; i<adj_matrix.vertex_count; i++) {
    for (int j=0; j<adj_matrix.vertex_count; j++) {
      fscanf(f, "%d", &adj_matrix.matrix[i][j]);
    }
  }

  if (fclose(f) == -1) {
    perror("close");
    exit(EXIT_FAILURE);
  }

#ifdef MPI
  }
  mpi_check(MPI_Bcast(&(adj_matrix.matrix[0][0]), adj_matrix.vertex_count * adj_matrix.vertex_count, MPI_INT, 0, MPI_COMM_WORLD));
#endif

  int* maximal_clqs_vrts_cnt = NULL;
  int* maximal_cliques = NULL;//[adj_matrix.vertex_count][adj_matrix.vertex_count];

#ifdef OMPI
  if (pid == 0) {
#endif
    maximal_cliques = calloc(adj_matrix.vertex_count, sizeof(int));
#ifdef OMPI
  }
#endif

#ifdef MPI
  int* sendbuf = calloc(recvlen, sizeof(int));
  for (unsigned i = 0; i < recvlen; i++) {
    sendbuf[i] = 1;
  }
  
  if (pid != 0) {
    mpi_check(MPI_Gatherv(sendbuf, recvlen, MPI_INT, NULL, NULL, NULL, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD));
  } else {
    mpi_check(MPI_Gatherv(sendbuf, recvlen, MPI_INT, maximal_cliques, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD));
    for (int i=0; i<adj_matrix.vertex_count; i++) {
      printf("%d ", maximal_cliques[i]);
    }
  }

  mpi_check(MPI_Finalize());
#else
  for (int i = 0; i < adj_matrix.vertex_count; i++) {
    /* procedure_1(clique, &clique_size, &adj_matrix); */
  }
#endif

  free2dint(&adj_matrix.matrix);
  
  return 0;
}
