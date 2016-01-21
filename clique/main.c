#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef MPI
#include <mpi.h>
#elif OPENMP
#include <omp.h>
#endif

typedef int bool;
#define true 1
#define false 0

typedef struct {
  size_t vertex_count;
  int** matrix;
} AdjacencyMatrix;

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

#ifdef LINEAR
      for (int j=0; j<(*adj_matrix).vertex_count; j++) {
        if (!belongs_clique(j, clique, clique_size) &&
            is_adjoinable_vertex(clique,clique_size, j, adj_matrix)) {
          current_rho++;
        }
      }
#elif MPI
      int nprocs;
      mpi_check(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
      int pid;
      mpi_check(MPI_Comm_rank(MPI_COMM_WORLD, &pid));
      
      int* displs;
      if (pid == 0) {
        displs = calloc(nprocs, sizeof(int));
      }
      int* recvcounts = calloc(nprocs, sizeof(int));
      UMPI_Scatterv_lens(INTS_SIZE, recvcounts, displs, 0, MPI_COMM_WORLD);

      int recvlen = recvcounts[pid];
      int* recvbuf = calloc(recvlen, sizeof(int));

      if (pid == 0) {
        srand(time(NULL));

        int msg[INTS_SIZE];
        for (int i = 0; i < INTS_SIZE; i++) {
          msg[i] = rand();
        }

        MPI_Scatterv(msg, recvcounts, displs, MPI_INT, recvbuf, recvlen, MPI_INT, 0, MPI_COMM_WORLD);
      } else {
        MPI_Scatterv(NULL, NULL, NULL, MPI_DATATYPE_NULL, recvbuf, recvlen, MPI_INT, 0, MPI_COMM_WORLD);
      }
      for (int j=0; j<recvlen; j++) {
        if (!belongs_clique(j, clique, clique_size) &&
            is_adjoinable_vertex(clique,clique_size, j, adj_matrix)) {
          current_rho++;
        }
      }
      for (int i = 0; i < recvlen; i++) {
        recvbuf[i] += 1;
      }
      
      if (pid != 0) {
        MPI_Gatherv(recvbuf, recvlen, MPI_INT, NULL, NULL, NULL, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
      } else {
        int msg[INTS_SIZE];
        MPI_Gatherv(recvbuf, recvlen, MPI_INT, msg, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
      }
#endif
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

int main(int argc, char** argv) {
  FILE *f = fopen("adjacency_matrix.txt", "r");
  if (f == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  AdjacencyMatrix adj_matrix;
  fscanf(f, "%zu", &adj_matrix.vertex_count);
  adj_matrix.matrix = malloc(adj_matrix.vertex_count * sizeof(int*));

  for (int i=0; i<adj_matrix.vertex_count; i++) {
    adj_matrix.matrix[i] = malloc(adj_matrix.vertex_count * sizeof(int));
    for (int j=0; j<adj_matrix.vertex_count; j++) {
      fscanf(f, "%d", &adj_matrix.matrix[i][j]);
    }
  }

#ifdef MPI
  MPI_Init(&argc, &argv);

  // Set a sane error handler (return errors and don't kill our process)
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  
  MPI_Datatype mpi_adjacency_matrix_type;
  {
    int blocklengths[] = { 1, adj_matrix.vertex_count * adj_matrix.vertex_count };
    MPI_Datatype types[] = { MPI_UNSIGNED, MPI_INT };
    MPI_Aint offsets[] = { offsetof(struct border, na), offsetof(struct border, a), offsetof(struct border, b) };

    mpi_check(MPI_Type_create_struct(sizeof(blocklengths) / sizeof(int), blocklengths, offsets, types, &mpi_border_type));
    mpi_check(MPI_Type_commit(&mpi_adjacency_matrix_type));
  }

  int pid;
  mpi_check(MPI_Comm_rank(MPI_COMM_WORLD, &pid));
  int nprocs;
  mpi_check(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
  int* displs = NULL;
  if (pid == 0) {
    displs = calloc(nprocs, sizeof(int));
  }
  int* recvcounts = calloc(nprocs, sizeof(int));

  if (pid == 0) {
#endif

  int clique[adj_matrix.vertex_count];
  int clique_size = 0;

  int maximal_cliques[adj_matrix.vertex_count];
  

  int max_clique[adj_matrix.vertex_count];
  int max_clique_size = 0;

#ifdef MPI
  mpi_check(MPI_Bcast(&adj_matrix, 1, , pid, MPI_COMM_WORLD));
  }
#endif

  for (int i = 0; i < adj_matrix.vertex_count; i++) {
    clique[0] = i;
    clique_size = 1;

    procedure_1(clique, &clique_size, &adj_matrix);

    if (clique_size > max_clique_size) {
      max_clique_size = clique_size;
      memcpy(max_clique, clique, sizeof(int) * adj_matrix.vertex_count);
    }
  }

  for (int v_i = 0; v_i < max_clique_size; v_i++) {
    printf("%d\n", max_clique[v_i]);
  }

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
