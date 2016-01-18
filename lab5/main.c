#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_MPI
#include <mpi.h>
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

bool belongs_clique(int v_i, int* clique, int* clique_size) {
  for (int i=0; i<*clique_size; i++) {
    if (v_i == clique[i])
      return true;
  }
  return false;
}

bool is_adjoinable_vertex(int* clique, int* clique_size,
                          int v_i, AdjacencyMatrix* adj_matrix) {
  for (int i=0; i<*clique_size; i++) {
    if (!(*adj_matrix).matrix[clique[i]][v_i]) {
      return false;
    }
  }
  return true;
}

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

int main(int argc, char** argv) {
  #ifdef USE_MPI
  MPI_Init(&argc, &argv);

  int nprocs;
  mpi_check(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
  int pid;
  mpi_check(MPI_Comm_rank(MPI_COMM_WORLD, &pid));
  #endif

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

  int clique[adj_matrix.vertex_count];
  int clique_size = 0;

  int max_clique[adj_matrix.vertex_count];
  int max_clique_size = 0;

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

  #ifdef USE_MPI
  MPI_Finalize();
  #endif USE_MPI
  
  return 0;
}
