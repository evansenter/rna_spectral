#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include "constants.h"
#include "functions.h"
#include "initializers.h"

double* convert_structures_to_transition_matrix(SOLUTION* all_structures, int num_structures, int use_min) {
  int i, j;
  double col_sum;
  double* transition_matrix = malloc(num_structures * num_structures * sizeof(double));

  for (i = 0; i < num_structures; ++i) {
    col_sum = 0;

    for (j = 0; j < num_structures; ++j) {
      if (i != j) {
        if (use_min) {
          transition_matrix[i + num_structures * j] = MIN(1, exp(-((double)all_structures[j].energy - (double)all_structures[i].energy) / RT));
        } else {
          transition_matrix[i + num_structures * j] = exp(-((double)all_structures[j].energy - (double)all_structures[i].energy) / RT);
        }

        col_sum += transition_matrix[i + num_structures * j];
        #ifdef INSANE_DEBUG
        printf("%d\t%d\t%.4e\n", i, j, transition_matrix[i + num_structures * j]);
        #endif
      }
    }

    #ifdef INSANE_DEBUG
    printf("%d col_sum:\t%.4e\n\n", i, col_sum);
    #endif
    transition_matrix[i * num_structures + i] = -col_sum;
  }

  return transition_matrix;
}

EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double* transition_matrix, int num_structures) {
  int i, j;
  EIGENSYSTEM eigensystem;
  eigensystem = init_eigensystem(num_structures);
  gsl_matrix_view matrix_view      = gsl_matrix_view_array(transition_matrix, eigensystem.length, eigensystem.length);
  gsl_vector_complex* eigenvalues  = gsl_vector_complex_alloc(eigensystem.length);
  gsl_matrix_complex* eigenvectors = gsl_matrix_complex_alloc(eigensystem.length, eigensystem.length);
  gsl_eigen_nonsymmv_workspace* workspace = gsl_eigen_nonsymmv_alloc(eigensystem.length);
  gsl_eigen_nonsymmv_params(1, workspace);
  gsl_eigen_nonsymmv(&matrix_view.matrix, eigenvalues, eigenvectors, workspace);
  gsl_eigen_nonsymmv_free(workspace);

  for (i = 0; i < eigensystem.length; ++i) {
    eigensystem.values[i]               = GSL_REAL(gsl_vector_complex_get(eigenvalues, i));
    gsl_vector_complex_view eigenvector = gsl_matrix_complex_column(eigenvectors, i);

    for (j = 0; j < eigensystem.length; ++j) {
      eigensystem.vectors[i + eigensystem.length * j] = GSL_REAL(gsl_vector_complex_get(&eigenvector.vector, j));
    }
  }

  free(transition_matrix);
  gsl_vector_complex_free(eigenvalues);
  gsl_matrix_complex_free(eigenvectors);
  return eigensystem;
}

void invert_matrix(EIGENSYSTEM eigensystem) {
  int i, j, signum;
  gsl_matrix* matrix_to_invert = gsl_matrix_alloc(eigensystem.length, eigensystem.length);
  gsl_matrix* inversion_matrix = gsl_matrix_alloc(eigensystem.length, eigensystem.length);
  gsl_permutation* permutation = gsl_permutation_alloc(eigensystem.length);

  for (i = 0; i < eigensystem.length; ++i) {
    for (j = 0; j < eigensystem.length; ++j) {
      gsl_matrix_set(matrix_to_invert, i, j, eigensystem.vectors[i * eigensystem.length + j]);
    }
  }

  gsl_linalg_LU_decomp(matrix_to_invert, permutation, &signum);
  gsl_linalg_LU_invert(matrix_to_invert, permutation, inversion_matrix);

  for (i = 0; i < eigensystem.length; ++i) {
    for (j = 0; j < eigensystem.length; ++j) {
      eigensystem.inverse_vectors[i * eigensystem.length + j] = gsl_matrix_get(inversion_matrix, i, j);
    }
  }

  gsl_matrix_free(matrix_to_invert);
  gsl_matrix_free(inversion_matrix);
  gsl_permutation_free(permutation);
}

double probability_at_time(EIGENSYSTEM eigensystem, double timepoint, int start_index, int target_index) {
  // This function is hard-wired to only consider the kinetics for folding from a distribution where p_{0}(start_index) == 1.
  int i;
  double cumulative_probability = 0;

  for (i = 0; i < eigensystem.length; ++i) {
    cumulative_probability +=
      eigensystem.vectors[target_index * eigensystem.length + i] *
      eigensystem.inverse_vectors[i * eigensystem.length + start_index] *
      exp(eigensystem.values[i] * timepoint);
  }

  return cumulative_probability;
}

void print_array(char* title, double* matrix, int length) {
  int i;
  printf("%s\n", title);

  for (i = 0; i < length; ++i) {
    printf("%+.4f\n", matrix[i]);
  }

  printf("\n");
}


void print_matrix(char* title, double* matrix, int length) {
  int i, j;
  printf("%s\n", title);

  for (i = 0; i < length; ++i) {
    for (j = 0; j < length; ++j) {
      printf("%+.4f\t", matrix[i * length + j]);
    }

    printf("\n");
  }

  printf("\n");
}
