#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "initializers.h"
#include "functions.h"

EIGENSYSTEM init_eigensystem(int length) {
  EIGENSYSTEM eigensystem = {
    .values          = malloc(length * sizeof(double)),
    .vectors         = malloc(length* length * sizeof(double)),
    .inverse_vectors = malloc(length* length * sizeof(double)),
    .length          = length
  };
  return eigensystem;
}

void free_eigensystem(EIGENSYSTEM eigensystem) {
  free(eigensystem.values);
  free(eigensystem.vectors);
  free(eigensystem.inverse_vectors);
}

void print_eigensystem(EIGENSYSTEM eigensystem) {
  print_array("eigensystem.values", eigensystem.values, eigensystem.length);
  print_matrix("eigensystem.vectors", eigensystem.vectors, eigensystem.length);
  print_matrix("eigensystem.inverse_vectors", eigensystem.inverse_vectors, eigensystem.length);
}
