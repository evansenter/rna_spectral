#ifndef SPECTRAL_FUNCTIONS_H
#define SPECTRAL_FUNCTIONS_H

#include "data_structures.h"

double* convert_structures_to_transition_matrix(SOLUTION*, int, int);
EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double*, int);
void invert_matrix(EIGENSYSTEM);
double probability_at_time(EIGENSYSTEM, double, int, int);
void print_array(char*, double*, int);
void print_matrix(char*, double*, int);

#endif