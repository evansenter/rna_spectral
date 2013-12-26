#ifndef SPECTRAL_DATA_STRUCTURES_H
#define SPECTRAL_DATA_STRUCTURES_H

#include "vienna/data_structures.h"

typedef struct {
  short verbose;
  char* sequence;
  char* energy_grid_file;
  char* start_structure;
  char* end_structure;
  double temperature;
  double start_time;
  double end_time;
  double step_size;
  short lonely_bp;
  short energy_cap;
  short eigen_only;
  short benchmark;
  short use_min;
} SPECTRAL_PARAMS;

typedef struct {
  double* values;
  double* vectors;
  double* inverse_vectors;
  int length;
} EIGENSYSTEM;

#endif
