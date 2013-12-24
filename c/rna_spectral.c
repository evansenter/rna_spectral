#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include "constants.h"
#include "vienna/functions.h"
#include "params.h"
#include "functions.h"
#include "initializers.h"

#define TIMING(start, stop, task) printf("Time in ms for %s: %.2f\n", task, (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0));

extern double temperature;
SPECTRAL_PARAMS parameters;
  
int main(int argc, char* argv[]) {
  #ifdef DEBUG
    struct timeval fullStart, fullStop, start, stop;
    gettimeofday(&fullStart, NULL);
    gettimeofday(&start, NULL);
  #endif
  
  paramT*        vienna_params;
  model_detailsT vienna_details;
  
  set_model_details(&vienna_details);
  parameters          = parse_spectral_args(argc, argv);
  vienna_details.noLP = !parameters.lonely_bp;
  vienna_params       = get_scaled_parameters(temperature, vienna_details);
  
  int i, seq_length, energy_cap, from_index = -1, to_index = -1, num_structures = 0;
  double mfe_energy, step_counter;
  double* transition_matrix;
  EIGENSYSTEM eigensystem;
  
  char* sequence  = parameters.sequence;
  seq_length      = strlen(sequence);
  char* empty_str = malloc((seq_length + 1) * sizeof(char));
  char* mfe_str   = malloc((seq_length + 1) * sizeof(char));
  
  for (i = 0; i < seq_length; ++i) {
    empty_str[i] = '.';
  }
  empty_str[i] = '\0';
  
  mfe_energy = (double)fold_par(sequence, mfe_str, vienna_params, 0, 0);
  energy_cap = parameters.energy_cap ? (int)(2 * abs(mfe_energy * 100)) : 1000000;
  
  SOLUTION* all_structures = subopt_par(sequence, empty_str, vienna_params, energy_cap, 0, 0, NULL);
  while (all_structures[num_structures].structure != NULL) {
    num_structures++;
  }
  
  for (i = 0; i < num_structures; ++i) {
    if (parameters.start_structure != NULL) {
      if (!strcmp(parameters.start_structure, all_structures[i].structure)) {
        from_index = i;
      }
    } else {
      if (!strcmp(empty_str, all_structures[i].structure)) {
        parameters.start_structure = all_structures[i].structure;
        from_index                 = i;
      }
    }
    
    if (parameters.end_structure != NULL) {
      if (!strcmp(parameters.end_structure, all_structures[i].structure)) {
        to_index = i;
      }
    } else {
      if (!strcmp(mfe_str, all_structures[i].structure)) {
        parameters.end_structure = all_structures[i].structure;
        to_index                 = i;
      }
    }
  }
  
  if (parameters.verbose) {
    printf("sequence:\t%s\n", sequence);
    printf("start:\t\t%s\t%+.2f kcal/mol\t(%d)\n", parameters.start_structure, all_structures[from_index].energy, from_index);
    printf("stop:\t\t%s\t%+.2f kcal/mol\t(%d)\n", parameters.end_structure, all_structures[to_index].energy, to_index);
    printf("energy cap:\t%+.2f kcal/mol above MFE\n", energy_cap / 100.);
    printf("num str:\t%d\n", num_structures);
  }
  
  #ifdef SUPER_HEAVY_DEBUG
    for (i = 0; i < num_structures; ++i) {
      printf("%s\t%+.2f\t%d\n", all_structures[i].structure, all_structures[i].energy, i);
    }
  
    printf("\n");
  #endif
    
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "initialization")
    gettimeofday(&start, NULL);
  #endif
  
  // if (parameters.energy_grid_file != NULL) {
  //   // Not yet implemented.
  // } else {
    transition_matrix = convert_structures_to_transition_matrix(all_structures, num_structures, parameters.use_min);
  // }
  
  #ifdef INSANE_DEBUG
    print_matrix("transition_matrix", transition_matrix, num_structures);
  #endif
  
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "convert_structures_to_transition_matrix")
    gettimeofday(&start, NULL);
  #endif
  
  eigensystem = convert_transition_matrix_to_eigenvectors(transition_matrix, num_structures);
  
  if (parameters.eigen_only) {
    for (i = 0; i < seq_length; ++i) {
      printf("%+.8f\n", eigensystem.values[i]);
    }
  } else {
    #ifdef DEBUG
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "convert_transition_matrix_to_eigenvectors")
      gettimeofday(&start, NULL);
    #endif

    invert_matrix(eigensystem);

    #ifdef DEBUG
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "invert_matrix")
      gettimeofday(&start, NULL);
    #endif

    #ifdef INSANE_DEBUG
      print_eigensystem(eigensystem);
    #endif

    for (step_counter = parameters.start_time; step_counter <= parameters.end_time; step_counter += parameters.step_size) {
      printf("%f\t%+.8f\n", step_counter, probability_at_time(eigensystem, pow(10, step_counter), from_index, to_index));
    }
  }
  
  free_eigensystem(eigensystem);
  
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "generate_population_points")
    gettimeofday(&fullStop, NULL);
    TIMING(fullStart, fullStop, "total")
  #endif
  
  return 0;
}
