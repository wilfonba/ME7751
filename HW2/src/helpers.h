/**
 * @author      : Ben Wilfong (bwilfong3@gatech.edu)
 * @file        : helpers
 * @created     : Wednesday Sep 04, 2024 11:51:29 EDT
 */

#ifndef HELPERS_H

#define HELPERS_H

#include "derived_types.h"

void left_cycle(int* idx);

void read_user_input(int argc, char** argv, struct ps* problem);

void initialize_variables(double ****PT, double ***Q, struct ps* problem);

void estimate_error(double*** T, struct ps* problem, int* idx, double error, double* z);

#endif /* end of include guard HELPERS_H */

