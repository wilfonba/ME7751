/**
 * @author      : Ben Wilfong (bwilfong3@gatech.edu)
 * @file        : iterative_solvers
 * @created     : Wednesday Sep 04, 2024 11:51:04 EDT
 */

#ifndef ITERATIVE_SOLVERS_H

#define ITERATIVE_SOLVERS_H

#include "derived_types.h"

void jacobi_iteration(double*** PT, double** Q, struct ps* problem, int* idx);

void gauss_seidel_iteration(double **T, double** Q, struct ps* problem, int* idx);

void SOR_iteration(double** T, double**Q, double omega, struct ps* problem, int* idx);

#endif /* end of include guard ITERATIVE_SOLVERS_H */

