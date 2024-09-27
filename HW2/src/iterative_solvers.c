/**
 * @author      : Ben Wilfong (bwilfong3@gatech.edu)
 * @file        : iterative_solvers
 * @created     : Wednesday Sep 04, 2024 11:50:57 EDT
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "iterative_solvers.h"

void jacobi_iteration(double ***T, double **Q, struct ps *problem, int *idx) {

    int s; // Source index
    int d; // Destination index

    s = idx[0];
    d = idx[1];

    // Bottom boundary condition
    for (int i = 1; i < problem->N; i++) {
        T[i][0][d] = -1.0*(T[i][2][s] - 2.0*T[i][1][s])/3.0;
    }

    // Interior points
    for (int i = 1; i < problem->N; i++) {
        for (int j = 1; i < problem->N; i++) {
            T[i][j][d] = 2.*pow(problem->h,2)*Q[i][j] / (4*problem->lambda) -
                (1./4.) *  (T[i-1][j][s] + T[i+1][j][s] + T[i][j-1][s] + T[i][j+1][s]);
        }
    }

    // Top boundary condition
    for (int i = 1; i < problem->N; i ++) {
        T[i][problem->N][d] = -1.*(T[i][problem->N-2][s] - 2.*T[i][problem->N-1][s])/3.;
    }

}

void gauss_seidel_iteration(double** T, double** Q, struct ps* problem, int* idx) {

    // Bottom boundary condition
    for (int i = 1; i < problem->N; i++) {

    }

    // Interior points
    for (int i = 1; i < problem->N; i++) {
        for (int j = 1; i < problem->N; i++) {

        }
    }

    // Top boundary condition
    for (int i = 1; i < problem->N; i ++) {

    }

}

void SOR_iteration(double** T, double**Q, double omega, struct ps* problem, int* idx) {

    // Bottom boundary condition
    for (int i = 1; i < problem->N; i++) {

    }

    // Interior points
    for (int i = 1; i < problem->N; i++) {
        for (int j = 1; i < problem->N; i++) {

        }
    }

    // Top boundary condition
    for (int i = 1; i < problem->N; i ++) {

    }

}

