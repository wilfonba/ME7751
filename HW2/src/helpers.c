/**
 * @author      : Ben Wilfong (bwilfong3@gatech.edu)
 * @file        : helpers
 * @created     : Wednesday Sep 04, 2024 11:51:23 EDT
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "helpers.h"

void left_cycle(int* idx) {
    int temp = idx[0];
    idx[0] = idx[1];
    idx[1] = idx[2];
    idx[2] = temp;
}

void read_user_input(int argc, char** argv, struct ps* problem) {

    uint64_t i;
    uint64_t dummyi;
    double dummyf;

    problem->tol = 1e-3; // Default convergence tolerance
    problem->solver = -100; // Default solver

    for (i=1;i<argc;i++) {
        if (sscanf(argv[i],"-N=%lli",&dummyi)==1) {
            problem->N = dummyi;
        } else if (sscanf(argv[i],"--solver=%lli",&dummyi)==1) {
            problem->solver = dummyi;
        } else if (sscanf(argv[i],"--solver=%lli",&dummyi)==1) {
            problem->tol = dummyf;
        }
    }

    problem->lambda = 0.; // Default from assignment

    printf("Computing a solution to the 2D heat equation with the following parameters:\n");
    printf("N = %llu\n", problem->N);
    printf("solver = %hhu\n", problem->solver);
    fflush(stdout);
}

void initialize_variables(double ****PT, double ***PQ, struct ps* problem) {

    double ***T;
    double **Q;

    // Calculating grid points
    problem->h = 1/problem->N;

    problem->xs = (double *) malloc((problem->N + 1) * sizeof(double));
    problem->ys = (double *) malloc((problem->N + 1) * sizeof(double));

    for (int i = 0; i <= problem->N; i++) {
        problem->xs[i] = i*problem->h;
        problem->ys[i] = i*problem->h;
    }

    // Allocate array to store temperature
    T = (double ***) malloc((problem->N + 1) * sizeof(double **));

    for (int i = 0; i <= problem->N; i++) {
        T[i] = (double **) malloc((problem->N + 1) * sizeof(double*));
    }

    for (int i = 0; i <= problem->N; i++) {
        for (int j = 0; j <= problem->N; j ++) {
            T[i][j] = (double *) malloc(3 * sizeof(double));
        }
    }

    // Allocate array to store Q
    Q = (double **) malloc((problem->N + 1) * sizeof(double *));

    for (int i = 0; i <= problem->N; i++) {
        Q[i] = (double *) malloc((problem->N + 1) * sizeof(double));
    }

    // Apply initial guess and compute Q
    for (int i = 0; i <= problem->N; i ++) {
        for (int j = 0; j <= problem->N; j ++) {
            Q[i][j] = 4*pow(problem->ys[i],3) - 6*pow(problem->ys[i],2) + 2
                - 6*(1 - pow(problem->xs[j],2))*(3*problem->ys[i] - 1);
            T[i][j][0] = Q[i][j];
        }
    }

    *PT = T;
    *PQ = Q;

}

void estimate_error(double*** T, struct ps* problem, int* idx, double error, double* z) {

    int s;
    int d;
    int o;

    s = idx[0];
    d = idx[1];
    o = idx[2];

    double l;
    double delta;

    for (int i = 0; i <= problem->N; i++) {
        for (int j = 0; j <= problem->N; j++) {
            z[d] += T[i][j][o] * T[i][j][d] - T[i][j][s] * T[i][j][s];
            delta += pow(T[i][j][d] - T[i][j][s], 2);
        }
    }

    delta = sqrt(delta);

    l = sqrt(z[d]/z[s]);

    error = delta/sqrt(pow(l,2) + 1);

}
