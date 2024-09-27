/**
 * @author      : Ben Wilfong (bwilfong3@gatech.edu)
 * @file        : main
 * @created     : Wednesday Sep 04, 2024 11:50:47 EDT
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "helpers.h"
#include "derived_types.h"
#include "iterative_solvers.h"

int main(int argc, char **argv)
{

    // Get user input
    struct ps problem;

    read_user_input(argc, argv, &problem);

    // Initialize arrays

    double*** T; // Temperature array (N + 1 x N + 1 x 2)
    double** Q; // Source term array (N + 1 x N + 1)
    printf("T = %p\n", (void *)T);
    int idx[] = {0, 1, 2};
    int n_iter = 0;
    int max_iter = 100;
    double z[] = {0, 0, 0};

    double error = 100*problem.tol;

    initialize_variables(&T, &Q, &problem);

    printf("%f, %f, %f\n", T[0][0][0], T[0][1][0], T[1][1][0]);

    printf("T = %p\n", (void *)T);

    switch (problem.solver) {
        case 0:

            while (1) {

                jacobi_iteration(T, Q, &problem, idx);

                estimate_error(T, &problem, idx, error, z);

                left_cycle(idx); // Swap source and destination indices

                n_iter += 1;

                if (n_iter > max_iter) {
                    printf("\033[0;31m"); //Set the text to the color red.
                    printf("Iterative solver failed to converge.\n");
                    printf("\033[0m");
                    exit(136);
                }
            }

            break;
        case 1:

            while (1) {

                left_cycle(idx); // Swap source and destination indices

                n_iter += 1;

                if (n_iter > max_iter) {
                    printf("\033[0;31m"); //Set the text to the color red.
                    printf("Iterative solver failed to converge.\n");
                    printf("\033[0m");
                    exit(136);
                }
            }

            break;
        case 2:

            while (1) {

                left_cycle(idx); // Swap source and destination indices

                n_iter += 1;

                if (n_iter > max_iter) {
                    printf("\033[0;31m"); //Set the text to the color red.
                    printf("Iterative solver failed to converge.\n");
                    printf("\033[0m");
                    exit(136);
                }
            }

            break;
        default:
            printf("\033[0;31m"); //Set the text to the color red.
            printf("Invalid solver provided.\n");
            printf("\033[0m");
            exit(136);
    }


    return 0;
}

