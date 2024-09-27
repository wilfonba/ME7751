/**
 * @author      : Ben Wilfong (bwilfong3@gatech.edu)
 * @file        : derived_types
 * @created     : Thursday Sep 05, 2024 13:18:50 EDT
 */

#ifndef DERIVED_TYPES_H

#define DERIVED_TYPES_H

struct ps{
    uint64_t N;
    uint8_t solver;
    double* xs;
    double* ys;
    double h;
    double lambda;
    double tol;
};

#endif /* end of include guard DERIVED_TYPES_H */

