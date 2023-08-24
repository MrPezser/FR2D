//
// Created by Tsail on 8/24/2023.
//

#include <cstdio>
#include "indexing.h"
#include "EulerFlux.h"

#ifndef FR2D_OUTPUT_H
#define FR2D_OUTPUT_H

void printgrid(const char *title, int *inpoel, int nelem, int npoin, double *x, double *y);
void printscalar(const char *title, const char *varname, const char *varname2, const char *varname3,
                 const char *varname4, int *inpoel, int nelem, int npoin, double *x, double* y,
                 double *unkel);

#endif //FR2D_OUTPUT_H
