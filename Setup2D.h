//
// Created by Tsail on 7/31/2023.
//

#ifndef FR2D_SETUP2D_H
#define FR2D_SETUP2D_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "indexing.h"

void LoadStruct2Unstruct(int imx, int jmx, int nelem, int nface, int** inpoel, int** inpofa, int** inelfa);
void CalcCoordJacobian(int order, int npoin, int nelem, const int* inpoel, const double* x_in, const double* y_in,\
                       const double* xi, const double* eta, double* eldrdxi, double* eldxidr, double* eljac);
void FacePoint2PointMap(int ndegr, int* facpts);


#endif //FR2D_SETUP2D_H
