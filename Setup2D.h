//
// Created by Tsail on 7/31/2023.
//

#ifndef FR2D_SETUP2D_H
#define FR2D_SETUP2D_H

#include <cmath>
#include <cstdlib>
#include "indexing.h"

void LoadStruct2Unstruct(int imx, int jmx, int nelem, int nface, int** inpoel, int** inpofa, int** inelfa);
void CalcCoordJacobian(int ndegr, int npoin, int nelem, int* inpoel, double* x_in, double* y_in, double* xi,\
                        double* eta, double* eldrdxi);
void FacePoint2PointMap(int ndegr, int* facpts);


#endif //FR2D_SETUP2D_H
