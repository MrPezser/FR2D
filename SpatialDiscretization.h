//
// Created by Tsail on 7/26/2023.
//

#include "indexing.h"
#include "EulerFlux.h"
#include <cmath>

//Oh, this stuff here? just don't worry about it. Used to make sense to be here I guess...
#define GAM 1.4
#define FACEFLUX(a, b, c, d) LeerFlux(GAM, a, b, c, d)
#define FLUX(a,b,c)  EulerFlux(GAM, a, b, c)
//...on second thought probably not... should just stop being lazy and make an input deck

#ifndef FR1D_SPATIALDISCRETIZATION_H
#define FR1D_SPATIALDISCRETIZATION_H

void CalcDudt(int nelem, int ndegr, int nface, int nvar,\
              const int* inelfa, const int* facpts, const double* u, const double* Dmatrix, const double* Dradau,\
              const double* eldrdxi, const double* eldxidr, const double* eljac,\
              double* dudt );
#endif //FR1D_SPATIALDISCRETIZATION_H
