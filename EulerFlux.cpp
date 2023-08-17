//
// Created by Tsail on 6/24/2023.
//

#include <cmath>
#include <stdexcept>

#include "EulerFlux.h"

using namespace std;

void getPrimativesPN(const double gam, const double *unkel, double *rho, double *vx, double *vy, double *p,\
                     double *c, double *M) {
    // Gets the primative variables at a specific node
    rho[0]        = unkel[0];
    double rhovx  = unkel[1];
    double rhovy  = unkel[2];
    double rhoe   = unkel[3];

    //Density Limiter
    //rho[0] = fmax(1e-8, rho[0]);

    //break down into primatives
    vx[0] = rhovx / rho[0];
    vy[0] = rhovy / rho[0];
    double v2 = (vx[0]*vx[0]) + (vy[0]*vy[0]);
    p[0] = (gam - 1) * (rhoe - 0.5 * rho[0] * v2);

    //Pressure Limit
    p[0] = fmax(1e-8, p[0]);

    c[0] = sqrt(gam * p[0] / fmax(1e-8, rho[0]));
    M[0] = sqrt(v2) / c[0]; ///=========== might need to sign this??? ================= (maybe not)

    if (isnan(M[0])) {throw overflow_error("Nonpositive Density\n");}
}

void EulerFlux(const int dim, const double gam, const double *u, double* flux){
    flux[0] = 0.0;
    flux[1] = 0.0;
    flux[2] = 0.0;
    flux[3] = 0.0;

    double rho, vx, vy, p, c, M;
    getPrimativesPN(gam, u, &rho, &vx, &vy, &p, &c, &M);

    if (dim == 0) {
        //X flux component
        flux[0] = u[1];
        flux[1] = (u[1] * vx) + p;
        flux[2] = (u[1] * vy);
        flux[3] = vx * (u[3] + p);
    } else {
        // Y flux component
        flux[0] = u[2];
        flux[1] = (u[2] * vx);
        flux[2] = (u[2] * vy) + p;
        flux[3] = vx * (u[3] + p);
    }
}

void RoeFDS(double gam, const double* uL, const double *uR, double* roeFlux){
    /// NOT UPDATED FOR 2D
    double rhoL, vL, pL, cL, ML;
    getPrimativesPN(gam, uL, &rhoL, &vL, &pL, &cL, &ML);
    double hL = (uL[2] + pL) / uL[0];

    double rhoR, vR, pR, cR, MR;
    getPrimativesPN(gam, uR, &rhoR, &vR, &pR, &cR, &MR);
    double hR = (uR[2] + pR) / uR[0];

    //find Roe-averaged state vars
    double R, rho, v, h0, c, v2;
    R = sqrt(rhoR / rhoL);
    rho = R * rhoL;
    v  = (R*vR + vL) / (R+1);
    h0 = (R*hR + hL) / (R+1);
    v2 = v*v;
    c  = sqrt((gam-1.0)*(h0 - 0.5*v2));

    double deltaU = vR - vL;
    double deltaP = pR - pL;

    //find the approx flux jacobian eigenvalues with the Harten correction
    double lam[3], phi[3];
    lam[0] = v;
    lam[1] = v + c;
    lam[2] = v - c;

    double eps = 1e-6;//deltaU;

    for (int i=0; i<3; i++){
        if (fabs(lam[i]) < eps){
            phi[i] = 0.5*(lam[i]*lam[i]/eps + eps);//(lam[i]*lam[i] + eps*eps) / (2*eps);
        } else {
            phi[i] = fabs(lam[i]);
        }
    }

    //transformed conservative vars
    double dw1, dw2, dw3;
    dw1 = (rhoR-rhoL) - (deltaP / (c*c));
    dw2 = deltaU + (deltaP / (rho*c));
    dw3 = deltaU - (deltaP / (rho*c));

    //eigenvectors of approx jac matrix
    double e1[3], e2[3], e3[3];

    e1[0] = 1;
    e1[1] = v;
    e1[2] = 0.5*v2;

    double coeff = 0.5*rho/c;
    e2[0] = coeff;
    e2[1] = coeff*(v + c);
    e2[2] = coeff*(h0 + (v*c));

    coeff = -coeff;
    e3[0] = coeff;
    e3[1] = coeff * (v - c);
    e3[2] = coeff * (h0 - (v*c));

    //Calculate left and right flux states
    /*auto*fluxLeft = (double*)malloc(3*sizeof(double));
    auto*fluxRight = (double*)malloc(3*sizeof(double));
    EulerFlux(gam, uL, fluxLeft);
    EulerFlux(gam, uR, fluxRight);*/
    double fluxLeft[3], fluxRight[3];
    fluxLeft[0]  = uL[1];
    fluxRight[0] = uR[1];

    fluxLeft[1]  = (((gam - 1) / gam) * rhoL * hL)  +  ((0.5 * (gam + 1) / gam) * uL[1] * vL);
    fluxRight[1] = (((gam - 1) / gam) * rhoR * hR)  +  ((0.5 * (gam + 1) / gam) * uR[1] * vR);

    fluxLeft[2]  = uL[1]*hL;
    fluxRight[2] = uR[1]*hR;

    //Compile the final Roe flux


    roeFlux[0] = 0.5 * (fluxRight[0] + fluxLeft[0] - (e1[0]*phi[0]*dw1) - (e2[0]*phi[1]*dw2) - (e3[0]*phi[2]*dw3));
    roeFlux[1] = 0.5 * (fluxRight[1] + fluxLeft[1] - (e1[1]*phi[0]*dw1) - (e2[1]*phi[1]*dw2) - (e3[1]*phi[2]*dw3));
    roeFlux[2] = 0.5 * (fluxRight[2] + fluxLeft[2] - (e1[2]*phi[0]*dw1) - (e2[2]*phi[1]*dw2) - (e3[2]*phi[2]*dw3));

}

double F1pm(const int isPlus, const double M, const double rho, const double c){

    if (isPlus==1){
        if (M<=-1.0) {
            //F1+
            return 0.0;
        }
        if (M>=1.0) {
            return rho * M * c;
        }
        return 0.25*rho*c*(M+1)*(M+1);
    }

    if (isPlus==0) {
        if (M <= -1.0) {
            //F1-
            return rho * M * c;// *-1
        }
        if (M >= 1.0) {
            return 0.0;
        }
        return -0.25 * rho * c * (M - 1) * (M - 1);
    }

    return NAN;
}

void LeerFlux(const double gam, const double* uL, const double* uR, const int nvec, double *flux){
    /*
     * gam - ratio of specific heats
     * uL - Left state variables
     * uR - Right state variables
     * nvec - flux direction (0=x+, 1=y+)
     * flux - output common flux in the given direction
     */
    double F1L, F1R, fPlus[4], fMins[4], nx, ny;
    double rhoL, vxL, vyL, pL, cL, ML, vnL, MnL;
    double rhoR, vxR, vyR, pR, cR, MR, vnR, MnR;

    nx = (1.0 - (double)nvec);
    ny = (double)nvec;

    flux[0] = 0.0;
    flux[1] = 0.0;
    flux[2] = 0.0;
    flux[4] = 0.0;

    getPrimativesPN(gam, uL, &rhoL, &vxL, &vyL, &pL, &cL, &ML);
    getPrimativesPN(gam, uR, &rhoR, &vxR, &vyR, &pR, &cR, &MR);

    ///Find normal Mach Number, Mn, and pressure direction

    vnL = nx*vxL + ny*vyL;
    MnL = vnL / cL;

    vnR = nx*vxR + ny*vyR;
    MnR = vnR / cR;

    F1L = F1pm(1, MnL, rhoL, cL);
    F1R = F1pm(0, MnR, rhoR, cR);

    if(isnan(F1L+F1R)){
        throw overflow_error("getting NAN mass flux!");
    }

    //Compute the positive fluxes (left of face)
    if (MnL > 1.0) {
        fPlus[0] = rhoL*vnL;
        fPlus[1] = rhoL*vnL*vxL + pL*nx;
        fPlus[2] = rhoL*vnL*vxL + pL*ny;
        fPlus[3] = vnL*(uL[3] + pL);
    } else if (MnL < -1.0 ){
        fPlus[0] = 0;
        fPlus[1] = 0;
        fPlus[2] = 0;
        fPlus[3] = 0;
    } else {
        fPlus[0] = F1L;
        fPlus[1] = F1L * (vxL + nx*(2*cL-vnL)/gam);
        fPlus[2] = F1L * (vyL + ny*(2*cL-vnL)/gam);
        double A = ((gam-1) * vnL) + 2*cL;
        fPlus[3] = F1L * ( 0.5*(vxL*vxL + vyL*vyL - vnL*vnL) + A*A*0.5/(gam*gam - 1.0)) ;
    }

    //Compute the negative fluxes (right of face)
    if (MnR > 1.0) {
        fMins[0] = 0;
        fMins[1] = 0;
        fMins[2] = 0;
        fMins[3] = 0;
    } else if (MnR < -1.0 ){
        fMins[0] = rhoR*vnR;
        fMins[1] = rhoR*vnR*vxR + pR*nx;
        fMins[2] = rhoR*vnR*vxR + pR*ny;
        fMins[3] = vnR*(uR[3] + pR);
    } else {
        fMins[0] = F1R;
        fMins[1] = F1R * (vxR + nx*(-2*cR-vnR)/gam);
        fMins[2] = F1R * (vyR + ny*(-2*cR-vnR)/gam);
        double A = ((gam-1) * vnR) - 2*cL;
        fMins[3] = F1R * ( 0.5*(vxR*vxR + vyR*vyR - vnR*vnR) + A*A*0.5/(gam*gam - 1.0)) ;
    }


    flux[0] = fPlus[0] + fMins[0];
    flux[1] = fPlus[1] + fMins[1];
    flux[2] = fPlus[2] + fMins[2];
    flux[3] = fPlus[3] + fMins[3];
}