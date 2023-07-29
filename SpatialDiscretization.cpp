//
// Created by Tsail on 7/26/2023.
//

#include <stdexcept>
#include "SpatialDiscretization.h"

#define MU (5e-4)

void FluxFaceCorrection(int nx, int ndegr, int nvar, const double* u, const double* Dradau, double* fcorr_xi){
    double fL[nvar], fR[nvar], common_flux[nvar];

    //Compute face contributions to the corrected flux gradient
    for (int iface=0; iface<nx; iface++){
        int ielem, iep1;
        ielem = iface;
        if (iface == nx-1) {
            //Periodic boundary condition
            iep1 = iface;//0;
        } else {
            //Interior Cell
            iep1 = iface+1;
        }


        //Calculate the common flux at the face
        FACEFLUX(&u[iu3(ielem, ndegr-1, 0, ndegr)], &u[iu3(iep1, 0, 0, ndegr)], common_flux);
        //common_flux[0] = FACEFLUX(a, u[iu3(iem1, ndegr-1, 0, ndegr)], u[iu3(ielem, 0, 0, ndegr)]);

        //Calculate the element's local fluxes at the face
        FLUX(&u[iu3(ielem , ndegr-1 , 0, ndegr)], &fL[0]);
        FLUX(&u[iu3(iep1, 0       , 0, ndegr)], &fR[0]);
        //fL[0] = FLUX(a, u[iu3(iem1, ndegr-1, 0, ndegr)]);
        //fR[0] = FLUX(a, u[iu3(ielem, 0, 0, ndegr)]);


        for (int inode=0; inode<ndegr; inode++){
            for (int kvar = 0; kvar<nvar; kvar++) {
                fcorr_xi[iu3(ielem,ndegr-1-inode, kvar, ndegr)] -= (common_flux[kvar] - fL[kvar]) * Dradau[inode];   //ndegr-1-
                fcorr_xi[iu3(iep1, inode        , kvar, ndegr)] += (common_flux[kvar] - fR[kvar]) * Dradau[inode];
            }
        }
    }
}

void DiscontinuousFlux(const int nx, const int np, const int ndegr, const int nvar, const double* u, const double* Dmatrix, double* fcorr_xi) {
    double flux_node[np];
    double fxi[nvar];

    for (int inode=0; inode<ndegr; inode++){
        for (int kvar=0; kvar<nvar;kvar++){
            flux_node[iup(inode, kvar, nvar)] = 0.0;
        }
    }

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nx; ielem++) {

        for (int inode = 0; inode < ndegr; inode++) {
            double flux[3];
            FLUX(&u[iu3(ielem, inode, 0, ndegr)], &flux[0]);
            //flux_node[inode][0] = FLUX(a, u[iu3(ielem, inode , 0, ndegr)]);

            flux_node[iup(inode, 0, nvar)] = flux[0];
            flux_node[iup(inode, 1, nvar)] = flux[1];
            flux_node[iup(inode, 2, nvar)] = flux[2];
        }

        //compute discontinuous flux slope (in reference element of width 2)
        for (int inode = 0; inode < ndegr; inode++) {
            for (int kvar = 0; kvar < nvar; kvar++) {
                fxi[kvar] = 0.0;
            }

            //Use the Lagrange polynomial derivative matrix to perform the MatVec -> interior contribution to flux slope
            for (int jnode = 0; jnode < ndegr; jnode++) {
                //if (jnode == inode) continue;
                for (int kvar = 0; kvar < nvar; kvar++) {
                    fxi[kvar] += flux_node[iup(jnode, kvar, nvar)] * Dmatrix[iup(inode, jnode, ndegr)];
                }
            }

            // add the interior contribution to the flux slope
            for (int kvar = 0; kvar < nvar; kvar++) {
                fcorr_xi[iu3(ielem, inode, kvar, ndegr)] += fxi[kvar];
            }
        }
    }

}


void StateFaceCorrection(int nx, int ndegr, int nvar, const double* u, const double* Dradau, double* ucorr_xi){
    double uL[nvar], uR[nvar], common_state[nvar];

    //Compute face contributions to the corrected flux gradient
    for (int iface=0; iface<nx; iface++){
        int ielem, iep1;
        ielem = iface;
        if (iface == nx-1) {
            //Periodic boundary condition
            iep1 = iface;//0;
        } else {
            //Interior Cell
            iep1 = iface+1;
        }

        uL[0] = u[iu3(ielem, ndegr-1, 0, ndegr)];
        uL[1] = u[iu3(ielem, ndegr-1, 1, ndegr)];
        uL[2] = u[iu3(ielem, ndegr-1, 2, ndegr)];

        uR[0] = u[iu3(iep1, 0, 0, ndegr)];
        uR[1] = u[iu3(iep1, 0, 1, ndegr)];
        uR[2] = u[iu3(iep1, 0, 2, ndegr)];

        //Calculate the common state at the face - use centered / Bassi-Rebay
        common_state[0] = 0.5 * (uL[0] + uR[0]);
        common_state[1] = 0.5 * (uL[1] + uR[1]);
        common_state[2] = 0.5 * (uL[2] + uR[2]);

        for (int inode=0; inode<ndegr; inode++){
            for (int kvar = 0; kvar<nvar; kvar++) {
                ucorr_xi[iu3(ielem,ndegr-1-inode, kvar, ndegr)] -= (common_state[kvar] - uL[kvar]) * Dradau[inode];
                ucorr_xi[iu3(iep1, inode        , kvar, ndegr)] += (common_state[kvar] - uR[kvar]) * Dradau[inode];
            }
        }
    }

}

void DiscontinuousState(const int nx, const int np, const int ndegr, const int nvar, const double* u, const double* Dmatrix, double* ucorr_xi){
    double state_node[np];
    double uxi[nvar];

    for (int inode=0; inode<ndegr; inode++){
        for (int kvar=0; kvar<nvar;kvar++){
            state_node[iup(inode, kvar, nvar)] = 0.0;
        }
    }

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nx; ielem++) {

        for (int inode = 0; inode < ndegr; inode++) {
            state_node[iup(inode, 0, nvar)] = u[iu3(ielem, inode, 0, ndegr)];
            state_node[iup(inode, 1, nvar)] = u[iu3(ielem, inode, 1, ndegr)];
            state_node[iup(inode, 2, nvar)] = u[iu3(ielem, inode, 2, ndegr)];
        }

        //compute discontinuous flux slope (in reference element of width 2)
        for (int inode = 0; inode < ndegr; inode++) {
            for (int kvar = 0; kvar < nvar; kvar++) {
                uxi[kvar] = 0.0;
            }

            //Use the Lagrange polynomial derivative matrix to perform the MatVec -> interior contribution to flux slope
            for (int jnode = 0; jnode < ndegr; jnode++) {
                //if (jnode == inode) continue;
                for (int kvar = 0; kvar < nvar; kvar++) {
                    uxi[kvar] += state_node[iup(jnode, kvar, nvar)] * Dmatrix[iup(inode, jnode, ndegr)];
                }
            }

            // add the interior contribution to the flux slope
            for (int kvar = 0; kvar < nvar; kvar++) {
                ucorr_xi[iu3(ielem, inode, kvar, ndegr)] += uxi[kvar];
            }
        }
    }
}


void StateDerivFaceCorrection(int nx, int ndegr, int nvar, const double* ucx, const double* Dradau, double* ucx_xi){
    double uL[nvar], uR[nvar], common_state[nvar];

    //Compute face contributions to the corrected flux gradient
    for (int iface=0; iface<nx; iface++){
        int ielem, iep1;
        ielem = iface;
        if (iface == nx-1) {
            //Periodic boundary condition
            iep1 = iface;//0;
        } else {
            //Interior Cell
            iep1 = iface+1;
        }

        uL[0] = ucx[iu3(ielem, ndegr-1, 0, ndegr)];
        uL[1] = ucx[iu3(ielem, ndegr-1, 1, ndegr)];
        uL[2] = ucx[iu3(ielem, ndegr-1, 2, ndegr)];

        uR[0] = ucx[iu3(iep1, 0, 0, ndegr)];
        uR[1] = ucx[iu3(iep1, 0, 1, ndegr)];
        uR[2] = ucx[iu3(iep1, 0, 2, ndegr)];

        //Calculate the common state at the face - use centered / Bassi-Rebay
        common_state[0] = 0.5 * (uL[0] + uR[0]);
        common_state[1] = 0.5 * (uL[1] + uR[1]);
        common_state[2] = 0.5 * (uL[2] + uR[2]);

        for (int inode=0; inode<ndegr; inode++){
            for (int kvar = 0; kvar<nvar; kvar++) {
                ucx_xi[iu3(ielem,ndegr-1-inode, kvar, ndegr)] -= (common_state[kvar] - uL[kvar]) * Dradau[inode];
                ucx_xi[iu3(iep1, inode        , kvar, ndegr)] += (common_state[kvar] - uR[kvar]) * Dradau[inode];
            }
        }
    }
}

void DiscontinuousStateDeriv(const int nx, const int np, const int ndegr, const int nvar, const double* ucx, const double* Dmatrix, double* ucx_xi){
    double state_node[np];
    double uxi[nvar];

    for (int inode=0; inode<ndegr; inode++){
        for (int kvar=0; kvar<nvar;kvar++){
            state_node[iup(inode, kvar, nvar)] = 0.0;
        }
    }

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nx; ielem++) {

        for (int inode = 0; inode < ndegr; inode++) {
            state_node[iup(inode, 0, nvar)] = ucx[iu3(ielem, inode, 0, ndegr)];
            state_node[iup(inode, 1, nvar)] = ucx[iu3(ielem, inode, 1, ndegr)];
            state_node[iup(inode, 2, nvar)] = ucx[iu3(ielem, inode, 2, ndegr)];
        }

        //compute discontinuous flux slope (in reference element of width 2)
        for (int inode = 0; inode < ndegr; inode++) {
            for (int kvar = 0; kvar < nvar; kvar++) {
                uxi[kvar] = 0.0;
            }

            //Use the Lagrange polynomial derivative matrix to perform the MatVec -> interior contribution to flux slope
            for (int jnode = 0; jnode < ndegr; jnode++) {
                //if (jnode == inode) continue;
                for (int kvar = 0; kvar < nvar; kvar++) {
                    uxi[kvar] += state_node[iup(jnode, kvar, nvar)] * Dmatrix[iup(inode, jnode, ndegr)];
                }
            }

            // add the interior contribution to the flux slope
            for (int kvar = 0; kvar < nvar; kvar++) {
                ucx_xi[iu3(ielem, inode, kvar, ndegr)] += uxi[kvar];
            }
        }
    }
}

void CalcDudt(const int nelem, const int ndegr, const int nvar, const double a, const double dx, const double* u, const double* Dmatrix, const double* Dradau, double* dudt ){
    ///Calculates the solution update given a function to find flux
    int nu    = nelem * ndegr * nvar;
    int tdegr = ndegr * ndegr;
    int np    = tdegr * nvar;
    double fcorr_xi[nu];//, ucorr_xi[nu], ucorr_x[nu], ucx_xi[nu];

    //Initialize dudt & corrected flux slope
    for (int i=0;i<nu; i++){
        dudt[i] = 0.0;
        fcorr_xi[i] = 0.0;
        //ucorr_xi[i] = 0.0;
        //ucx_xi[i] = 0.0;
    }

    ///Need Boundary Conditions!!!

    //Flux Reconstruction (P_ndegr)
    FluxFaceCorrection(nelem, ndegr, nvar, u, Dradau, fcorr_xi);
    DiscontinuousFlux(nelem, np, ndegr, nvar, u, Dmatrix, fcorr_xi);

    /*
    StateFaceCorrection(nx, ndegr, nvar, u, Dradau, ucorr_xi);
    DiscontinuousState(nx, np, ndegr, nvar, u, Dmatrix, ucorr_xi);

    for (int iu=0; iu<nu; iu++) {
        ucorr_x[iu] = ucorr_xi[iu]*(2/dx);
    }

    StateDerivFaceCorrection(nx, ndegr, nvar, ucorr_x, Dradau, ucx_xi);
    DiscontinuousStateDeriv(nx, np, ndegr, nvar, ucorr_x, Dmatrix, ucx_xi);
     */


    for (int ielem=0; ielem<nelem; ielem++) {
        for (int inode=0; inode<ndegr ;inode++) {
            for (int jnode=0; jnode<ndegr; jnode++) {
                for (int kvar = 0; kvar < nvar; kvar++) {

                    int knode = iup(inode, jnode, ndegr);
                    //convert from the reference element to the real element and negate to put make it dudt (= -f_x)
                    dudt[iu3(ielem, knode, kvar, ndegr)] = -(2 / dx) * (fcorr_xi[iu3(ielem, knode, kvar, ndegr)]);
                    // - MU*ucx_xi[iu3(ielem,  inode, kvar, ndegr)]);

                    if (_isnan(dudt[iu3(ielem, knode, kvar, ndegr)])) {
                        throw std::overflow_error("dudt NAN\n");
                    }
                }
            }
        }
    }

}