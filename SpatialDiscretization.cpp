//
// Created by Tsail on 7/26/2023.
//

#include <stdexcept>
#include "SpatialDiscretization.h"

#define MU (5e-4)

void FluxFaceCorrection(const int nface, const int ndegr, const int nvar,\
                        const double* u, const double* Dradau, const int* inelfa, const int* facpts,\
                        double* fcorr_xi){
    ///Vertical and horizontal need to be separated into f anf g
    /*
     * Calculates the contribution of the face discontinuities to the flux slope
     */

    double fL[nvar], fR[nvar], common_flux[nvar];

    //Compute face contributions to the corrected flux gradient
    for (int iface=0; iface<nface; iface++){
        int Lelem, Relem, Lpoin, Rpoin, id, id2;
        Lelem = inelfa[iu(iface,0,nface)];
        Relem = inelfa[iu(iface,1,nface)];

        //Neumann-0 boundary conditions
        if (Lelem == -1) Lelem = Relem;
        if (Relem == -1) Relem = Lelem;

        if (Relem) throw "somthang messed up";


        int facori = inelfa[iu(iface,2,nface)]; //Face Orientation

        //Loop over the fact points
        for (int ifpt=0; ifpt<ndegr; ifpt++){
            if (facori == 0){
                //Vertical face
                Lpoin = facpts[iu(ifpt, 0, ndegr)];
                Rpoin = facpts[iu(ifpt, 1, ndegr)];

                id  = ndegr * iu(ifpt, 0, ndegr);
                id2 = ndegr * iu(ifpt, 1, ndegr);
            } else {
                //Horizontal face
                Lpoin = facpts[iu(ifpt, 2, ndegr)];
                Rpoin = facpts[iu(ifpt, 3, ndegr)];

                id  = ndegr * iu(ifpt, 2, ndegr);
                id2 = ndegr * iu(ifpt, 3, ndegr);
            }

            //Calculate the common flux at the face
            FACEFLUX(&u[iu3(Lelem, Lpoin, 0, ndegr)], &u[iu3(Relem, Rpoin, 0, ndegr)], common_flux);

            //Calculate the element's local fluxes at the face
            FLUX(&u[iu3(Lelem, Lpoin, 0, ndegr)], &fL[0]);
            FLUX(&u[iu3(Relem, Rpoin, 0, ndegr)], &fR[0]);

            for (int inode=0; inode<ndegr; inode++){
                for (int kvar = 0; kvar<nvar; kvar++) {
                    ///Figure out the new node mapping... omg this is horrendous
                    //each face point contributes to the perpendicular set of solution points
                    //need information on those points -

                    int Lpoin2 = id + inode;
                    int Rpoin2 = id2 + inode;
                    int tdegr = ndegr*ndegr;
                    fcorr_xi[iu3(Lelem,Lpoin2, kvar, tdegr)] -= (common_flux[kvar] - fL[kvar]) * Dradau[inode];
                    fcorr_xi[iu3(Relem,Rpoin2, kvar, tdegr)] += (common_flux[kvar] - fR[kvar]) * Dradau[inode];
                }
            }

        }

    }
}

void DiscontinuousFlux(const int nelem, const int ndegr, const double* u, const double* Dmatrix, double* fcorr_xi) {
    ///Vertical and horizontal need to be separated into f anf g
    double fxi[NVAR];
    int tdegr= ndegr*ndegr;
    double flux_node[NVAR*tdegr];

    for (int inode=0; inode<tdegr; inode++){
        for (int kvar=0; kvar<NVAR;kvar++){
            flux_node[iu(inode, kvar, tdegr)] = 0.0;
        }
    }

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nelem; ielem++) {

        for (int inode = 0; inode < tdegr; inode++) {
            double flux[NVAR];
            FLUX(&u[iu3(ielem, inode, 0, ndegr)], &flux[0]);

            flux_node[iu(inode, 0, tdegr)] = flux[0];
            flux_node[iu(inode, 1, tdegr)] = flux[1];
            flux_node[iu(inode, 2, tdegr)] = flux[2];
            flux_node[iu(inode, 3, tdegr)] = flux[3];
        }

        //compute discontinuous flux slope (in reference element of width 2)
        for (int i = 0; i < ndegr; i++) {
            for (int j = 0; j < ndegr; j++) {
                int inode = iu(i,j,ndegr);

                for (double &kvar: fxi) {   // No idea what this means but Clang-Tidy told me to do it
                    kvar = 0.0;
                }

                //Use the Lagrange polynomial derivative matrix to perform the MatVec -> interior contribution to flux slope
                for (int k = 0; k < ndegr; k++) {
                    //Each derivative affects the xi=const and the eta = const values, jnode is iterating both simultaneously


                    for (int kvar = 0; kvar < NVAR; kvar++) {
                        int vnode = iu(i,k,ndegr);
                        int hnode = iu(k,j,ndegr);

                        //Vertical
                        fxi[kvar] += flux_node[iu(vnode, kvar, tdegr)] * Dmatrix[iu(i, k, tdegr)];
                        //Horizontal
                        fxi[kvar] += flux_node[iu(hnode, kvar, tdegr)] * Dmatrix[iu(j, k, tdegr)];

                        //might need to remove duplicate?? maybe not
                    }
                }

                // add the interior contribution to the flux slope
                for (int kvar = 0; kvar < NVAR; kvar++) {
                    fcorr_xi[iu3(ielem, inode, kvar, tdegr)] += fxi[kvar];
                }
            }
        }
    }

}

void CalcDudt(const int nelem, const int ndegr, const int nface, const int nvar, const double dx,\
              const int* inelfa, const int* facpts, const double* u, const double* Dmatrix, const double* Dradau,\
              double* dudt ){
    ///Calculates the solution update given a function to find flux
    int nu    = nelem * ndegr * nvar;
    int tdegr = ndegr * ndegr;
    double fcorr_xi[nu];//, ucorr_xi[nu], ucorr_x[nu], ucx_xi[nu];

    //Initialize dudt & corrected flux slope
    for (int i=0;i<nu; i++){
        dudt[i] = 0.0;
        fcorr_xi[i] = 0.0;
    }
    ///Need to add the second direction flux derivative, gcorr

    ///Need Boundary Conditions!!!

    //Flux Reconstruction (P_ndegr)
    FluxFaceCorrection(nface, ndegr, NVAR, u, Dradau, inelfa, facpts, fcorr_xi);
    DiscontinuousFlux(nelem, ndegr, u, Dmatrix, fcorr_xi);


    for (int ielem=0; ielem<nelem; ielem++) {
        for (int inode=0; inode<tdegr ;inode++) {
                for (int kvar = 0; kvar < nvar; kvar++) {
                    //convert from the reference element to the real element and negate to put make it dudt (= -f_x)
                    /// need a little linear algebra to get from f_xi, g_eta to f_x and g_y
                    dudt[iu3(ielem, inode, kvar, ndegr)] = -(2 / dx) * (fcorr_xi[iu3(ielem, inode, kvar, ndegr)]);\
                    throw std::overflow_error("dude you forgot to impliment the coordinate transform out of reference coords\n");

                    if (_isnan(dudt[iu3(ielem, inode, kvar, ndegr)])) {
                        throw std::overflow_error("try again loser\n dudt NAN\n");
                    }
                }
        }
    }

}