//
// Created by Tsail on 7/26/2023.
//

#include <stdexcept>
#include "SpatialDiscretization.h"

#define MU (5e-4)

void FluxFaceCorrection(const int nface, const int ndegr, const int nelem,\
                        const double* u, const double* Dradau, const int* inelfa, const int* facpts,\
                        const double* eldxidr, const double* eljac, double* fcorr_xi){
    ///Vertical and horizontal need to be separated into f anf g
    /*
     * Calculates the contribution of the face discontinuities to the flux slope (i.e. the correction term)
     */

    double fL[NVAR], fR[NVAR], gL[NVAR], gR[NVAR], local_L[NVAR], local_R[NVAR], flux_comm_L[NVAR], flux_comm_R[NVAR];
    int tdegr = ndegr*ndegr;

    for (int iface=0; iface<nface; iface++){
        int Lelem, Relem, Lpoin, Rpoin, Lid, Rid;
        Lelem = inelfa[iu(iface,0,nface)];
        Relem = inelfa[iu(iface,1,nface)];

        //Neumann=0 boundary conditions
        int bcside = 0; //says which side of face is boundary
        if (Lelem < 0) {
            Lelem = Relem;
            bcside = -1;
        }
        if (Relem < 0) {
            Relem = Lelem;
            bcside = 1;
        }

        if (Relem) throw "somthang messed up";


        int facori = inelfa[iu(iface,2,nface)]; //Face Orientation
        int offset, ind;
        double ddxL, ddyL, ddxR, ddyR, JmagL, JmagR;
        double nvec[4];

        //Loop over the fact points
        for (int ifpt=0; ifpt<ndegr; ifpt++){
            if (facori == 0){
                //"Vertical" face - dxid() derivs
                Lpoin = facpts[iu(ifpt, 0, ndegr)];
                Rpoin = facpts[iu(ifpt, 1, ndegr)];

                //Find the transform derivatives
                ind = 2 * iu(Lelem, ifpt + 0, nelem);
                ddxL = eldxidr[ind + 0];
                ddyL = eldxidr[ind + 1];
                ind = 2 * iu(Relem, ifpt + 0, nelem);
                ddxR = eldxidr[ind + 0];
                ddyR = eldxidr[ind + 1];


            } else {
                //"Horizontal" face - detad() derivs
                Lpoin = facpts[iu(ifpt, 2, ndegr)];
                Rpoin = facpts[iu(ifpt, 3, ndegr)];

                //Find the transform derivatives
                ind = 2 * iu(Lelem, ifpt + ndegr, nelem);
                ddxL = eldxidr[ind + 0];
                ddyL = eldxidr[ind + 1];
                ind = 2 * iu(Relem, ifpt + ndegr, nelem);
                ddxR = eldxidr[ind + 0];
                ddyR = eldxidr[ind + 1];

            }


            JmagL = eljac[iu(Lelem, Lpoin, nelem)];
            JmagR = eljac[iu(Relem, Rpoin, nelem)];

            //direction of the left reference direction in the real plane
            double len = sqrt(ddxL*ddxL + ddyL*ddyL);
            nvec[0] = ddxL / len;
            nvec[1] = ddyL / len;
            //dir of the Right element's reference lines
            len = sqrt(ddxR*ddxR + ddyR*ddyR);
            nvec[2] = ddxR / len;
            nvec[3] = ddyR / len;

            //Calculate the flux in the respective grid directions
            FACEFLUX(&u[iu3(Lelem, Lpoin, 0, ndegr)], &u[iu3(Relem, Rpoin, 0, ndegr)], &nvec[0], flux_comm_L);
            FACEFLUX(&u[iu3(Lelem, Lpoin, 0, ndegr)], &u[iu3(Relem, Rpoin, 0, ndegr)], &nvec[2], flux_comm_R);



            //Calculate the element's local fluxes at the face
            FLUX(0, &u[iu3(Lelem, Lpoin, 0, ndegr)], &fL[0]);
            FLUX(1, &u[iu3(Lelem, Lpoin, 0, ndegr)], &gL[0]);

            FLUX(0, &u[iu3(Relem, Rpoin, 0, ndegr)], &fR[0]);
            FLUX(1, &u[iu3(Relem, Rpoin, 0, ndegr)], &gR[0]);

            //Rotate the Local Fluxes to be in line with the reference directions
            for (int ivar=0; ivar<NVAR; ivar++){
                local_L[ivar] = JmagL * (ddxL*fL[ivar] + ddyL*gL[ivar]);
                local_R[ivar] = JmagR * (ddxR*fR[ivar] + ddyR*gR[ivar]);
            }


            //Add contributions to the corrected flux
            for (int inode=0; inode<ndegr; inode++){
                for (int kvar = 0; kvar<NVAR; kvar++) {
                    ///Figure out the new node mapping... omg this is horrendous
                    //each face point contributes to the perpendicular set of solution points
                    //need information on those points -
                    //Points indexed; ipoin = iu(i_xi, j_eta, ni)
                    //For Vertical Face
                    //  Lfacpt: (xi_max,  eta=ifpt)
                    //  Rfacpt: (xi=0,    eta=ifpt)
                    //  Lpoin:  (xi=nnode - inode,eta=ifpt)
                    //  Rpoin:  (xi=inode, eta=ifpt)
                    int Lpoin2, Rpoin2;
                    if (facori == 0) {
                        Lpoin2 = iu(ndegr - inode, ifpt, ndegr);
                        Rpoin2 = iu(inode, ifpt, ndegr);
                    } else {
                        // Horizontal Face
                        Lpoin2 = iu(ifpt, inode, ndegr);
                        Rpoin2 = iu(ifpt, ndegr - inode, ndegr);
                    }

                    if (bcside == 0) {
                        fcorr_xi[iu3(Lelem, Lpoin2, kvar, tdegr)] -=
                                (local_L[kvar] - flux_comm_L[kvar]) * Dradau[inode];
                        fcorr_xi[iu3(Relem, Rpoin2, kvar, tdegr)] +=
                                (local_R[kvar] - flux_comm_R[kvar]) * Dradau[inode];
                    } else if (bcside == -1) {
                        fcorr_xi[iu3(Relem, Rpoin2, kvar, tdegr)] +=
                                (local_R[kvar] - flux_comm_R[kvar]) * Dradau[inode];
                    } else {
                        fcorr_xi[iu3(Lelem, Lpoin2, kvar, tdegr)] -=
                                (local_L[kvar] - flux_comm_L[kvar]) * Dradau[inode];
                    }
                }
            }

        }

    }
}

void DiscontinuousFlux(const int nelem, const int ndegr, const double* u, const double* Dmatrix, double* fcorr_xi) {
    ///Vertical and horizontal need to be separated into f anf g
    double fxi[NVAR], geta[NVAR];
    int tdegr= ndegr*ndegr;
    double f_node[NVAR*tdegr], g_node[NVAR*tdegr];

    for (int inode=0; inode<tdegr; inode++){
        for (int kvar=0; kvar<NVAR;kvar++){
            f_node[iu(inode, kvar, tdegr)] = 0.0; //X fluxes
            g_node[iu(inode, kvar, tdegr)] = 0.0; //Y fluxes
        }
    }

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nelem; ielem++) {
        for (int inode = 0; inode < tdegr; inode++) {
            double flux[2*NVAR];
            FLUX(0,&u[iu3(ielem, inode, 0, tdegr)], &flux[0]);
            FLUX(1,&u[iu3(ielem, inode, 0, tdegr)], &flux[NVAR]);

            f_node[iu(inode, 0, tdegr)] = flux[0];
            f_node[iu(inode, 1, tdegr)] = flux[1];
            f_node[iu(inode, 2, tdegr)] = flux[2];
            f_node[iu(inode, 3, tdegr)] = flux[3];

            g_node[iu(inode, 0, tdegr)] = flux[4];
            g_node[iu(inode, 1, tdegr)] = flux[5];
            g_node[iu(inode, 2, tdegr)] = flux[6];
            g_node[iu(inode, 3, tdegr)] = flux[7];
        }

        //compute discontinuous flux slope (in reference element of width 2x2)
        for (int i = 0; i < ndegr; i++) {
            for (int j = 0; j < ndegr; j++) {
                int inode = iu(i,j,ndegr);

                for (int kvar = 0; kvar<NVAR; kvar++) {
                    fxi[kvar] = 0.0;
                    geta[kvar] = 0.0;
                }

                //Use the Lagrange polynomial derivative matrix to perform the MatVec -> interior contribution to flux slope
                for (int k = 0; k < ndegr; k++) {
                    //Each derivative affects the xi=const and the eta = const values, jnode is iterating both simultaneously
                    for (int kvar = 0; kvar < NVAR; kvar++) {
                        int vnode = iu(i,k,ndegr);
                        int hnode = iu(k,j,ndegr);

                        //Vertical
                        fxi[kvar] += g_node[iu(vnode, kvar, tdegr)] * Dmatrix[iu(i, k, tdegr)];
                        //Horizontal
                        fxi[kvar] += f_node[iu(hnode, kvar, tdegr)] * Dmatrix[iu(j, k, tdegr)];

                        //might need to remove duplicate?? maybe not    probably not
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

void CalcDudt(const int nelem, const int ndegr, const int nface, const int nvar,\
              const int* inelfa, const int* facpts, const double* u, const double* Dmatrix, const double* Dradau,\
              const double* eldxidr, const double* eljac,\
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

    ///Need Boundary Conditions!!! (currently neumann = 0 everywhere)

    //Flux Reconstruction (P_ndegr)
    FluxFaceCorrection(nface, ndegr, nelem, u, Dradau, inelfa, facpts, eldxidr, eljac, fcorr_xi);
    DiscontinuousFlux(nelem, ndegr, u, Dmatrix, fcorr_xi);


    for (int ielem=0; ielem<nelem; ielem++) {
        for (int inode=0; inode<tdegr ;inode++) {
                for (int kvar = 0; kvar < nvar; kvar++) {
                    //convert from the reference element to the real element and negate to put make it dudt (= -f_x)

                    double Jacobi = eljac[iu(ielem, inode, nelem)]; //coordinate transform jacobian
                    dudt[iu3(ielem, inode, kvar, ndegr)] = -(1 / Jacobi) * (fcorr_xi[iu3(ielem, inode, kvar, ndegr)]);\
                    throw std::overflow_error("dude you forgot to delete this error\n");

                    if (_isnan(dudt[iu3(ielem, inode, kvar, ndegr)])) {
                        throw std::overflow_error("try again loser\n dudt NAN\n");
                    }
                }
        }
    }

}