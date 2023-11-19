//
// Created by Tsail on 7/26/2023.
//

#include <stdexcept>
#include "SpatialDiscretization.h"


void FluxFaceCorrection(const int nface, const int ndegr, const int nelem,\
                        const double* u, const double* Dradau, const int* inelfa, const int* facpts,\
                        const double* eldrdxi, const double* eldxidr, const double* eljac, double* fcorr_xi){
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


        int facori = inelfa[iu(iface,2,nface)]; //Face Orientation
        int ind;
        double dxdL, dydL, dxdR, dydR, JmagL, JmagR, ddxL, ddxR, ddyL, ddyR;
        double nvec[4];

        //Loop over the fact points
        for (int ifpt=0; ifpt<ndegr; ifpt++){
            if (facori != 0){
                //"Horizontal" face - detad() derivs
                Lpoin = facpts[iu(ifpt, 2, ndegr)];
                Rpoin = facpts[iu(ifpt, 3, ndegr)];
                //Find the transform derivatives
                ind = 2 * iu(Lelem, ifpt + 0, nelem);
                dxdL = eldrdxi[ind + 0];
                dydL = eldrdxi[ind + 1];
                ddxL = eldxidr[ind + 0];
                ddyL = eldxidr[ind + 1];

                ind = 2 * iu(Relem, ifpt + 0, nelem);
                dxdR = eldrdxi[ind + 0];
                dydR = eldrdxi[ind + 1];
                ddxR = eldxidr[ind + 0];
                ddyR = eldxidr[ind + 1];

            } else {
                //"Vertical" face - dxid() derivs
                Lpoin = facpts[iu(ifpt, 0, ndegr)];/////////////
                Rpoin = facpts[iu(ifpt, 1, ndegr)];
                //Find the transform derivatives
                ind = 2 * iu(Lelem, ifpt + ndegr, nelem);
                dxdL = eldrdxi[ind + 0];
                dydL = eldrdxi[ind + 1];
                ddxL = eldxidr[ind + 0];
                ddyL = eldxidr[ind + 1];

                ind = 2 * iu(Relem, ifpt + ndegr, nelem);
                dxdR = eldrdxi[ind + 0];
                dydR = eldrdxi[ind + 1];
                ddxR = eldxidr[ind + 0];
                ddyR = eldxidr[ind + 1];
            }


            JmagL = eljac[iu(Lelem, Lpoin, nelem)];
            JmagR = eljac[iu(Relem, Rpoin, nelem)];


            //direction of the left reference direction in the real plane
            double len = sqrt(dxdL*dxdL + dydL*dydL);
            nvec[0] = dxdL / len;
            nvec[1] = dydL / len;
            //dir of the Right element's reference lines
            len = sqrt(dxdR*dxdR + dydR*dydR);
            nvec[2] = dxdR / len;
            nvec[3] = dydR / len;



            //Calculate the flux in the respective grid directions
            FACEFLUX(&u[iu3(Lelem, Lpoin, 0, tdegr)], &u[iu3(Relem, Rpoin, 0, tdegr)], &nvec[0], &flux_comm_L[0]);
            FACEFLUX(&u[iu3(Lelem, Lpoin, 0, tdegr)], &u[iu3(Relem, Rpoin, 0, tdegr)], &nvec[2], &flux_comm_R[0]);

            for (int kvar=0; kvar< NVAR; kvar++){
                ///THIS IS PROBABLY WRONG AND I NEED TO CHANGE IT BUT I AM GOING TO SEE IF IT WORKS FOR THIS SIMPLIFIED CASE
                flux_comm_L[kvar] *= sqrt(dxdL*dxdL + dydL*dydL);
                flux_comm_R[kvar] *= sqrt(dxdR*dxdR + dydR*dydR);
            }


            //Calculate the element's local fluxes at the face
            FLUX(0, &u[iu3(Lelem, Lpoin, 0, tdegr)], &fL[0]);
            FLUX(1, &u[iu3(Lelem, Lpoin, 0, tdegr)], &gL[0]);

            FLUX(0, &u[iu3(Relem, Rpoin, 0, tdegr)], &fR[0]);
            FLUX(1, &u[iu3(Relem, Rpoin, 0, tdegr)], &gR[0]);

            //Rotate the Local Fluxes to be in line with the reference directions
            for (int ivar=0; ivar<NVAR; ivar++){
                local_L[ivar] = JmagL * (ddxL*fL[ivar] + ddyL*gL[ivar]);
                local_R[ivar] = JmagR * (ddxR*fR[ivar] + ddyR*gR[ivar]);

                ASSERT(!_isnan(local_L[ivar]) && !_isnan(local_R[ivar]), "NAN local face fluxes")
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
                        Lpoin2 = iu(ndegr-1-inode, ifpt, ndegr);
                        Rpoin2 = iu(inode, ifpt, ndegr);
                    } else {
                        // Horizontal Face
                        Lpoin2 = iu(ifpt, inode, ndegr);
                        Rpoin2 = iu(ifpt, ndegr-1-inode, ndegr);
                    }

                    ASSERT(Lpoin2 < tdegr && Rpoin2 < tdegr, "you done messed up fam")

                    if (kvar == 2 && (fabs(local_L[kvar] - flux_comm_L[kvar]) > 1e-5 || fabs(local_R[kvar] - flux_comm_R[kvar]) > 1e-5)) {
                    //    printf("thar she blows");
                    }
                    if (bcside == 0) {
                        fcorr_xi[iu3(Lelem, Lpoin2, kvar, tdegr)] += (local_L[kvar] - flux_comm_L[kvar]) * Dradau[inode];
                        fcorr_xi[iu3(Relem, Rpoin2, kvar, tdegr)] -= (local_R[kvar] - flux_comm_R[kvar]) * Dradau[inode];

                    } else if (bcside == -1) {
                        fcorr_xi[iu3(Relem, Rpoin2, kvar, tdegr)] -= (local_R[kvar] - flux_comm_R[kvar]) * Dradau[inode];
                    } else {
                        fcorr_xi[iu3(Lelem, Lpoin2, kvar, tdegr)] += (local_L[kvar] - flux_comm_L[kvar]) * Dradau[inode];
                    }

                }

            }

        }

    }
}

void DiscontinuousFlux(const int nelem, const int ndegr, const double* u, const double* Dmatrix, const double* eldxidr, const double* eljac, double* fcorr_xi) {
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
        for (int i = 0; i < ndegr; i++) {
            for (int j = 0; j < ndegr; j++) {
                int inode = iu(i,j,ndegr);
                double dxidx, dxidy, detady, detadx;
                //xi derivatives
                int ind = 2 * iu(ielem, j + 0, nelem);
                dxidx = eldxidr[ind + 0];
                dxidy = eldxidr[ind + 1];
                //eta derivatives
                ind = 2 * iu(ielem, i + ndegr, nelem);
                detadx = eldxidr[ind + 0];
                detady = eldxidr[ind + 1];

                double Jac = eljac[iu(ielem,inode,nelem)];

                double f_flux[NVAR], g_flux[NVAR];
                FLUX(0, &u[iu3(ielem, inode, 0, tdegr)], &f_flux[0]);
                FLUX(1, &u[iu3(ielem, inode, 0, tdegr)], &g_flux[0]);

                f_node[iu(inode, 0, tdegr)] = Jac * (dxidx*f_flux[0] + dxidy*g_flux[0]);
                f_node[iu(inode, 1, tdegr)] = Jac * (dxidx*f_flux[1] + dxidy*g_flux[1]);
                f_node[iu(inode, 2, tdegr)] = Jac * (dxidx*f_flux[2] + dxidy*g_flux[2]);
                f_node[iu(inode, 3, tdegr)] = Jac * (dxidx*f_flux[3] + dxidy*g_flux[3]);

                g_node[iu(inode, 0, tdegr)] = Jac * (detadx*f_flux[0] + detady*g_flux[0]);
                g_node[iu(inode, 1, tdegr)] = Jac * (detadx*f_flux[1] + detady*g_flux[1]);
                g_node[iu(inode, 2, tdegr)] = Jac * (detadx*f_flux[2] + detady*g_flux[2]);
                g_node[iu(inode, 3, tdegr)] = Jac * (detadx*f_flux[3] + detady*g_flux[3]);
            }
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
                    //Each derivative affects the xi=const and the eta = const values, k is iterating both simultaneously
                    for (int kvar = 0; kvar < NVAR; kvar++) {
                        int vnode = iu(i,k,ndegr);
                        int hnode = iu(k,j,ndegr);

                        //Vertical
                        fxi[kvar] += g_node[iu(vnode, kvar, tdegr)] * Dmatrix[iu(i, k, ndegr)];
                        //Horizontal
                        fxi[kvar] += f_node[iu(hnode, kvar, tdegr)] * Dmatrix[iu(j, k, ndegr)];
                    }
                }

                if (fabs(fxi[2]) > 1e-50){
                    //printf("Thar she blows!");
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
              const int* inelfa, const int* facpts, const double* u, const double* Dmatrix, const double* Dradau,
              const double* eldrdxi, const double* eldxidr, const double* eljac,\
              double* dudt ){
    ///Calculates the solution update given a function to find flux
    int tdegr = ndegr * ndegr;
    int nu    = nelem * tdegr * nvar;

    //double fcorr_xi[nu];//, ucorr_xi[nu], ucorr_x[nu], ucx_xi[nu];
    auto* fcorr_xi = (double*)malloc(nu*sizeof(double));

    //Initialize dudt & corrected flux slope
    for (int i=0;i<nu; i++){
        dudt[i] = 0.0;
        fcorr_xi[i] = 0.0;
    }

    ///Need Boundary Conditions!!! (currently neumann = 0 everywhere)

    //Flux Reconstruction (P_ndegr)
    FluxFaceCorrection(nface, ndegr, nelem, u, Dradau, inelfa, facpts, eldrdxi, eldxidr, eljac, fcorr_xi);
    DiscontinuousFlux(nelem, ndegr, u, Dmatrix, eldxidr, eljac, fcorr_xi);


    for (int ielem=0; ielem<nelem; ielem++) {
        for (int inode=0; inode<tdegr ;inode++) {
                for (int kvar = 0; kvar < nvar; kvar++) {
                    //convert from the reference element to the real element and negate to put make it dudt (= -f_x)

                    double Jacobi = eljac[iu(ielem, inode, nelem)]; //coordinate transform jacobian
                    dudt[iu3(ielem, inode, kvar, tdegr)] = -(1 / Jacobi) * (fcorr_xi[iu3(ielem, inode, kvar, tdegr)]);\
                    //throw std::overflow_error("dude you forgot to delete this error\n");

                    if (_isnan(dudt[iu3(ielem, inode, kvar, tdegr)])) {
                        throw std::overflow_error("try again loser\n dudt NAN\n");
                    }
                }
        }
    }
    free(fcorr_xi);
}