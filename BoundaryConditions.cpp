//
// Created by Tsail on 8/25/2023.
//
/*

#include "BoundaryConditions.h"

void SubsonInflo(double gam, double vxint, double vyint, double cint, const double *unkel0, double nx, double ny,
                 double *rhogst, double *rhovxgst, double *rhovygst, double *rhoegst) {
    //INPUT = int and fs values
    //OUTPUT = ghost values

    //Get the free stream in primative variables
    double rhoinfty, vxinfty, vyinfty, pinfty, cinfty, Minfty;
    getPrimativesP0(gam, unkel0, 0, &rhoinfty, &vxinfty, &vyinfty, &pinfty, &cinfty, &Minfty);

    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (vxinfty * nx) + (vyinfty * ny);
    double vinftyTANx = vxinfty - vinftyDOTn * nx;
    double vinftyTANy = vyinfty - vinftyDOTn * ny;

    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (gam - 1)) * (cint - cinfty);

    //c of ghost cell
    double cgst = cinfty + 0.5 * (gam - 1) * (vgstDOTn - vinftyDOTn);

    //rho of ghost cell
    rhogst[0] = rhoinfty * pow(cgst * cgst / (cinfty * cinfty), (gam - 1));

    //p  of ghost cell
    double pgst = pinfty * pow(rhogst[0] / rhoinfty, gam);

    //finish by finding the conserved variables
    rhovxgst[0] = rhoinfty * (vinftyTANx + vgstDOTn * nx);
    rhovygst[0] = rhoinfty * (vinftyTANy + vgstDOTn * ny);
    rhoegst[0] = (pgst / (gam - 1)) + 0.5 * (rhovxgst[0] * rhovxgst[0] + rhovygst[0] * rhovygst[0]) / rhogst[0];
}

void SubsonOutfl(double gam, double rhoint, double pint, double vxint, double vyint, double cint, const double *unkel0,
                 double nx, double ny,
                 double *rhogst, double *rhovxgst, double *rhovygst, double *rhoegst) {
    //INPUT = int and fs values
    //OUTPUT = ghost values

    //Get the free stream in primative variables
    double rhoinfty, vxinfty, vyinfty, pinfty, cinfty, Minfty;
    getPrimativesP0(gam, unkel0, 0, &rhoinfty, &vxinfty, &vyinfty, &pinfty, &cinfty, &Minfty);

    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (vxinfty * nx) + (vyinfty * ny);
    double vintTANx = vxint - vintDOTn * nx;
    double vintTANy = vyint - vintDOTn * ny;

    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (gam - 1)) * (cint - cinfty);

    //c of ghost cell
    double cgst = cinfty + 0.5 * (gam - 1) * (vgstDOTn - vinftyDOTn);

    //rho of ghost cell
    rhogst[0] = rhoint * pow(cgst * cgst / (cint * cint), (gam - 1));

    //p  of ghost cell
    double pgst = pint * pow(rhogst[0] / rhoint, gam);

    //finish by finding the rest of the conserved variables
    rhovxgst[0] = rhoinfty * (vintTANx + vgstDOTn * nx);
    rhovygst[0] = rhoinfty * (vintTANy + vgstDOTn * ny);
    rhoegst[0] = (pgst / (gam - 1)) + 0.5 * (rhovxgst[0] * rhovxgst[0] + rhovygst[0] * rhovygst[0]) / rhogst[0];
}

void BoundaryUpdate(const double gam, const int nbfac, const int nface, const int *bface, const int *inelfa,
                    const double *unkel0, const double *geofa, double *unkel) {
    // Give ghost cells values
    // loop over faces
    for (int iface = 0; iface < nface; iface++) {

        int Lelem = inelfa[iu(iface, 0, nface)];
        int Relem = inelfa[iu(iface, 1, nface)];

        if (Relem && Lelem >= 0) {
            //skip non-boundary points
            continue;
        }
        // ghost cell number
        int iegst = inelfa[iu(iface, 1, nface)];

        // interior cell number
        int ieint = inelfa[iu(iface, 0, nface)];

        // ghost cell boundary type
        int btype = bface[(iu(ibfac, 2, nbfac))];
        ASSERT(btype == 2 || btype == 4, "unknown boundary condition")


        //=============================CALCULATE NORMAL VELOCITY=========================================
        double rhoin = unkel[I3(ieint, 0, 0, 1)];

        //Find the length of the face
        double Lnx = geofa[ID(0, iface, nface)];
        double Lny = geofa[ID(1, iface, nface)];
        double length = sqrt(Lnx * Lnx + Lny * Lny);

        //Fine the unit normal and velocity components
        double nx = Lnx / length;
        double ny = Lny / length;
        double vxint = unkel[I3(ieint, 0, 1, 1)] / rhoin;
        double vyint = unkel[I3(ieint, 0, 2, 1)] / rhoin;

        //calculate Vn
        double vintDOTn = vxint * nx + vyint * ny;

        //=============================APPLY BOUNDARY CONDITIONS=========================================

        //Wall boundary
        if (btype == 2) {
            unkel[I3(iegst, 0, 0, 1)] = rhoin;
            unkel[I3(iegst, 0, 3, 1)] = unkel[I3(ieint, 0, 3, 1)]; //rho E stays the same
            //find ghost velocity to enforce wall || vg = vie - 2(vieDOTn) * n
            double vxgst = vxint - 2 * vintDOTn * nx;
            double vygst = vyint - 2 * vintDOTn * ny;

            unkel[I3(iegst, 0, 1, 1)] = vxgst * rhoin;
            unkel[I3(iegst, 0, 2, 1)] = vygst * rhoin;
        }

        //Freestream Boundary
        if (btype == 4) {
            //Get the full set of primitive variables+ on the element
            double pint, cint, Mint;
            getPrimativesP0(gam, unkel, ieint, &rhoin, &vxint, &vyint, &pint, &cint, &Mint);
            //Based on characteristics, firt step is to find the normal Mach number, Vn/c
            double MachN = vintDOTn / cint;

            if (MachN <= -1) {
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Supersonic Inflow - fully determined by free stream
                unkel[I3(iegst, 0, 0, 1)] = unkel0[0];
                unkel[I3(iegst, 0, 1, 1)] = unkel0[1];
                unkel[I3(iegst, 0, 2, 1)] = unkel0[2];
                unkel[I3(iegst, 0, 3, 1)] = unkel0[3];


            }
            if (MachN <= 0 && MachN > -1) {
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Subsonic Inflow
                double rhogst, rhovxgst, rhovygst, rhoegst;
                SubsonInflo(gam, vxint, vyint, cint, unkel0, nx, ny, &rhogst, &rhovxgst, &rhovygst, &rhoegst);

                unkel[I3(iegst, 0, 0, 1)] = rhogst;
                unkel[I3(iegst, 0, 1, 1)] = rhovxgst;
                unkel[I3(iegst, 0, 2, 1)] = rhovygst;
                unkel[I3(iegst, 0, 3, 1)] = rhoegst;


            }
            if (MachN >= 1) {
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Supersonic Outflow - fully determined by interior
                unkel[I3(iegst, 0, 0, 1)] = unkel[I3(ieint, 0, 0, 1)];
                unkel[I3(iegst, 0, 1, 1)] = unkel[I3(ieint, 0, 1, 1)];
                unkel[I3(iegst, 0, 2, 1)] = unkel[I3(ieint, 0, 2, 1)];
                unkel[I3(iegst, 0, 3, 1)] = unkel[I3(ieint, 0, 3, 1)];
            }


            if (MachN > 0 && MachN < 1) {
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Subsonic Outflow
                double rhogst, rhovxgst, rhovygst, rhoegst;
                SubsonOutfl(gam, rhoin, pint, vxint, vyint, cint, unkel0, nx, ny, &rhogst, &rhovxgst, &rhovygst,
                            &rhoegst);

                unkel[I3(iegst, 0, 0, 1)] = rhogst;
                unkel[I3(iegst, 0, 1, 1)] = rhovxgst;
                unkel[I3(iegst, 0, 2, 1)] = rhovygst;
                unkel[I3(iegst, 0, 3, 1)] = rhoegst;

            }
        }
        //printf("ielint:%d \tu1:%f\tu2:%f\tu3:%f\tu4:%f\n", ieint, unkel[I3(iegst,0,0,1)], unkel[I3(iegst,0,1,1)], unkel[I3(iegst,0,2,1)], unkel[I3(iegst,0,3,1)]);
    }
}

*/