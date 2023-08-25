#include <iostream>
#include <cmath>
#include <cstdlib>

#include "indexing.h"
#include "basis.h"
#include "EulerFlux.h"
#include "SpatialDiscretization.h"
#include "Setup2D.h"
#include "output.h"


using namespace std;

void veccopy(double* a, const double* b, size_t n){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

double Initialize(double x){
    ///This defines the intial state of the solution space

    /*
    //Gaussian bump and step combo
    if (x < 0.6) {
        double beta = 0.01;
        return 1 + exp(-(x-0.3)*(x-0.3) / beta);
    } else {
        if (x < 0.8) {
            return 2.0;
        } else {
            return 1.0;
        }
    }*/

    return x + 1.0;

    //return 2.0 + sin(2.0*M_PI*x);

    //return 1.0 + exp(-40*(x-0.5)*(x-0.5));
}

void InitializeEuler(double x, double y, double* u){
    double rho = 1.0;
    double vx = 1.0;
    double vy = 0.0;
    double p = 1.0;

    rho = Initialize(x);

    /*  //shock problem
    if (x < 0.5 && y < 0.5){
        //rho = 1.0;
        vx = 0.0;
        //vy = 0.0;
        //p = 1.0;
    } else {
        rho = 0.125;
        vx = 0.0;
        //vy = 0.0;
        p = 0.1;
    }*/

    u[0] = rho;                             //rho
    u[1] = rho * vx;                        //rho Vx
    u[2] = rho * vy;                        //rho Vy
    double vmags = vx*vx + vy*vy;
    u[3] = 0.5*rho*vmags + (p/(GAM-1.0));     //rho e
}

int main() {
    ///hardcoded inputs
    //Input grid informaiton
    int imx = 31;
    int jmx = 31;
    int nelem = (imx-1) * (jmx-1);
    int nface = (2*nelem + (imx-1) + (jmx-1));
    int npoin = imx*jmx;

    int ndegr = 1;             //Degrees of freedom per element in one dimension
    int tdegr = ndegr*ndegr;   //Total degrees of freedom per element
    //int nvar = NVAR;              //Number of variables
    int nu = nelem * tdegr * NVAR;

    double cfl = 0.01 / (ndegr*ndegr);          //CFL Number

    double tmax = 0.05;

    //Find the solution points in one reference dimension
    auto* xi = (double*)malloc(ndegr*sizeof(double));
    GenerateLobattoPoints(ndegr, xi);

    //Find the derivative matrix of the associated lagrange polynomials
    //Derivatie of L_j at position x_i
    auto* Dmatrix = (double*)malloc(ndegr*ndegr*sizeof(double));
    GenerateLagrangeDMatrix(ndegr, xi, Dmatrix);

    //Find the derivatives of the Radau polynomial of appropriate degree
    auto* Dradau = (double*)malloc((1+ndegr)*sizeof(double));
    GenerateRadauDerivatives(ndegr, xi, Dradau);

    //Setup for 2D based on a structured Quad grid (structure input, processing done unstructred)
    int *inpoel, *inpofa, *inelfa;
    LoadStruct2Unstruct(imx, jmx, nelem, nface, &inpoel, &inpofa, &inelfa);


    //test grid
    double dx = 1.0 / (double)(imx-1);
    double dy = 1.0 / (double)(jmx-1);


    //auto* x = (double*)malloc(npoin*sizeof(double));
    //auto* y = (double*)malloc(npoin*sizeof(double));
    double x[npoin];
    double y[npoin];

    for (int j = 0; j < jmx; j++) {
        for (int i = 0; i < imx; i++) {
            int ipoin = iu(i,j,imx);
            //defining x & y position of cell corners
            x[ipoin] = i * dx;
            y[ipoin] = j * dy;
        }
    }

    auto* eldrdxi = (double*)malloc(nelem*ndegr*4*sizeof(double));
    auto* eldxidr = (double*)malloc(nelem*ndegr*4*sizeof(double));
    auto* eljac   = (double*)malloc(nelem*ndegr*ndegr*sizeof(double));
    CalcCoordJacobian(ndegr, npoin, nelem, inpoel, x, y, xi, xi, eldrdxi, eldxidr, eljac);

    auto* facpts = (int*)malloc(4*ndegr*sizeof(int));
    FacePoint2PointMap(ndegr,facpts);

    printgrid("Title", inpoel, nelem, npoin, x, y);

    //This should be recalculated each time step
    double dt = (cfl * fmin(dx, dy)); ///Remember to do this

    //Aprox number of iterations required to get to the given tmax
    int niter = ceil(tmax/dt);

    //Allocate Arrays
    auto* u = (double*)malloc(2*nu*sizeof(double));
    auto* u0 = (double*)malloc(2*nu*sizeof(double));
    auto* dudt = (double*)malloc(2*nu*sizeof(double));


    //Generate Grid (currently uniform 2D) & initialize solution
    for (int ielem=0; ielem<nelem; ielem++){
        for (int j=0; j<ndegr; j++) {
            for (int k=0; k<ndegr; k++) {
                int jnode = iu(k,j,ndegr);
                int ipoin = inpoel[iu(ielem, 0, nelem)]; //BL corner

                double xnode = x[ipoin] + (xi[k] + 1.0) * (0.5 * dx);  ///assuming cartesean grid
                double ynode = y[ipoin] + (xi[j] + 1.0) * (0.5 * dy);
                InitializeEuler(xnode, ynode, &u[iu3(ielem, jnode, 0, tdegr)]);
            }
        }
    }

    veccopy(u0,u,nu);

    //printf("elem:%d \t rho*u: %f\n", 99, u[iu3(99, 0, 1, tdegr)]);

    // Begin Time Marching (3 stage TVD RK)
    auto* u_tmp = (double*)malloc(nu*sizeof(double));

    for (int iter=0; iter<niter; iter++){
        veccopy(u_tmp, u, nu);
        //1st stage
        CalcDudt(nelem, ndegr, nface, NVAR, inelfa, facpts, u, Dmatrix, Dradau, eldrdxi, eldxidr, eljac, dudt);
        for (int i=0; i<nu; i++){
            //u_tmp[i] += dt * dudt[i];
            u[i] += dt * dudt[i];

            if (isnan(u[i])){
                throw overflow_error("Kaboom!\n");
            }
        }
        /*
        //2nd stage
        CalcDudt(nx, ndegr, nvar, a, dx, u_tmp, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            u_tmp[i] = 0.75*u[i] + 0.25*( u_tmp[i] + dt*dudt[i]);
        }

        //3rd stage
        CalcDudt(nx, ndegr, nvar, a, dx, u_tmp, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            u[i] = (1.0/3.0)*u[i] + (2.0/3.0)*(u_tmp[i] + dt*dudt[i]);
        }*/

        if (iter % 10 == 0){printf("iter:%10d\t%7.2f%% Complete\n",iter, 100.0*(double)iter/(double)niter);}
    }

    printf("iter=%d\tdt=%f\n", niter, dt);

    printscalar("FV Output", "rho", "rho_u", "rho_v", "rho_e", inpoel, nelem, npoin, x, y, u);

    /*
    //Printout Solution
    FILE* fout = fopen("waveout.tec", "w");
    fprintf(fout, "x\tu\tu0\n");

    for (int i=0;i<nx;i++) {
        for (int j=0; j<ndegr; j++) {
            for (int k=0; k<ndegr; k++) {

                int idegr = iup(k,j,ndegr);
                double xk = x[k] + xi[k] * (0.5 * dx);
                double yj = y[j] + xi[j] * (0.5 * dy);
                fprintf(fout, "%f\t%f\t%f\t%f\t%f\t", xk, yj, u[iu3(i, j, 0, ndegr)], u[iu3(i, j, 1, ndegr)],
                        u[iu3(i, j, 2, ndegr)]);
            }
        }
    }
    fclose(fout);
    */

    /*
    free(xi);
    free(u);
    free(u0);
    free(dudt);
    free(Dmatrix);
    free(Dradau); */
}
