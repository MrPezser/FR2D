#include <iostream>
#include <cmath>

#include "indexing.h"
#include "basis.h"
#include "EulerFlux.h"
#include "SpatialDiscretization.h"
#include "Setup2D.h"



using namespace std;

void veccopy(double* a, const double* b, size_t n){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

double Initialize(double x){
    ///This defines the intial state of the solution space


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
    }

    //return 2.0 + sin(2.0*M_PI*x);

    //return 1.0 + exp(-40*(x-0.5)*(x-0.5));
}

void InitializeEuler(double x, double y, double* u){
    double rho = 1.0;
    double v = 1.0;
    double p = 1.0;

    //rho = Initialize(x);

    if (x < 0.5 && y < 0.5){
        //rho = 1.0;
        v = 0.0;
        //p = 1.0;
    } else {
        rho = 0.125;
        v = 0.0;
        p = 0.1;
    }

    u[0] = rho;                             //rho
    u[1] = rho * v;                         //rho V
    u[2] = 0.5*rho*v*v + (p/(GAM-1.0));       //rho e
}

int main() {
    ///hardcoded inputs
    //Input grid informaiton
    int imx = 10;
    int jmx = 10;
    int nelem = (imx-1) * (jmx-1);
    int nface = (2*nelem + (imx-1) + (jmx-1));
    int npoin = imx*jmx;

    int ndegr = 4;             //Degrees of freedom per element in one dimension
    int tdegr = ndegr*ndegr;   //Total degrees of freedom per element
    int nvar = 3;              //Number of variables
    int nu = nelem * ndegr * nvar;

    double cfl = 0.01 / (ndegr*ndegr);          //CFL Number
    double a = 1.0;             //Wave Speed

    double tmax = 0.2;

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
    double dx = 1.0 / (double)imx;
    double dy = 1.0 / (double)jmx;

    auto* x = (double*)malloc(npoin*sizeof(double));
    auto* y = (double*)malloc(npoin*sizeof(double));
    for (int j = 0; j < jmx; j++) {
        for (int i = 0; i < imx; i++) {
            int ipoin = iu(i,j,imx);
            //defining x & y position of cell corners
            x[ipoin] = i * dx;
            y[ipoin] = j * dy;
        }
    }

    double* eldrdxi;
    CalcCoordJacobian(ndegr, npoin, nelem, inpoel, x, y, xi, xi, eldrdxi);

    int* facpts;
    FacePoint2PointMap(ndegr,facpts);

    //This should be recalculated each time step
    double dt = (cfl * fmin(dx, dy)); ///Remember to do this

    //Aprox number of iterations required to get to the given tmax
    int niter = ceil(tmax/dt);

    //Allocate Arrays
    auto* u = (double*)malloc(nu*sizeof(double));
    auto* u0 = (double*)malloc(nu*sizeof(double));
    auto* dudt = (double*)malloc(nu*sizeof(double));


    //Generate Grid (currently uniform 2D) & initialize solution
    double dx = 1.0 / (double)imx;
    double dy = 1.0 / (double)jmx;

    for (int i=0; i<nelem; i++){
        for (int j=0; j<ndegr; j++) {
            for (int k=0; k<ndegr; k++) {
                //defining x & y position of cell centers
                x[k] = (k+0.5) * dx;
                y[j] = (j+0.5) * dy;

                int jnode = iu(k,j,ndegr);
                InitializeEuler(x[k] + xi[k] * (0.5 * dx), y[j] + xi[j] * (0.5 * dx), &u[iu3(i, jnode, 0, tdegr)]);
            }
        }
    }

    veccopy(u0,u,nu);

    // Begin Time Marching (3 stage TVD RK)
    auto* u_tmp = (double*)malloc(nu*sizeof(double));

    for (int iter=0; iter<niter; iter++){
        veccopy(u_tmp, u, nu);
        //1st stage
        CalcDudt(nelem, ndegr, nface, nvar, dx, inelfa, facpts, u, Dmatrix, Dradau, dudt );
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

        if (iter % 100 == 0){printf("iter:%10d\t%7.2f%% Complete\n",iter, 100.0*(double)iter/(double)niter);}
    }

    printf("iter=%d\tdt=%f\n", niter, dt);

    //Printout Final Solution
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


    free(xi);
    free(u);
    free(u0);
    free(dudt);
    free(Dmatrix);
    free(Dradau);
}
