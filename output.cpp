//
// Created by Tsail on 8/24/2023.
//

#include "output.h"

void printgrid(const char *title, int *inpoel, const int nelem, const int npoin, double *x, double *y) {
    FILE* fout = fopen("grid.tec", "w");

    //printf("\nDisplaying Grid Header\n");

    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\"\n");
    fprintf(fout, "ZONE T=\"fezone\", N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n", npoin, nelem);

    //printf("Printing Coordinate Information\n");

    for (int i=0; i<npoin; i++) {
        fprintf(fout, "%lf %lf\n", x[i], y[i]);
    }
    //printf("Printing Connectivity Information\n");
    for (int i=0; i<nelem; i++) {
        fprintf(fout, "%d, %d, %d, %d\n", inpoel[iu(i,0,nelem)]+1, inpoel[iu(i,1,nelem)]+1, inpoel[iu(i,2,nelem)]+1, inpoel[iu(i,3,nelem)]+1);
    }
    fclose(fout);
}

void printscalar(const char *title, const char *varname, const char *varname2, const char *varname3,
                 const char *varname4, int *inpoel, int nelem, int npoin, double *x, double* y,
                 double *unkel) {
    FILE* fout = fopen("plotP0.tec", "w");

    fprintf(fout, "TITLE = %s\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\", \"%s\", \"%s\", \"%s\", \"%s\", \"Mach\"\n", varname, varname2, varname3, varname4);
    fprintf(fout, "ZONE T=\"fezone\", N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL, varlocation=([1,2]=nodal,[3,4,5,6,7]=cellcentered)\n", npoin, nelem);

    //x node
    for (int i = 0; i<npoin; i++) {
        fprintf(fout, "%lf, ", x[i]);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }
    fprintf(fout, "\n\n");

    //y node
    for (int i = 0; i<npoin; i++) {
        fprintf(fout, "%lf,", y[i]);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }
    fprintf(fout, "\n\n");

    //var 1
    for (int i = 0; i<nelem; i++) {
        fprintf(fout, "%lf, ", unkel[iu3(i, 0, 0, 1)]);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }
    fprintf(fout, "\n\n");

    //var 2
    for (int i = 0; i<nelem; i++) {
        fprintf(fout, "%lf, ", unkel[iu3(i, 0, 1, 1)]);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }
    fprintf(fout, "\n\n");

    //var 3
    for (int i = 0; i<nelem; i++) {
        fprintf(fout, "%lf, ", unkel[iu3(i, 0, 2, 1)]);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }
    fprintf(fout, "\n\n");

    //var 4
    for (int i = 0; i<nelem; i++) {
        fprintf(fout, "%lf, ", unkel[iu3(i, 0, 3, 1)]);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }
    fprintf(fout, "\n\n");


    //Mach number
    for (int i = 0; i<nelem; i++) {
        double rho, vx, vy, p, c, M;
        getPrimativesPN(1.4, &unkel[iu3(i, 0, 0, 1)], &rho, &vx, &vy, &p, &c, &M);
        fprintf(fout, "%lf, ", M);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }
    fprintf(fout, "\n\n");



    for (int i=0; i<nelem; i++) {
        fprintf(fout, "%d, %d, %d, %d\n", inpoel[iu(i,0,nelem)]+1, inpoel[iu(i,1,nelem)]+1, inpoel[iu(i,2,nelem)]+1, inpoel[iu(i,3,nelem)]+1);
        if ((i+1) % 200 == 0) {fprintf(fout,"\n");}
    }

    fclose(fout);
}