#include "Setup2D.h"

void LoadStruct2Unstruct(const int imx, const int jmx, int nelem, int nface,
                         int** inpoel, int** inpofa, int** inelfa) {
    /*Processes given structured grid points and outputs unstructured data format
     * Inputs:
     * imx, jmx - max i and j values
     *
     * Outputs:
     * nelem - number of elements
     * nface - number of faces
     * inpoel- interconnectivity points from elements
     * inpofa- interconnectivity points from faces
     * inelfa- interconnectivity elements from faces
    */

    int npoin = imx*jmx;

    //
    *inpoel = (int*)malloc(4*nelem*sizeof(int));
    *inpofa = (int*)malloc(2*nface*sizeof(int));
    *inelfa = (int*)malloc(2*nface*sizeof(int));


    //loop through every element
    int iface;
    for (int j=0; j<jmx-1; j++){
        for (int i=0; i<imx-1; i++){
            ///assign points using i-major order p(i,j) = j*imx + i = iup(i,j,imx)
            //point i,j, is the bottom left corner of element if i,j ~ x,y
            //points assigned counterclockwise
            int ielem = iu(i,j,imx-1);
            (*inpoel)[iu(ielem,0,nelem)] = iu(i  ,j  ,imx);   //Bot Lef
            (*inpoel)[iu(ielem,1,nelem)] = iu(i+1,j  ,imx);   //Bot Rig
            (*inpoel)[iu(ielem,2,nelem)] = iu(i+1,j+1,imx);   //Top Rig
            (*inpoel)[iu(ielem,3,nelem)] = iu(i  ,j+1,imx);   //Tol Lef


            ///assign points of faces and record elements

            //Left Face
            iface = 2 * ielem;
            (*inpofa)[iu(iface,0,nface)] = iu(i,j  ,imx);    //Bot Lef
            (*inpofa)[iu(iface,1,nface)] = iu(i,j+1,imx);    //Top Lef

            (*inelfa)[iu(iface,0,nface)] = ielem-1;    //elem left
            (*inelfa)[iu(iface,1,nface)] = ielem;      //elem right
            (*inelfa)[iu(iface,2,nface)] = 0; //Vertical face
            if (i==0){ //boundary
                (*inelfa)[iu(iface,0,nface)] = -1;
            }

            //Bottom Face
            iface = 2 * ielem + 1;
            (*inpofa)[iu(iface,0,nface)] = iu(i  ,j,imx);     //Bot Lef
            (*inpofa)[iu(iface,1,nface)] = iu(i+1,j,imx);     //Bot Rig

            (*inelfa)[iu(iface,0,nface)] = ielem;              //elem left of face
            (*inelfa)[iu(iface,1,nface)] = ielem - (imx-1);    //elem right of face
            (*inelfa)[iu(iface,2,nface)] = 1; //Horizontal face
            if (j==0){ //boundary
                (*inelfa)[iu(iface,1,nface)] = -1;
            }

            if (i==imx-2){ //very far right
                iface = 2*nelem + 1 + j;
                (*inpofa)[iu(iface,0,nface)] = iu(i+1,j  ,imx);    //Bot Rig
                (*inpofa)[iu(iface,1,nface)] = iu(i+1,j+1,imx);    //Top Rig

                (*inelfa)[iu(iface,0,nface)] = -1;     //elem left of face
                (*inelfa)[iu(iface,1,nface)] = ielem;  //elem right of face
                (*inelfa)[iu(iface,2,nface)] = 0; //Vertical face
            }
            if (j==jmx-2){ //very far top
                iface = 2*nelem + 1 + jmx-1 + i;
                (*inpofa)[iu(iface,0,nface)] = iu(i  ,j+1,imx);     //Top Lef
                (*inpofa)[iu(iface,1,nface)] = iu(i+1,j+1,imx);     //Top Rig

                (*inelfa)[iu(iface,0,nface)] = ielem; //elem left of face
                (*inelfa)[iu(iface,1,nface)] = -1;    //elem right of face
                (*inelfa)[iu(iface,2,nface)] =  1; //Horizontal face
            }


        }
    }
}

void CalcCoordJacobian(int order, int npoin, int nelem, const int* inpoel, const double* x_in, const double* y_in,\
                       const double* xi, const double* eta, double* eldrdxi, double* eldxidr, double* eljac){
    /*Calculates the coordinate transformation matrix d(x,y)/d(xi,eta) for each column and row in each quad element
     *Input:
     *  ndegr  - degrees of freedom per dimension of the cell
     *  npoin  - number of points
     *  inpoel - points corresponding to each element
     *  coords - x and y coordinated of point
     *  xi     - solution points in reference domain xi direction
     *  eta    - solution points in ref. dir. eta
     *
     *Output:
     * eldrdxi - derivatives of the coordinate transform for each row/col of each element ( assuming each row/col has its own linear transformation => general quadrilaterals)
     *          d(x,y)/d(eta)= f(xi)  d(x,y)/d(xi) = f(eta)
     *          indexed: (ielem*ndegr*2 + idegr)*2 + ixy
     *              ielem = element index  |  ndegr = above  |  idegr = current degree of freedom  |  ixy = dx or dy
     *              idegr: [0,1,...,(2*ndegr)-1] = [d(xi) @ eta = 0,1,...ndegr-1,d(eta) @ xi=0,1,...ndegr-1]
     *              ixy: [0,1] = [dx/d{} , dy/d{}]
     * eldxidr - inverse of the above derivatives. At the same index that above has d$/d# , this has d#/d$
     * eljac    - magnitude of the dr/dxi transform jacobian matrix
     */

    //On each element we need to figure out the slope of each (xi = const.) and (eta = const.) line.

    //Points indexed; ipoin = iu(i_xi, j_eta, ni)
    int ndegr = order*order;

    for (int ielem=0; ielem < nelem; ielem++){
        double x[4], y[4];
        int ipoin;

        //Get the coordinates of each corner point
        for (int icrnr = 0; icrnr<4; icrnr++) {
            ipoin = inpoel[iu(ielem, icrnr, nelem)];
            x[icrnr] = x_in[ipoin];
            y[icrnr] = y_in[ipoin];
        }

        //Interpolate along the edges of the quadrilateral to find the endpoints for the constant xi/eta slices
        double etamax_x, etamin_x, ximax_x, ximin_x;
        double etamax_y, etamin_y, ximax_y, ximin_y;
        // {*}min/max contains the x or y values at the maximum or minimum eat/xi values
        // i.e. the endpoints of the constant xi/eta lines
        //  xi ~~ i dir ||  eta ~~ j dir

        for (int islice=0; islice<ndegr; islice++) {
            double etai = eta[islice];
            double xii = xi[islice];

            double alpha = (xii + 1.0) * 0.5;
            double beta = (etai + 1.0) * 0.5;

            etamax_x = alpha * x[2] + (1.0 - alpha) * x[3]; //alpha going from TL corner to TR
            etamin_x = alpha * x[1] + (1.0 - alpha) * x[0]; //alpha going from BL corner to BR
            ximax_x = beta * x[2] + (1.0 - beta) * x[1];  //beta going from BR corner to TR
            ximin_x = beta * x[3] + (1.0 - beta) * x[0];  //beta going from BL corner to TL

            etamax_y = alpha * y[2] + (1.0 - alpha) * y[3]; //alpha going from TL corner to TR
            etamin_y = alpha * y[1] + (1.0 - alpha) * y[0]; //alpha going from BL corner to BR
            ximax_y = beta * y[2] + (1.0 - beta) * y[1];  //beta going from BR corner to TR
            ximin_y = beta * y[3] + (1.0 - beta) * y[0];  //beta going from BL corner to TL

            //Define the slope of the lines along which the solution points lie
            double dxdxi  = (ximax_x - ximin_x) / 2.0;
            double dydxi  = (ximax_y - ximin_y) / 2.0;
            double dxdeta = (etamax_x - etamin_x) / 2.0;
            double dydeta = (etamax_y - etamin_y) / 2.0;

            double dxidx  = 1.0 / dxdxi;
            double dxidy  = 1.0 / dydxi;
            double detadx = 1.0 / dxdeta;
            double detady = 1.0 / dydeta;

            //xi derivatives
            int ind = 2 * iu(ielem, islice + 0, nelem);
            eldrdxi[ind + 0] = dxdxi;
            eldrdxi[ind + 1] = dydxi;

            //eta derivatives
            ind = 2 * iu(ielem, islice + ndegr, nelem);
            eldrdxi[ind + 0] = dxdeta;
            eldrdxi[ind + 1] = dydeta;
        }
        for (int i_xi=0; i_xi<order; i_xi++){
            for (int j_eta=0; j_eta<order; j_eta++) {
                ipoin = iu(i_xi, j_eta, order);

                //xi derivatives
                int ind = 2 * iu(ielem, j_eta + 0, nelem);
                double dxdxi = eldrdxi[ind + 0];
                double dydxi = eldrdxi[ind + 1];

                //eta derivatives
                ind = 2 * iu(ielem, i_xi + ndegr, nelem);
                double dxdeta = eldrdxi[ind + 0];
                double dydeta = eldrdxi[ind + 1];

                //Jacobian
                double jac = (dxdxi * dydeta) - (dxdeta * dydxi);
                ind = iu(ielem, ipoin, nelem);
                eljac[ind] = jac;
            }
        }
    }

}

void FacePoint2PointMap(int ndegr, int* facpts) {
    /*
     *      Records the indexing for points on the left and right of a given face assuming 1-to-1 2D quads
     */

    //Points indexed; ipoin = iu(i_xi, j_eta, ni)
    //Need interior points perpendicular from the face points

    for (int i=0; i<ndegr; i++){

        int id  = ndegr * iu(i, 0, ndegr);
        int id2 = ndegr * iu(i, 1, ndegr);
        int id3 = ndegr * iu(i, 2, ndegr);
        int id4 = ndegr * iu(i, 3, ndegr);
        for (int j=0; j<ndegr; j++){
            //j=0 = xi_max set of points which will lie on the left side of the face
            //j++ incrementing towards the interior until arriving at the other face

            //Vertical Face
            //Left
            facpts[id +j] = iu(ndegr-1-j, i, ndegr); // Left side of vertical face
            //Right
            facpts[id2+j] = iu(j        , i, ndegr); // Right side of vertical face

            //Horizontal Face
            //Left|Top
            facpts[id3+j] = iu(i, j        , ndegr);    // Left|Top  side of vertical face
            //Right|Upper
            facpts[id4+j] = iu(i, ndegr-1-j, ndegr);    // Right|Bot side of vertical face
        }

    }
}