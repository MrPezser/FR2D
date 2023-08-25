//
// Created by Tsail on 6/14/2023.
//

#ifndef FR1D_INDEXING_H
#define FR1D_INDEXING_H

#define ASSERT(cond, msg) if(!(cond)){printf("Failed Assert: %s:%u %s\nYou idiot! How could you not realize that %s\n", __FILE__, __LINE__, #cond, msg); exit(0);}

//indexing into state variable array u of degree nj for the jth node on the ith element
//column major i guess, maybe not
#define iu(i, j, ni)  (((j)*(ni)) + (i))
#define NVAR 4
#define iu3(ielem, jdegr, kvar, ndegr) ((((ielem)*(ndegr)) + (jdegr))*NVAR + (kvar))
#define ASSERT(cond, msg) if(!(cond)){printf("Failed Assert: %s:%u %s\nYou idiot! How could you not realize that %s\n", __FILE__, __LINE__, #cond, msg); exit(0);}

#endif //FR1D_INDEXING_H
