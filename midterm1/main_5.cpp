#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include "MAT.h"
#include "VEC.h"

using namespace std;

int main(int argc, char* argv[])
{
    int len, dim, vol;
    double g = 0;
    dim = 14;
    MAT A(dim);
    VEC b(dim);
    /*
    A[0][0] = 3;
    A[0][1] = -1;
    A[0][4] = -1;
    A[1][1] = 3;
    A[1][2] = -1;
    A[1][5] = -1;
    A[2][1] = -1;
    A[2][2] = 2;
    A[2][6] = -1;
    A[3][3] = 3;
    A[3][4] = -1;
    A[3][7] = -1;
    A[4][4] = 4;
    A[4][3] = -1;
    A[4][5] = -1;
    A[4][8] = -1;
    A[5][1] = -1;
    A[5][4] = -1;
    A[5][5] = 4;
    A[5][6] = -1;
    A[5][9] = -1;
    A[6][2] = -1;
    A[6][5] = -1;
    A[6][6] = 3;
    A[6][10] = -1;
    A[7][3] = -1;
    A[7][7] = 3;
    A[7][8] = -1;
    A[7][11] = -1;
    A[8][4] = -1;
    A[8][7] = -1;
    A[8][8] = 4;
    A[8][9] = -1;
    A[8][12] = -1;
    A[9][5] = -1;
    A[9][8] = -1;
    A[9][9] = 3;
    A[10][6] = -1;
    A[10][9] = -1;
    A[10][10] = 3;
    A[11][7] = -1;
    A[11][11] = 2;
    A[11][12] = -1;
    A[12][8] = -1;
    A[12][11] = -1;
    A[12][12] = 3;
    A[12][13] = -1;
    A[13][9] = -1;
    A[13][12] = -1;
    A[13][13]= 3;
    */
    fstream fin(argv[1]);
    fin >> dim;
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            fin >> A[i][j];
        }
    }

    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
    
    double V = 1.0;
    double alpha = 0.003;
    //b[0] = V * alpha;
    //A[0][0] *= alpha;


    int iter;
    int maxIter = 3000;
    double lambda;
    VEC lambda_v1(maxIter);
    VEC lambda_v2(maxIter);
    VEC lambda_v3(maxIter);
    VEC lambda_v4(maxIter);
    VEC q1(dim);
    VEC q2(dim);
    VEC q3(dim);
    VEC q4(dim);
    VEC all_err1(maxIter);
    q1[0] = 0.5;
    q2[0] = 1;
    q3[0] = 1;
    q4[0] = 1;

    iter = EVpwr_all_err(A, q1, lambda_v1, maxIter, 1, all_err1);
    printf("largest lambda = %g\n", lambda_v1[iter]);
    iter = EVipwr(A, q1, lambda, maxIter);
    printf("smallest lambda = %g\n", lambda);
    printf("k = %g\n", lambda_v1[iter]/lambda);

    MAT B(A*A);
    iter = EVpwr_all_err(B, q1, lambda_v2, maxIter, 1, all_err1);
    printf("largest lambda = %g\n", lambda_v2[iter]);
    iter = EVipwr(B, q1, lambda, maxIter);
    printf("smallest lambda = %g\n", lambda);
    printf("k = %g\n", lambda_v2[iter]/lambda);
/*
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            cout<<B[i][j]<<" ";
        }
        cout<<endl;
    }
*/

    MAT C(A + A.tpose());
    iter = EVpwr_all_err(C, q1, lambda_v3, maxIter, 1, all_err1);
    printf("largest lambda = %g\n", lambda_v3[iter]);
    iter = EVipwr(C, q1, lambda, maxIter);
    printf("smallest lambda = %g\n", lambda);
    printf("k = %g\n", lambda_v3[iter]/lambda);
 /* 
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            cout<<C[i][j]<<" ";
        }
        cout<<endl;
    }
  */
   
    MAT I(dim);
    for(int i=0; i<dim; i++){
        I[i][i] = 1;
    }
    MAT D(A+I); 
    iter = EVpwr_all_err(D, q1, lambda_v4, maxIter, 1, all_err1);
    printf("largest lambda = %g\n", lambda_v4[iter]);
    iter = EVipwr(D, q1, lambda, maxIter);
    printf("smallest lambda = %g\n", lambda);
    printf("k = %g\n", lambda_v4[iter]/lambda);
    
    
    VEC x(dim);
    VEC y(dim);
    A[0][0] = 1;
    b[0] = 1;
    A = luFact(A);
    y = fwdSubs(A, b);
    x = bckSubs(A, y);


    double req;
    req = V/((V-x[1]) + (V-x[3]));
    printf("R_eq = %g\n", req);
}


