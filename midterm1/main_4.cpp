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
    
    A[0][0] = 1;
    A[1][0] = -1;
    A[1][1] = 3;
    A[1][2] = -1;
    A[1][4] = -1;
    A[2][1] = -1;
    A[2][2] = 2;
    A[2][5] = -1;
    A[3][0] = -1;
    A[3][3] = 3;
    A[3][4] = -1;
    A[3][7] = -1;
    A[4][1] = -1;
    A[4][3] = -1;
    A[4][4] = 2;
    //A[4][5] = -1*g;
    //A[4][8] = -1*g;
    A[5][2] = -1;
    //A[5][4] = -1*g;
    A[5][5] = 2;
    A[5][6] = -1;
    //A[5][9] = -1*g;
    A[6][5] = -1;
    A[6][6] = 2;
    A[6][10] = -1;
    A[7][3] = -1;
    A[7][7] = 2;
    A[7][8] = -1;
    //A[8][4] = -1*g;
    A[8][7] = -1;
    A[8][8] = 2;
    //A[8][9] = -1*g;
    A[8][11] = -1;
    //A[9][5] = -1*g;
    //A[9][8] = -1*g;
    A[9][9] = 2;
    A[9][10] = -1;
    A[9][12] = -1;
    A[10][6] = -1;
    A[10][9] = -1;
    A[10][10] = 3;
    A[10][13] = -1;
    A[11][8] = -1;
    A[11][11] = 2;
    A[11][12] = -1;
    A[12][9] = -1;
    A[12][11] = -1;
    A[12][12] = 3;
    A[12][13] = -1;
    A[13][13] = 1;
    
    double V = 1.0;
    double alpha = 0.003;
    b[0] = V * alpha;
    A[0][0] *= alpha;
    A[13][13] *= alpha;


    int iter;
    int maxIter = 3000;
    double lambda;
    VEC lambda_v1(maxIter);
    VEC q1(dim);
    VEC q2(dim);
    VEC all_err1(maxIter);
    q1[0] = 1;
    q2[0] = 1;

    iter = EVpwr_all_err(A, q1, lambda_v1, maxIter, 1, all_err1);
    printf("largest lambda = %g\n", lambda_v1[iter]);
    iter = EVipwr(A, q2, lambda, maxIter);
    printf("smallest lambda = %g\n", lambda);
    printf("k = %g\n", lambda_v1[iter]/lambda);

    VEC x(dim);
    VEC y(dim);
    A[0][0] = 1;
    A[13][13] = 1;
    b[0] = 1;
    MAT G(A);
    A = luFact(A);
    y = fwdSubs(A, b);
    x = bckSubs(A, y);
/*
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            cout<<G[i][j]<<" ";
        }
        cout<<endl;
    }

    for(int i=0; i<dim; i++){
        cout<<x[i]<<endl;
    }
*/
    double req;
    req = V/((V-x[1]) + (V-x[3]));
    printf("R_eq = %g\n", req);
}


