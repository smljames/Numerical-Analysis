/*
EE407002 hw06. Matrix Condition Numbers
102061125 Chen Kuan-Chun
*/

/*
Usage:
$ g++ hw06.cpp MAT.cpp VEC.cpp
$ ./a.out 20 (for example: 20 resistors per side)
*/

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
    // ----------initialize----------

    int len, dim, vol; // len: # of nodes per side
                       // dim: total # of nodes
                       // vol:voltage
    double g; // g: conductance
    
    if (argv[1] == string("2")){
        len = 3;
        dim = 3*3;
        g = 0.001;
        vol = 1;
    }
    else if (argv[1] == string("4")){
        len = 5;
        dim = 5*5;
        g = 0.002;
        vol = 1;
    }
    else if (argv[1] == string("10")){
        len = 11;
        dim = 11*11;
        g = 0.005;
        vol = 1;
    }
    else if (argv[1] == string("20")){
        len = 21;
        dim = 21*21;
        g = 0.01;
        vol = 1;
    }
    else if (argv[1] == string("40")){
        len = 41;
        dim = 41*41;
        g = 0.02;
        vol = 1;
    }
    else if (argv[1] == string("50")){
        len = 51;
        dim = 51*51;
        g = 0.025;
        vol = 1;
    }
    else if (argv[1] == string("60")){
        len = 61;
        dim = 61*61;
        g = 0.025;
        vol = 1;
    }
    else if (argv[1] == string("80")){
        len = 81;
        dim = 81*81;
        g = 0.025;
        vol = 1;
    }
    else if (argv[1] == string("100")){
        len = 101;
        dim = 101*101;
        g = 0.025;
        vol = 1;
    }
    
    
       
    // ----------build the matrix and vector of the system----------

    MAT A(dim);
    VEC b(dim);
    
    for(int i=0; i<dim; i++){
        // each corner of resistor network
        if (i==0 || i==0+len-1 || i==dim-len || i==dim-1){
            A[i][i] = 2;
        }
        // resistors not on the side of the resistor network
        else if (i>len && i<dim-len-1 && (i%len!=(len-1)) && (i%len!=0)){
            A[i][i] = 4;
        }
        // rest of the resistors of the resistor network
        else{
            A[i][i] = 3;
        }
    }
    
    // for each resistor connecting node beside
    for(int i=0; i<dim; i++){
        if (i%len != 0 && (i-1) >= 0)    A[i][i-1] = -1;
        if ((i+1)%len != 0 && (i+1) <= (dim-1)) A[i][i+1] = -1;
        if (i-len >= 0)  A[i][i-len] = -1;
        if (i+len <= dim-1) A[i][i+len] = -1;
    }
    
    int vol_node = (len-1)/2; // the node connecting to the fixed voltage
    int gd_node = (dim-1)-(len-1)/2; // ground node
    double alpha;

    alpha = 0.003;

    for(int i=0; i<dim; i++){
        if (i == vol_node){
            A[vol_node][vol_node] = 1 * alpha;
        }
        else{
            A[vol_node][i] = 0;
        }
    }

    for(int i=0; i<dim; i++){
        if (i == gd_node){
            A[gd_node][gd_node] = 1;
        }
        else{
            A[gd_node][i] = 0;
        }
    }

    b[vol_node] = vol * alpha;
    
    MAT G(A);
    luFact(G);

    // output matrix file
    ofstream matrix_file;

    matrix_file.open("matrix.txt");
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            matrix_file << A[i][j] << " ";
        }
        matrix_file << "\n";
    }
    matrix_file.close();
    cout << "matrix is written in file" << endl;

    // output LU matrix file
    ofstream lumatrix_file;

    lumatrix_file.open("lu.txt");
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            lumatrix_file << G[i][j] << " ";
        }
        lumatrix_file << "\n";
    }
    lumatrix_file.close();
    cout << "LU matrix is written in file" << endl;


    
    


    // ----------calculate node voltages vector----------
    /* 
    VEC q1(dim);
    VEC q2(dim);
    VEC q3(dim);
    VEC q4(dim);
    int start, end; // record the clock
    int iter; // # of iteration that iterative solutions use
    int maxIter; // max # of iteration we allow
    double lambda = 0.5;
    double lambda1 = 0.5;
    double lambda2 = 0.5;
    double lambda3 = 0.5;
    double lambda4 = 0.5;

    cout << len-1 << " resistors per side" << endl; 

    // -----iterative method-----
    maxIter = 100;
    q1[0] = 1;
    q2[0] = 1;
    q3[0] = 1;
    q4[0] = 1;


    VEC iters(maxIter);
    VEC lambda_v1(maxIter);
    VEC lambda_v2(maxIter);
    VEC lambda_v3(maxIter);
    VEC lambda_v4(maxIter);
    VEC all_err1(maxIter);
    VEC all_err2(maxIter);
    VEC all_err3(maxIter);
    VEC all_err4(maxIter);

    lambda = 0.5;
    lambda_v1[0] = lambda;
    lambda_v2[0] = lambda;
    lambda_v3[0] = lambda;
    lambda_v4[0] = lambda;


    start = clock();
    iter = EVpwr_all_err(A, q1, lambda_v1, maxIter, 1, all_err1);
    iter = EVpwr_all_err(A, q2, lambda_v2, maxIter, 2, all_err2);
    iter = EVpwr_all_err(A, q3, lambda_v3, maxIter, 3, all_err3);
    iter = EVpwr_all_err(A, q4, lambda_v4, maxIter, 4, all_err4);
    
    cout << lambda_v1[iter] << endl;
    cout << lambda_v2[iter] << endl;
    cout << lambda_v3[iter] << endl;
    cout << lambda_v4[iter] << endl;
    end = clock();

    ofstream file1;
    ofstream file2;
    ofstream file3;
    ofstream file4;

    file1.open("err1.txt");
    for(int i=0; i<iter; i++){
        file1 << all_err1[i] << "\n";
    }
    file1.close();

    file2.open("err2.txt");
    for(int i=0; i<iter; i++){
        file2 << all_err2[i] << "\n";
    }
    file2.close();

    file3.open("err3.txt");
    for(int i=0; i<iter; i++){
        file3 << all_err3[i] << "\n";
    }
    file3.close();

    file4.open("err4.txt");
    for(int i=0; i<iter; i++){
        file4 << all_err4[i] << "\n";
    }
    file4.close();

    cout << "all error are written in file" << endl;


    ofstream lambda_file1;
    ofstream lambda_file2;
    ofstream lambda_file3;
    ofstream lambda_file4;

    lambda_file1.open("lambda1.txt");
    for(int i=0; i<iter; i++){
        lambda_file1 << lambda_v1[i] << "\n";
    }
    lambda_file1.close();

    lambda_file2.open("lambda2.txt");
    for(int i=0; i<iter; i++){
        lambda_file2 << lambda_v2[i] << "\n";
    }
    lambda_file2.close();

    lambda_file3.open("lambda3.txt");
    for(int i=0; i<iter; i++){
        lambda_file3 << lambda_v3[i] << "\n";
    }
    lambda_file3.close();

    lambda_file4.open("lambda4.txt");
    for(int i=0; i<iter; i++){
        lambda_file4 << lambda_v4[i] << "\n";
    }
    lambda_file4.close();

    cout << "all lambda are written in file" << endl;

    ofstream matrix_file;

    matrix_file.open("matrix.txt");
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            matrix_file << A[i][j] << " ";
        }
        matrix_file << "\n";
    }
    matrix_file.close();
    cout << "matrix_is written in file" << endl;
    */
}

