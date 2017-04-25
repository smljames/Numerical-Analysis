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
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("4")){
        len = 5;
        dim = 5*5;
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("10")){
        len = 11;
        dim = 11*11;
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("20")){
        len = 21;
        dim = 21*21;
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("40")){
        len = 41;
        dim = 41*41;
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("50")){
        len = 51;
        dim = 51*51;
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("60")){
        len = 61;
        dim = 61*61;
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("80")){
        len = 81;
        dim = 81*81;
        g = (len-1)/2000.0;
        vol = 1;
    }
    else if (argv[1] == string("100")){
        len = 101;
        dim = 101*101;
        g = (len-1)/2000.0;
        vol = 1;
    }
    
    
       
    // ----------build the matrix and vector of the system----------

    MAT A(dim);
    VEC b(dim);
    
    for(int i=0; i<dim; i++){
        // each corner of resistor network
        if (i==0 || i==0+len-1 || i==dim-len || i==dim-1){
            A[i][i] = 2*g;
        }
        // resistors not on the side of the resistor network
        else if (i>len && i<dim-len-1 && (i%len!=(len-1)) && (i%len!=0)){
            A[i][i] = 4*g;
        }
        // rest of the resistors of the resistor network
        else{
            A[i][i] = 3*g;
        }
    }
    
    // for each resistor connecting node beside
    for(int i=0; i<dim; i++){
        if (i%len != 0 && (i-1) >= 0)    A[i][i-1] = -1*g;
        if ((i+1)%len != 0 && (i+1) <= (dim-1)) A[i][i+1] = -1*g;
        if (i-len >= 0)  A[i][i-len] = -1*g;
        if (i+len <= dim-1) A[i][i+len] = -1*g;
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
            A[gd_node][gd_node] = 1 * alpha;
        }
        else{
            A[gd_node][i] = 0;
        }
    }

    b[vol_node] = vol * alpha;
    
    
    // ----------calculate node voltages vector----------
    
    VEC q1(dim);
    VEC q2(dim);
    VEC q3(dim);
    VEC q4(dim);
    int start1, end1; // record the clock
    int start2, end2; // record the clock
    int start3, end3; // record the clock
    int start4, end4; // record the clock
    int iter; // # of iteration that iterative solutions use
    int maxIter; // max # of iteration we allow
    double lambda = 0.5;
    double lambda1 = 0.5;
    double lambda2 = 0.5;
    double lambda3 = 0.5;
    double lambda4 = 0.5;

    cout << len-1 << " resistors per side" << endl; 

    // -----iterative method-----
    maxIter = 1000;
    for(int i=0; i<dim; i++) q1[i] = 1;
    for(int i=0; i<dim; i++) q2[i] = 1;
    for(int i=0; i<dim; i++) q3[i] = 1;
    for(int i=0; i<dim; i++) q4[i] = 1;


    VEC iters(maxIter);
    VEC lambda_v1(maxIter);
    VEC lambda_v2(maxIter);
    VEC lambda_v3(maxIter);
    VEC lambda_v4(maxIter);
    VEC all_err1(maxIter);
    VEC all_err2(maxIter);
    VEC all_err3(maxIter);
    VEC all_err4(maxIter);

    lambda = 0;
    lambda_v1[0] = lambda;
    lambda_v2[0] = lambda;
    lambda_v3[0] = lambda;
    lambda_v4[0] = lambda;
    

    start1 = clock();
    iter = EVpwr_all_err(A, q1, lambda_v1, maxIter, 1, all_err1);
    end1 = clock();
    
    start2 = clock();
    iter = EVpwr_all_err(A, q2, lambda_v2, maxIter, 2, all_err2);
    end2 = clock();
    
    start3 = clock();
    iter = EVpwr_all_err(A, q3, lambda_v3, maxIter, 3, all_err3);
    end3 = clock();
    
    start4 = clock();
    iter = EVpwr_all_err(A, q4, lambda_v4, maxIter, 4, all_err4);
    end4 = clock();
    
    cout << lambda_v1[iter] << endl;
    cout << lambda_v2[iter] << endl;
    cout << lambda_v3[iter] << endl;
    cout << lambda_v4[iter] << endl;

    cout<<((float)(end1-start1)/CLOCKS_PER_SEC)/maxIter<<endl;
    cout<<((float)(end2-start2)/CLOCKS_PER_SEC)/maxIter<<endl;
    cout<<((float)(end3-start3)/CLOCKS_PER_SEC)/maxIter<<endl;
    cout<<((float)(end4-start4)/CLOCKS_PER_SEC)/maxIter<<endl;

    ofstream file1;
    ofstream file2;
    ofstream file3;
    ofstream file4;

    file1.open("../output/err1.txt");
    for(int i=0; i<iter; i++){
        file1 << all_err1[i] << "\n";
    }
    file1.close();

    file2.open("../output/err2.txt");
    for(int i=0; i<iter; i++){
        file2 << all_err2[i] << "\n";
    }
    file2.close();

    file3.open("../output/err3.txt");
    for(int i=0; i<iter; i++){
        file3 << all_err3[i] << "\n";
    }
    file3.close();

    file4.open("../output/err4.txt");
    for(int i=0; i<iter; i++){
        file4 << all_err4[i] << "\n";
    }
    file4.close();

    cout << "all error are written in file" << endl;


    ofstream lambda_file1;
    ofstream lambda_file2;
    ofstream lambda_file3;
    ofstream lambda_file4;

    lambda_file1.open("../output/lambda1.txt");
    for(int i=0; i<iter; i++){
        lambda_file1 << lambda_v1[i] << "\n";
    }
    lambda_file1.close();

    lambda_file2.open("../output/lambda2.txt");
    for(int i=0; i<iter; i++){
        lambda_file2 << lambda_v2[i] << "\n";
    }
    lambda_file2.close();

    lambda_file3.open("../output/lambda3.txt");
    for(int i=0; i<iter; i++){
        lambda_file3 << lambda_v3[i] << "\n";
    }
    lambda_file3.close();

    lambda_file4.open("../output/lambda4.txt");
    for(int i=0; i<iter; i++){
        lambda_file4 << lambda_v4[i] << "\n";
    }
    lambda_file4.close();

    cout << "all lambda are written in file" << endl;

    ofstream matrix_file;

    matrix_file.open("../output/matrix.txt");
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            matrix_file << A[i][j] << " ";
        }
        matrix_file << "\n";
    }
    matrix_file.close();
    cout << "matrix_is written in file" << endl;
}

