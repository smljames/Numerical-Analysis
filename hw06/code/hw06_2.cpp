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
    
    // -----iterative method-----
    VEC q0(dim); // initial guess eigenvector
    int start, end; // record the clock
    int iter; // # of iteration that iterative solutions use
    int maxIter; // max # of iteration we allow
    double lambda_l; // largest non-1 lambda using Pwr
    double lambda_s; // smallest lambda using iPwr
    double lambda_sft_s; // smallest lambda using shifted iPwr

    cout << len-1 << " resistors per side" << endl; 

    maxIter = 10000;
    for(int i=0; i<dim; i++) q0[i] = 1; // initialize q0
    start = clock();
    iter = EVpwr(A, q0, lambda_l, 1e-09, maxIter);
    end = clock();
    cout<<"[ Largest ]"<<endl;
    cout<<"iter = "<<iter<<"\t"<<"non-1 lambda = "<<lambda_l<<endl;
    cout<<"time = "<<(float)(end-start)/CLOCKS_PER_SEC<<"\t"\
        <<"time/iter = "<<(float)(end-start)/CLOCKS_PER_SEC/iter<<"\n"<<endl;

    for(int i=0; i<dim; i++) q0[i] = 1;
    start = clock();
    iter = EVipwr(A, q0, lambda_s, 1e-09, maxIter);
    end = clock();
    cout<<"[ smallest ]"<<endl;
    cout<<"iter = "<<iter<<"\t"<<"lambda = "<<lambda_s<<endl;
    cout<<"time = "<<(float)(end-start)/CLOCKS_PER_SEC<<"\t"\
        <<"time/iter = "<<(float)(end-start)/CLOCKS_PER_SEC/iter<<"\n"<<endl;

    for(int i=0; i<dim; i++) q0[i] = 1;
    start = clock();
    iter = EVipwrShft(A, q0, lambda_sft_s, 0.00076, 1e-09, maxIter);
    end = clock();
    cout<<"[ shifted iPwr smallest ]"<<endl;
    cout<<"iter = "<<iter<<"\t"<<"lambda = "<<lambda_sft_s<<endl;
    cout<<"time = "<<(float)(end-start)/CLOCKS_PER_SEC<<"\t"\
        <<"time/iter = "<<(float)(end-start)/CLOCKS_PER_SEC/iter<<"\n"<<endl;
    
    cout<<"condition number = "<<lambda_l/lambda_s<<"\n"<<endl;
}

