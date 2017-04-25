/*
EE407002 hw04. Linear Iterative Methods
102061125 Chen Kuan-Chun
*/

/*
Usage:
$ g++ hw05.cpp MAT.cpp VEC.cpp
$ ./a.out 10 (for example: 10 resistors per side)
*/

#include <iostream>
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

    MAT A(dim), G(dim);
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
    
    for(int i=0; i<dim; i++){
        if (i == vol_node){
            A[vol_node][vol_node] = 1;
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

    b[vol_node] = vol;
    
    G = A;

    // ----------calculate node voltages vector----------
    
    VEC y(dim);
    VEC x(dim); // node voltage vector
    VEC x_iter(dim); // initail vector for iterative solutions 
    int start, end; // record the clock
    int iter; // # of iteration that iterative solutions use
    double total_i; // total current i
    double eq_r; // equivalent resistance
    double v_ne, v_ea, v_sw, v_ne1, v_ea1, v_sw1;

    cout << endl;
    cout << len-1 << " resistors per side" << endl;
    
    // -----LU decomposition method-----
    //cout << "---LU decomposition---" << endl;

    //G = luFact(G);// G is a copy of A
    //y = fwdSubs(G, b);
    //x = bckSubs(G, y);
    
    //v_ne = x[len-1];
    //v_ea = x[len*(len+1)/2-1];
    //v_sw = x[dim-len];
   
    //cout << v_ne <<endl; 
    //cout << v_ea <<endl; 
    //cout << v_sw <<endl; 

    // -----iterative method-----
    //cout << "---Conjugate Gradient Method---" << endl;

    start = clock();
    iter = cg(A, b, x_iter, 10000, 4.16091e-08);
    end = clock();

    v_ne1 = x_iter[len-1];
    v_ea1 = x_iter[len*(len+1)/2-1];
    v_sw1 = x_iter[dim-len];
     
    //cout << v_ne1 <<endl; 
    //cout << v_ea1 <<endl; 
    //cout << v_sw1 <<endl; 
    
    // print out the results
    cout << "iter: " << iter << endl;
    cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << endl;
    //cout << "v_ne_err: " << fabs(v_ne-v_ne1) << endl;
    //cout << "v_ea_err: " << fabs(v_ea-v_ea1) << endl;
    //cout << "v_sw_err: " << fabs(v_sw-v_sw1) << endl;    
}

