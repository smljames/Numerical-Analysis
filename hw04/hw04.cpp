/*
EE407002 hw04. Linear Iterative Methods
102061125 Chen Kuan-Chun
*/

/*
Usage:
$ g++ hw04.cpp MAT.cpp VEC.cpp
$ ./a.out 10 (for example: 10 resistors per side)
*/

#include <iostream>
#include <time.h>
#include <math.h>
#include "MAT.h"
#include "VEC.h"

using namespace std;

int main(int argc, char* argv[])
{
    /*----------initialize----------*/

    int len, dim, vol; //len: # of nodes per side
                       //dim: total # of nodes
                       //vol:voltage
    double g; //g: conductance
    
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
    
    
       
    /*----------build the matrix and vector of the system----------*/

    MAT A(dim);
    VEC b(dim);
    
    double start, end;

    for(int i=0; i<dim; i++){
        /*each corner of resistor network*/
        if (i==0 || i==0+len-1 || i==dim-len || i==dim-1){
            A[i][i] = 2;
        }
        /*resistors not on the side of the resistor network*/
        else if (i>len && i<dim-len-1 && (i%len!=(len-1)) && (i%len!=0)){
            A[i][i] = 4;
        }
        /*rest of the resistors of the resistor network*/
        else{
            A[i][i] = 3;
        }
    }
    
    /*for each resistor connecting node beside*/
    for(int i=0; i<dim; i++){
        if (i%len != 0 && (i-1) >= 0)    A[i][i-1] = -1;
        if ((i+1)%len != 0 && (i+1) <= (dim-1)) A[i][i+1] = -1;
        if (i-len >= 0)  A[i][i-len] = -1;
        if (i+len <= dim-1) A[i][i+len] = -1;
    }
    
    int vol_node = (len-1)/2; //the node connecting to the fixed voltage
    int gd_node = (dim-1)-(len-1)/2; //ground node
    
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
    


    /*----------calculate node voltages vector----------*/
    
    VEC y(dim);
    VEC x(dim); //node voltage vector
    VEC x_iter(dim); //initail vector for iterative solutions 
    int iter; //# of iteration that iterative solutions use
    double total_i; //total current i
    double eq_r; //equivalent resistance
    double v_ne, v_ea, v_sw;

    cout << len-1 << " resistors per side" << endl;
    
    /*-----iterative method-----*/
    
    cout << "---Jacobi Method---" << endl;
    iter = jacobi(A, b, x_iter, 15000, 0.0000001);
    
    /* 
    cout << "---Jacobi Over-Relaxation Method---" << endl;
    iter = jacobi_over_relaxation(A, b, x_iter, 10000, 0.95, 0.0000001);
    */
    /*
    cout << "---Gauss-Seidel Method---" << endl;
    iter = gaussSeidel(A, b, x_iter, 10000, 0.0000001);
    */
    /*
    cout << "---Successive Over-Relaxation Method---" << endl;
    iter = successive_over_relaxation(A, b, x_iter, 10000, 1.5, 0.0000001);
    */
    /*
    cout << "---Symmetric Guass-Seidel Method---" << endl;
    iter = sgs(A, b, x_iter, 10000, 0.0000001);
    */
    /*
    cout << "---Symmetric Successive Over-Relaxation Method---" << endl;
    iter = symmetric_successive_over_relaxation(A, b, x_iter, 10000, 1.5, 0.0000001);
    */

    v_ne = x_iter[len-1];
    v_ea = x_iter[len*(len+1)/2-1];
    v_sw = x_iter[dim-len];
    
    cout << "v_ne: " << v_ne << endl;
    cout << "v_ea: " << v_ea << endl;
    cout << "v_sw: " << v_sw << endl;
    cout << "iter: " << iter << endl;

    /*-----LU decomposition method-----*/    
    cout << "---LU decomposition---" << endl;

    A = luFact(A);
    y = fwdSubs(A, b);
    x = bckSubs(A, y);
    
    v_ne = x[len-1];
    v_ea = x[len*(len+1)/2-1];
    v_sw = x[dim-len];
    
    cout << "v_ne: " << v_ne << endl;
    cout << "v_ea: " << v_ea << endl;
    cout << "v_sw: " << v_sw << endl;    
}

