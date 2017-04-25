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
    MAT G(dim);
    VEC b(dim);
    
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

    G = A;

    /*print out the matrix and b vector*/     
    /* 
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            cout << A[i][j] << ' ';
        }
        cout << '\t' << '\t' << b[i] << endl;
    }
    */
    


    /*----------calculate node voltages vector----------*/
    
    VEC y(dim);
    VEC x(dim); //node voltage vector
    VEC x_iter1(dim), x_iter2(dim), x_iter3(dim); //initail vector for iterative solutions 
    int iter1, iter2, iter3; //# of iteration that iterative solutions use
    double tol1, tol2, tol3;
    double total_i; //total current i
    double eq_r; //equivalent resistance
    double v_ne, v_ea, v_sw;
    double v_ne1, v_ea1, v_sw1;
    double v_ne2, v_ea2, v_sw2;
    double v_ne3, v_ea3, v_sw3;
    int start1, end1; //record time
    int start2, end2;
    int start3, end3;

    G = luFact(G);
    y = fwdSubs(G, b);
    x = bckSubs(G, y);

    cout << len-1 << " resistors per side" << endl;
    cout << "[ 1 - norm ]" << endl;
    cout << "---cpu time---" << endl;
    //cout << "---tol---" << endl;

    /*-----iterative method-----*/
    
    /*---Jacobi Method---*/
    start1=clock();
    iter1 = jacobi(A, b, x_iter1, 100000, 5.50378e-08);
    end1=clock();
    cout << (float)(end1-start1)/CLOCKS_PER_SEC << endl;
    //tol1 = tol_jacobi(A, b, x_iter1, x);
    //cout << tol1 << endl;
    
    /* 
    cout << "---Jacobi Over-Relaxation Method---" << endl;
    iter = jacobi_over_relaxation(A, b, x_iter, 10000, 0.95, 0.0000001);
    */

    /*---Gauss-Seidel Method---*/
    start2=clock();
    iter2 = gaussSeidel(A, b, x_iter2, 100000, 1.08374e-07);
    end2=clock();
    cout << (float)(end2-start2)/CLOCKS_PER_SEC << endl;
    //tol2 = tol_gaussSeidel(A, b, x_iter2, x);
    //cout << tol2 << endl;

    /*
    cout << "---Successive Over-Relaxation Method---" << endl;
    iter = successive_over_relaxation(A, b, x_iter, 10000, 1.5, 0.0000001);
    */

    /*---Symmetric Guass-Seidel Method---*/
    start3=clock();
    iter3 = sgs(A, b, x_iter3, 100000, 2.15797e-07);
    end3=clock();
    cout << (float)(end3-start3)/CLOCKS_PER_SEC << endl;
    //tol3 = tol_sgs(A, b, x_iter3, x);
    //cout << tol3 << endl;
    
    
    cout << "---iteration---" << endl;
    cout << iter1 << endl;
    cout << iter2 << endl;
    cout << iter3 << endl;
    /*
    cout << "---Symmetric Successive Over-Relaxation Method---" << endl;
    iter = symmetric_successive_over_relaxation(A, b, x_iter, 10000, 1.5, 0.0000001);
    */
    /*
    v_ne = x_iter[len-1];
    v_ea = x_iter[len*(len+1)/2-1];
    v_sw = x_iter[dim-len];
    
    cout << "v_ne: " << v_ne << endl;
    cout << "v_ea: " << v_ea << endl;
    cout << "v_sw: " << v_sw << endl;
    cout << "iter: " << iter << endl;
    */

    /*-----LU decomposition method-----*/    
    
    //cout << "---LU decomposition---" << endl;
    /*
    A = luFact(A);
    y = fwdSubs(A, b);
    x = bckSubs(A, y);
    */

    v_ne = x[len-1];
    v_ea = x[len*(len+1)/2-1];
    v_sw = x[dim-len];
        
    v_ne1 = x_iter1[len-1];
    v_ea1 = x_iter1[len*(len+1)/2-1];
    v_sw1 = x_iter1[dim-len];
    
    v_ne2 = x_iter2[len-1];
    v_ea2 = x_iter2[len*(len+1)/2-1];
    v_sw2 = x_iter2[dim-len];
    
    v_ne3 = x_iter3[len-1];
    v_ea3 = x_iter3[len*(len+1)/2-1];
    v_sw3 = x_iter3[dim-len];
    

    cout << "---error---" << endl;

    cout << "jacobi error" << endl;
    cout << fabs(v_ne-v_ne1) << endl;
    cout << fabs(v_ea-v_ea1) << endl;
    cout << fabs(v_sw-v_sw1) << endl;

    cout << "gauss-seidel error" << endl;
    cout << fabs(v_ne-v_ne2) << endl;
    cout << fabs(v_ea-v_ea2) << endl;
    cout << fabs(v_sw-v_sw2) << endl;

    cout << "symmetric gauss-seidel error" << endl;
    cout << fabs(v_ne-v_ne3) << endl;
    cout << fabs(v_ea-v_ea3) << endl;
    cout << fabs(v_sw-v_sw3) << endl;

    //cout << "v_ne: " << v_ne << endl;
    //cout << "v_ea: " << v_ea << endl;
    //cout << "v_sw: " << v_sw << endl;    
        
}

