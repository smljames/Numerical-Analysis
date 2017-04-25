/*
EE407002 hw03. Resistor Networks
102061125 Chen Kuan-Chun
*/

/*
Usage:
$ g++ hw03.cpp MAT.cpp VEC.cpp
$ ./a.out 10 (for example: 10 resistors per side)
*/

#include <iostream>
#include <iomanip>
#include <time.h>
#include "MAT.h"
#include "VEC.h"

using namespace std;

int main(int argc, char* argv[])
{
    /*
    len: # of nodes per side
    dim: total # of nodes
    vol: voltage
    g  : conductance
    */
    int len, dim, vol;
    double g;
    
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
       
    /*------------------------------------------*/

    MAT A(dim);
    VEC b(dim);
    
    double start, end;

    /*build the matrix and vector of the system*/
    for(int i=0; i<dim; i++){
        /*each corner*/
        if (i==0 || i==0+len-1 || i==dim-len || i==dim-1){
            A[i][i] = 2;
        }
        /*resistors not on the side*/
        else if (i>len && i<dim-len-1 && (i%len!=(len-1)) && (i%len!=0)){
            A[i][i] = 4;
        }
        /*rest of the resistors*/
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

    //b = A.tpose() * b;
    //A = A.tpose() * A;

    /*print out the matrix and b vector*/     
    
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(i == 1 | i == 7) break;
            cout << setw(2) << A[i][j] << ' ';;
        }
        cout << endl;
        //cout << '\t' << b[i] << endl;
    }
    
    
    /*-----------------------------------------*/
    
    /*use LU decomposition to solve the system*/
    VEC y(dim);
    VEC x(dim); //node voltage vector
    MAT G(dim); //LU decomposition matrix of A

    start = clock();
    G = luFact(A);
    y = fwdSubs(G, b);
    x = bckSubs(G, y);
    end = clock();
    
    /*print out the node voltage vector*/
    /*
    cout << "\n";
    for(int i=0; i<dim; i++){
        cout << x[i] << endl;
    }
    */
    
    /*chech the answer, Ax = b ?*/ 
    /*
    double sum;
    
    cout << "\n";

    for(int i=0; i<dim; i++){
        sum = 0;
        for(int j=0; j<dim; j++){
            sum += A[i][j]*x[j];
        }
        cout << sum << endl;
    }
    */
    
    /*-----------------------------------------*/

    double total_i; //total current i
    double eq_r; //equivalent resistance
    double v_ne, v_ea, v_sw;
    
    total_i = ((vol-x[(len-1)/2+1])+(vol-x[(len-1)/2-1])+(vol-x[(len-1)/2+len]))*g;
    eq_r = vol/total_i;
    v_ne = x[len-1];
    v_ea = x[len*(len+1)/2-1];
    v_sw = x[dim-len];
    
    cout << len-1 << " resistors per side" << endl;
    cout << "equivalent resistance: " << eq_r << endl;
    cout << "v_ne: " << v_ne << endl;
    cout << "v_ea: " << v_ea << endl;
    cout << "v_sw: " << v_sw << endl;

    cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << "\n" << endl;
}

