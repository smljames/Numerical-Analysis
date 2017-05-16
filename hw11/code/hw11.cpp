/*
EE407002 hw11. Diode Networks
102061125 Chen Kuan-CHun
*/

#include <iostream>
#include <math.h>
#include "VEC.h"
#include "MAT.h"

using namespace std;


// -----global variables-----
double Is = 1; 
double phi0 = 0.026;
double RL = 0.01;   // resistence of RL
double V;           // input voltage
VEC func1(VEC &x);  // two KCL equations
VEC func2(VEC &x);



int main(int argc, char* argv[]){
    double h;               // delta of x for the difference approximation
    double tol;             // tolerence
    int p;                  // step for new calculation of Jacobian matrix
    int maxiter;
    double v1, v2;          // node voltage v1, v2
    double T1, T2, T3, T4;  // temperature variables for problem 2
    double phi1, phi2, phi3, phi4; // phi variables for problem 2
    double i1, i2, i3, i4;  // current for each diodes
    double iR;              // current for resistor RL
    VEC v(2);               // [v1, v2] vector



    // -----set up parameters-----
    h = 1e-05;
    p = 1;
    tol = 1e-09;
    maxiter = 1000;
   
   
    
    // -----input voltage: 0, 0.02, ..., 1-----
    V = 0;
    while(V <= 1.001){
        // Cyclic Jacobian Updates method 
        cj(func1, v, h, p, tol, maxiter);
    
        v1 = v[0];
        v2 = v[1];
        
        i1 = Is*(exp((V-v1)/phi0)-1);
        i2 = Is*(exp((-1)*v1/phi0)-1);
        i3 = Is*(exp((v2-V)/phi0)-1);
        i4 = Is*(exp(v2/phi0)-1);
        iR = (v1-v2)/RL;
 
        cout << V << "\t" << v1 << "\t" << v2 << "\t" 
             << i1 << "\t" << i2 << "\t" << i3 << "\t" 
             << i4 << "\t" << iR << endl;

        V += 0.02;
    }
    

    // -----input voltage: 0, -0.02, ..., -1-----
    V = 0;
    while(V >= -1.001){
        // Cyclic Jacobian Updates method 
        cj(func1, v, h, p, tol, maxiter);
    
        v1 = v[0];
        v2 = v[1];
        
        i1 = Is*(exp((V-v1)/phi0)-1);
        i2 = Is*(exp((-1)*v1/phi0)-1);
        i3 = Is*(exp((v2-V)/phi0)-1);
        i4 = Is*(exp(v2/phi0)-1);
        iR = (v1-v2)/RL;
 
        cout << V << "\t" << v1 << "\t" << v2 << "\t" 
             << i1 << "\t" << i2 << "\t" << i3 << "\t" 
             << i4 << "\t" << iR << endl;
        
        V -= 0.02;
    }
    
    return 0;
}

VEC func1(VEC &x){
    double v1 = x[0];
    double v2 = x[1];
    
    double i1 = Is*(exp((V-v1)/phi0)-1);
    double i2 = Is*(exp(((-1)*v1)/phi0)-1);
    double i3 = Is*(exp((v2-V)/phi0)-1);
    double i4 = Is*(exp(v2/phi0)-1);
    double iR = (v1-v2)/RL;
    
    VEC result(len(x));
     
    result[0] = i1 + i2 - iR;
    result[1] = i3 + i4 - iR;

    return result;
}

VEC func2(VEC &x){
    double v1 = x[0];
    double v2 = x[1];
    
    double T1 = x[2];
    double T2 = x[3];
    double T3 = x[4];
    double T4 = x[5];
    
    double phi1 = phi0 * T1 / 300;
    double phi2 = phi0 * T2 / 300;
    double phi3 = phi0 * T3 / 300;
    double phi4 = phi0 * T4 / 300;
    
    double i1 = Is*(exp((V-v1)/phi1)-1);
    double i2 = Is*(exp(((-1)*v1)/phi2)-1);
    double i3 = Is*(exp((v2-V)/phi3)-1);
    double i4 = Is*(exp(v2/phi4)-1);
    double iR = (v1-v2)/RL;

    VEC result(len(x));

    result[0] = i1 + i2 - iR;
    result[1] = i3 + i4 - iR;

    result[2] = 300 + 2 * i1 * (V-v1) - T1;
    result[3] = 300 + 2 * i2 * ((-1)*v1) - T2;
    result[4] = 300 + 2 * i3 * (v2-V) - T3;
    result[5] = 300 + 2 * i4 * (v2) - T4;

    return result;
}
