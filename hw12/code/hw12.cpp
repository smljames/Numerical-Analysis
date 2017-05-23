/*
EE407002 hw12. RLC Circuit.
102061125 Chen Kuan-CHun
*/

#include <iostream>
#include "VEC.h"
#include "MAT.h"

using namespace std;

double V = 1;
double R = 1;
double L = 1;
double C = 1;
double h = 0.01;

VEC func(VEC &x, int method);

int main(int argc, char* argv[]){
    double t;
    VEC v(3);

    v[0] = 0; // iL
    v[1] = 1; // v1
    v[2] = 0; // v2

    t = 0;
    while(t <= 10){
        v = func(v, 2); 
        cout<<v[2]<<endl;
        t += h;
    }

    return 0;
}

VEC func(VEC &x, int method){
    VEC result(3);
    
    if(method == 1){
        result[0] = x[0] + h*(x[1]-x[2])/L;
        result[1] = -1*result[0]*R + V;
        result[2] = x[2] + h*x[0]/C;
    }
    else if(method == 2){
        result[0] = (x[0] + h*(V-x[2])/L) / (1+ h*R/L + h*h/(L*C));
        result[1] = -1*result[0]*R + V;
        result[2] = x[2] + h*result[0]/C;
    }
    else if(method == 3){
        result[0] = ((1 - h*h/(4*L*C)) * x[0] + h * (V-2*x[2]+x[1]) / (2*L)) \
                    / (1 + h*R/(2*L) + h*h/(4*L*C));
        result[1] = -1*result[0] * R + V;
        result[2] = x[2] + h * (result[0] + x[0]) / (2*C);
    }

    return result;
}
