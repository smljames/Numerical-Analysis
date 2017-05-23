#include <iostream>
#include "MAT.h"
#include "VEC.h"

using namespace std;

double V = 1;
double R = 1;
double C = 1;
double RC = R*C;
double h = RC/5;
VEC func(VEC &x);

int main(){
    int iter, maxiter;
    VEC v(1);

    maxiter = 300;

    iter = FwdEuler(func, v, h, maxiter);

    return 0;
}

VEC func(VEC &x){
    VEC result(len(x));

    result[0] = x[0] + h*(V-x[0])/(RC);

    return result;
}
