#include <iostream>
#include <math.h>
#include "MAT.h"
#include "VEC.h"
using namespace std;

double func(double x);
double func_(double x);

int main(int argc, char*argv[]){
    double x;
    double y;
    double step;
    double tol;
    int maxiter;

    maxiter = 100;
    x = 2;
    step = 0.5;
    tol = 1e-09;

    y = newton(func, func_, x, step, tol, maxiter);
    cout<<y<<endl;
    
    return 0;
}

double func(double x){
    double result;
    double x2 = x*x;

    result = x2 + 2*x - 3;
    
    return result;
}

double func_(double x){
    double result;

    result = 2*x + 2;

    return result;
}
