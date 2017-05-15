/*
EE407002 hw10. Numerical Integration
102061125 Chen Kuan-Chun
*/

/*
usage:
$ g++ hw10.cpp MAT.cpp VEC.cpp
$ ./a.out min max segment order

min: min of integral range
max: max of integral range
segment: divide whole range to several regions
order: n-th order
*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "MAT.h"
#include "VEC.h"
using namespace std;

double func(double x);

int main(int argc, char* argv[]){
    double result;
    double h;
    double x_max;
    double x_min;
    int i;
    int segment;
    int order;

    x_min = atof(argv[1]);
    x_max = atof(argv[2]);
    segment = atoi(argv[3]);
    order = atoi(argv[4]);
    
    VEC x(segment+1);
    VEC y(segment+1);
    
    h = (x_max-x_min) / segment;
    for(i=0; i<=segment; i++){
        x[i] = x_min + i*h;
        y[i] = func(x[i]);
    }

    result = integral(func, 0, 2, segment, order);
    printf("The integral of exp(x) from %g to %g is %g\n", x_min, x_max, result);
    printf("Absolute error is %g\n", fabs(exp(2)-exp(0)-result));
    
    return 0;
}

double func(double x){
    return exp(x);
}
