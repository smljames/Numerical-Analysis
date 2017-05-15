#include <iostream>
#include <math.h>
#include "MAT.h"
#include "VEC.h"
using namespace std;

double func(double x);


int main(int argc, char* argv[]){
    double result;
    result = integral(func, 0, 2, 12, 2);
    cout<<"The integral of exp(x) from 0 ot 2 is "<<result<<endl;
    cout<<"error = "<<fabs(func(2)-1.0-result)<<endl;

    return 0;
}

double func(double x){
    return exp(x);
}
