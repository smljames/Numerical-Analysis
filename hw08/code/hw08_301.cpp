/*
EE407002 hw08. Polynomial Interpolations
102061125 Chen Kuan-Chun
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include "MAT.h"
#include "VEC.h"

using namespace std;

int main(int argc, char* argv[])
{
    int i, j;
    int len, dim;

    //-----read the data file-----
    fstream file;
    string str;

    file.open(argv[1]);
    file.seekg(3, ios::beg);
    len = 0;
    while(file>>str){
        len++;
    }
    dim = len/2;
    file.close();


    VEC XDATA(dim);
    VEC YDATA(dim);

    file.open(argv[1]);
    file.seekg(3, ios::beg);
    for(i=0; i<dim; i++){
        file>>XDATA[i]>>YDATA[i];
    }
    YDATA.print();

    //-----lagrange interpolations-----
//    VEC X(301);
//    VEC Y(301);
    
//    for(i=0; i<301; i++){
//        X[i] = i+475;
//        Y[i] = Lagrange(X[i], XDATA, YDATA);
//    }

    return 0;
}
