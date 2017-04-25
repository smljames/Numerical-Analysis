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
    
    fstream file, f301;
    string str;

    //calculate the dimension
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
    VEC Y301(301);
    
    //f?.dat
    file.open(argv[1]);
    file.seekg(3, ios::beg);
    for(i=0; i<dim; i++){
        file>>XDATA[i]>>YDATA[i];
    }
    
    //f301.dat
    f301.open("f301.dat");
    f301.seekg(3, ios::beg);
    for(i=0; i<301; i++){
        f301>>str>>Y301[i];
    }


    //-----lagrange interpolations-----

    VEC X(301);
    VEC Y(301);
    
    for(i=0; i<301; i++){
        X[i] = i+475;
        Y[i] = Lagrange(X[i], XDATA, YDATA);
    }


    //-----output-----

    VEC Y_151(151);
    VEC Y301_151(151);

    for(i=75; i<=225; i++){
        Y_151[i-75] = Y[i];
        Y301_151[i-75] = Y301[i];
    }

    cout<<max(Y-Y301)<<endl;
    cout<<max(Y_151-Y301_151)<<endl;
    
    return 0;
}
