/*
EE407002 hw09. Spline Interpolations
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
    file.close();
    dim = len/2;


    VEC XDATA(dim); // given x data
    VEC YDATA(dim); // given y data
    VEC Y301(301); // y coordinate of simulated waveform
    
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


    //-----spline interpolations-----

    VEC M(dim); // moments on each given x
    VEC Y(301); // output of spline interpolation
    splineM(dim, XDATA, YDATA, M); // calculate the moments
    
    for(i=475; i<=775; i++){
         Y[i-475] = spline(i, dim, XDATA, YDATA, M); // spline interpolation
    }


    //-----output-----

    VEC err(301);

    for(i=0; i<301; i++){
        err[i] = fabs(Y301[i]-Y[i]);
    }

    cout<<"Interpolated Values: "<<endl;
    print(Y);
    cout<<"Maximum Absolute Error: "<<endl;
    cout<<max(err)<<endl;
    
    return 0;
}
