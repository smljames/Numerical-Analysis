/*
EE407002 Project. Polynomial Roots Finder.
102061125 Kuan-Chun Chen
*/

/*
Usage:
$ g++ proj.cpp Complex.cpp VEC.cpp
$ ./a.out "polynomial file"

Example:
$ ./a.out ta1.dat

"ta1.dat"
3
1 x^3
-6 x^2
11 x^1
-6 x^0
*/

#include <stdio.h>
#include <stdlib.h>
#include "Complex.h"
#include "VEC.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]){
    // -----read the file-----
    fstream file;
    string str;
    int len;
    int dim;
    double tmp;

    file.open(argv[1]);
    file.seekg(2, ios::beg);
    len = 0;
    while(file>>str){
        len++;
    }
    file.close();
    dim = len/2;



    VEC p(dim);
    double max=0;

    file.open(argv[1]);
    file.seekg(2, ios::beg);
    for(int i=dim-1; i>=0; i--){
        file>>tmp>>str;
        p[i] = rtoc(tmp);
    }

    // -----find roots-----
    VEC root(dim-1);
    Complex x(1,1);
    double tol = 1e-10;
    double maxiter = 10000;

    
    root = NewtonHorner(p, x, tol, maxiter);
    print(root);
    
    

    
    root = LinsQuadratic(p, tol, maxiter);
    print(root);
    
    

    
    root = Bairstow(p, tol, maxiter);
    print(root);
    

    VEC ans(dim-1);
    
    for(int i=0; i<dim; i++){
         ans[i] = polyval(p, root[i]);
         ans[i] = rtoc(cfabs(ans[i]));
    }
    //print(ans);
    
    return 0;
}
