/*
EE407002 hw07. Matrix Eigenvalues
102061125 Chen Kuan-Chun
*/

/*
Usage:
$ g++ hw07.cpp MAT.cpp VEC.cpp
$ ./a.out < m3.dat 
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
    int dim;
    int maxiter = 10000;
    int iter1, iter2;
    int start1, end1;
    int start2, end2;


    cin>>dim;
    MAT A(dim);
    MAT B(dim);
    MAT C(dim);
    
    // -----fill the components in the matrix A-----
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            cin>>A[i][j];
        }
    }

    VEC eig(dim);
    VEC eig_shift(dim);
    VEC largest(3);
    VEC smallest(3);
    VEC largest_shift(3);
    VEC smallest_shift(3);

    B = A;
    C = A;
    
    start1 = clock(); 
    iter1 = EVqr(B, 1e-09, maxiter);
    end1 = clock();

    start2 = clock();
    iter2 = EVqrShifted(C, 0.5, 1e-09, maxiter);
    end2 = clock();


    // -----record eigenvalues-----

    for(i=0; i<dim; i++){
        eig[i] = B[i][i];
    }

    for(i=0; i<dim; i++){
        eig_shift[i] = C[i][i];
    }

    eig = eig.sort();
    eig_shift = eig_shift.sort();

    for(i=0; i<3; i++){
        largest[i] = eig[dim-1-i];
        largest_shift[i] = eig_shift[dim-1-i];
        smallest[i] = eig[i];
        smallest_shift[i] = eig_shift[i];
    }


    // -----outputs-----

    cout<<"[ dim ]"<<endl;
    cout<<dim<<endl;
    cout<<"[ time ]"<<endl;
    cout<<(float)(end1-start1)/CLOCKS_PER_SEC<<endl;
    cout<<(float)(end2-start2)/CLOCKS_PER_SEC<<endl;
    cout<<"[ iter ]"<<endl;
    cout<<iter1<<endl;
    cout<<iter2<<endl;
    cout<<"[ 3 largest eigenvalues ]"<<endl;
    largest.print();
    largest_shift.print();
    cout<<"[ 3 smallest eigenvalues ]"<<endl;
    smallest.print();
    smallest_shift.print();

    cout<<endl;    
    return 0;
}
