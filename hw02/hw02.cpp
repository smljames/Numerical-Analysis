/*
EE407002 hw02 LU decomposition, forward substitution & backward substitution
102061125 Chen Kuan-Chun
*/

#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <time.h>
#include "MAT.h"
#include "VEC.h"

using namespace std;

/*
README

compile:
    $ g++ hw02.cpp MAT.cpp VEC.cpp

excute:
    $ ./a.out "the matrix to solve"
    (ex. $ ./a.out m3.dat)

You can just: $ bash run_all_dat.sh
It will feed in the matrices from m3.dat to m10.dat
*/

int main(int argc, char* argv[])
{
    int dim; //dimension

    fstream fin(argv[1]); //open the file
    fin >> dim;

    MAT A(dim);
    VEC b(dim);

    /*read the matrix*/
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            fin >> A[i][j];
        }
    }

    /*read the vector*/
    for(int i=0; i<dim; i++){
        fin >> b[i];
    }
    
    double lu_start, lu_end; //variables for recording the execution time

    lu_start = clock();
    A = luFact(A);
    lu_end = clock();

    /*print out the orthogonal matrix*/
    /*
    for(int i=0; i<dim; i++){
        for (int j=0; j<dim; j++){
            cout << setprecision(5)  << A[i][j] << "\t";
        }
        cout << endl;
    }
    for(int i=0; i<dim; i++){
        cout << setprecision(5) <<  b[i] << endl;
    }
    */

    VEC y(dim);
    VEC x(dim);
    double subs_start, subs_end; //variables for recording the execution time

    subs_start = clock();
    y = fwdSubs(A, b);
    x = bckSubs(A, y);
    subs_end = clock();

    //cout << "LU decomposition time: " << (float)(lu_end - lu_start)/CLOCKS_PER_SEC << " s" << endl;
    //cout << "forward & backward substitution time: " << (float)(subs_end - subs_start)/CLOCKS_PER_SEC << " s"  << endl;
    
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    for(int i=0; i<dim; i++){
        cout << x[i] << endl;
    }
} 
