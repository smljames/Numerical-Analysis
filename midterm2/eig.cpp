#include <iostream>
#include "MAT.h"
#include "VEC.h"

using namespace std;

int main(int argc, char* argv[]){
    int dim;

    cin>>dim;

    MAT A(dim);

    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            cin>>A[i][j];
        }
    }

    int iter, maxiter;
    VEC eig(dim);

    maxiter = 1000;
    iter = EVqr(A, 1e-09, maxiter);
    for(int i=0; i<dim; i++) eig[i] = A[i][i];

    print(eig);

    return 0;
}
