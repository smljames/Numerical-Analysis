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

    print(A);
    cout<<"det(A) = "<< det(A)<<endl;
    return 0;
}
