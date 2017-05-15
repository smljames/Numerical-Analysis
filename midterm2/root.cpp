#include <iostream>
#include <math.h>
#include "MAT.h"
#include "VEC.h"

using namespace std;

int main(int argc, char* argv[]){
    int dim;
    
    cin>>dim;

    VEC p(dim);
    
    for(int i=0; i<dim; i++){
        cin>>p[i];
    }
    print(p);

    VEC root(dim-1);
    root = polyroot(p, 4, 1e-09, 1000);
    print(root);
    
}
