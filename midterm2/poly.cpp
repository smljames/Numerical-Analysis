#include <iostream>
#include <math.h>
#include "MAT.h"
#include "VEC.h"

using namespace std;

double func1(double x);
double func1_(double x);

int main(int argc, char* argv[]){
    int dim;
    
    cin>>dim;

    VEC p(dim);
    VEC p1(dim-1);
    
    for(int i=0; i<dim; i++){
        cin>>p[i];
    }
    print(p);

    double x, y;
    
    x = 2;
    y = polyval(p, x);
    cout<<"x = "<<x<<endl;
    cout<<"f(x) = ";
    for(int i=dim-1; i>=0; i--) cout<<p[i]<<" ";
    cout<<endl;
    p1 = polydiff(p);
    cout<<"f_(x) = ";
    for(int i=dim-2; i>=0; i--) cout<<p1[i]<<" ";
    cout<<endl;

    cout<<"f1(x) = "<<func1(x)<<endl;
    cout<<"f1'(x) = "<<func1_(x)<<endl;

    return 0;
}

double func1(double x){
    double y;
    y = x*x + 2*x + log10(x);
    return y;
}

double func1_(double x){
    double y;
    y = 2*x + 2 + 1/x;
    return y;
}
