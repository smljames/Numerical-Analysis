/*
EE407002 102061125 Chen Kuan-Chun
*/

// VEC class functions
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "VEC.h"

using namespace std;

VEC::VEC(int n) // uninit constructor
{
    dim=n;
    val=(double *)calloc(n,sizeof(double));
}
VEC::VEC(const VEC &v1) // copy constructor
{
    dim=v1.dim;
    val=(double *)calloc(dim,sizeof(double));
    for (int i=0; i<dim; i++) {
        val[i]=v1.val[i];
    }
}
VEC::VEC(int n,double *v) // init constructor
{
    dim=n;
    val=(double *)calloc(n,sizeof(double));
    for (int i=0; i<n; i++) val[i]=v[i];
}
VEC::~VEC() // destructor
{
    free(val);
}
int VEC::len() // return dimension of the vector
{
    return dim;
}
VEC &VEC::operator-() // unary operator - : negative value
{
    for (int i=0; i<dim; i++) val[i]=-val[i];
    return *this;
}
VEC &VEC::operator=(const VEC v1) // assignment
{
    dim=v1.dim;
    for (int i=0; i<dim; i++) {
        val[i]=v1.val[i];
    }
    return *this;
}
VEC &VEC::operator+=(const VEC v1) // V += v1
{
    for (int i=0; i<dim; i++) {
        val[i]+=v1.val[i];
    }
    return *this;
}
VEC &VEC::operator-=(const VEC v1) // V -= v1
{
    for (int i=0; i<dim; i++) {
        val[i]-=v1.val[i];
    }
    return *this;
}
VEC &VEC::operator*=(double a) // V *= dbl
{
    for (int i=0; i<dim; i++) {
        val[i]*=a;
    }
    return *this;
}
VEC &VEC::operator/=(double a) // V /= dbl
{
    for (int i=0; i<dim; i++) {
        val[i]/=a;
    }
    return *this;
}
VEC VEC::operator+(const VEC v1) // V + v1
{
    VEC s(*this);
    for (int i=0; i<dim; i++) s.val[i]+=v1.val[i];
    return s;
}
VEC VEC::operator-(const VEC v1) // V - v1
{
    VEC s(*this);
    for (int i=0; i<dim; i++) s.val[i]-=v1.val[i];
    return s;
}
double VEC::operator*(VEC v1) // inner product
{
    double a=0;
    VEC s(*this);
    for (int i=0; i<dim; i++) 
        a += s.val[i] * v1.val[i];
    return a;
}
VEC VEC::operator*(double a) // V * dbl
{
    VEC s(*this);
    for (int i=0; i<dim; i++) s.val[i]*=a;
    return s;
}
VEC VEC::operator/(double a) // V / dbl
{
    VEC s(*this);
    for (int i=0; i<dim; i++) s.val[i]/=a;
    return s;
}
double &VEC::operator[](int n) // indexing
{
    if (n<0) n=0;
    else if (n>=dim) n=dim-1;
    return val[n];
}
VEC operator*(double a, const VEC v1) // dbl x V
{
    VEC s(v1);
    for (int i=0; i<v1.dim; i++) s.val[i]=v1.val[i]*a;
    return s;
}
VEC *newVEC(int n) // allocate a dynamic VEC
{
    VEC *vptr;
    vptr=(VEC *)malloc(sizeof(VEC));
    vptr->dim=n;
    vptr->val=(double*)calloc(n,sizeof(double));
    return vptr;
}

double one_norm(VEC v1)
{   
    double output=0;
    VEC s(v1);

    for(int i=0; i<s.dim; i++){
        output += fabs(s[i]);
    }
    return output;
}

double two_norm(VEC v1)
{
    double output=0;
    VEC s(v1);

    output = s * s;
    output = sqrt(output);

    return output;
}

double inf_norm(VEC v1)
{
    double output=0;
    VEC s(v1);

    output = fabs(s[0]);
    for(int i=1; i<s.dim; i++){
        if(fabs(s[i]) > output) output = fabs(s[i]);
    }

    return output;
}
void print(VEC v1)
{
    VEC s(v1);
    for(int i=0; i<s.dim; i++){
        cout<<s[i]<<"\n";
    }
    cout<<endl;
}
VEC sort(VEC v1)
{
    VEC s(v1);
    int i, j, k;
    double a = 0;

    for(i=0; i<s.dim; i++){
        j = i;
        for(k=i; k<s.dim; k++){
            if(s[j] > s[k]){
                j = k;
            }
        }
        a = s[i];
        s[i] = s[j];
        s[j] = a;
    }
    return s;
}
double max(VEC v1)
{
    VEC s(v1);
    double max = fabs(s[0]);

    for(int i=1; i<s.dim; i++){
        if(fabs(s[i]) > max) max = fabs(s[i]);
    }
    return max;
}
int ceil_ind(double x, VEC &X)
{
    for(int i=1; i<X.dim; i++){
        if(x <= X[i] and x >= X[i-1]) return i;
    }
}
VEC inverse(VEC &v1)
{
    int dim = v1.dim;
    VEC s(dim);

    for(int i=0; i<dim; i++) s[i] = v1[dim-1-i];
    
    return s;
}
int len(VEC &v1)
{
    return v1.dim;
}
double integral(double(*f)(double), double min, double max, int N, int order){
    double i;
    double j;
    double k;
    double result = 0;
    double segment = (max - min) / N; // width of each regions
    double delta = segment / order; // width of each subregions
    VEC w(order+1); // weight vector

    switch(order){
        case 1:
            w[0] = 1/2.;
            w[1] = w[0];
            break;
        case 2:
            w[0] = 1/3.;
            w[1] = 4/3.;
            w[2] = w[0];
            break;
        case 3:
            w[0] = 3/8.;
            w[1] = 9/8.;
            w[2] = w[1];
            w[3] = w[0];
            break;
        case 4:
            w[0] = 14/45.;
            w[1] = 64/45.;
            w[2] = 24/45.;
            w[3] = w[1];
            w[4] = w[0];
            break;
        case 5:
            w[0] = 95/288.;
            w[1] = 375/288.;
            w[2] = 250/288.;
            w[3] = w[2];
            w[4] = w[1];
            w[5] = w[0];
            break;
        case 6:
            w[0] = 41/140.;
            w[1] = 216/140.;
            w[2] = 27/140.;
            w[3] = 272/140.;
            w[4] = w[2];
            w[5] = w[1];
            w[6] = w[0];
    }

    for(i=0; i<N; i++){
        for(j=min+i*segment, k=0; k<=order; j+=delta, k++){
            result += f(j) * w[k];
        }
    }
    result = delta * result;
    
    return result;
}
