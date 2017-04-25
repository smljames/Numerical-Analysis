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

double VEC::one_norm()
{   
    double output=0;
    VEC s(*this);

    for(int i=0; i<s.dim; i++){
        output += fabs(s[i]);
    }
    return output;
}

double VEC::two_norm()
{
    double output=0;
    VEC s(*this);

    output = s * s;
    output = sqrt(output);

    return output;
}

double VEC::inf_norm()
{
    double output=0;
    VEC s(*this);

    output = fabs(s[0]);
    for(int i=1; i<s.dim; i++){
        if(fabs(s[i]) > output) output = fabs(s[i]);
    }

    return output;
}
void VEC::print()
{
    VEC s(*this);
    for(int i=0; i<s.dim; i++){
        cout<<s[i]<<"\t";
    }
    cout<<endl;
}
VEC VEC::sort()
{
    VEC s(*this);
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


