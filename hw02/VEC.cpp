// VEC class functions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VEC.h"
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
double one_norm(){
    double a = 0;
    VEC s(*this);
    for(int i=0; i<s.dim; i++){
        a += fabs(s[i]);
    }
    return a;
}
double two_norm(){
    double a = 0;
    VEC s(*this);
    for(int i=0; i<s.dim; i++){
        a += (s[i]) * (s[i]);
    }
    a = sqrt(a);
    return a;
}
double inf_norm(){
    double a = 0;
    VEC s(*this);
    for(int i=0; i<s.dim; i++){
        if(fabs(s[i]) > a) a = fabs(s[i]);
    }
    return a;
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