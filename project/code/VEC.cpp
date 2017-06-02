// VEC class functions
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "VEC.h"
#include "Complex.h"

using namespace std;

VEC::VEC(int n) // uninit constructor
{
    dim=n;
    val=(Complex *)calloc(n,sizeof(Complex));
}
VEC::VEC(const VEC &v1) // copy constructor
{
    dim=v1.dim;
    val=(Complex *)calloc(dim,sizeof(Complex));
    for (int i=0; i<dim; i++) {
        val[i]=v1.val[i];
    }
}
VEC::~VEC() // destructor
{
    free(val);
}
int VEC::len() // return dimension of the vector
{
    return dim;
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
Complex VEC::operator*(VEC v1) // inner product
{
    Complex a(0,0);
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
Complex &VEC::operator[](int n) // indexing
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
    vptr->val=(Complex*)calloc(n,sizeof(Complex));
    return vptr;
}
void print(VEC v1)
{
    VEC s(v1);
    for(int i=0; i<s.dim; i++){
        print(s[i]);
    }
    cout<<endl;
}
int len(VEC &v1)
{
    return v1.dim;
}
Complex polyval(VEC &p, Complex x)
{
    int i, j;
    int dim = len(p);
    Complex ans(0,0);
    Complex tmp1(1,0);
    Complex tmp2(1,0);
 
    for(i=0; i<dim; i++){
        tmp1 = tmp2;
        for(j=0; j<i; j++){
            tmp1 *= x;
        }
        ans += p[i]*tmp1;
    }

    return ans;
}
VEC NewtonHorner(VEC p, Complex x, double tol, int maxiter)
{
    int j;
    int n = len(p)-1;
    int iter;
    double err;
    Complex f(0,0);
    Complex f_(0,0);
    Complex fx(0,0);
    VEC b(n+1);
    VEC c(n+1);
    VEC z(n);

    while(n>=1){
        err = 1+tol;
        iter = 0;
        while((err>=tol) and (iter<maxiter)){
            b[n] = p[n];
            c[n] = b[n];
            for(j=n-1; j>=0; j--)   b[j] = p[j] + x*b[j+1];
            for(j=n-1; j>=1; j--)   c[j] = b[j] + x*c[j+1];
            f = b[0];
            f_ = c[1];
            x = x - f/f_;
            err = cfabs(f);
            fx = polyval(p, x);
            iter++;
        }
        z[n-1] = x;
        for(j=0; j<n; j++)  p[j] = b[j+1];
        x = z[n-1];
        n--;
    }
    return z;
}
VEC LinsQuadratic(VEC a, double tol, int maxiter)
{
    int i, j;
    int dim = len(a);
    int iter;
    double err;
    VEC b(dim);
    VEC vzero(dim);
    VEC root(dim-1);
    Complex cinit(1,0);
    Complex p(0,0);
    Complex q(0,0);
    Complex R(0,0);
    Complex S(0,0);

    j = 0;
    while(dim > 0){
        p = a[dim-2]/a[dim-1];
        q = a[dim-3]/a[dim-1];
        if(dim > 3){
            b = vzero;
            b[0] = a[2];
            R = a[1];
            S = a[0];
            iter = 0;
            err = 1+tol;
            while(err > 1e-09 and iter < maxiter){
                p = p + R/b[0];
                q = q + S/b[0];
                for(i=dim-3; i>=0; i--){
                    b[i] = a[i+2] - p*b[i+1] - q*b[i+2];
                }
                R = a[1] - p*b[0] - q*b[1];
                S = a[0] - q*b[0];
                err = max(cfabs(R), cfabs(S));
                iter++;
            }
            if(p.r()*p.r()-4*q.r() < 0){
                root[j].x = -1*p.r()/2;
                root[j].y = sqrt(4*q.r()-p.r()*p.r())/2;
                j++;
                root[j].x = -1*p.r()/2;
                root[j].y = -1 * sqrt(4*q.r()-p.r()*p.r())/2;
                j++;
            }
            else{
                root[j].x = (-1*p.r() + sqrt(p.r()*p.r()-4*q.r()))/2;
                root[j].y = 0;
                j++;
                root[j].x = (-1*p.r() - sqrt(p.r()*p.r()-4*q.r()))/2;
                root[j].y = 0;
                j++;
            }
            a = b;
            dim -= 2;
        }
        else if(dim == 3){
            if(b[1].r()*b[1].r()-4*b[2].r()*b[0].r() < 0){
                root[j].x = (-1*b[1].r()) / (2*b[2].r());
                root[j].y = sqrt(4*b[2].r()*b[0].r()-b[1].r()*b[1].r()) / (2*b[2].r());
                j++;
                root[j].x = (-1*b[1].r()) / (2*b[2].r());
                root[j].y = -1 * sqrt(4*b[2].r()*b[0].r()-b[1].r()*b[1].r()) / (2*b[2].r());
                j++;
            }
            else{
                root[j].x = (-1*b[1].r() + sqrt(b[1].r()*b[1].r()-4*b[2].r()*b[0].r())) / (2*b[2].r());
                root[j].y = 0;
                j++;
                root[j].x = (-1*b[1].r() - sqrt(b[1].r()*b[1].r()-4*b[2].r()*b[0].r())) / (2*b[2].r());
                root[j].y = 0;
                j++;
            }
            dim -= 3;
        }
        else if(dim == 2){
            root[j].x = -1*b[0].r()/b[1].r();
            root[j].y = 0;
            dim -= 2;
        }
    }

    return root;
}
VEC Bairstow(VEC a, double tol, int maxiter)
{
    int i, j;
    int n = len(a);
    int iter = 0;
    double err = 1+tol;
    VEC vzero(n);
    VEC b(n);
    VEC c(n);
    VEC d(n);
    VEC root(n-1);
    Complex czero(0,0);
    Complex cinit(0.5,0);
    Complex p(0,0);
    Complex q(0,0);
    Complex R(0,0);
    Complex S(0,0);
    Complex dR_dp(0,0);
    Complex dR_dq(0,0);
    Complex dS_dp(0,0);
    Complex dS_dq(0,0);
    Complex det(0,0);

    j = 0;
    while(n > 0){
        p = cinit;
        q = cinit;
        if(n > 3){
            b = vzero;
            c = vzero;
            d = vzero;
            iter = 0;
            err = 1+tol;
            while(err >= tol and iter < maxiter){
                b[n-3] = a[n-1];
                b[n-4] = a[n-2] - p*b[n-3];
                for(i=n-5; i>=0; i--){
                    b[i] = a[i+2] - p*b[i+1] - q*b[i+2];
                }
                R = a[1] - p*b[0] - q*b[1];
                S = a[0] - q*b[0];
        
                c[n-3] = czero;
                c[n-4] = -1*b[n-3] - p*c[n-3];
                for(i=n-5; i>=0; i--){
                    c[i] = -1*b[i+1] - p*c[i+1] - q*c[i+2];
                }
                dR_dp = -1*b[0] - p*c[0] - q*c[1];
                dS_dp = -1*q*c[0];
        
                d[n-3] = czero;
                d[n-4] = -1*p*d[n-3];
                for(i=n-5; i>=0; i--){
                    d[i] = -1*p*d[i+1] - b[i+2] - q*d[i+2];
                }
                dR_dq = -1*p*d[0] - b[1] - q*d[1];
                dS_dq = -1*b[0] - q*d[0];

                det = dR_dp*dS_dq - dR_dq*dS_dp;

                p = p - (dS_dq*R - dR_dq*S)/det;
                q = q - (-1*dS_dp*R + dR_dp*S)/det;
        
                err = max(cfabs(R), cfabs(S));
                iter++;
            }
            if(p.r()*p.r()-4*q.r() < 0){
                root[j].x = -1*p.r()/2;
                root[j].y = sqrt(4*q.r()-p.r()*p.r())/2;
                j++;
                root[j].x = -1*p.r()/2;
                root[j].y = -1 * sqrt(4*q.r()-p.r()*p.r())/2;
                j++;
            }
            else{
                root[j].x = (-1*p.r() + sqrt(p.r()*p.r()-4*q.r()))/2;
                root[j].y = 0;
                j++;
                root[j].x = (-1*p.r() - sqrt(p.r()*p.r()-4*q.r()))/2;
                root[j].y = 0;
                j++;
            }
            a = b;
            n -= 2;
        }
        else if(n == 3){
            if(b[1].r()*b[1].r()-4*b[2].r()*b[0].r() < 0){
                root[j].x = (-1*b[1].r()) / (2*b[2].r());
                root[j].y = sqrt(4*b[2].r()*b[0].r()-b[1].r()*b[1].r()) / (2*b[2].r());
                j++;
                root[j].x = (-1*b[1].r()) / (2*b[2].r());
                root[j].y = -1 * sqrt(4*b[2].r()*b[0].r()-b[1].r()*b[1].r()) / (2*b[2].r());
                j++;
            }
            else{
                root[j].x = (-1*b[1].r() + sqrt(b[1].r()*b[1].r()-4*b[2].r()*b[0].r())) / (2*b[2].r());
                root[j].y = 0;
                j++;
                root[j].x = (-1*b[1].r() - sqrt(b[1].r()*b[1].r()-4*b[2].r()*b[0].r())) / (2*b[2].r());
                root[j].y = 0;
                j++;
            }
            n -= 3;
        }
        else if(n == 2){
            root[j].x = -1*b[0].r()/b[1].r();
            root[j].y = 0;
            n -= 2;
        }

    }
    return root;
}
