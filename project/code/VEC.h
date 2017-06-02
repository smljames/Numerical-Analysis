/*
EE407002 102061125 Chen Kuan-Chun
*/

// vector class
#ifndef VEC_H
#define VEC_H
#include "Complex.h"
class VEC {
    private:
    int dim; // vector length
    Complex *val; // array to store vector
    public:
    VEC(int n); // uninit constructor, val set to 0
    VEC(const VEC &v1); // copy constructor
    ~VEC(); // destructor
    int len(); // dimension of the vector
    VEC &operator=(const VEC v1); // assignment
    VEC &operator+=(const VEC v1); // V += v1;
    VEC &operator-=(const VEC v1); // V -= v1;
    VEC &operator*=(double a); // V *= dbl;
    VEC &operator/=(double a); // V /= dbl;
    VEC operator+(const VEC v1); // V + v1
    VEC operator-(const VEC v1); // V - v1
    Complex operator*(VEC v1); // inner product
    VEC operator*(double a); // V * dbl
    VEC operator/(double a); // V / dbl
    Complex &operator[](int n); // indexing
    friend void print(VEC v1);
    friend VEC operator*(double a,const VEC v1); // dbl x V
    friend VEC *newVEC(int n); // create dynamic VEC
    friend int len(VEC &v1); 
};
VEC operator*(double a,const VEC v1);
VEC *newVEC(int n); // create dynamic VEC
void print(VEC v1);
int len(VEC &v1); // return the dimension of v1
Complex polyval(VEC &p, Complex x); // calculate f(x)
VEC NewtonHorner(VEC p, Complex x, double tol, int maxiter);
VEC LinsQuadratic(VEC p, double tol, int maxiter);
VEC Bairstow(VEC p, double tol, int maxiter);
#endif
