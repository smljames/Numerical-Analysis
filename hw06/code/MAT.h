/*
EE407002 102061125 Chen Kuan-Chun
*/

// matrix class
#ifndef MAT_H
#define MAT_H
#include "VEC.h"
class MAT {
    private:
    int n; // define nxn matrix
    VEC **va; // array of n pointers to vectors
    public:
    MAT(int dim); // uninit constructor
    MAT(const MAT &m1); // copy constructor
    MAT(int dim,double *v); // init constructor
    ~MAT(); // destructor
    int dim(); // return dimension of the matrix
    MAT tpose(); // transpose
    MAT &operator-(); // unary operator, negative value
    MAT &operator=(MAT m1); // assignment
    MAT &operator+=(MAT &m1); // m += m1;
    MAT &operator-=(MAT &m1); // m -= m1;
    MAT &operator*=(double a); // m *= dbl;
    MAT &operator/=(double a); // m /= dbl;
    MAT operator+(MAT m1); // m1 + m2
    MAT operator-(MAT m1); // m1 - m2
    MAT operator*(MAT m1); // m1 * m2
    VEC &operator[](int m); // m'th row
    VEC operator*(VEC v1); // m x v1
    MAT operator*(double a); // m * dbl
    MAT operator/(double a); // m / dbl
    friend MAT operator*(double a,MAT &m1); // dbl x m
    friend VEC operator*(VEC &v1,MAT &m1); // vT x m
    friend void GE(MAT A, VEC b);
    friend void GE_PP(MAT A, VEC b); 
    friend MAT &luFact(MAT &m1); // LU decomposition
    friend VEC fwdSubs(MAT &m1, VEC b); // forward substitution
    friend VEC bckSubs(MAT &m1, VEC b); // backward substitution
    friend int jacobi(MAT &A, VEC b, VEC &x, int maxIter, double tol);
    friend double tol_jacobi(MAT &A, VEC b, VEC &x, VEC x_lu);
    friend int jacobi_over_relaxation(MAT &A, VEC b, VEC &x, int maxIter, double w, double tol);
    friend int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol);
    friend double tol_gaussSeidel(MAT &A,VEC b,VEC &x, VEC x_lu);
    friend int successive_over_relaxation(MAT &A, VEC b, VEC &x, int maxIter, double w, double tol);
    friend int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol); 
    friend double tol_sgs(MAT &A,VEC b,VEC &x, VEC x_lu);
    friend int symmetric_successive_over_relaxation(MAT &A, VEC b, VEC &x, int maxIter, double w, double tol);
    friend int cg(MAT &A, VEC b, VEC &x, int maxIter, double tol);
    friend double tol_cg(MAT &A, VEC b, VEC &x, VEC x_lu);
    friend int EVpwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);
    friend int EVpwr_all_err(MAT &A, VEC &q0, VEC &lambda, int maxiter, int method, VEC &all_err);
    friend int EVipwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);
    friend int EVipwrShft(MAT &A, VEC &q0, double &lambda, double mu, double tol, int maxiter);
};
MAT operator*(double a,MAT &m1); // dbl x m
VEC operator*(VEC &v1,MAT &m1); // vT x m
void GE(MAT A, VEC b);
void GE_PP(MAT A, VEC b); 
MAT &luFact(MAT &m1); // LU decomposition
VEC fwdSubs(MAT &m1, VEC b); // forward substitution
VEC bckSubs(MAT &m1, VEC b); // backward substitution
int jacobi(MAT &A, VEC b, VEC &x, int maxIter, double tol);
double tol_jacobi(MAT &A, VEC b, VEC &x, VEC x_lu);
int jacobi_over_relaxation(MAT &A, VEC b, VEC &x, int maxIter, double w, double tol);
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol);
double tol_gaussSeidel(MAT &A,VEC b,VEC &x, VEC x_lu);
int successive_over_relaxation(MAT &A, VEC b, VEC &x, int maxIter, double w, double tol);
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol);
double tol_sgs(MAT &A,VEC b,VEC &x, VEC x_lu);
int symmetric_successive_over_relaxation(MAT &A, VEC b, VEC &x, int maxIter, double w, double tol);
int cg(MAT &A, VEC b, VEC &x, int maxIter, double tol);
double tol_cg(MAT &A, VEC b, VEC &x, VEC x_lu);
int EVpwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);
int EVpwr_all_err(MAT &A, VEC &q0, VEC &lambda, int maxiter, int method, VEC &all_err);
int EVipwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);
int EVipwrShft(MAT &A, VEC &q0, double &lambda, double mu, double tol, int maxiter);
#endif
