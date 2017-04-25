// MAT class fuctions
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"

using namespace std;

MAT::MAT(int dim) // uninit constructor
{
    n=dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++) {
        va[i]=newVEC(n);
    }
}

MAT::MAT(const MAT &m1) // copy constructor
{
    VEC **vsrc=m1.va; // to get around not indexing const MAT
    n=m1.n;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++) {
        va[i]=newVEC(n);
        (*va[i])=(*vsrc[i]); // VEC assignment
    }
}

MAT::MAT(int dim,double *v) // init constructor
{
    n=dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++) {
        va[i]=newVEC(n);
        for (int j=0; j<n; j++) {
            (*va[i])[j]=*(v++); // array indexing + VEC indexing
        }
    }
}

MAT::~MAT() // destructor
{
    for (int i=n-1; i>=0; i--)
        (*va[i]).~VEC();
        free(va);
}

int MAT::dim() // return dimension of the matrix
{
    return n;
}

MAT MAT::tpose() // matrix transpose
{
    MAT mnew(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            mnew[i][j]=(*va[j])[i];
        }
    }
    return mnew;
}

MAT &MAT::operator-() // unary operator - : negative value
{
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            (*va[i])[j]=-(*va[i])[j];
    return *this;
}

MAT &MAT::operator=(MAT m1) // assignment
{
    for (int i=0; i<n; i++)
        (*va[i])=m1[i]; // VEC assignment
    return *this;
}

MAT &MAT::operator+=(MAT &m1) // m += m1
{
    for (int i=0; i<n; i++)
        (*va[i])+=m1[i]; // VEC += operation
    return *this;
}

MAT &MAT::operator-=(MAT &m1) // m -= m1
{
    for (int i=0; i<n; i++)
        (*va[i])-=m1[i]; // VEC -= operation
    return *this;
}

MAT &MAT::operator*=(double a) // m *= dbl
{
    for (int i=0; i<n; i++)
        (*va[i])*=a; 
    return *this;
}

MAT &MAT::operator/=(double a) // m /= dbl
{
    for (int i=0; i<n; i++)
        (*va[i])/=a;
    return *this;
}

MAT MAT::operator+(MAT m1) // addition
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])+m1[i]; // VEC addition and assignment
    return s;
}

MAT MAT::operator-(MAT m1) // subtraction
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])-m1[i]; // VEC subtraction and assignment
    return s;
}

MAT MAT::operator*(MAT m1) // matrix-matrix product
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=0;
            for (int k=0; k<n; k++)
                z[i][j]+=((*va[i])[k]*m1[k][j]);
        }
    }
    return z;
}

VEC &MAT::operator[](int m) // m'th row
{
    if (m<0) m=0;
    else if (m>=n) m=n-1;
    return *va[m];
}

VEC MAT::operator*(VEC v1) // M * v
{
    VEC s(n);
    for (int i=0; i<n; i++) {
        s[i]=(*va[i])*v1; // VEC inner product
    }
    return s;
}

MAT MAT::operator*(double a) // m * dbl
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])*a;
    return s;
}

MAT MAT::operator/(double a) // m / dbl
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])/a;
    return s;
}

MAT operator*(double a, MAT &m1)
{
    MAT s(m1.n);
    for (int i=0; i<m1.n; i++)
        s[i]=m1[i]*a;
    return s;
}

VEC operator*(VEC &v1, MAT &m1) // vT x M
{
    VEC v2(m1.n);
    for (int i=0; i<m1.n; i++) {
        v2[i]=0;
        for (int j=0; j<m1.n; j++) {
            v2[i] += v1[j]*m1[j][i];
        }
    }
    return v2;
}

MAT &luFact(MAT &m1)
{
    int i, j, k;

    for(i=0; i<m1.n; i++) {
        for (j=i+1; j<m1.n; j++) {
            m1[j][i] /= m1[i][i];
        }
        for(j=i+1; j<m1.n; j++) {
            for(k=i+1; k<m1.n; k++) {
                m1[j][k] -= m1[j][i]*m1[i][k];
            }
        }
    }
    return m1;
}

VEC fwdSubs(MAT &m1, VEC b)
{
    int i, j;
    VEC y(b);

    for (i=0; i<m1.n-1; i++) {
        for (j=i+1; j<m1.n; j++) y[j] -= m1[j][i]*y[i];
    }
    return y;
}

VEC bckSubs(MAT &m1, VEC b)
{
    int i, j, k;
    VEC x(b);
    for (i=m1.n-1; i>=0; i--) {
        x[i] /= m1[i][i];
        for (j=i-1; j>=0; j--) 
            x[j] -= m1[j][i]*x[i];
    }
    return x;
}

int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    diff: the difference between new x and old x
    */
    VEC y(x), z(A.n);
    int iter=0, row, col;
    double sigma, diff=1.0;
    
    while(iter <= maxIter and diff > tol){ 
        iter++;
        
        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*x[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row]; //divided by g?
        }

        z = y-x; //difference between new vector and old vector
        diff = z.one_norm();
        x = y; //update the vector
    }
    return iter;
}

double tol_jacobi(MAT &A,VEC b,VEC &x,VEC x_lu)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    tol: the difference between new x and old x
    */
    VEC y(x), z(A.n), err_vec(A.n);
    int iter=1, row, col;
    double sigma, err, tol;
    
    err_vec = x_lu-x;
    err = err_vec.inf_norm();

    while(err >= 1e-7){
        x = y;

        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*x[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row]; //divided by g?
        }
        err_vec = x_lu-y;
        err = err_vec.inf_norm();
    }
    
    z = y-x;
    tol = z.one_norm();
     
    return tol;
}

int jacobi_over_relaxation(MAT &A,VEC b,VEC &x, int maxIter, double w, double tol)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    err: we use a p-norm to calculate the error
    */
    VEC y(x), z(A.n);
    int iter=1, row, col;
    double sigma, err=1.0;
    
    while(iter <= maxIter and err > tol){ 
        err = 0; //initialize the error
        
        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*x[col];
                }
            }
            y[row] = w*(b[row]-sigma)/A[row][row] + (1-w)*x[row];
        }

        z = y-x; //difference between new vector and old vector
 
        /*1-norm*/
        /*        
        for(row=0; row<A.n; row++){
            err += fabs(z[row]);
        }
        */
        
        /*2-norm*/
        
        for(row=0; row<A.n; row++){
            err += (z[row]*z[row]);
        }
        err = sqrt(err);
        

        /*infinity-norm*/
        /*        
        err = abs(z[0]);
        for(row=1; row<A.n; row++){
            if(abs(z[row]) > err) err = abs(z[row]);
        }
        */
        
        x = y; //update the vector
        iter++;
    }
    return iter-1;
}

int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    diff: the difference between new x and old x
    */
    VEC y(x), z(A.n);
    int iter=0, row, col;
    double sigma, diff=1.0;
    
    while(iter <= maxIter and diff > tol){ 
        iter++;
        
        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row];
        }

        z = y-x; //difference between new vector and old vector
        diff = z.one_norm(); 
        x = y; //update the vector
    }
    return iter;
}

double tol_gaussSeidel(MAT &A,VEC b,VEC &x, VEC x_lu)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    err: we use a p-norm to calculate the error
    */
    VEC y(x), z(A.n), err_vec(A.n);
    int iter=1, row, col;
    double sigma, err, tol;
    
    err_vec = x_lu-x;
    err = err_vec.inf_norm();

    while(err >= 1e-7){
        x = y;
        
        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row];
        }

        err_vec = x_lu-y;
        err = err_vec.inf_norm();
    }
    
    z = y-x;
    tol = z.one_norm();

    return tol;
}

int successive_over_relaxation(MAT &A,VEC b,VEC &x,int maxIter, double w, double tol)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    err: we use a p-norm to calculate the error
    */
    VEC y(x), z(A.n);
    int iter=1, row, col;
    double sigma, err=1.0;
    
    while(iter <= maxIter and err > tol){ 
        err = 0; //initialize the error
        
        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = w*(b[row]-sigma)/A[row][row] + (1-w)*x[row];
        }

        z = y-x; //difference between new vector and old vector
 
        /*1-norm*/
        /*        
        for(row=0; row<A.n; row++){
            err += fabs(z[row]);
        }
        */
        
        /*2-norm*/
        
        for(row=0; row<A.n; row++){
            err += (z[row]*z[row]);
        }
        err = sqrt(err);
        

        /*infinity-norm*/
        /*        
        err = abs(z[0]);
        for(row=1; row<A.n; row++){
            if(abs(z[row]) > err) err = abs(z[row]);
        }
        */
        
        x = y; //update the vector
        iter++;
    }
    return iter-1;
}

int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    diff: the difference between new x and old x
    */
    VEC y(x), z(A.n);
    int iter=0, row, col;
    double sigma, diff=1.0;
    
    while(iter <= maxIter and diff > tol){ 
        iter++;

        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row];
        }

        for(row=A.n-1; row>=0; row--){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row];
        }

        z = y-x; //difference between new vector and old vector
        diff = z.one_norm(); 
        x = y; //update the vector
    }
    return iter;
}

double tol_sgs(MAT &A,VEC b,VEC &x, VEC x_lu)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    err: we use a p-norm to calculate the error
    */
    VEC y(x), z(A.n), err_vec(A.n);
    int iter=1, row, col;
    double sigma, err, tol;
    
    err_vec = x_lu-x;
    err = err_vec.inf_norm();

    while(err >= 1e-7){
        x = y;

        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row];
        }

        for(row=A.n-1; row>=0; row--){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = (b[row]-sigma)/A[row][row];
        }

        err_vec = x_lu-y;
        err = err_vec.inf_norm();
    }
    
    z = y-x;
    tol = z.one_norm();
    
    return tol;
}

int symmetric_successive_over_relaxation(MAT &A,VEC b,VEC &x,int maxIter, double w, double tol)
{
    /*
    y: record the updated vector in new iteration
    
    z: a vector to record the difference 
       between new vector and old vector
    
    iter: iteration in while-loop
    
    err: we use a p-norm to calculate the error
    */
    VEC y(x), z(A.n);
    int iter=1, row, col;
    double sigma, err=1.0;
    
    while(iter <= maxIter and err > tol){ 
        err = 0; //initialize the error
        
        for(row=0; row<A.n; row++){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = w*(b[row]-sigma)/A[row][row] + (1-w)*x[row];
        }

        for(row=A.n-1; row>=0; row--){
            sigma = 0;
            for(col=0; col<A.n; col++){
                if(row != col){
                    sigma += A[row][col]*y[col];
                }
            }
            y[row] = w*(b[row]-sigma)/A[row][row] + (1-w)*y[row];
        }

        z = y-x; //difference between new vector and old vector
 
        /*1-norm*/
        /*        
        for(row=0; row<A.n; row++){
            err += fabs(z[row]);
        }
        */
        
        /*2-norm*/
        
        for(row=0; row<A.n; row++){
            err += (z[row]*z[row]);
        }
        err = sqrt(err);
        

        /*infinity-norm*/
        /*        
        err = abs(z[0]);
        for(row=1; row<A.n; row++){
            if(abs(z[row]) > err) err = abs(z[row]);
        }
        */
        
        x = y; //update the vector
        iter++;
    }
    return iter-1;
}
