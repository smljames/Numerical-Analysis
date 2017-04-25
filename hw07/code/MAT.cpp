/*
EE407002 102061125 Chen Kuan-Chun
*/

// MAT class fuctions
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
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
                if((*va[i])[k]!=0 and m1[k][j]!=0){
                    z[i][j]+=((*va[i])[k]*m1[k][j]);
                }
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

void MAT::print()
{
    MAT s(*this);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<fixed<<setprecision(4)<<s[i][j]<<"\t";
        }
        printf("\n");
    }
}

void GE(MAT A, VEC b)
{
    double y;

    for(int i=0; i<A.n-1; i++){
        for(int j=0; j<A.n-1; j++){
            y = A[j][i]/A[i][i];
            for(int k=i; k<=A.n-1; k++){
                A[j][k] -= y * A[i][k];
            }
            b[j] -= y * b[i];
        }
    }
}

void GE_PP(MAT A, VEC b)
{
    int i, j, k;
    double y;

    for(i=0; i<A.n-1; i++){
        y = fabs(A[i][i]);
        for(k=i, j=i+1; j<=A.n-1; j++){
            if(fabs(A[i][j]) > y){
                y = fabs(A[j][i]);
                k = j;
            }
        }
        if(i != k){
            for(j=i; j<A.n; j++){
                y = A[i][j];
                A[i][j] = A[k][j];
                A[k][j] = y;
            }
            y = b[i];
            b[i] = b[k];
            b[k] = y;
        }
        for(j=i+1; j<=A.n-1; j++){
            y = A[j][i]/A[i][i];
            for(k=i; k<=A.n-1; k++){
                A[j][k] -= y * A[i][k];
                b[j] -= y * b[i];
            }
        }
    }
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
            y[row] = (b[row]-sigma)/A[row][row]; //divided by g?
        }

        z = y-x; //difference between new vector and old vector
 
        // 1-norm
        //err = z.one_norm();
        
        // 2-norm
        //err = z.two_norm();

        // infinity-nrom
        err = z.inf_norm();
        
        
        x = y; //update the vector
        iter++;
    }
    return iter-1;
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
 
        // 1-norm
        //err = z.one_norm();
        
        // 2-norm
        //err = z.two_norm();

        // infinity-nrom
        err = z.inf_norm();
        
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
            y[row] = (b[row]-sigma)/A[row][row];
        }

        z = y-x; //difference between new vector and old vector
 
        // 1-norm
        //err = z.one_norm();
        
        // 2-norm
        //err = z.two_norm();

        // infinity-nrom
        err = z.inf_norm();
        
        
        x = y; //update the vector
        iter++;
    }
    return iter-1;
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
 
        // 1-norm
        //err = z.one_norm();
        
        // 2-norm
        //err = z.two_norm();

        // infinity-nrom
        err = z.inf_norm();
        
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
 
        // 1-norm
        //err = z.one_norm();
        
        // 2-norm
        //err = z.two_norm();

        // infinity-nrom
        err = z.inf_norm();
        
        
        x = y; //update the vector
        iter++;
    }
    return iter-1;
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
 
        // 1-norm
        //err = z.one_norm();
        
        // 2-norm
        //err = z.two_norm();

        // infinity-nrom
        err = z.inf_norm();
        
        x = y; //update the vector
        iter++;
    }
    return iter-1;
}

int cg(MAT &A, VEC b, VEC &x, int maxIter, double tol)
{
    // r: define r=b-Ax, the residue of x
    // p: search path direction
    // Ap: the vector of A*p
    VEC r(b-A*x), p(r), Ap(A.n);
    int iter=0;
    double alpha, beta, err=1.0;

    while(iter <= maxIter and err > tol){
        iter++;
        Ap = A*p;
        alpha = (p*r)/(p*Ap);
        x += alpha*p;
        r -= alpha*(Ap);
        beta = (p*(A*r))/(p*Ap);
        p = r - beta*p;
        err = sqrt(r*r/A.n);
    }

    return iter;
}

double tol_cg(MAT &A, VEC b, VEC &x, VEC x_lu)
{
    // r: define r=b-Ax, the residue of x
    // p: search path direction
    // Ap: the vector of A*p
    VEC r(b-A*x), p(r), Ap(A.n);
    double alpha, beta, err=1.0, tol;

    while(err > 1e-07){
        Ap = A*p;
        alpha = (p*r)/(p*Ap);
        x += alpha*p;
        r -= alpha*(Ap);
        beta = (p*(A*r))/(p*Ap);
        p = r - beta*p;
        err = (x_lu - x).inf_norm();
    }
    tol = sqrt(r*r/A.n);

    return tol;
}

int EVpwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter)
{
    int iter = 0;
    double err = 1.0;
    VEC Aq(A.n); 
    VEC r(A.n);

    while(iter < maxiter and err > tol){
        Aq = A * q0;
        r = Aq - lambda * q0;
        lambda = q0 * Aq;
        q0 = Aq / Aq.two_norm();                
        err = r.two_norm();
        iter++;
    }
    return iter;
}

int EVpwr_all_err(MAT &A, VEC &q0, VEC &lambda_v, int maxiter, int method, VEC &all_err)
{
    int iter = 0;
    double lambda = lambda_v[0];
    double lambda_old = lambda_v[0];
    VEC q0_old(q0);
    VEC Aq(A.n); 
    VEC r(A.n);
    VEC u(A.n);
    VEC w(A.n);

    while(iter < maxiter){
        Aq = A * q0;
        lambda = q0 * Aq;
        q0 = Aq / Aq.two_norm();
        
        switch(method){
            case 1:
                all_err[iter] = fabs(lambda - lambda_old);
                break;

            case 2:
                all_err[iter] = (q0 - q0_old).two_norm();
                break;

            case 3:
                r = Aq - lambda * q0;
                all_err[iter] = r.two_norm();
                break;

            case 4:
                r = Aq - lambda_old * q0_old;
                u = q0_old * A;
                w = u / u.two_norm();
                all_err[iter] = r.two_norm() / fabs(w * q0_old);
        }
        lambda_v[iter] = lambda;
        lambda_old = lambda;        
        q0_old = q0;
        iter++;
    }
    return iter;
}

int EVipwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter)
{
    int iter=0;
    double err = 1.0;
    VEC y(A.n);
    VEC z(A.n);
    VEC r(A.n);
    MAT G(A);
    G = luFact(G);
    while(iter < maxiter and err > tol){
        y = fwdSubs(G, q0);
        z = bckSubs(G, y);
        r = A * q0 - lambda * q0;
        q0 = z/z.two_norm();
        lambda = q0 * (A * q0);
        err = r.two_norm();
        
        iter++;
    }
    return iter;
}      

int EVipwrShft(MAT &A, VEC &q0, double &lambda, double mu, double tol, int maxiter)
{
    int iter=0;
    double err = 1.0;
    VEC y(A.n);
    VEC z(A.n);
    VEC r(A.n);
    MAT G(A);
    MAT I(A.n);

    for(int i=0; i<A.n; i++){
        I[i][i] = 1;
    }

    G = G - mu * I;
    G = luFact(G);

    while(iter < maxiter and err > tol){
        y = fwdSubs(G, q0);
        z = bckSubs(G, y);
        r = A * q0 - lambda * q0;
        q0 = z/z.two_norm();
        lambda = q0 * (A * q0);
        err = r.two_norm();

        iter++;
    }
    return iter;
}   

void qrFact(MAT &A, MAT &Q, MAT &R)
{
    int i, j;
    MAT At(A.tpose());

    A = A.tpose(); // col <-> row
    R[0][0] = A[0].two_norm();
    Q[0] = A[0] / R[0][0];
    for(j=1; j<A.n; j++){
        Q[j] = A[j];
        for(i=0; i<j; i++){
            R[i][j] = Q[i] * Q[j];
            Q[j] -= R[i][j] * Q[i];
        }
        R[j][j] = Q[j].two_norm();
        Q[j] /= R[j][j];
    }
    Q = Q.tpose(); // col <-> row
}

int EVqr(MAT &A, double tol, int maxiter)
{
    MAT Q(A.n);
    MAT R(A.n);
    int i, j;
    int iter = 0;
    double err = 1.0;

    while(err > tol and iter < maxiter){
        // -----QR decomposition and QR iteration-----
        qrFact(A, Q, R);
        A = R * Q;

        err = fabs(A[1][0]);
        for(i=2; i<A.n; i++){
            if(A[i][i-1] > err) err = fabs(A[i][i-1]);
        }

        iter++;
    }
    return iter;
}

int EVqrShifted(MAT &A, double mu, double tol, int maxiter)
{
    MAT R(A.n);
    MAT Q(A.n);
    MAT I(A.n);
    int i, j;
    int iter = 0;
    double err = 1.0;
    
    for(i=0; i<A.n; i++) I[i][i] = 1;

    while(err > tol and iter < maxiter){
        // -----QR decomposition and QR iteration-----
        A = A - mu * I;
        qrFact(A, Q, R);
        A = R * Q + mu * I; // A = R * Q

        err = fabs(A[1][0]);
        for(i=2; i<A.n; i++){
            if(A[i][i-1] > err) err = fabs(A[i][i-1]);
        }

        iter++;
    }
    return iter;
}

