/*
EE407002 102061125 Chen Kuan-Chun
*/

// vector class
#ifndef VEC_H
#define VEC_H
class VEC {
    private:
    int dim; // vector length
    double *val; // array to store vector
    public:
    VEC(int n); // uninit constructor, val set to 0
    VEC(const VEC &v1); // copy constructor
    VEC(int n,double *v); // init constructor
    ~VEC(); // destructor
    int len(); // dimension of the vector
    VEC &operator-(); // unary operator, negative value
    VEC &operator=(const VEC v1); // assignment
    VEC &operator+=(const VEC v1); // V += v1;
    VEC &operator-=(const VEC v1); // V -= v1;
    VEC &operator*=(double a); // V *= dbl;
    VEC &operator/=(double a); // V /= dbl;
    VEC operator+(const VEC v1); // V + v1
    VEC operator-(const VEC v1); // V - v1
    double operator*(VEC v1); // inner product
    VEC operator*(double a); // V * dbl
    VEC operator/(double a); // V / dbl
    double &operator[](int n); // indexing
    double one_norm();
    double two_norm();
    double inf_norm();
    friend VEC operator*(double a,const VEC v1); // dbl x V
    friend VEC *newVEC(int n); // create dynamic VEC
};
VEC operator*(double a,const VEC v1);
VEC *newVEC(int n); // create dynamic VEC
#endif
