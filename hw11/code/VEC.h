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
    friend double one_norm(VEC v1);
    friend double two_norm(VEC v1);
    friend double inf_norm(VEC v1);
    friend void print(VEC v1);
    friend VEC sort(VEC v1);
    friend VEC operator*(double a,const VEC v1); // dbl x V
    friend VEC *newVEC(int n); // create dynamic VEC
    friend double max(VEC v1);
    friend int ceil_ind(double x, VEC &X); // find the index of ceil component
    friend VEC inverse(VEC &v1);
    friend int len(VEC &v1);
};
VEC operator*(double a,const VEC v1);
VEC *newVEC(int n); // create dynamic VEC
double max(VEC v1);
double one_norm(VEC v1);
double two_norm(VEC v1);
double inf_norm(VEC v1);
void print(VEC v1);
VEC sort(VEC v1);
int ceil_ind(double x, VEC &X);
VEC inverse(VEC &v1);
int len(VEC &v1);
double integral(double(*f)(double), double min, double max, int N, int order);
#endif