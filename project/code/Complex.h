// complex number class
#ifndef COMPLEX_H
#define COMPLEX_H
#include <stdio.h>
#include <stdlib.h>

class Complex {
    public:
        double x, y;
        Complex(double r, double i);        // constructor
        Complex(const Complex &C);          // copy constructor
        double r() const;                   // get real part
        double i() const;                   // get imaginary part
        Complex operator+(const Complex z1);      // Z + z1
        Complex operator+(double db1);      // Z + db1
        Complex operator-(const Complex z1);      // Z - z1
        Complex operator-(double db1);      // Z - db1
        Complex operator*(const Complex z1);      // Z * z1
        Complex operator*(double db1);      // Z * db1
        Complex operator/(const Complex z1);      // Z / z1
        Complex operator/(double db1);      // Z / db1
        int operator==(Complex z1);         // testing for equality
        int operator!=(Complex z1);         // testing for inequality
        Complex operator+=(const Complex z1);  // Z += z1
        Complex operator+=(double db1);   // Z += db1
        Complex operator-=(const Complex z1);  // Z -= z1
        Complex operator-=(double db1);   // Z -= db1
        Complex operator*=(const Complex z1);  // Z *= z1
        Complex operator*=(double db1);   // Z *= db1
        Complex operator/=(const Complex z1);  // Z /= z1
        Complex operator/=(double db1);   // Z /= db1
        friend Complex operator+(double db1, Complex z1); // db1 + z1
        friend Complex operator-(double db1, Complex z1); // db1 - z1
        friend Complex operator*(double db1, Complex z1); // db1 * z1
        friend Complex operator/(double db1, Complex z1); // db1 / z1
        friend double cfabs(Complex z1);
        friend void print(Complex z1);
};
double cfabs(Complex z1);                   // complex absolute value |z1|
void print(Complex z1);                     // print out z1
Complex rtoc(double i);                     // real to complex
#endif
