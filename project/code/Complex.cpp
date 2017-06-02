// Function definitions
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Complex.h"
using namespace std;

Complex::Complex(double r, double i)
{
    x=r; y=i;
}
Complex::Complex(const Complex &z)
{
    x=z.x; y=z.y;
}
double Complex::r() const
{
    return x;
}
double Complex::i() const
{
    return y;
}
Complex Complex::operator+=(const Complex z)
{
    x+=z.x;
    y+=z.y;
    return *this;
}
Complex Complex::operator+=(double db1)
{
    x+=db1;
    return *this;
}
Complex Complex::operator-=(const Complex z)
{
    x-=z.x;
    y-=z.y;
    return *this;
}
Complex Complex::operator-=(double db1)
{
    x-=db1;
    return *this;
}
Complex Complex::operator*=(const Complex z)
{
    double x_old = x;
    double y_old = y;
    x=x_old*(z.x)-y_old*(z.y);
    y=x_old*(z.y)+y_old*(z.x);
    return *this;
}
Complex Complex::operator*=(double db)
{
    x*=db;
    y*=db;
    return *this;
}
Complex Complex::operator/=(const Complex z)
{
    double x_old = x;
    double y_old = y;
    double tmp=z.x*z.x+z.y*z.y;
    x=(x_old*z.x+y_old*z.y)/tmp;
    y=(y_old*z.x-x_old*z.y)/tmp;
    return *this;
}
Complex Complex::operator/=(double db)
{
    x/=db;
    y/=db;
    return *this;
}
Complex Complex::operator+(const Complex z1)
{
    Complex z(*this);
    return z+=z1;
}
Complex Complex::operator+(double db)
{
    Complex z(*this);
    return z+=db;
}
Complex operator+(double db, Complex z1)
{
    Complex z(z1);
    return z+=db;
}
Complex Complex::operator-(const Complex z1)
{
    Complex z(*this);
    return z-=z1;
}
Complex Complex::operator-(double db)
{
    Complex z(*this);
    return z-=db;
}
Complex operator-(double db, Complex z1)
{
    Complex z(z1);
    return z-=db;
}
Complex Complex::operator*(const Complex z1)
{
    Complex z(*this);
    return z*=z1;
}
Complex Complex::operator*(double db)
{
    Complex z(*this);
    return z*=db;
}
Complex operator*(double db, Complex z1)
{
    Complex z(z1);
    return z*=db;
}
Complex Complex::operator/(const Complex z1)
{
    Complex z(*this);
    return z/=z1;
}
Complex Complex::operator/(double db)
{
    Complex z(*this);
    return z/=db;
}
Complex operator/(double db, Complex z1)
{
    Complex z(z1);
    return z/=db;
}
int Complex::operator==(Complex z1)
{
    Complex z(*this);
    if (z.r()==z1.r() && z.i()==z1.i()) return 1;
    else return 0;
}
int Complex::operator!=(Complex z1)
{
    Complex z(*this);
    if (z.r()!=z1.r() || z.i()!=z1.i()) return 1;
    else return 0;
}
double cfabs(Complex z)
{
    return sqrt(z.r()*z.r()+z.i()*z.i());
}
void print(Complex z)
{
   if(z.y>0) cout<<z.x<<"+"<<z.y<<"i"<<endl;
   else if(z.y<0) cout<<z.x<<z.y<<"i"<<endl;
   else cout<<z.x<<endl;
}
Complex rtoc(double i)
{
    Complex z(i,0);
    return z;
}
