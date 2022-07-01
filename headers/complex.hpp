#ifndef COMPLEX_H
#define COMPLEX_H



template <class T>
class Complex {
    public:
        T re;
        T im;


        Complex(T re, T im);
        Complex operator + (const Complex<T>& rhs);
        Complex operator - (const Complex<T>& rhs);
        Complex operator * (const Complex<T>& rhs);
        Complex operator / (const Complex<T>& rhs);

};




#endif