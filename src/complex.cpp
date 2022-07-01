#include <iostream>

#include "../headers/complex.hpp"
#include "../headers/vector_ops.hpp"

/**
 * @brief Main constructor of the Complex class
 * 
 * 
 */
template <class T>
Complex<T>::Complex(T re, T im) {
    this->re = re;
    this->im = im;
}