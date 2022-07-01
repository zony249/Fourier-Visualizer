#include <iostream>


#include "../headers/complex.hpp"

// Testing Purposes
#include "../headers/vector_ops.hpp"

int main(int argc, char* argv[]) {

    int size = 18;
    float a[size];
    float b[size];


    for (int i = 0; i < size; i++) {
        a[i] = (float)i;
        b[i] = (float)i;
    }

    float * c = vector_add(a, b, size);

    for (int i = 0; i < size; ++i) {
        std::cout << c[i] << " ";
    } std::cout << std::endl;
    
    return 0;
}