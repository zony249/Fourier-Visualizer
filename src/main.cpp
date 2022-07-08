#include <iostream>
#include <chrono>

#include "../headers/complex.hpp"

// Testing Purposes
#include "../headers/vector_ops.hpp"


using namespace std::chrono;

int main(int argc, char* argv[]) {

    long size = 100000000;
    float* a = (float*)malloc(size * sizeof(float));
    float* b = (float*)malloc(size * sizeof(float));


    for (long i = 0; i < size ; i++) {
        a[i] = (float)(i);
        b[i] = (float)(i);
    }


    // Vector Operation (TIME THIS)
    auto start = high_resolution_clock::now();

    float * c = vector_add(a, b, size);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "vector_add Execution Time: " << duration.count() << std::endl;
    std::cout << "c[100]: " << c[100] << std::endl; 





    
    float * d = (float*)malloc(size * sizeof(float));
    auto start1 = high_resolution_clock::now();

    for (long i = 0; i < size; ++i) {
        d[i] = a[i] + b[i];
    };

    auto stop1 = high_resolution_clock::now();
    auto duration1 = duration_cast<nanoseconds>(stop1 - start1);
    std::cout << "scalar_add Execution Time: " << duration1.count() << std::endl;
    std::cout << "d[100]: " << d[100] << std::endl; 






    // Print out
    for (long i = 0; i < size; ++i) {
        // std::cout << c[i] << " ";
    } std::cout << std::endl;


    // Vector Sum (TIME THIS)
    auto start2 = high_resolution_clock::now();

    float sum = vector_sum(a, size);

    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<nanoseconds>(stop2 - start2);

    std::cout << "vector_sum Execution Time: " << duration2.count() << std::endl;


    std::cout << "Sum of A (Vector Implementation): " << sum << std::endl;






    sum = 0.;
    // Scalar sum (TIME THIS)
    auto start3 = high_resolution_clock::now();

    for (long i = 0; i < size; i++) {
        sum += a[i];
    };

    auto stop3 = high_resolution_clock::now();
    auto duration3 = duration_cast<nanoseconds>(stop3 - start3);
    std::cout << "scalar_sum Execution Time: " << duration3.count() << std::endl;



    std::cout << "Sum of A (Scalar Implementation): " << sum << std::endl << std::endl;








    auto start4 = high_resolution_clock::now();

    float e = vector_dot(a, b, size);

    auto stop4 = high_resolution_clock::now();
    auto duration4 = duration_cast<nanoseconds>(stop4 - start4);
    std::cout << "vector_dot Execution Time: " << duration4.count() << std::endl;
    std::cout << "Dot Product of A and B: " << e << std::endl; 






    auto start5 = high_resolution_clock::now();

    float f = 0.0;
    for (long i = 0; i < size; ++i) {
        f += a[i] * b[i];
    }

    auto stop5 = high_resolution_clock::now();
    auto duration5 = duration_cast<nanoseconds>(stop5 - start5);
    std::cout << "scalar_dot Execution Time: " << duration5.count() << std::endl;
    std::cout << "Dot Product of A and B: " << f << std::endl; 


    return 0;
}