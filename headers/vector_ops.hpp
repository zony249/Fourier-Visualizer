#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

float* vector_add(float* a, float* b, long size);
float* vector_sub(float* a, float* b, long size);
float* vector_mul(float* a, float* b, long size);
float* vector_div(float* a, float* b, long size);

float vector_sum(float* a, long size);
float vector_dot(float* a, float* b, long size);


double* vector_add(double* a, double* b, long size);
double* vector_sub(double* a, double* b, long size);
double* vector_mul(double* a, double* b, long size);
double* vector_div(double* a, double* b, long size);

double vector_sum(double* a, long size);
double vector_dot(double* a, double* b, long size);

#endif