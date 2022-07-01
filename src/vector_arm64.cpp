#include <arm_neon.h>
#include <stdlib.h>
#include "../headers/vector_ops.hpp"

float* vector_add(float* a, float* b, long size) {
    long remainder = size % 16;

    float* c = (float *) malloc(size * sizeof(float)); 

    for (long i = 0; i + 16 <= size; i+=16) {
        float32x4_t r1 = vld1q_f32(a + i);
        float32x4_t r2 = vld1q_f32(b + i);

        float32x4_t r3 = vld1q_f32(a + i + 4);
        float32x4_t r4 = vld1q_f32(b + i + 4);
        
        float32x4_t r5 = vld1q_f32(a + i + 8);
        float32x4_t r6 = vld1q_f32(b + i + 8);
        
        float32x4_t r7 = vld1q_f32(a + i + 12);
        float32x4_t r8 = vld1q_f32(b + i + 12);


        float32x4_t rd1 = vaddq_f32(r1, r2);
        float32x4_t rd2 = vaddq_f32(r3, r4);
        float32x4_t rd3 = vaddq_f32(r5, r6);
        float32x4_t rd4 = vaddq_f32(r7, r8);

        vst1q_f32(c + i, rd1);
        vst1q_f32(c + i + 4, rd2);
        vst1q_f32(c + i + 8, rd3);
        vst1q_f32(c + i + 12, rd4);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] + b[i];
    }

    return c;
}