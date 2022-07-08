#include <arm_neon.h>
#include <stdlib.h>
#include "../headers/vector_ops.hpp"



//////////////////////////////
//                          //
//     FLOAT OPERATIONS     //
//                          //
//////////////////////////////


float* vector_add(float* a, float* b, long size) {
    long remainder = size % 16;

    float* c = (float *) malloc(size * sizeof(float)); 

    for (long i = 0; i + 16 <= size; i+=16) {
        float32x4x4_t r1 = vld4q_f32(a + i);
        float32x4x4_t r2 = vld4q_f32(b + i);
        float32x4x4_t rd;

        rd.val[0] = vaddq_f32(r1.val[0], r2.val[0]);
        rd.val[1] = vaddq_f32(r1.val[1], r2.val[1]);
        rd.val[2] = vaddq_f32(r1.val[2], r2.val[2]);
        rd.val[3] = vaddq_f32(r1.val[3], r2.val[3]);

        vst4q_f32(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] + b[i];
    }

    return c;
}




float* vector_sub(float* a, float* b, long size) {
    long remainder = size % 16;

    float* c = (float *) malloc(size * sizeof(float)); 

    for (long i = 0; i + 16 <= size; i+=16) {
        float32x4x4_t r1 = vld4q_f32(a + i);
        float32x4x4_t r2 = vld4q_f32(b + i);
        float32x4x4_t rd;

        rd.val[0] = vsubq_f32(r1.val[0], r2.val[0]);
        rd.val[1] = vsubq_f32(r1.val[1], r2.val[1]);
        rd.val[2] = vsubq_f32(r1.val[2], r2.val[2]);
        rd.val[3] = vsubq_f32(r1.val[3], r2.val[3]);

        vst4q_f32(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] - b[i];
    }

    return c;
}


float* vector_mul(float* a, float* b, long size) {
    long remainder = size % 16;

    float* c = (float *) malloc(size * sizeof(float)); 

    for (long i = 0; i + 16 <= size; i+=16) {
        float32x4x4_t r1 = vld4q_f32(a + i);
        float32x4x4_t r2 = vld4q_f32(b + i);
        float32x4x4_t rd;

        rd.val[0] = vmulq_f32(r1.val[0], r2.val[0]);
        rd.val[1] = vmulq_f32(r1.val[1], r2.val[1]);
        rd.val[2] = vmulq_f32(r1.val[2], r2.val[2]);
        rd.val[3] = vmulq_f32(r1.val[3], r2.val[3]);

        vst4q_f32(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] * b[i];
    }

    return c;
}


float* vector_div(float* a, float* b, long size) {
    long remainder = size % 16;

    float* c = (float *) malloc(size * sizeof(float)); 

    for (long i = 0; i + 16 <= size; i+=16) {
        float32x4x4_t r1 = vld4q_f32(a + i);
        float32x4x4_t r2 = vld4q_f32(b + i);
        float32x4x4_t rd;

        rd.val[0] = vdivq_f32(r1.val[0], r2.val[0]);
        rd.val[1] = vdivq_f32(r1.val[1], r2.val[1]);
        rd.val[2] = vdivq_f32(r1.val[2], r2.val[2]);
        rd.val[3] = vdivq_f32(r1.val[3], r2.val[3]);

        vst4q_f32(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] / b[i];
    }

    return c;
}

float vector_sum(float* a, long size) {
    if (size == 1) {
        return a[0];
    } else {
        long remainder = size % 16;
        float* b = (float*) malloc((size / 16 + 1) * sizeof(float));
        b[size/16] = 0.0;

        long b_idx = 0;
        for (long i = 0; i + 16 <= size; i+=16) {
            float32x4_t r1 = vld1q_f32(a + i);
            float32x4_t r2 = vld1q_f32(a + i + 4);
            float32x4_t r3 = vld1q_f32(a + i + 8);
            float32x4_t r4 = vld1q_f32(a + i + 12);

            float32x4_t d1 = {vaddvq_f32(r1), vaddvq_f32(r2), vaddvq_f32(r3), vaddvq_f32(r4)};

            b[b_idx] = vaddvq_f32(d1);


            b_idx ++;
        }

        for (long i = size - remainder; i < size; ++i) {
            b[size/16] += a[i];
        }

        return vector_sum(b, size / 16 + 1);

    }
}


float vector_dot(float* a, float* b, long size) {

    long remainder = size % 16;
    float32x4_t acc = vdupq_n_f32(0);
    for (long i = 0; i + 16 <= size; i+=16){
        float32x4x4_t r1 = vld4q_f32(a + i);
        float32x4x4_t r2 = vld4q_f32(b + i);

        for (long j = 0; j < 4; j++) {
            acc = vmlaq_f32(acc, r1.val[j], r2.val[j]);
        }
    }

    // Deal with remainders
    float sumr = 0.0;
    for (long i = size - remainder; i < size; ++i) {
        sumr += a[i] * b[i];
    }

    // Deal with accumulation
    float suma = 0.0;
    for (long i = 0; i < 4; ++i) {
        suma += acc[i];
    }

    return suma + sumr;

}






























/////////////////////////////
//                         //
//    DOUBLE OPERATIONS    //
//                         //
/////////////////////////////

double* vector_add(double* a, double* b, long size) {
    long remainder = size % 8;

    double* c = (double *) malloc(size * sizeof(double)); 

    for (long i = 0; i + 8 <= size; i+=8) {
        float64x2x4_t r1 = vld4q_f64(a + i);
        float64x2x4_t r2 = vld4q_f64(b + i);
        float64x2x4_t rd;
        
        for (long j = 0; j < 4; ++j) {
            rd.val[j] = vaddq_f64(r1.val[j], r2.val[j]);
        }
        vst4q_f64(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] + b[i];
    }

    return c;
}


double* vector_sub(double* a, double* b, long size) {
    long remainder = size % 8;

    double* c = (double *) malloc(size * sizeof(double)); 

    for (long i = 0; i + 8 <= size; i+=8) {
        float64x2x4_t r1 = vld4q_f64(a + i);
        float64x2x4_t r2 = vld4q_f64(b + i);
        float64x2x4_t rd;
        
        for (long j = 0; j < 4; ++j) {
            rd.val[j] = vsubq_f64(r1.val[j], r2.val[j]);
        }
        vst4q_f64(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] - b[i];
    }

    return c;
}


double* vector_mul(double* a, double* b, long size) {
    long remainder = size % 8;

    double* c = (double *) malloc(size * sizeof(double)); 

    for (long i = 0; i + 8 <= size; i+=8) {
        float64x2x4_t r1 = vld4q_f64(a + i);
        float64x2x4_t r2 = vld4q_f64(b + i);
        float64x2x4_t rd;
        
        for (long j = 0; j < 4; ++j) {
            rd.val[j] = vmulq_f64(r1.val[j], r2.val[j]);
        }
        vst4q_f64(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] * b[i];
    }

    return c;
}


double* vector_div(double* a, double* b, long size) {
    long remainder = size % 8;

    double* c = (double *) malloc(size * sizeof(double)); 

    for (long i = 0; i + 8 <= size; i+=8) {
        float64x2x4_t r1 = vld4q_f64(a + i);
        float64x2x4_t r2 = vld4q_f64(b + i);
        float64x2x4_t rd;
        
        for (long j = 0; j < 4; ++j) {
            rd.val[j] = vdivq_f64(r1.val[j], r2.val[j]);
        }
        vst4q_f64(c + i, rd);
    }

    for (long i = size - remainder; i < size; ++i) {
        c[i] = a[i] / b[i];
    }

    return c;
}



double vector_sum(double* a, long size) {
    if (size == 1) {
        return a[0];
    } else {
        long remainder = size % 8;
        double* b = (double*) malloc((size / 8 + 1) * sizeof(double));
        b[size/8] = 0.0;

        long b_idx = 0;
        for (long i = 0; i + 8 <= size; i+=8) {
            float64x2_t r1 = vld1q_f64(a + i);
            float64x2_t r2 = vld1q_f64(a + i + 2);
            float64x2_t r3 = vld1q_f64(a + i + 4);
            float64x2_t r4 = vld1q_f64(a + i + 6);

            float64x2_t d1 = {vaddvq_f64(r1), vaddvq_f64(r2)};
            float64x2_t d2 = {vaddvq_f64(r3), vaddvq_f64(r4)};


            b[b_idx] = vaddvq_f64(d1);
            b[b_idx] += vaddvq_f64(d2);

            b_idx ++;
        }

        for (long i = size - remainder; i < size; ++i) {
            b[size/8] += a[i];
        }

        return vector_sum(b, size / 8 + 1);

    }
}


double vector_dot(double* a, double* b, long size) {
    long remainder = size % 8;
    float64x2_t acc = vdupq_n_f64(0);
    for (long i = 0; i + 8 <= size; i+=8){
        float64x2x4_t r1 = vld4q_f64(a + i);
        float64x2x4_t r2 = vld4q_f64(b + i);

        for (long j = 0; j < 4; j++) {
            acc = vmlaq_f64(acc, r1.val[j], r2.val[j]);
        }
    }

    // Deal with remainders
    double sumr = 0.0;
    for (long i = size - remainder; i < size; ++i) {
        sumr += a[i] * b[i];
    }

    // Deal with accumulation
    double suma = 0.0;
    for (long i = 0; i < 2; ++i) {
        suma += acc[i];
    }

    return suma + sumr;
}