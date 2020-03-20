
#pragma once

#include <iostream>
//#include <zmmintrin.h>
#include <immintrin.h>

/*
				__declspec(align(64)) Real S1[4];
				__declspec(align(64)) Real S2[4];

				__m256d _InvBlock = _mm256_set_pd(0.0, InverseBlock(2, j), InverseBlock(1, j), InverseBlock(0, j));
				__m256d _R1Block = _mm256_set_pd(0.0, R1_block(k_index, ei, 2), R1_block(k_index, ei, 1), R1_block(k_index, ei, 0));				__m256d _R2Block = _mm256_set_pd(0.0, R2_block(k_index, ei, 2), R2_block(k_index, ei, 1), R2_block(k_index, ei, 0));				__m256d _Prod1 = _mm256_mul_pd(_R1Block, _InvBlock);
				__m256d _Prod2 = _mm256_mul_pd(_R2Block, _InvBlock);

				_mm256_store_pd(S1, _Prod1);
				_mm256_store_pd(S2, _Prod2);

				Sum1 = S1[0] + S1[1] + S1[2] + S1[3];
				Sum2 = S2[0] + S2[1] + S2[2] + S2[3];


*/


// !!!!!!!!!!!!!!!!!!!
// https://github.com/JGU-HPC/parallelprogrammingbook/blob/master/chapter3/AVX/matrix_matrix_mult.cpp

/*

v#include <random>       // prng
#include <cstdint>      // uint32_t
#include <iostream>     // std::cout
#include <immintrin.h>  // AVX intrinsics

// timers distributed with this book
#include "../include/hpc_helpers.hpp"

void init(float * data, uint64_t length) {

    std::mt19937 engine(42);
    std::uniform_real_distribution<float> density(-1, 1);

    for (uint64_t i = 0; i < length; i++)
        data[i] = density(engine);
}

inline float hsum_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 maxs = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, maxs); // high half -> low half
    maxs        = _mm_add_ss(maxs, shuf);
    return        _mm_cvtss_f32(maxs);
}

inline float hsum_avx(__m256 v) {
    __m128 lo = _mm256_castps256_ps128(v);   // low 128
    __m128 hi = _mm256_extractf128_ps(v, 1); // high 128
           lo = _mm_add_ps(lo, hi);          // max the low 128
    return hsum_sse3(lo);                    // and inline the sse3 version
}

void plain_dmm(float * A,
               float * B,
               float * C,
               uint64_t M,
               uint64_t L,
               uint64_t N,
               bool parallel) {

    #pragma omp parallel for collapse(2) if(parallel)
    for (uint64_t i = 0; i < M; i++)
        for (uint64_t j = 0; j < N; j++) {
            float accum = float(0);
            for (uint64_t k = 0; k < L; k++)
                accum += A[i*L+k]*B[j*L+k];
            C[i*N+j] = accum;
       }
}

void avx_dmm(float * A,
             float * B,
             float * C,
             uint64_t M,
             uint64_t L,
             uint64_t N,
             bool parallel) {

    #pragma omp parallel for collapse(2) if(parallel)
    for (uint64_t i = 0; i < M; i++)
        for (uint64_t j = 0; j < N; j++) {

            __m256 X = _mm256_setzero_ps();
            for (uint64_t k = 0; k < L; k += 8) {
                const __m256 AV = _mm256_load_ps(A+i*L+k);
                const __m256 BV = _mm256_load_ps(B+j*L+k);
                X = _mm256_add_ps(X, _mm256_mul_ps(AV, BV));
            }

            C[i*N+j] = hsum_avx(X);
       }
}

void avx_dmm_unroll_2(float * A,
                      float * B,
                      float * C,
                      uint64_t M,
                      uint64_t L,
                      uint64_t N,
                      bool parallel) {

    #pragma omp parallel for collapse(2) if(parallel)
    for (uint64_t i = 0; i < M; i++)
        for (uint64_t j = 0; j < N; j++) {

            __m256 X = _mm256_setzero_ps();
            __m256 Y = _mm256_setzero_ps();
            for (uint64_t k = 0; k < L; k += 16) {
                const __m256 AVX = _mm256_load_ps(A+i*L+k+0);
                const __m256 BVX = _mm256_load_ps(B+j*L+k+0);
                const __m256 AVY = _mm256_load_ps(A+i*L+k+8);
                const __m256 BVY = _mm256_load_ps(B+j*L+k+8);
                X = _mm256_add_ps(X, _mm256_mul_ps(AVX, BVX));
                Y = _mm256_add_ps(X, _mm256_mul_ps(AVY, BVY));
            }

            C[i*N+j] = hsum_avx(X)+hsum_avx(Y);
       }
}

int main () {

    const uint64_t M = 1UL <<  10;
    const uint64_t L = 1UL <<  11;
    const uint64_t N = 1UL <<  12;

    TIMERSTART(alloc_memory)
    auto A = static_cast<float*>(_mm_malloc(M*L*sizeof(float) , 32));
    auto B = static_cast<float*>(_mm_malloc(N*L*sizeof(float) , 32));
    auto C = static_cast<float*>(_mm_malloc(M*N*sizeof(float) , 32));
    TIMERSTOP(alloc_memory)

    TIMERSTART(init)
    init(A, M*L);
    init(B, N*L);
    TIMERSTOP(init)

    TIMERSTART(plain_dmm_single)
    plain_dmm(A, B, C, M, L, N, false);
    TIMERSTOP(plain_dmm_single)

    TIMERSTART(plain_dmm_multi)
    plain_dmm(A, B, C, M, L, N, true);
    TIMERSTOP(plain_dmm_multi)

    TIMERSTART(avx_dmm_single)
    avx_dmm(A, B, C, M, L, N, false);
    TIMERSTOP(avx_dmm_single)

    TIMERSTART(avx_dmm_multi)
    avx_dmm(A, B, C, M, L, N, true);
    TIMERSTOP(avx_dmm_multi)

    TIMERSTART(avx_dmm_unroll_2_single)
    avx_dmm_unroll_2(A, B, C, M, L, N, false);
    TIMERSTOP(avx_dmm_unroll_2_single)

    TIMERSTART(avx_dmm_unroll_2_multi)
    avx_dmm_unroll_2(A, B, C, M, L, N, true);
    TIMERSTOP(avx_dmm_unroll_2_multi)

    TIMERSTART(free_memory)
    _mm_free(A);
    _mm_free(B);
    _mm_free(C);
    TIMERSTOP(free_memory)

C++ functions for the computation of transposed matrix multiplication : 
-plain_tmm is a straightforward sequential non - vectorized implementation and 
-avx2_tmm is an AVX2 - vectorized implementation.

We first need to include the immintrin.h header file to build an application that uses intrinsics.

The inner loop computes the scalar product of row i of A and row j of B.
The vectorized function operates on eight values from A and B simultaneously while 
its non - vectorized counterpart operates only on single elements.

The commands _mm256_load_ps(A + i * L + k) and 
_mm256_load_ps(B + j * L + k) load eight consecutive single - precision floating - point numbers 
from the matrices A and B^T into the 256 - bit registers AV and BV.

Note that this intrinsic only loads data from aligned memory addresses.
Thus, we need to guarantee that the two matrices have actually been allocated on 32 - byte boundaries.
This can for example be realized with the _mm_malloc command :

auto A = static_cast<float*>(_mm_malloc(M*L * sizeof(float), 32));
auto B = static_cast<float*>(_mm_malloc(N*L * sizeof(float), 32));

The intrinsic function _mm256_fmadd_ps(AV, BV, X) then multiplies the eight floats of AV 
by the eight floats of BV and adds each of them to the values stored in the vector X

At the end of the inner loop X stores eight partial scalar products.
These eight values are then added horizontally with the call to a user - defined function hsum_avx(X) 
(the implementation of which is left as an exercise in Section 3.3) to produce the full scalar product.







inline float hsum_sse3(__m128 v) {
	__m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
	__m128 maxs = _mm_add_ps(v, shuf);
	shuf        = _mm_movehl_ps(shuf, maxs); // high half -> low half
	maxs        = _mm_add_ss(maxs, shuf);
	return        _mm_cvtss_f32(maxs);
}

inline float hsum_avx(__m256 v) {
	__m128 lo = _mm256_castps256_ps128(v);   // low 128
	__m128 hi = _mm256_extractf128_ps(v, 1); // high 128
		   lo = _mm_add_ps(lo, hi);          // max the low 128
	return hsum_sse3(lo);                    // and inline the sse3 version
}


*/


