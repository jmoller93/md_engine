#pragma once
#ifndef EVALUATOR_LJREPUL
#define EVALUATOR_LJREPUL

#include "cutils_math.h"

class EvaluatorLJRepul {
    public:
        inline __device__ float3 force(float3 dr, float params[3], float lenSqr, float multiplier) {
            if (multiplier)  {
                if(lenSqr < params[2]/0.89/0.89)
                {
                    float epstimes24 = params[1];
                    float sig6 = params[2]*params[2]*params[2];
                    float p1 = epstimes24*2.0f*sig6*sig6;
                    float p2 = epstimes24*sig6;
                    float r2inv = 1.0f/lenSqr;
                    float r6inv = r2inv*r2inv*r2inv;
                    float forceScalar = r6inv * r2inv * (p1 * r6inv - p2) * multiplier;
                    printf("dist is %f\t sigm is %f\t forc is %f\n", sqrtf(lenSqr), sqrtf(params[2])/0.89, forceScalar);
                    return dr * forceScalar;
                }
                return make_float3(0,0,0);
            }
            return make_float3(0,0,0);
        }
        inline __device__ float energy(float params[3], float lenSqr, float multiplier) {
            float epstimes24 = params[1];
            float sig6 = params[2]*params[2]*params[2];
            float r2inv = 1/lenSqr;
            float r6inv = r2inv*r2inv*r2inv;
            float sig6r6inv = sig6 * r6inv;
            return 0.5f * (epstimes24 / 24)*sig6r6inv*(sig6r6inv-2.0f) * multiplier; //0.5 b/c we need to half-count energy b/c pairs are redundant
        }

};

#endif
