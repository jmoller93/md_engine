#pragma once
#ifndef EVALUATOR_DH
#define EVALUATOR_DH

#include "cutils_math.h"

class ChargeEvaluatorDH {
    public:
        float lambdai;
        float epsi;
        inline __device__ float3 force(float3 dr, float lenSqr, float qi, float qj, float multiplier) {
            float r2inv = 1.0f/lenSqr;
            float rinv = sqrtf(r2inv);
            float len = sqrtf(lenSqr);
            float forceScalar = qi*qj*epsi*expf(-len*lambdai)*(rinv*(rinv+lambdai)) * multiplier;
            return dr * forceScalar;
        }
        inline __device__ float energy(float lenSqr, float qi, float qj, float multiplier) {
            printf("DSF engs not implemented\n");
            return 0;
        }
        ChargeEvaluatorDH(float lambdai_, float epsi_) : lambdai(lambdai_), epsi(epsi_)  {};

};

#endif
