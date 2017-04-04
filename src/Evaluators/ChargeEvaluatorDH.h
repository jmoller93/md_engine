#pragma once
#ifndef EVALUATOR_DH
#define EVALUATOR_DH

#include "cutils_math.h"

class ChargeEvaluatorDH {
    public:
        float lambdai;
        float epsi;
        float qqr_to_eng;
        inline __device__ float3 force(float3 dr, float lenSqr, float qi, float qj, float multiplier, float molParam) {
            float r2inv = 1.0f/lenSqr;
            float rinv = sqrtf(r2inv);
            float len = sqrtf(lenSqr);
            float forceScalar = qqr_to_eng*qi*qj*epsi*expf(-len*lambdai)*r2inv*(rinv+lambdai) * multiplier * molParam;
            //if(forceScalar != 0) { printf("force is this %f\n", forceScalar);}
            return dr * forceScalar;
        }
        inline __device__ float energy(float lenSqr, float qi, float qj, float multiplier) {
            float len = sqrtf(lenSqr);
            return qqr_to_eng*qi*qj*epsi*expf(-len*lambdai) * multiplier / len;
        }
        ChargeEvaluatorDH(float lambdai_, float epsi_, float qqr_to_eng_) : lambdai(lambdai_), epsi(epsi_), qqr_to_eng(qqr_to_eng_)  {};

};

#endif
