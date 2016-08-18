#pragma once
#ifndef BONDEVALUATORGOLIKE_H
#define BONDEVALUATORGOLIKE_H

#include "Bond.h"

class BondEvaluatorGoLike{
public:
    inline __device__ float3 force(float3 bondVec, float rSqr, BondGoLikeType bondType) {
        float eps = bondType.eps;
        float sig = bondType.sig;
        float sr2 = sig*sig/rSqr;
        float sr6 = sr2*sr2*sr2;
        float sr12 = sr6*sr6;
        float sr10 = sr6*sr2*sr2;
        
        float fbond = (240.0f * eps) / rSqr * (sr12 - sr10);

        float3 force = bondVec * fbond;
        return force;
    }
    inline __device__ float energy(float3 bondVec, float rSqr, BondGoLikeType bondType) {
        float eps = bondType.eps;
        float sig = bondType.sig;
        float sr2 = sig*sig/rSqr;
        float sr6 = sr2*sr2*sr2;
        float sr12 = sr6*sr6;
        float sr10 = sr6*sr2*sr2;
        float eng = 4.0f * eps * (5.0f * sr12 - 60.f * sr10);
        return 0.5f * eng; //0.5 for splitting between atoms
    }
};
#endif
