#pragma once
#ifndef BONDEVALUATORHARMONICEXTEND_H
#define BONDEVALUATORHARMONICEXTEND_H

#include "Bond.h"

class BondEvaluatorHarmonicExtend {
public:
    inline __device__ float3 force(float3 bondVec, float rSqr, BondHarmonicExtendType bondType) {
        float r = sqrtf(rSqr);
        float dr = r - bondType.r0;
        float rk = 2.0f *bondType.k * dr;
        float rk3 = rk * 200.0f * dr * dr;
        //float rk3 = rk * dr * dr;
        if (r > 0) {//MAKE SURE ALL THIS WORKS, I JUST BORROWED FROM LAMMPS
            float fBond = -(rk+rk3)/r;
            //printf("f %f\n", fBond);
            return bondVec * fBond;
        } 
        return make_float3(0, 0, 0);
    }
    inline __device__ float energy(float3 bondVec, float rSqr, BondHarmonicExtendType bondType) {
        float r = sqrtf(rSqr);
        float dr = r - bondType.r0;
        float dr2 = dr * dr;
        float eng = bondType.k * dr2 *( 1.0f + 100.0f * dr2);
        return 0.5f * eng; //0.5 for splitting between atoms
    }
};
#endif
