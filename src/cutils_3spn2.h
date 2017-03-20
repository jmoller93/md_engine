#pragma once
#ifndef CUTILS_3SPN2_H
#define CUTILS_3SPN2_H

#include "globalDefs.h"
#include "cutils_math.h"

inline __device__ float morsAttrEnrgy(const float drsi, const float alfa, const float epsi, const float sigma) {
    float dist = 1.0f / drsi;
    float ener;
    if (dist > sigma) {
        float argu = alfa * (dist - sigma);
        float mors = (1.0f - expf(-argu));
        ener = epsi * mors * mors - epsi;
    }
    else {
        ener = -epsi;
    }
    return ener;
}

inline __device__ float morsAttrForc(const float drsi, const float alfa, const float epsi, const float sigma) {
    float dist = 1.0f / drsi;
    float forc;
    if (dist > sigma) {
        float argu = alfa * (dist - sigma);
        float mors = expf(-argu);
        forc = -2.0f * alfa * epsi * drsi * mors * (1.0f - mors);
    }
    else {
        forc = 0.0f;
    }
    return forc;
}

inline __device__ float morsRepEnrgy(const float drsi, const float alfa, const float epsi, const float sigma) {
    float dist = 1.0f / drsi;
    float ener;
    if (dist < sigma) {
        float argu = alfa * (dist - sigma);
        float mors = (1.0f - expf(-argu));
        ener = epsi * mors * mors;
    }
    else {
        ener = 0.0f;
    }
    return ener;
}

inline __device__ float morsRepForc(const float drsi, const float alfa, const float epsi, const float sigma) {
    float dist = 1.0f / drsi;
    float forc;
    if (dist < sigma) {
        float argu = alfa * (dist - sigma);
        float mors = expf(-argu);
        forc = -2.0f * alfa * epsi * drsi * mors * (1.0f - mors);
    }
    else {
        forc = 0.0f;
    }
    return forc;
}

#endif
