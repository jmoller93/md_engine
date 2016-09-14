#pragma once
#ifndef EVALUATOR_BASE_PAIR_3SPN2
#define EVALUATOR_BASE_PAIR_3SPN2

#include "cutils_math.h"
#include "cutils_3spn2.h"
#include "BasePair.h"
class BasePairEvaluator3SPN2{
public:
    float alpha;
    float range;
    BasePairEvaluator3SPN2() {}
    BasePairEvaluator3SPN2(float alpha_, float range_) {
        alpha = alpha_;
        range = range_;
    }
    //evaluator.force(theta, basepairType, s, distSqrs, directors, invDotProd);
    /*inline __device__ float3 force(BasePair3SPN2Type basepairType, float phi, float alpha, float range, float thetas[2], float3 cVector, float3 dVector, float scValues[3], float invMagProds[3], float invLens[3], float invLenSqrs[3], float c12Mags[3], float c, float3 directors[3], int MyIdxInBasePair) {
        float3 myForce = {0,0,0};
        float dPhi = phi - basepairType.phi0;
        float dTheta1 = thetas[0] - basepairType.theta1;
        float dTheta2 = thetas[1] - basepairType.theta2;
        float cone = range * 2.0f;
        
        float phiFactor = 0.5f * (1.0f + cosf(dPhi));
        float ftor = 0.5f * sinf(dPhi);
        float energyMors;
        float forceMors;

        if(lengthSqr(directors[1]) < basepairType.sigma) {
            //Use a purely repulsive Morse potential
            energyMors = morsRepEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
            forceMors = morsRepForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

            if (MyIdxInBasePair == 1) {
                myForce.x -= forceMors * directors[1].x;
                myForce.y -= forceMors * directors[1].y;
                myForce.z -= forceMors * directors[1].z;
            }
            else if (MyIdxInBasePair == 3) { 
                myForce.x += forceMors * directors[1].x;
                myForce.y += forceMors * directors[1].y;
                myForce.z += forceMors * directors[1].z;
            }
        }

        if ((dTheta1 >= -M_PI/(cone)) && (dTheta1 <= M_PI/(cone))) {
            if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {

                energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                if (MyIdxInBasePair == 0) {
                    float derivofPoten = -ftor * invLens[0] * scValues[0] * energyMors;
                    myForce.x += derivofPoten * cVector.x;
                    myForce.y += derivofPoten * cVector.y;
                    myForce.z += derivofPoten * cVector.z;
                }
                else if (MyIdxInBasePair == 1) {
                    //I may want to store the lengths instead of continuously calculating them
                    float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * energyMors;
                    float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * energyMors;
                    myForce.x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                    myForce.y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                    myForce.z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;
                }
                else if (MyIdxInBasePair == 2) {
                    float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * energyMors;
                    float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * energyMors;
                    myForce.x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                    myForce.y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                    myForce.z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                }
                else if (MyIdxInBasePair == 3) {
                    float derivofPoten = -ftor * invLens[2] * scValues[1] * energyMors;
                    myForce.x += derivofPoten * dVector.x;
                    myForce.y += derivofPoten * dVector.y;
                    myForce.z += derivofPoten * dVector.z;
                }
            }
            else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || ((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range)))) {
                float cosine2 = cosf(range*dTheta2);
                float sine2 = sinf(range*dTheta2);
                float hbon_cosine = 1.0f - cosine2 * cosine2;

                energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                float prefactor2 = cone * cosine2 * sine2 * 1.0f/ sqrtf(1.0f - thetas[1] * thetas[1]);

                if (MyIdxInBasePair == 1) {
                    float3 fr1;
                    fr1.x = prefactor2 * (invLens[1] * (thetas[1] * directors[1].x * invLens[1] - directors[2].x * invLens[2])) * energyMors - hbon_cosine * directors[2].x * forceMors;
                    fr1.y = prefactor2 * (invLens[1] * (thetas[1] * directors[1].y * invLens[1] - directors[2].y * invLens[2])) * energyMors - hbon_cosine * directors[2].y * forceMors;
                    fr1.z = prefactor2 * (invLens[1] * (thetas[1] * directors[1].z * invLens[1] - directors[2].z * invLens[2])) * energyMors - hbon_cosine * directors[2].z * forceMors;
                    myForce.x += fr1.x * phiFactor;
                    myForce.y += fr1.y * phiFactor;
                    myForce.z += fr1.z * phiFactor;
                }
                else if (MyIdxInBasePair == 2) {
                    float3 fr2;
                    fr2.x = prefactor2 * (invLens[2] * (directors[1].x * invLens[1] - thetas[1] * directors[2].x * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].x * invLens[1])) * energyMors + hbon_cosine * directors[1].x * forceMors;
                    fr2.y = prefactor2 * (invLens[2] * (directors[1].y * invLens[1] - thetas[1] * directors[2].y * invLens[2]) + invLens[1] * (directors[2].y * invLens[2] - thetas[1] * directors[1].y * invLens[1])) * energyMors + hbon_cosine * directors[1].y * forceMors;
                    fr2.z = prefactor2 * (invLens[2] * (directors[1].z * invLens[1] - thetas[1] * directors[2].z * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].z * invLens[1])) * energyMors + hbon_cosine * directors[1].z * forceMors;
                    myForce.x += fr2.x * phiFactor;
                    myForce.y += fr2.y * phiFactor;
                    myForce.z += fr2.z * phiFactor;
                }
                else if (MyIdxInBasePair == 3) {
                    float3 fr3;
                    fr3.x = prefactor2 * (invLens[2] * (directors[2].x * invLens[2] * thetas[1] - directors[1].x * invLens[1])) * energyMors;
                    fr3.y = prefactor2 * (invLens[2] * (directors[2].y * invLens[2] * thetas[1] - directors[1].y * invLens[1])) * energyMors;
                    fr3.z = prefactor2 * (invLens[2] * (directors[2].z * invLens[2] * thetas[1] - directors[1].z * invLens[1])) * energyMors;
                    myForce.x += fr3.x * phiFactor;
                    myForce.y += fr3.y * phiFactor;
                    myForce.z += fr3.z * phiFactor;
                }

                //This is mostly copied from Dan Hinckley
                if (MyIdxInBasePair == 0) {
                    float derivofPoten = -ftor * invLens[0] * scValues[0] * energyMors;
                    myForce.x += derivofPoten * cVector.x;
                    myForce.y += derivofPoten * cVector.y;
                    myForce.z += derivofPoten * cVector.z;
                }
                else if (MyIdxInBasePair == 1) {
                    //I may want to store the lengths instead of continuously calculating them
                    float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * energyMors;
                    float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * energyMors;
                    myForce.x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                    myForce.y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                    myForce.z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;
                }
                else if (MyIdxInBasePair == 2) {
                    float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * energyMors;
                    float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * energyMors;
                    myForce.x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                    myForce.y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                    myForce.z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                }
                else if (MyIdxInBasePair == 3) {
                    float derivofPoten = -ftor * invLens[2] * scValues[1] * energyMors;
                    myForce.x += derivofPoten * dVector.x;
                    myForce.y += derivofPoten * dVector.y;
                    myForce.z += derivofPoten * dVector.z;
                }

            }
            else if (((dTheta1 >= M_PI/(cone)) && (dTheta1 <= M_PI/(range))) || ((dTheta1 <= -M_PI/(cone)) && (dTheta1 >= -M_PI/(range)))) {
                float cosine = cosf(range*dTheta1);
                float sine = sinf(range*dTheta1);
                float hbon_cosine = 1.0f - cosine * cosine;

                if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {
                    energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                    forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                    float prefactor = cone * cosine * sine * 1.0f/ sqrtf(1.0f - thetas[0] * thetas[0]);

                    if (MyIdxInBasePair == 0) {
                        float3 fr1;
                        fr1.x = prefactor * (invLens[0] * (thetas[0] * directors[0].x * invLens[0] + directors[1].x * invLens[1])) * energyMors;
                        fr1.y = prefactor * (invLens[0] * (thetas[0] * directors[0].y * invLens[0] + directors[1].y * invLens[1])) * energyMors;
                        fr1.z = prefactor * (invLens[0] * (thetas[0] * directors[0].z * invLens[0] + directors[1].z * invLens[1])) * energyMors;
                        myForce.x += fr1.x * phiFactor;
                        myForce.y += fr1.y * phiFactor;
                        myForce.z += fr1.z * phiFactor;
                    }
                    else if (MyIdxInBasePair == 1) {
                        float3 fr2;
                        fr2.x = prefactor * (invLens[0] * (-directors[1].x * invLens[1] - thetas[0] * directors[0].x * invLens[0])+ invLens[1] * (directors[0].x * invLens[0] + thetas[0] * directors[1].x * invLens[1])) * energyMors - hbon_cosine * directors[1].x * forceMors;
                        fr2.y = prefactor * (invLens[0] * (-directors[1].y * invLens[1] - thetas[0] * directors[0].y * invLens[0]) + invLens[1] * (directors[0].y * invLens[0] + thetas[0] * directors[1].y * invLens[1])) * energyMors - hbon_cosine * directors[1].y * forceMors;
                        fr2.z = prefactor * (invLens[0] * (-directors[1].z * invLens[1] - thetas[0] * directors[0].z * invLens[0]) + invLens[1] * (directors[0].z * invLens[0] + thetas[0] * directors[1].z * invLens[1])) * energyMors - hbon_cosine * directors[1].z * forceMors;
                        myForce.x += fr2.x * phiFactor;
                        myForce.y += fr2.y * phiFactor;
                        myForce.z += fr2.z * phiFactor;
                    }
                    else if (MyIdxInBasePair == 3) {
                        float3 fr3;
                        fr3.x = prefactor * (invLens[1] * (-directors[1].x * invLens[1] * thetas[0] - directors[0].x * invLens[0])) * energyMors + hbon_cosine * directors[1].x * forceMors;
                        fr3.y = prefactor * (invLens[1] * (-directors[1].y * invLens[1] * thetas[0] - directors[0].y * invLens[0])) * energyMors + hbon_cosine * directors[1].y * forceMors;
                        fr3.z = prefactor * (invLens[1] * (-directors[1].z * invLens[1] * thetas[0] - directors[0].z * invLens[0])) * energyMors + hbon_cosine * directors[1].z * forceMors;
                        myForce.x += fr3.x * phiFactor;
                        myForce.y += fr3.y * phiFactor;
                        myForce.z += fr3.z * phiFactor;
                    }

                    float first_term = energyMors * hbon_cosine;

                    //This is mostly copied from Dan Hinckley
                    if (MyIdxInBasePair == 0) {
                        float derivofPoten = -ftor * invLens[0] * scValues[0] * first_term;
                        myForce.x += derivofPoten * cVector.x;
                        myForce.y += derivofPoten * cVector.y;
                        myForce.z += derivofPoten * cVector.z;
                    }
                    else if (MyIdxInBasePair == 1) {
                        //I may want to store the lengths instead of continuously calculating them
                        float frb1 = ftor * (length(directors[1]) - length(directors[0])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * first_term;
                        float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * first_term;
                        myForce.x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                        myForce.y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                        myForce.z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;
                    }
                    else if (MyIdxInBasePair == 2) {
                        float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[1] * first_term;
                        float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * first_term;
                        myForce.x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                        myForce.y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                        myForce.z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                    }
                    else if (MyIdxInBasePair == 3) {
                        float derivofPoten = -ftor * invLens[2] * scValues[1] * first_term;
                        myForce.x += derivofPoten * dVector.x;
                        myForce.y += derivofPoten * dVector.y;
                        myForce.z += derivofPoten * dVector.z;
                    }
                }
                else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || ((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range)))) {
                    float cosine2 = cosf(range*dTheta2);
                    float sine2 = sinf(range*dTheta2);
                    float hbon_cosine2 = 1.0f - cosine2 * cosine2;

                    energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                    forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                    float prefactor2 = cone * cosine2 * sine2 * 1.0f / sqrtf(1.0f - thetas[1] * thetas[1]);
                    float prefactor  = cone * cosine  * sine  * 1.0f / sqrtf(1.0f - thetas[0] * thetas[0]);

                    if (MyIdxInBasePair == 0) {
                        float3 fr1;
                        fr1.x = prefactor * (invLens[0] * (thetas[0] * directors[0].x * invLens[0] + directors[1].x * invLens[1])) * energyMors * hbon_cosine2;
                        fr1.y = prefactor * (invLens[0] * (thetas[0] * directors[0].y * invLens[0] + directors[1].y * invLens[1])) * energyMors * hbon_cosine2;
                        fr1.z = prefactor * (invLens[0] * (thetas[0] * directors[0].z * invLens[0] + directors[1].z * invLens[1])) * energyMors * hbon_cosine2;
                        myForce.x += fr1.x * phiFactor;
                        myForce.y += fr1.y * phiFactor;
                        myForce.z += fr1.z * phiFactor;
                    }
                    else if (MyIdxInBasePair == 1) {
                        float3 fr2;
                        fr2.x = prefactor * (invLens[0] * (-directors[1].x * invLens[1] - thetas[0] * directors[0].x * invLens[0]) + invLens[1] * (directors[0].x * invLens[0] - thetas[0] * directors[1].x * invLens[1])) * energyMors + prefactor2 * (invLens[1] * (thetas[1] * directors[1].x * invLens[1] - directors[2].x * invLens[2])) * hbon_cosine * energyMors - hbon_cosine * hbon_cosine2 * directors[1].x * forceMors; 
                        fr2.y = prefactor * (invLens[0] * (-directors[1].y * invLens[1] - thetas[0] * directors[0].y * invLens[0]) + invLens[1] * (directors[0].y * invLens[0] - thetas[0] * directors[1].y * invLens[1])) * energyMors + prefactor2 * (invLens[1] * (thetas[1] * directors[1].y * invLens[1] - directors[2].y * invLens[2])) * hbon_cosine * energyMors - hbon_cosine * hbon_cosine2 * directors[1].y * forceMors; 
                        fr2.z = prefactor * (invLens[0] * (-directors[1].z * invLens[1] - thetas[0] * directors[0].z * invLens[0]) + invLens[1] * (directors[0].z * invLens[0] - thetas[0] * directors[1].z * invLens[1])) * energyMors + prefactor2 * (invLens[1] * (thetas[1] * directors[1].z * invLens[1] - directors[2].z * invLens[2])) * hbon_cosine * energyMors - hbon_cosine * hbon_cosine2 * directors[1].z * forceMors; 
                        myForce.x += fr2.x * phiFactor;
                        myForce.y += fr2.y * phiFactor;
                        myForce.z += fr2.z * phiFactor;
                    }
                    else if (MyIdxInBasePair == 2) {
                        float3 fr3;
                        fr3.x = prefactor2 * (invLens[2] * (directors[2].x * invLens[2] * thetas[1] - directors[1].x * invLens[1])) * energyMors * hbon_cosine;
                        fr3.y = prefactor2 * (invLens[2] * (directors[2].y * invLens[2] * thetas[1] - directors[1].y * invLens[1])) * energyMors * hbon_cosine;
                        fr3.z = prefactor2 * (invLens[2] * (directors[2].z * invLens[2] * thetas[1] - directors[1].z * invLens[1])) * energyMors * hbon_cosine;
                        myForce.x += fr3.x * phiFactor;
                        myForce.y += fr3.y * phiFactor;
                        myForce.z += fr3.z * phiFactor;
                    }
                    else if (MyIdxInBasePair == 3) {
                        float3 fr4;
                        fr4.x = prefactor2 * (invLens[2] * (directors[1].x * invLens[1] - thetas[1] * directors[2].x * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].x * invLens[1])) * energyMors * hbon_cosine + prefactor * (invLens[1] * (- thetas[0] * directors[1].x * invLens[1] - directors[0].x * invLens[0])) * hbon_cosine2 * energyMors + hbon_cosine * hbon_cosine2 * directors[1].x * forceMors; 
                        fr4.y = prefactor2 * (invLens[2] * (directors[1].y * invLens[1] - thetas[1] * directors[2].y * invLens[2]) + invLens[1] * (directors[2].y * invLens[2] - thetas[1] * directors[1].y * invLens[1])) * energyMors * hbon_cosine + prefactor * (invLens[1] * (- thetas[0] * directors[1].y * invLens[1] - directors[0].y * invLens[0])) * hbon_cosine2 * energyMors + hbon_cosine * hbon_cosine2 * directors[1].y * forceMors; 
                        fr4.z = prefactor2 * (invLens[2] * (directors[1].z * invLens[1] - thetas[1] * directors[2].z * invLens[2]) + invLens[1] * (directors[2].z * invLens[2] - thetas[1] * directors[1].z * invLens[1])) * energyMors * hbon_cosine + prefactor * (invLens[1] * (- thetas[0] * directors[1].z * invLens[1] - directors[0].z * invLens[0])) * hbon_cosine2 * energyMors + hbon_cosine * hbon_cosine2 * directors[1].z * forceMors; 
                        myForce.x += fr4.x * phiFactor;
                        myForce.y += fr4.y * phiFactor;
                        myForce.z += fr4.z * phiFactor;
                    }

                    float first_term = energyMors * hbon_cosine * hbon_cosine2;
                    //This is mostly copied from Dan Hinckley
                    if (MyIdxInBasePair == 0) {
                        float derivofPoten = -ftor * invLens[0] * scValues[0] * first_term;
                        myForce.x += derivofPoten * cVector.x;
                        myForce.y += derivofPoten * cVector.y;
                        myForce.z += derivofPoten * cVector.z;
                    }
                    else if (MyIdxInBasePair == 1) {
                        //I may want to store the lengths instead of continuously calculating them
                        float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * first_term;
                        float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * first_term;
                        myForce.x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                        myForce.y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                        myForce.z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;
                    }
                    else if (MyIdxInBasePair == 2) {
                        float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * first_term;
                        float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * first_term;
                        myForce.x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                        myForce.y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                        myForce.z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                    }
                    else if (MyIdxInBasePair == 3) {
                        float derivofPoten = -ftor * invLens[2] * scValues[1] * first_term;
                        myForce.x += derivofPoten * dVector.x;
                        myForce.y += derivofPoten * dVector.y;
                        myForce.z += derivofPoten * dVector.z;
                    }

                }

            }
        }
        return myForce;


    }*/


    inline __device__ void forces(BasePair3SPN2Type basepairType, float phi, float alpha, float range, float thetas[2], float3 cVector, float3 dVector, float scValues[3], float invMagProds[3], float invLens[3], float invLenSqrs[3], float c12Mags[3], float c, float3 directors[3], float3 forces[4]) {
        float dPhi = phi - basepairType.phi0;
        float dTheta1 = thetas[0] - basepairType.theta1;
        float dTheta2 = thetas[1] - basepairType.theta2;
        
        float phiFactor = 0.5f * (1.0f + cosf(dPhi));
        float ftor = 0.5f * sinf(dPhi);
        float cone = range * 2.0f;
        float energyMors;
        float forceMors;

        if(lengthSqr(directors[1]) < basepairType.sigma) {
            //Use a purely repulsive Morse potential
            energyMors = morsRepEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
            forceMors = morsRepForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

            forces[1].x -= forceMors * directors[1].x;
            forces[1].y -= forceMors * directors[1].y;
            forces[1].z -= forceMors * directors[1].z;
            forces[3].x += forceMors * directors[1].x;
            forces[3].y += forceMors * directors[1].y;
            forces[3].z += forceMors * directors[1].z;
        }

        if ((dTheta1 >= -M_PI/(cone)) && (dTheta1 <= M_PI/(cone))) {
            if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {

                energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                float fra1 = -ftor * invLens[0] * scValues[0] * energyMors;
                forces[0].x += fra1 * cVector.x;
                forces[0].y += fra1 * cVector.y;
                forces[0].z += fra1 * cVector.z;

                //I may want to store the lengths instead of continuously calculating them
                float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * energyMors;
                float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * energyMors;
                forces[1].x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                forces[1].y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                forces[1].z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;

                float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * energyMors;
                float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * energyMors;
                forces[2].x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                forces[2].y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                forces[2].z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;

                float frd1 = -ftor * invLens[2] * scValues[1] * energyMors;
                forces[3].x += frd1 * dVector.x;
                forces[3].y += frd1 * dVector.y;
                forces[3].z += frd1 * dVector.z;
            }
            else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || ((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range)))) {
                float cosine2 = cosf(range*dTheta2);
                float sine2 = sinf(range*dTheta2);
                float hbon_cosine = 1.0f - cosine2 * cosine2;

                energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                float prefactor2 = cone * cosine2 * sine2 * 1.0f/ sqrtf(1.0f - thetas[1] * thetas[1]);

                float3 fr1;
                fr1.x = prefactor2 * (invLens[1] * (thetas[1] * directors[1].x * invLens[1] - directors[2].x * invLens[2])) * energyMors - hbon_cosine * directors[2].x * forceMors;
                fr1.y = prefactor2 * (invLens[1] * (thetas[1] * directors[1].y * invLens[1] - directors[2].y * invLens[2])) * energyMors - hbon_cosine * directors[2].y * forceMors;
                fr1.z = prefactor2 * (invLens[1] * (thetas[1] * directors[1].z * invLens[1] - directors[2].z * invLens[2])) * energyMors - hbon_cosine * directors[2].z * forceMors;
                forces[1].x += fr1.x * phiFactor;
                forces[1].y += fr1.y * phiFactor;
                forces[1].z += fr1.z * phiFactor;

                float3 fr2;
                fr2.x = prefactor2 * (invLens[2] * (directors[1].x * invLens[1] - thetas[1] * directors[2].x * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].x * invLens[1])) * energyMors + hbon_cosine * directors[1].x * forceMors;
                fr2.y = prefactor2 * (invLens[2] * (directors[1].y * invLens[1] - thetas[1] * directors[2].y * invLens[2]) + invLens[1] * (directors[2].y * invLens[2] - thetas[1] * directors[1].y * invLens[1])) * energyMors + hbon_cosine * directors[1].y * forceMors;
                fr2.z = prefactor2 * (invLens[2] * (directors[1].z * invLens[1] - thetas[1] * directors[2].z * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].z * invLens[1])) * energyMors + hbon_cosine * directors[1].z * forceMors;
                forces[2].x += fr2.x * phiFactor;
                forces[2].y += fr2.y * phiFactor;
                forces[2].z += fr2.z * phiFactor;

                float3 fr3;
                fr3.x = prefactor2 * (invLens[2] * (directors[2].x * invLens[2] * thetas[1] - directors[1].x * invLens[1])) * energyMors;
                fr3.y = prefactor2 * (invLens[2] * (directors[2].y * invLens[2] * thetas[1] - directors[1].y * invLens[1])) * energyMors;
                fr3.z = prefactor2 * (invLens[2] * (directors[2].z * invLens[2] * thetas[1] - directors[1].z * invLens[1])) * energyMors;
                forces[3].x += fr3.x * phiFactor;
                forces[3].y += fr3.y * phiFactor;
                forces[3].z += fr3.z * phiFactor;

                float derivofPoten = -ftor * invLens[0] * scValues[0] * energyMors;
                forces[0].x += derivofPoten * cVector.x;
                forces[0].y += derivofPoten * cVector.y;
                forces[0].z += derivofPoten * cVector.z;

                //I may want to store the lengths instead of continuously calculating them
                float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * energyMors;
                float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * energyMors;
                forces[1].x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                forces[1].y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                forces[1].z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;

                float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * energyMors;
                float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * energyMors;
                forces[2].x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                forces[2].y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                forces[2].z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;

                float frd1 = -ftor * invLens[2] * scValues[1] * energyMors;
                forces[3].x += frd1 * dVector.x;
                forces[3].y += frd1 * dVector.y;
                forces[3].z += frd1 * dVector.z;

            }
            else if (((dTheta1 >= M_PI/(cone)) && (dTheta1 <= M_PI/(range))) || ((dTheta1 <= -M_PI/(cone)) && (dTheta1 >= -M_PI/(range)))) {
                float cosine = cosf(range*dTheta1);
                float sine = sinf(range*dTheta1);
                float hbon_cosine = 1.0f - cosine * cosine;

                if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {
                    energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                    forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                    float prefactor = cone * cosine * sine * 1.0f/ sqrtf(1.0f - thetas[0] * thetas[0]);

                    float3 fr1;
                    fr1.x = prefactor * (invLens[0] * (thetas[0] * directors[0].x * invLens[0] + directors[1].x * invLens[1])) * energyMors;
                    fr1.y = prefactor * (invLens[0] * (thetas[0] * directors[0].y * invLens[0] + directors[1].y * invLens[1])) * energyMors;
                    fr1.z = prefactor * (invLens[0] * (thetas[0] * directors[0].z * invLens[0] + directors[1].z * invLens[1])) * energyMors;
                    forces[0].x += fr1.x * phiFactor;
                    forces[0].y += fr1.y * phiFactor;
                    forces[0].z += fr1.z * phiFactor;

                    float3 fr2;
                    fr2.x = prefactor * (invLens[0] * (-directors[1].x * invLens[1] - thetas[0] * directors[0].x * invLens[0]) + invLens[1] * (directors[0].x * invLens[0]+ thetas[0] * directors[1].x * invLens[1])) * energyMors - hbon_cosine * directors[1].x * forceMors;
                    fr2.y = prefactor * (invLens[0] * (-directors[1].y * invLens[1] - thetas[0] * directors[0].y * invLens[0]) + invLens[1] * (directors[0].y * invLens[0]+ thetas[0] * directors[1].y * invLens[1])) * energyMors - hbon_cosine * directors[1].y * forceMors;
                    fr2.z = prefactor * (invLens[0] * (-directors[1].z * invLens[1] - thetas[0] * directors[0].z * invLens[0]) + invLens[1] * (directors[0].z * invLens[0]+ thetas[0] * directors[1].z * invLens[1])) * energyMors - hbon_cosine * directors[1].z * forceMors;
                    forces[1].x += fr2.x * phiFactor;
                    forces[1].y += fr2.y * phiFactor;
                    forces[1].z += fr2.z * phiFactor;

                    float3 fr3;
                    fr3.x = prefactor * (invLens[1] * (-directors[1].x * invLens[1] * thetas[0] - directors[0].x * invLens[0])) * energyMors + hbon_cosine * directors[1].x * forceMors;
                    fr3.y = prefactor * (invLens[1] * (-directors[1].y * invLens[1] * thetas[0] - directors[0].y * invLens[0])) * energyMors + hbon_cosine * directors[1].y * forceMors;
                    fr3.z = prefactor * (invLens[1] * (-directors[1].z * invLens[1] * thetas[0] - directors[0].z * invLens[0])) * energyMors + hbon_cosine * directors[1].z * forceMors;
                    forces[3].x += fr3.x * phiFactor;
                    forces[3].y += fr3.y * phiFactor;
                    forces[3].z += fr3.z * phiFactor;

                    float first_term = energyMors * hbon_cosine;

                    //This is mostly copied from Dan Hinckley
                    float fra1 = -ftor * invLens[0] * scValues[0] * first_term;
                    forces[0].x += fra1 * cVector.x;
                    forces[0].y += fra1 * cVector.y;
                    forces[0].z += fra1 * cVector.z;

                    //I may want to store the lengths instead of continuously calculating them
                    float frb1 = ftor * (length(directors[1]) - length(directors[0])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * first_term;
                    float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * first_term;
                    forces[1].x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                    forces[1].y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                    forces[1].z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;

                    float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[1] * first_term;
                    float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * first_term;
                    forces[2].x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                    forces[2].y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                    forces[2].z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                    
                    float frd1 = -ftor * invLens[2] * scValues[1] * first_term;
                    forces[3].x += frd1 * dVector.x;
                    forces[3].y += frd1 * dVector.y;
                    forces[3].z += frd1 * dVector.z;
                }
                else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || (((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range))))) {
                    float cosine2 = cosf(range*dTheta2);
                    float sine2 = sinf(range*dTheta2);
                    float hbon_cosine2 = 1.0f - cosine2 * cosine2;

                    energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                    forceMors = morsAttrForc(invLens[1], alpha, basepairType.epsi, basepairType.sigma);

                    float prefactor2 = cone * cosine2 * sine2 * 1.0f / sqrtf(1.0f - thetas[1] * thetas[1]);
                    float prefactor  = cone * cosine  * sine  * 1.0f / sqrtf(1.0f - thetas[0] * thetas[0]);

                    float3 fr1;
                    fr1.x = prefactor * (invLens[0] * (thetas[0] * directors[0].x * invLens[0] + directors[1].x * invLens[1])) * energyMors * hbon_cosine2;
                    fr1.y = prefactor * (invLens[0] * (thetas[0] * directors[0].y * invLens[0] + directors[1].y * invLens[1])) * energyMors * hbon_cosine2;
                    fr1.z = prefactor * (invLens[0] * (thetas[0] * directors[0].z * invLens[0] + directors[1].z * invLens[1])) * energyMors * hbon_cosine2;
                    forces[0].x += fr1.x * phiFactor;
                    forces[0].y += fr1.y * phiFactor;
                    forces[0].z += fr1.z * phiFactor;
                    
                    float3 fr2;
                    fr2.x = prefactor * (invLens[0] * (-directors[1].x * invLens[1] - thetas[0] * directors[0].x * invLens[0]) + invLens[1] * (directors[0].x * invLens[0] - thetas[0] * directors[1].x * invLens[1])) * energyMors + prefactor2 * (invLens[1] * (thetas[1] * directors[1].x * invLens[1] - directors[2].x * invLens[2])) * hbon_cosine * energyMors - hbon_cosine * hbon_cosine2 * directors[1].x * forceMors; 
                    fr2.y = prefactor * (invLens[0] * (-directors[1].y * invLens[1] - thetas[0] * directors[0].y * invLens[0]) + invLens[1] * (directors[0].y * invLens[0] - thetas[0] * directors[1].y * invLens[1])) * energyMors + prefactor2 * (invLens[1] * (thetas[1] * directors[1].y * invLens[1] - directors[2].y * invLens[2])) * hbon_cosine * energyMors - hbon_cosine * hbon_cosine2 * directors[1].y * forceMors; 
                    fr2.z = prefactor * (invLens[0] * (-directors[1].z * invLens[1] - thetas[0] * directors[0].z * invLens[0]) + invLens[1] * (directors[0].z * invLens[0] - thetas[0] * directors[1].z * invLens[1])) * energyMors + prefactor2 * (invLens[1] * (thetas[1] * directors[1].z * invLens[1] - directors[2].z * invLens[2])) * hbon_cosine * energyMors - hbon_cosine * hbon_cosine2 * directors[1].z * forceMors; 
                    forces[1].x += fr2.x * phiFactor;
                    forces[1].y += fr2.y * phiFactor;
                    forces[1].z += fr2.z * phiFactor;
                
                    float3 fr3;
                    fr3.x = prefactor2 * (invLens[2] * (directors[2].x * invLens[2] * thetas[1] - directors[1].x * invLens[1])) * energyMors * hbon_cosine;
                    fr3.y = prefactor2 * (invLens[2] * (directors[2].y * invLens[2] * thetas[1] - directors[1].y * invLens[1])) * energyMors * hbon_cosine;
                    fr3.z = prefactor2 * (invLens[2] * (directors[2].z * invLens[2] * thetas[1] - directors[1].z * invLens[1])) * energyMors * hbon_cosine;
                    forces[2].x += fr3.x * phiFactor;
                    forces[2].y += fr3.y * phiFactor;
                    forces[2].z += fr3.z * phiFactor;

                    float3 fr4;
                    fr4.x = prefactor2 * (invLens[2] * (directors[1].x * invLens[1] - thetas[1] * directors[2].x * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].x * invLens[1])) * energyMors * hbon_cosine + prefactor * (invLens[1] * (- thetas[0] * directors[1].x * invLens[1] - directors[0].x * invLens[0])) * hbon_cosine2 * energyMors + hbon_cosine * hbon_cosine2 * directors[1].x * forceMors; 
                    fr4.y = prefactor2 * (invLens[2] * (directors[1].y * invLens[1] - thetas[1] * directors[2].y * invLens[2]) + invLens[1] * (directors[2].y * invLens[2] - thetas[1] * directors[1].y * invLens[1])) * energyMors * hbon_cosine + prefactor * (invLens[1] * (- thetas[0] * directors[1].y * invLens[1] - directors[0].y * invLens[0])) * hbon_cosine2 * energyMors + hbon_cosine * hbon_cosine2 * directors[1].y * forceMors; 
                    fr4.z = prefactor2 * (invLens[2] * (directors[1].z * invLens[1] - thetas[1] * directors[2].z * invLens[2]) + invLens[1] * (directors[2].z * invLens[2] - thetas[1] * directors[1].z * invLens[1])) * energyMors * hbon_cosine + prefactor * (invLens[1] * (- thetas[0] * directors[1].z * invLens[1] - directors[0].z * invLens[0])) * hbon_cosine2 * energyMors + hbon_cosine * hbon_cosine2 * directors[1].z * forceMors; 
                    forces[3].x += fr4.x * phiFactor;
                    forces[3].y += fr4.y * phiFactor;
                    forces[3].z += fr4.z * phiFactor;

                    float first_term = energyMors * hbon_cosine * hbon_cosine2;
                    //This is mostly copied from Dan Hinckley
                    float fra1 = -ftor * invLens[0] * scValues[0] * first_term;
                    forces[0].x += fra1 * cVector.x;
                    forces[0].y += fra1 * cVector.y;
                    forces[0].z += fra1 * cVector.z;

                    //I may want to store the lengths instead of continuously calculating them
                    float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * first_term;
                    float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * first_term;
                    forces[1].x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                    forces[1].y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                    forces[1].z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;

                    float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * first_term;
                    float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * first_term;
                    forces[2].x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                    forces[2].y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                    forces[2].z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                    
                    float frd1 = -ftor * invLens[2] * scValues[1] * first_term;
                    forces[3].x += frd1 * dVector.x;
                    forces[3].y += frd1 * dVector.y;
                    forces[3].z += frd1 * dVector.z;

                }

            }
        }

    }


    inline __device__ float energy(BasePair3SPN2Type basepairType, float phi, float alpha, float range, float thetas[2], float invLens[3], float invLenSqrs[3], float3 directors[3]) {
        float dPhi = phi - basepairType.phi0;
        float dTheta1 = thetas[0] - basepairType.theta1;
        float dTheta2 = thetas[1] - basepairType.theta2;
        float cone = 2.0f * range;
        
        float phiFactor = 0.5f * (1.0f + cosf(dPhi));
        float ftor = 0.5f * sinf(dPhi);
        float energyMors;
        float eBasePair = 0.0f;

        if(lengthSqr(directors[1]) < basepairType.sigma) {
            //Use a purely repulsive Morse potential
            energyMors = morsRepEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
            eBasePair += energyMors * phiFactor;

        }

        if ((dTheta1 >= -M_PI/(cone)) && (dTheta1 <= M_PI/(cone))) {
            if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {

                energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                eBasePair += phiFactor * energyMors;

            }
            else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || ((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range)))) {
                float cosine2 = cosf(range*dTheta2);
                float sine2 = sinf(range*dTheta2);
                float hbon_cosine = 1.0f - cosine2 * cosine2;

                energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                eBasePair += phiFactor * hbon_cosine * energyMors;

            }
            else if (((dTheta1 >= M_PI/(cone)) && (dTheta1 <= M_PI/(range))) || ((dTheta1 <= -M_PI/(cone)) && (dTheta1 >= -M_PI/(range)))) {
                float cosine = cosf(range*dTheta1);
                float sine = sinf(range*dTheta1);
                float hbon_cosine = 1.0f - cosine * cosine;

                if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {
                    energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                    eBasePair += phiFactor * energyMors;

                }
                else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || ((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range)))) {
                    float cosine2 = cosf(range*dTheta2);
                    float sine2 = sinf(range*dTheta2);
                    float hbon_cosine2 = 1.0f - cosine2 * cosine2;

                    energyMors = morsAttrEnrgy(invLens[1], alpha, basepairType.epsi, basepairType.sigma);
                    eBasePair += phiFactor * hbon_cosine2 * energyMors;


                }
            }
        }
        return eBasePair;
    }
};

#endif

