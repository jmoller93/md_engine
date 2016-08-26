#pragma once
#ifndef EVALUATOR_BASE_PAIR
#define EVALUATOR_BASE_PAIR

#include "cutils_math.h"
#include "cutils_3spn2.h"
#include "BasePair.h"
#define EPSILON 0.00001f
class BasePairEvaluator{
public:

    //evaluator.force(theta, basepairType, s, distSqrs, directors, invDotProd);
    inline __device__ float3 force(BasePairType basepairType, float phi, float thetas[2], float3 cVector, float3 dVector, float scValues[3], float invMagProds[3], float invLens[3], float invLenSqrs[3], float c12Mags[3], float c, float3 directors[2], int myIdxInBasePair) {
        float3 myForce;
        float dPhi = phi - basepairType.phi0;
        float dTheta1 = thetas[0] - basepairType.theta1;
        float dTheta2 = thetas[1] - basepairType.theta2;
        
        float phiFactor = 0.5f * (1.0f + cosf(dPhi));
        float ftor = 0.5f * sinf(dPhi);
        float energyMors
        float forceMors

        if(lengthSqr(directors[1]) < basepairType.sigma)
        {
            //Use a purely repulsive Morse potential
            energyMors = morsRepEnrgy(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma);
            forceMors = morsRepForc(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma);

            if (myIdxInBasePair == 1) {
                myForce.x -= forceRep * directors[1].x;
                myForce.y -= forceRep * directors[1].y;
                myForce.z -= forceRep * directors[1].z;
            }
            else if (myIdxInBasePair == 3) { 
                myForce.x += forceRep * directors[1].x;
                myForce.y += forceRep * directors[1].y;
                myForce.z += forceRep * directors[1].z;
            }
        }

        if ((dTheta1 >= -M_PI/(basepairType.k*2.0f)) && (dTheta1 <= M_PI/(basepairType.k*2.0f))) {
            if ((dTheta2 >= -M_PI/(basepairType.k*2.0f)) && (dTheta2 <= M_PI/(basepairType.k*2.0f))) {

                energyMors = morsAttrEnrgy(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma)
                forceMors = morsAttrForc(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma)

                if (myIdxInBasePair == 0) {
                    float derivofPoten = -ftor * invLens[0] * scValues[0] * energyMors;
                    myForce.x += derivofPoten * cVector.x;
                    myForce.y += derivofPoten * cVector.y;
                    myForce.z += derivofPoten * cVector.z;
                }
                else if (myIdxinBasePair == 1) {
                    //I may want to store the lengths instead of continuously calculating them
                    float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * energyMors;
                    float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * energyMors;
                    myForce.x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                    myForce.y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                    myForce.z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;
                }
                else if (myIdxInBasePair == 2) {
                    float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * energyMors;
                    float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * energyMors;
                    myForce.x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                    myForce.y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                    myForce.z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                }
                else if (myIdxInBasePair == 3) {
                    float derivofPoten = -ftor * invLens[2] * scValues[1] * energyMors;
                    myForce.x += derivofPoten * dVector.x;
                    myForce.y += derivofPoten * dVector.y;
                    myForce.z += derivofPoten * dVector.z;
                }
            }
            else if ((dTheta2 >= M_PI/(basepairType.k*2.0f)) && (dTheta2 <= M_PI/(basepairType.k)) || ((dTheta2 <= -M_PI/(basepairType.k*2.0f)) && (dTheta2 >= -M_PI/(basepairType.k))) {
                float cosine2 = cosf(basepairType.k*dTheta2);
                float sine2 = sinf(basepairType.k*dTheta2);
                float hbon_cosine = 1.0f - cosine2 * consine2;

                energyMors = morsAttrEnrgy(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma)
                forceMors = morsAttrForc(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma)

                float prefactor2 = 2.0f * basepairType.k * cosine2 * sine2 * 1.0f/ sqrtf(1.0f - thetas[1] * thetas[1]);

                if (myIdxInBasePair == 1) {
                    float3 fr1;
                    fr1.x = prefactor2 * (invLens[1] * (thetas[1] * directors[1].x * invLens[1] - directors[2].x * invLens[2])) * energyMors - hbon_cosine * directors[2].x * forceMors;
                    fr1.y = prefactor2 * (invLens[1] * (thetas[1] * directors[1].y * invLens[1] - directors[2].y * invLens[2])) * energyMors - hbon_cosine * directors[2].y * forceMors;
                    fr1.z = prefactor2 * (invLens[1] * (thetas[1] * directors[1].z * invLens[1] - directors[2].z * invLens[2])) * energyMors - hbon_cosine * directors[2].z * forceMors;
                    myForce.x += fr1.x * phiFactor;
                    myForce.y += fr1.y * phiFactor;
                    myForce.z += fr1.z * phiFactor;
                }
                else if (myIdxInBasePair == 2) {
                    float3 fr2;
                    fr2.x = prefactor2 * (invLens[2] * (directors[1].x * invLens[1] - thetas[1] * directors[2].x * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].x * invLens[1])) * energyMors + hbon_cosine * directors[1].x * forceMors;
                    fr2.y = prefactor2 * (invLens[2] * (directors[1].y * invLens[1] - thetas[1] * directors[2].y * invLens[2]) + invLens[1] * (directors[2].y * invLens[2] - thetas[1] * directors[1].y * invLens[1])) * energyMors + hbon_cosine * directors[1].y * forceMors;
                    fr2.z = prefactor2 * (invLens[2] * (directors[1].z * invLens[1] - thetas[1] * directors[2].z * invLens[2]) + invLens[1] * (directors[2].x * invLens[2] - thetas[1] * directors[1].z * invLens[1])) * energyMors + hbon_cosine * directors[1].z * forceMors;
                    myForce.x += fr2.x * phiFactor;
                    myForce.y += fr2.y * phiFactor;
                    myForce.z += fr2.z * phiFactor;
                }
                else if (myIdxInBasePair == 3) {
                    float3 fr3;
                    fr3.x = prefactor2 * (invLens[2] * (directors[2].x * invLens[2] * thetas[1] - directors[1].x * invLens[1])) * energMors;
                    fr3.y = prefactor2 * (invLens[2] * (directors[2].y * invLens[2] * thetas[1] - directors[1].y * invLens[1])) * energMors;
                    fr3.z = prefactor2 * (invLens[2] * (directors[2].z * invLens[2] * thetas[1] - directors[1].z * invLens[1])) * energMors;
                    myForce.x += fr2.x * phiFactor;
                    myForce.y += fr2.y * phiFactor;
                    myForce.z += fr2.z * phiFactor;
                }

                //This is mostly copied from Dan Hinckley
                if (myIdxInBasePair == 0) {
                    float derivofPoten = -ftor * invLens[0] * scValues[0] * energyMors;
                    myForce.x += derivofPoten * cVector.x;
                    myForce.y += derivofPoten * cVector.y;
                    myForce.z += derivofPoten * cVector.z;
                }
                else if (myIdxinBasePair == 1) {
                    //I may want to store the lengths instead of continuously calculating them
                    float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * energyMors;
                    float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * energyMors;
                    myForce.x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                    myForce.y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                    myForce.z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;
                }
                else if (myIdxInBasePair == 2) {
                    float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * energyMors;
                    float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * energyMors;
                    myForce.x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                    myForce.y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                    myForce.z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                }
                else if (myIdxInBasePair == 3) {
                    float derivofPoten = -ftor * invLens[2] * scValues[1] * energyMors;
                    myForce.x += derivofPoten * dVector.x;
                    myForce.y += derivofPoten * dVector.y;
                    myForce.z += derivofPoten * dVector.z;
                }

            }
            else if ((dTheta1 >= M_PI/(basepairType.k*2.0f)) && (dTheta1 <= M_PI/(basepairType.k)) || ((dTheta1 <= -M_PI/(basepairType.k*2.0f)) && (dTheta1 >= -M_PI/(basepairType.k))) {
                float cosine = cosf(basepairType.k*dTheta1);
                float sine = sinf(basepairType.k*dTheta1);
                float hbon_cosine = 1.0f - cosine * consine;

                if ((dTheta2 >= -M_PI/(basepairType.k*2.0f)) && (dTheta2 <= M_PI/(basepairType.k*2.0f))) {
                    energyMors = morsAttrEnrgy(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma)
                    forceMors = morsAttrForc(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma)

                    float prefactor = 2.0f * basepairType.k * cosine * sine * 1.0f/ sqrtf(1.0f - thetas[0] * thetas[0]);

                    if (myIdxInBasePair == 0) {
                        float3 fr1;
                        fr1.x = prefactor * (invLens[0] * (thetas[0] * directors[0].x * invLens[0] + directors[1].x * invLens[1])) * energyMors;
                        fr1.y = prefactor * (invLens[0] * (thetas[0] * directors[0].y * invLens[0] + directors[1].y * invLens[1])) * energyMors;
                        fr1.z = prefactor * (invLens[0] * (thetas[0] * directors[0].z * invLens[0] + directors[1].z * invLens[1])) * energyMors;
                        myForce.x += fr1.x * phiFactor;
                        myForce.y += fr1.y * phiFactor;
                        myForce.z += fr1.z * phiFactor;
                    }
                    else if (myIdxInBasePair == 1) {
                        float3 fr2;
                        fr2.x = prefactor * (invLens[0] * (-directors[1].x * invLens[1] - thetas[0] * directors[0].x * invLens[0]) + invLens[1] * (directors[0].x * invLens[0].x + thetas[0] * directors[1].x * invLens[1].x) * energyMors - hbon_cosine * directors[1].x * forceMors;
                        fr2.y = prefactor * (invLens[0] * (-directors[1].y * invLens[1] - thetas[0] * directors[0].y * invLens[0]) + invLens[1] * (directors[0].y * invLens[0].y + thetas[0] * directors[1].y * invLens[1].y) * energyMors - hbon_cosine * directors[1].y * forceMors;
                        fr2.z = prefactor * (invLens[0] * (-directors[1].z * invLens[1] - thetas[0] * directors[0].z * invLens[0]) + invLens[1] * (directors[0].z * invLens[0].z + thetas[0] * directors[1].z * invLens[1].z) * energyMors - hbon_cosine * directors[1].z * forceMors;
                        myForce.x += fr2.x * phiFactor;
                        myForce.y += fr2.y * phiFactor;
                        myForce.z += fr2.z * phiFactor;
                    }
                    else if (myIdxInBasePair == 3) {
                        float3 fr3;
                        fr3.x = prefactor * (invLens[1] * (-directors[1].x * invLens[1] * thetas[0] - directors[0].x * invLens[0])) * energMors + hbon_cosine * direcotrs[1].x * forceMors;
                        fr3.y = prefactor * (invLens[1] * (-directors[1].y * invLens[1] * thetas[0] - directors[0].y * invLens[0])) * energMors + hbon_cosine * direcotrs[1].y * forceMors;
                        fr3.z = prefactor * (invLens[1] * (-directors[1].z * invLens[1] * thetas[0] - directors[0].z * invLens[0])) * energMors + hbon_cosine * direcotrs[1].z * forceMors;
                        myForce.x += fr3.x * phiFactor;
                        myForce.y += fr3.y * phiFactor;
                        myForce.z += fr3.z * phiFactor;
                    }

                    //This is mostly copied from Dan Hinckley
                    if (myIdxInBasePair == 0) {
                        float derivofPoten = -ftor * invLens[0] * scValues[0] * energyMors;
                        myForce.x += derivofPoten * cVector.x;
                        myForce.y += derivofPoten * cVector.y;
                        myForce.z += derivofPoten * cVector.z;
                    }
                    else if (myIdxinBasePair == 1) {
                        //I may want to store the lengths instead of continuously calculating them
                        float frb1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[0]) * invLens[0] * invLens[1] * scValues[0] * energyMors;
                        float frb2 = ftor * c12Mags[1] * invLens[1] * scValues[1] * energyMors;
                        myForce.x += (frb1 * cVector.x + frb2 * dVector.x) - forceMors * phiFactor;
                        myForce.y += (frb1 * cVector.y + frb2 * dVector.y) - forceMors * phiFactor;
                        myForce.z += (frb1 * cVector.z + frb2 * dVector.z) - forceMors * phiFactor;
                    }
                    else if (myIdxInBasePair == 2) {
                        float frc1 = ftor * (length(directors[1]) - length(directors[2])*c12Mags[1]) * invLens[2] * invLens[1] * scValues[0] * energyMors;
                        float frc2 = ftor * c12Mags[0] * invLens[1] * scValues[0] * energyMors;
                        myForce.x += (frc1 * dVector.x + frc2 * cVector.x) - forceMors * phiFactor;
                        myForce.y += (frc1 * dVector.y + frc2 * cVector.y) - forceMors * phiFactor;
                        myForce.z += (frc1 * dVector.z + frc2 * cVector.z) - forceMors * phiFactor;
                    }
                    else if (myIdxInBasePair == 3) {
                        float derivofPoten = -ftor * invLens[2] * scValues[1] * energyMors;
                        myForce.x += derivofPoten * dVector.x;
                        myForce.y += derivofPoten * dVector.y;
                        myForce.z += derivofPoten * dVector.z;
                    }
                }

        }


    }
    inline __device__ void forcesAll(BasePairType basepairType, float theta, float s, float c, float distSqrs[2], float3 directors[2], float invDistProd, float3 forces[3]) {
        float cot = c / s;
        float dTheta = theta - basepairType.theta0;
        
        float r2 = sqrtf(distSqrs[1]);
        float fmorse = 0.0f;
        float emorse = 0.0f;
        float3 frepul = make_float3(0.0f, 0.0f, 0.0f);

        if(r2 < basepairType.sigma)
        {
            //Use a purely repulsive Morse poten
            float argu = basepairType.alpha * (r2 - basepairType.sigma);
            fmorse = -2.0f * basepairType.alpha * basepairType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
            forces[1] -= fmorse * directors[1];
            forces[2] += fmorse * directors[1];
        }

        if ((dTheta >= -M_PI/(basepairType.k*2.0f)) && (dTheta <= M_PI/(basepairType.k*2.0f))) {
            if (r2 >= basepairType.sigma) {
                float argu = basepairType.alpha * (r2 - basepairType.sigma);
                fmorse = -2.0f * basepairType.alpha * basepairType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = basepairType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - basepairType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -basepairType.epsi;
            }
            forces[1] -= fmorse * directors[1];
            forces[2] += fmorse * directors[1];
        }

        else if (((dTheta >= M_PI/(basepairType.k*2.0f)) && (dTheta <= M_PI/basepairType.k))
            || ((dTheta <= -M_PI/(basepairType.k*2.0f)) && (dTheta >= - M_PI/basepairType.k))) {
            //Calculate attractive-only Morse
            if (r2 >= basepairType.sigma) {
                float argu = basepairType.alpha * (r2 - basepairType.sigma);
                fmorse = -2.0f * basepairType.alpha * basepairType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = basepairType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - basepairType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -basepairType.epsi;
            }

            float cosine = cosf(basepairType.k * dTheta);
            float cosine_term  = 1.0f - cosine * cosine;
            float sine = sinf(basepairType.k * dTheta);
            float prefactor = 2.0f * basepairType.k * cosine * sine * 1.0/sqrtf(1.0-c*c);
            float a = -prefactor * emorse;

            float a11 = a*c/distSqrs[0];
            float a12 = -a*invDistProd;
            float a22 = a*c/distSqrs[1];

            float b11 = -a*c*cot/distSqrs[0];
            float b12 = a*cot*invDistProd;
            float b22 = -a*c*cot/distSqrs[1];

            forces[0] += (directors[0] * a11 + directors[1] * a12) + (directors[0] * b11 + directors[1] * b12);
            forces[1] -=
                    (directors[0] * a11 + directors[1] * a12) + (directors[0] * b11 + directors[1] * b12) 
                    +
                    (directors[1] * a22 + directors[0] * a12) + (directors[1] * b22 + directors[0] * b12) 
                    +
                    (cosine_term * directors[1] * fmorse)
                    ;

            forces[2] += 
                    (directors[1] * a22 + directors[0] * a12) + (directors[1] * b22 + directors[0] * b12)
                    +
                    (cosine_term * directors[1] * fmorse)
                    ;
        }


    }
    inline __device__ float energy(BasePairType basepairType, float theta, float3 directors[2]) {

        //Need sep distance for energy here (even if bonded interaction)
        float r2 = sqrtf(lengthSqr(directors[1]));
        float dTheta = theta - basepairType.theta0;
        float emorse = 0.0f;

        if(r2 < basepairType.sigma)
        {
            //Use a purely repulsive Morse poten
            float argu = basepairType.alpha * (r2 - basepairType.sigma);
            emorse += basepairType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu));
        }

        if ((dTheta >= -M_PI/(basepairType.k*2.0f)) && (dTheta <= M_PI/(basepairType.k*2.0f))) {
            if (r2 >= basepairType.sigma) {
                float argu = basepairType.alpha * (r2 - basepairType.sigma);
                emorse += basepairType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - basepairType.epsi;
            }
            else {
                emorse += -basepairType.epsi;
            }
        }

        //The cone of shame
        else if (((dTheta >= M_PI/(basepairType.k*2.0f)) && (dTheta <= M_PI/basepairType.k))
            || ((dTheta <= -M_PI/(basepairType.k*2.0f)) && (dTheta >= - M_PI/basepairType.k))) {
            //Calculate attractive-only Morse
            float estck = 0.0f;
            if (r2 >= basepairType.sigma) {
                float argu = basepairType.alpha * (r2 - basepairType.sigma);
                estck = basepairType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - basepairType.epsi;
            }
            else {
                estck += -basepairType.epsi;
            }
            float cosine = cosf(dTheta * basepairType.k);
            float cosine_term = 1.0f - cosine * cosine;
            estck *= cosine_term;
            emorse += estck;
        }
    
        return emorse;
    }
};

#endif

