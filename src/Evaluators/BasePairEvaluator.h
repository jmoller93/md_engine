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
    inline __device__ float3 force(BasePairType basepairType, float phi, float thetas[2], float3 cVector, float3 dVector, float scValues[3], float invMagProds[3], float invLens[3], float c12Mags[3], float c, float3 directors[2], int myIdxInBasePair) {
        float3 myForce;
        float dPhi = phi - basepairType.phi0;
        float dTheta1 = thetas[0] - basepairType.theta1;
        float dTheta2 = thetas[1] - basepairType.theta2;
        
        float phiFactor = 0.5f * (1.0f + cosf(dPhi));
        float ftor = 0.5f * sinf(dPhi);

        if(length(directors[1]) < basepairType.sigma)
        {
            //Use a purely repulsive Morse potential
            float energyRep; //Repulsive energy contribution (comes into play later)
            float forceRep; //Repulsive force contribution
            energyRep = mors_rp_ener(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma);
            forceRep = mors_rp_frce(invLens[1], basepairType.alpha, basepairType.epsi, basepairType.sigma);

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
            if (myIdxInBasePair == 1) {
                frepul -= fmorse * directors[1];
            }
            else if (myIdxInBasePair == 2) {
                frepul += fmorse * directors[1];
            }
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

            if (myIdxInBasePair==0) {
                frepul += (directors[0] * a11 + directors[1] * a12) + (directors[0] * b11 + directors[1] * b12);
            } else if (myIdxInBasePair==1) {
                frepul -=
                    (directors[0] * a11 + directors[1] * a12) + (directors[0] * b11 + directors[1] * b12) 
                    +
                    (directors[1] * a22 + directors[0] * a12) + (directors[1] * b22 + directors[0] * b12) 
                    +
                    (cosine_term * directors[1] * fmorse)
                    ;

            } else {
                frepul += 
                    (directors[1] * a22 + directors[0] * a12) + (directors[1] * b22 + directors[0] * b12)
                    +
                    (cosine_term * directors[1] * fmorse)
                    ;
            }
        }

        return frepul;

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

