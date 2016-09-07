#pragma once
#ifndef EVALUATOR_ANGLE_COSINE_DELTA_
#define EVALUATOR_ANGLE_COSINE_DELTA

#include "cutils_math.h"
#include "Angle.h"
#define EPSILON 0.00001f
class AngleEvaluatorBaseStacking{
public:

    //evaluator.force(theta, angleType, s, distSqrs, directors, invDotProd);
    inline __device__ float3 force(AngleBaseStackingType angleType, float theta, float s, float c, float distSqrs[2], float3 directors[2], float invDistProd, int myIdxInAngle) {
        float cot = c / s;
        float dTheta = theta - angleType.theta0;
        
        float r2 = sqrtf(distSqrs[1]);
        float fmorse;
        float emorse;
        float3 frepul = make_float3(0.0f, 0.0f, 0.0f);

        if(r2 < angleType.sigma)
        {
            //Use a purely repulsive Morse potential
            float argu = angleType.alpha * (r2 - angleType.sigma);
            fmorse = -2.0f * angleType.alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
            if (myIdxInAngle == 1) {
                frepul -= fmorse * directors[1];
            }
            else if (myIdxInAngle == 2) {
                frepul += fmorse * directors[1];
            }
        }

        if ((dTheta >= -M_PI/(angleType.k*2.0f)) && (dTheta <= M_PI/(angleType.k*2.0f))) {
            if (r2 >= angleType.sigma) {
                float argu = angleType.alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * angleType.alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -angleType.epsi;
            }
            if (myIdxInAngle == 1) {
                frepul -= fmorse * directors[1];
            }
            else if (myIdxInAngle == 2) {
                frepul += fmorse * directors[1];
            }
        }

        else if (((dTheta >= M_PI/(angleType.k*2.0f)) && (dTheta <= M_PI/angleType.k))
            || ((dTheta <= -M_PI/(angleType.k*2.0f)) && (dTheta >= - M_PI/angleType.k))) {
            //Calculate attractive-only Morse
            if (r2 >= angleType.sigma) {
                float argu = angleType.alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * angleType.alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -angleType.epsi;
            }

            float cosine = cosf(angleType.k * dTheta);
            float cosine_term  = 1.0f - cosine * cosine;
            float sine = sinf(angleType.k * dTheta);
            float prefactor = 2.0f * angleType.k * cosine * sine * 1.0/sqrtf(1.0-c*c);
            float a = -prefactor * emorse;

            float a11 = a*c/distSqrs[0];
            float a12 = -a*invDistProd;
            float a22 = a*c/distSqrs[1];

            float b11 = -a*c*cot/distSqrs[0];
            float b12 = a*cot*invDistProd;
            float b22 = -a*c*cot/distSqrs[1];

            if (myIdxInAngle==0) {
                frepul += (directors[0] * a11 + directors[1] * a12) + (directors[0] * b11 + directors[1] * b12);
            } else if (myIdxInAngle==1) {
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
    inline __device__ void forcesAll(AngleBaseStackingType angleType, float theta, float s, float c, float distSqrs[2], float3 directors[2], float invDistProd, float3 forces[3]) {
        float cot = c / s;
        float dTheta = theta - angleType.theta0;
        
        float r2 = sqrtf(distSqrs[1]);
        float fmorse = 0.0f;
        float emorse = 0.0f;
        float3 frepul = make_float3(0.0f, 0.0f, 0.0f);

        if(r2 < angleType.sigma)
        {
            //Use a purely repulsive Morse poten
            float argu = angleType.alpha * (r2 - angleType.sigma);
            fmorse = -2.0f * angleType.alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
            forces[1] -= fmorse * directors[1];
            forces[2] += fmorse * directors[1];
        }

        if ((dTheta >= -M_PI/(angleType.k*2.0f)) && (dTheta <= M_PI/(angleType.k*2.0f))) {
            if (r2 >= angleType.sigma) {
                float argu = angleType.alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * angleType.alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -angleType.epsi;
            }
            forces[1] -= fmorse * directors[1];
            forces[2] += fmorse * directors[1];
        }

        else if (((dTheta >= M_PI/(angleType.k*2.0f)) && (dTheta <= M_PI/angleType.k))
            || ((dTheta <= -M_PI/(angleType.k*2.0f)) && (dTheta >= - M_PI/angleType.k))) {
            //Calculate attractive-only Morse
            if (r2 >= angleType.sigma) {
                float argu = angleType.alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * angleType.alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -angleType.epsi;
            }

            float cosine = cosf(angleType.k * dTheta);
            float cosine_term  = 1.0f - cosine * cosine;
            float sine = sinf(angleType.k * dTheta);
            float prefactor = 2.0f * angleType.k * cosine * sine * 1.0/sqrtf(1.0-c*c);
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
    inline __device__ float energy(AngleBaseStackingType angleType, float theta, float3 directors[2]) {

        //Need sep distance for energy here (even if bonded interaction)
        float r2 = sqrtf(lengthSqr(directors[1]));
        float dTheta = theta - angleType.theta0;
        float emorse = 0.0f;

        if(r2 < angleType.sigma)
        {
            //Use a purely repulsive Morse poten
            float argu = angleType.alpha * (r2 - angleType.sigma);
            emorse += angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu));
        }

        if ((dTheta >= -M_PI/(angleType.k*2.0f)) && (dTheta <= M_PI/(angleType.k*2.0f))) {
            if (r2 >= angleType.sigma) {
                float argu = angleType.alpha * (r2 - angleType.sigma);
                emorse += angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                emorse += -angleType.epsi;
            }
        }

        //The cone of shame
        else if (((dTheta >= M_PI/(angleType.k*2.0f)) && (dTheta <= M_PI/angleType.k))
            || ((dTheta <= -M_PI/(angleType.k*2.0f)) && (dTheta >= - M_PI/angleType.k))) {
            //Calculate attractive-only Morse
            float estck = 0.0f;
            if (r2 >= angleType.sigma) {
                float argu = angleType.alpha * (r2 - angleType.sigma);
                estck = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                estck += -angleType.epsi;
            }
            float cosine = cosf(dTheta * angleType.k);
            float cosine_term = 1.0f - cosine * cosine;
            estck *= cosine_term;
            emorse += estck;
        }
    
        return emorse;
    }
};

#endif

