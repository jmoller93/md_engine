#pragma once
#ifndef EVALUATOR_ANGLE_BASE_STACKING
#define EVALUATOR_ANGLE_BASE_STACKING

#include "cutils_math.h"
#include "Angle.h"
#define EPSILON 0.00001f
class AngleEvaluatorBaseStacking{
public:
    float range;
    float alpha;
    AngleEvaluatorBaseStacking() {};
    AngleEvaluatorBaseStacking(float alpha_, float range_) {
        range = range_;
        alpha = alpha_;
    }

    //evaluator.force(theta, angleType, s, distSqrs, directors, invDotProd);
    inline __device__ float3 force(AngleBaseStackingType angleType, float theta, float s, float c, float distSqrs[2], float3 directors[2], float invDistProd, int myIdxInAngle) {
        float dTheta = theta - angleType.theta0;
        
        float r2 = sqrtf(distSqrs[1]);
        float fmorse = 0;
        float emorse = 0;
        float invRange = 1.0f / range;
        float cone = M_PI * invRange;
        float coneHalf = cone * 0.5f;
        float3 frepul = make_float3(0.0f, 0.0f, 0.0f);

        if(r2 < angleType.sigma)
        {
            //Use a purely repulsive Morse potential
            float argu = alpha * (r2 - angleType.sigma);
            if (r2 > 0) {
                fmorse = -2.0f * alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
            }
            //printf("Epsi is %f, R2 is %f, sigma is %f, fmorse is %f\n", angleType.epsi, r2, angleType.sigma, fmorse);
            if (myIdxInAngle == 1) {
                frepul -= fmorse * directors[1];
            }
            else if (myIdxInAngle == 2) {
                frepul += fmorse * directors[1];
            }
        }

        if ((dTheta >= -coneHalf) && (dTheta <= coneHalf)) {
            if (r2 >= angleType.sigma) {
                float argu = alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
                //printf("Epsi is %f, R2 is %f, sigma is %f, fmorse is %f\n", angleType.epsi, r2, angleType.sigma, fmorse);
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

        else if (((dTheta >= coneHalf) && (dTheta <= cone))
            || ((dTheta <= -coneHalf) && (dTheta >= - cone))) {
            //Calculate attractive-only Morse
            if (r2 >= angleType.sigma) {
                float argu = alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
                //printf("Epsi is %f, R2 is %f, sigma is %f, fmorse is %f\n", angleType.epsi, r2, angleType.sigma, fmorse);
            }
            else {
                fmorse = 0.0f;
                emorse = -angleType.epsi;
            }

            float cosine = cosf(range * dTheta);
            float cosine_term  = 1.0f - cosine * cosine;
            float sine = sinf(range * dTheta);
            float prefactor = 2.0f * range * cosine * sine * 1.0f/sqrtf(1.0f-c*c);
            float a = -prefactor * emorse;

            float a11 = a*c/distSqrs[0];
            float a12 = -a*invDistProd;
            float a22 = a*c/distSqrs[1];

            if (myIdxInAngle==0) {
                frepul.x += (directors[0].x * a11 + directors[1].x * a12);
                frepul.y += (directors[0].y * a11 + directors[1].y * a12);
                frepul.z += (directors[0].z * a11 + directors[1].z * a12);
            } else if (myIdxInAngle==1) {
                frepul.x -= (directors[0].x * a11 + directors[1].x * a12) + (directors[1].x * a22 + directors[0].x * a12 + cosine_term * directors[1].x * fmorse);
                frepul.y -= (directors[0].y * a11 + directors[1].y * a12) + (directors[1].y * a22 + directors[0].y * a12 + cosine_term * directors[1].y * fmorse); 
                frepul.z -= (directors[0].z * a11 + directors[1].z * a12) + (directors[1].z * a22 + directors[0].z * a12 + cosine_term * directors[1].z * fmorse);
            } else {
                frepul.x += directors[1].x * a22 + directors[0].x * a12 + (cosine_term * directors[1].x * fmorse);
                frepul.y += directors[1].y * a22 + directors[0].y * a12 + (cosine_term * directors[1].y * fmorse);
                frepul.z += directors[1].z * a22 + directors[0].z * a12 + (cosine_term * directors[1].z * fmorse);
            }
        }

        return frepul;

    }
    inline __device__ void forcesAll(AngleBaseStackingType angleType, float theta, float s, float c, float distSqrs[2], float3 directors[2], float invDistProd, float3 forces[3]) {
        float dTheta = theta - angleType.theta0;
        
        float r2 = sqrtf(distSqrs[1]);
        float fmorse = 0.0f;
        float emorse = 0.0f;
        float3 frepul = make_float3(0.0f, 0.0f, 0.0f);

        if(r2 < angleType.sigma)
        {
            //Use a purely repulsive Morse poten
            float argu = alpha * (r2 - angleType.sigma);
            fmorse = -2.0f * alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
            forces[1] -= fmorse * directors[1];
            forces[2] += fmorse * directors[1];
        }

        if ((dTheta >= -M_PI/(range*2.0f)) && (dTheta <= M_PI/(range*2.0f))) {
            if (r2 >= angleType.sigma) {
                float argu = alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -angleType.epsi;
            }
            forces[1] -= fmorse * directors[1];
            forces[2] += fmorse * directors[1];
        }

        else if (((dTheta >= M_PI/(range*2.0f)) && (dTheta <= M_PI/range))
            || ((dTheta <= -M_PI/(range*2.0f)) && (dTheta >= - M_PI/range))) {
            //Calculate attractive-only Morse
            if (r2 >= angleType.sigma) {
                float argu = alpha * (r2 - angleType.sigma);
                fmorse = -2.0f * alpha * angleType.epsi * expf(-argu) * (1.0f - expf(-argu)) / r2;
                emorse = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                fmorse = 0.0f;
                emorse = -angleType.epsi;
            }

            float cosine = cosf(range * dTheta);
            float cosine_term  = 1.0f - cosine * cosine;
            float sine = sinf(range * dTheta);
            float prefactor = 2.0f * range * cosine * sine * 1.0/sqrtf(1.0-c*c);
            float a = -prefactor * emorse;

            float a11 = a*c/distSqrs[0];
            float a12 = -a*invDistProd;
            float a22 = a*c/distSqrs[1];

            float3 f1 = directors[0] * a11 + directors[1] * a12;
            float3 f3 = directors[1] * a22 + directors[0] * a12;

            forces[0] += f1;
            forces[1] -= (f1 + f3 + (cosine_term * directors[1] * fmorse));
            forces[2] += (f3 + (cosine_term * directors[1] * fmorse));
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
            float argu = alpha * (r2 - angleType.sigma);
            emorse += angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu));
        }

        if ((dTheta >= -M_PI/(range*2.0f)) && (dTheta <= M_PI/(range*2.0f))) {
            if (r2 >= angleType.sigma) {
                float argu = alpha * (r2 - angleType.sigma);
                emorse += angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                emorse += -angleType.epsi;
            }
        }

        //The cone of shame
        else if (((dTheta >= M_PI/(range*2.0f)) && (dTheta <= M_PI/range))
            || ((dTheta <= -M_PI/(range*2.0f)) && (dTheta >= - M_PI/range))) {
            //Calculate attractive-only Morse
            float estck = 0.0f;
            if (r2 >= angleType.sigma) {
                float argu = alpha * (r2 - angleType.sigma);
                estck = angleType.epsi * (1.0f - expf(-argu)) * (1.0f - expf(-argu)) - angleType.epsi;
            }
            else {
                estck += -angleType.epsi;
            }
            float cosine = cosf(dTheta * range);
            float cosine_term = 1.0f - cosine * cosine;
            estck *= cosine_term;
            emorse += estck;
        }
    
        return emorse;
    }
};

#endif

