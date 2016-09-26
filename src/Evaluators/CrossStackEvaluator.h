#pragma once
#ifndef EVALUATOR_CROSS_STACK_3SPN2
#define EVALUATOR_CROSS_STACK_3SPN2

#include "cutils_math.h"
#include "cutils_3spn2.h"
#include "CrossStack.h"
class CrossStackEvaluator3SPN2{
public:
    float alpha;
    float range;
    CrossStackEvaluator3SPN2() {}
    CrossStackEvaluator3SPN2(float alpha_, float range_) {
        alpha = alpha_;
        range = range_;
    }

    inline __device__ void forces(CrossStack3SPN2Type crossstackType, float thetas[3], float invMagProds[3], float invLens[3], float invLenSqrs[3], float c12Mags[3], float3 directors[5], float3 forces[6]) {

        //Difference in angle from reference value
        float dTheta1 = thetas[1] - crossstackType.theta1;
        float dTheta2 = thetas[2] - crossstackType.theta2;
        float dTheta3 = thetas[0] - crossstackType.theta3;

        //printf("dPhi is %f, dTheta1 is %f, dTheta2 is %f\n", dPhi * 180.0f/M_PI, dTheta1 * 180.0f/M_PI, dTheta2 * 180.0f/M_PI);
        //printf("phi is %f, theta1 is %f, theta2 is %f\n", phi * 180.0f/M_PI, thetas[0]* 180.0f/M_PI, thetas[1] * 180.0f/M_PI);
        //printf("c12Magsd %f\t%f\n", c12Mags[0], c12Mags[1]);
        //printf("range and alpha %f\t%f\n", range, alpha);
        
        //Store the "cone of interaction" values
        float invRange  = 1.0f/range;
        float cone = M_PI * invRange;
        float coneHalf = cone * 0.5f;
        float energyMors;
        float forceMors;

        //Cross-stack 1 interactions first (a-b-e)
        if ((dTheta3 >= -(coneHalf)) && (dTheta3 <= coneHalf)) {
            if ((dTheta1 >= -(coneHalf)) && (dTheta1 <= coneHalf)) {

                //Evaluate the Morse potentials
                forceMors = morsAttrForc(invLenSqrs[3], alpha, crossstackType.epsi, crossstackType.sigma1);

                //Only applied to sites b and e
                forces[4] -= forceMors * directors[3];
                forces[1] += forceMors * directors[3];
            }
            //The modulated form of the morse potential only on theta 1
            else if (((dTheta1 >= coneHalf) && (dTheta1 <= cone)) || ((dTheta1 <= -coneHalf) && (dTheta1 >= -cone))) {
                //Angle factors
                float cosine2 = cosf(range*dTheta1);
                float sine2 = sinf(range*dTheta1);
                float cross_term = 1.0f - cosine2 * cosine2;
                float prefactor2 = 2.0f * range * cosine2 * sine2 * 1.0f/ sqrtf(1.0f - c12Mags[1] * c12Mags[1]);

                //Evaluate the morse potential and force
                energyMors = morsAttrEnrgy(invLenSqrs[3], alpha, crossstackType.epsi, crossstackType.sigma1);
                forceMors = morsAttrForc(invLenSqrs[3], alpha, crossstackType.epsi, crossstackType.sigma1);

                //Site a forces
                float3 fr;
                fr = prefactor2 * (invLens[0] * (c12Mags[1] * directors[0] * invLens[0] - directors[3] * invLens[3])) * energyMors;
                forces[0] += fr;

                //Site b forces;
                fr = prefactor2 * (invLens[0] * (directors[3] * invLens[3] - directors[0] * invLens[0] * c12Mags[1]) + invLens[3] * (directors[0] * invLens[0] - c12Mags[1] * directors[3] * invLens[3] )) * energyMors + cross_term * directors[3] * forceMors; 
                forces[1] += fr;

                //Site e forces
                fr = prefactor2 * (invLens[3] * (directors[3] * invLens[3] * c12Mags[1] - directors[0] * invLens[0])) * energyMors - cross_term * directors[3] * forceMors; 
                forces[4] += fr;
            }
            else {
                //Nothing else happens here. I've cheated you. I've cheated ALL of you
            }
        }

        //Modulated theta3
        else if (((dTheta3 >= coneHalf) && (dTheta3 <= cone)) || ((dTheta3 <= -coneHalf) && (dTheta3 >= -cone))) {

            //More angle parameters for theta3
            float cosine2 = cosf(range*dTheta3);
            float sine2 = sinf(range*dTheta3);
            float basepair_term = 1.0f - cosine2 * cosine2;
            float prefactor = 2.0f * range * cosine2 * sine2 * 1.0f/ sqrtf(1.0f - c12Mags[0] * c12Mags[0]);
    
            //Full potential for theta1 and modulated theta3
            if ((dTheta1 >= -coneHalf) && (dTheta1 <= coneHalf)) {

                //Calculate Morse potential and forces again
                energyMors = morsAttrEnrgy(invLenSqrs[3], alpha, crossstackType.epsi, crossstackType.sigma1);
                forceMors = morsAttrForc(invLenSqrs[3], alpha, crossstackType.epsi, crossstackType.sigma1);

                //Site a forces
                float3 fr;
                fr = prefactor * invLens[0] * (c12Mags[0] * directors[0] * invLens[0] - directors[2] * invLens[2]) * energyMors;
                forces[0] += fr;

                //Site b forces
                fr = prefactor * invLens[0] * (directors[2] * invLens[2] - directors[0] * invLens[0] * c12Mags[0]) * energyMors + basepair_term * directors[3] * forceMors;
                forces[1] += fr;

                //Site c forces
                fr = prefactor * invLens[2] * (directors[2] * invLens[2] * c12Mags[0] - directors[0] * invLens[0]) * energyMors;
                forces[2] += fr;

                //Site d forces
                fr = prefactor * invLens[2] * (directors[0] * invLens[0] - directors[2] * invLens[2] * c12Mags[0]) * energyMors;
                forces[3] += fr;

                //Site e forces
                fr = -basepair_term * directors[3] * forceMors;
                forces[4] += fr;
            }
            //Both potentials are modulated
            else if (((dTheta1 >= (coneHalf)) && (dTheta1 <= cone)) || (((dTheta1 <= (-coneHalf)) && (dTheta1 >= -cone)))) {
                float cosine2 = cosf(range*dTheta1);
                float sine2 = sinf(range*dTheta1);
                float cross_term = 1.0f - cosine2 * cosine2;
                float prefactor2 = 2.0f * range * cosine2 * sine2 * 1.0f / sqrtf(1.0f - c12Mags[1] * c12Mags[1]);

                //More evaluations of the Morse potential
                energyMors = morsAttrEnrgy(invLenSqrs[3], alpha, crossstackType.epsi, crossstackType.sigma1);
                forceMors = morsAttrForc(invLenSqrs[3], alpha, crossstackType.epsi, crossstackType.sigma1);

                //Site a
                float3 fr;
                fr = prefactor * invLens[0] * (directors[0] * invLens[0] * c12Mags[0] - directors[2] * invLens[2]) * cross_term * energyMors + prefactor2 * (invLens[0] * (directors[0] * invLens[0] * c12Mags[1] - directors[3] * invLens[3])) * basepair_term * energyMors;
                forces[0] += fr;

                //Site b
                fr = prefactor * invLens[0] * (directors[2] * invLens[2] - directors[0] * invLens[0] * c12Mags[0]) * cross_term * energyMors + prefactor2 * (invLens[0] * (directors[3] * invLens[3] - c12Mags[1] * directors[0] * invLens[0]) + invLens[3] * (directors[0] * invLens[0] - c12Mags[1] * directors[3] * invLens[3])) * basepair_term * energyMors + basepair_term * cross_term * directors[3] * forceMors;
                forces[1] += fr; 
          
                //Site c
                fr = prefactor * invLens[2] * (directors[2] * invLens[2] * c12Mags[0] - directors[0] * invLens[0]) * cross_term * energyMors;
                forces[2] += fr;

                //Site d
                fr = prefactor * invLens[2] * (directors[0] * invLens[0] - directors[2] * invLens[2] * c12Mags[0]) * cross_term * energyMors;
                forces[3] += fr;

                //Site e
                fr = prefactor2 * (invLens[3] * (directors[3] * invLens[3] * c12Mags[1] - directors[0] * invLens[0])) * basepair_term * energyMors - basepair_term * cross_term * directors[3] * forceMors;
                forces[4] += fr;

            }
            else {
                //More nothingness
            }
        }
        else {
            //If an angle is here, ya goofed
        }
 
        //Cross-stack 2 interactions (c-d-f)
        //In a later version this may be switched to a second kernel call so that the nanoparticle sims can be run
        if ((dTheta3 >= -(coneHalf)) && (dTheta3 <= coneHalf)) {
            if ((dTheta2 >= -(coneHalf)) && (dTheta2 <= coneHalf)) {

                //Evaluate the Morse potentials
                forceMors = morsAttrForc(invLenSqrs[4], alpha, crossstackType.epsi, crossstackType.sigma2);

                //Only applied to sites d and f
                forces[5] -= forceMors * directors[4];
                forces[3] += forceMors * directors[4];
            }
            //The modulated form of the morse potential only on theta 2
            else if (((dTheta2 >= coneHalf) && (dTheta2 <= cone)) || ((dTheta2 <= -coneHalf) && (dTheta2 >= -cone))) {
                //Angle factors
                float cosine2 = cosf(range*dTheta2);
                float sine2 = sinf(range*dTheta2);
                float cross_term = 1.0f - cosine2 * cosine2;
                float prefactor2 = 2.0f * range * cosine2 * sine2 * 1.0f/ sqrtf(1.0f - c12Mags[2] * c12Mags[2]);

                //Evaluate the morse potential and force
                energyMors = morsAttrEnrgy(invLenSqrs[4], alpha, crossstackType.epsi, crossstackType.sigma2);
                forceMors = morsAttrForc(invLenSqrs[4], alpha, crossstackType.epsi, crossstackType.sigma2);

                //Site c forces
                float3 fr1;
                fr1 = prefactor2 * (invLens[2] * (c12Mags[2] * directors[2] * invLens[2] - directors[4] * invLens[4])) * energyMors;
                forces[2] += fr1;

                //Site d forces
                float3 fr2;
                fr2 = prefactor2 * (invLens[2] * (directors[4] * invLens[4] - directors[2] * invLens[2] * c12Mags[2]) + invLens[4] * (directors[2] * invLens[2] - c12Mags[2] * directors[4] * invLens[4] )) * energyMors + cross_term * directors[4] * forceMors; 
                forces[3] += fr2;

                //Site f forces
                float3 fr3;
                fr3 = prefactor2 * (invLens[4] * (directors[4] * invLens[4] * c12Mags[2] - directors[2] * invLens[2])) * energyMors - cross_term * directors[4] * forceMors; 
                forces[5] += fr3;
            }
            else {
                //Nothing else happens here. I've cheated you. I've cheated ALL of you
            }
        }

        //Modulated theta3
        else if (((dTheta3 >= coneHalf) && (dTheta3 <= cone)) || ((dTheta3 <= -coneHalf) && (dTheta3 >= -cone))) {

            //More angle parameters for theta3
            float cosine2 = cosf(range*dTheta3);
            float sine2 = sinf(range*dTheta3);
            float basepair_term = 1.0f - cosine2 * cosine2;
            float prefactor = 2.0f * range * cosine2 * sine2 * 1.0f/ sqrtf(1.0f - c12Mags[0] * c12Mags[0]);
    
            //Full potential for theta1 and modulated theta3
            if ((dTheta2 >= -coneHalf) && (dTheta2 <= coneHalf)) {

                //Calculate Morse potential and forces again
                energyMors = morsAttrEnrgy(invLenSqrs[4], alpha, crossstackType.epsi, crossstackType.sigma2);
                forceMors = morsAttrForc(invLenSqrs[4], alpha, crossstackType.epsi, crossstackType.sigma2);

                //Site a forces
                float3 fr;
                fr = prefactor * invLens[0] * (c12Mags[0] * directors[0] * invLens[0] - directors[2] * invLens[2]) * energyMors;
                forces[0] += fr;

                //Site b forces
                fr = prefactor * invLens[0] * (directors[2] * invLens[2] - directors[0] * invLens[0] * c12Mags[0]) * energyMors;
                forces[1] += fr;

                //Site c forces
                fr = prefactor * invLens[2] * (directors[2] * invLens[2] * c12Mags[0] - directors[0] * invLens[0]) * energyMors;
                forces[2] += fr;

                //Site d forces
                fr = prefactor * invLens[2] * (directors[0] * invLens[0] - directors[2] * invLens[2] * c12Mags[0]) * energyMors + basepair_term * directors[4] * forceMors;
                forces[3] += fr;

                //Site 5 forces
                fr = -basepair_term * directors[4] * forceMors;
                forces[5] += fr;
            }
            //Both potentials are modulated
            else if (((dTheta2 >= (coneHalf)) && (dTheta2 <= cone)) || (((dTheta2 <= (-coneHalf)) && (dTheta2 >= -cone)))) {
                float cosine2 = cosf(range*dTheta2);
                float sine2 = sinf(range*dTheta2);
                float cross_term = 1.0f - cosine2 * cosine2;
                float prefactor2 = 2.0f * range * cosine2 * sine2 * 1.0f / sqrtf(1.0f - c12Mags[2] * c12Mags[2]);

                //More evaluations of the Morse potential
                energyMors = morsAttrEnrgy(invLenSqrs[4], alpha, crossstackType.epsi, crossstackType.sigma2);
                forceMors = morsAttrForc(invLenSqrs[4], alpha, crossstackType.epsi, crossstackType.sigma2);

                //Site a
                float3 fr;
                fr = prefactor * invLens[0] * (directors[0] * invLens[0] * c12Mags[0] - directors[2] * invLens[2]) * cross_term * energyMors;
                forces[0] += fr;

                //Site b
                fr = prefactor * invLens[0] * (directors[2] * invLens[2] - directors[0] * invLens[0] * c12Mags[0]) * cross_term * energyMors;
                forces[1] += fr; 
          
                //Site c
                fr = prefactor * invLens[2] * (directors[2] * invLens[2] * c12Mags[0] - directors[0] * invLens[0]) * cross_term * energyMors + prefactor2 * (invLens[2] * (directors[2] * invLens[2] * c12Mags[2] - directors[4] * invLens[4])) * basepair_term * energyMors;
                forces[2] += fr;

                //Site d
                fr = prefactor * invLens[2] * (directors[0] * invLens[0] - directors[2] * invLens[2] * c12Mags[0]) * cross_term * energyMors + prefactor2 * (invLens[2] * (directors[4] * invLens[4] - c12Mags[2] * directors[2] * invLens[2]) + invLens[4] * (directors[2] * invLens[2] - c12Mags[2] * directors[4] * invLens[4])) * basepair_term * energyMors + basepair_term * cross_term * directors[4] * forceMors;
                forces[3] += fr;

                //Site f
                fr = prefactor2 * (invLens[4] * (directors[4] * invLens[4] * c12Mags[2] - directors[2] * invLens[2])) * basepair_term * energyMors - basepair_term * cross_term * directors[4] * forceMors;
                forces[4] += fr;

            }
            else {
                //More nothingness
            }
        }
        else {
            //If an angle is here, ya goofed
        }

    }


    inline __device__ float energy(CrossStack3SPN2Type crossstackType, float thetas[3], float invLens[5], float invLenSqrs[5], float3 directors[5]) {
        /*float dPhi = phi - crossstackType.phi0;
        float dTheta1 = thetas[0] - crossstackType.theta1;
        float dTheta2 = thetas[1] - crossstackType.theta2;
        float cone = 2.0f * range;
        
        float phiFactor = 0.5f * (1.0f + cosf(dPhi));
        float ftor = 0.5f * sinf(dPhi);
        float energyMors;
        float eCrossStack = 0.0f;

        if(lengthSqr(directors[1]) < crossstackType.sigma) {
            //Use a purely repulsive Morse potential
            energyMors = morsRepEnrgy(invLens[1], alpha, crossstackType.epsi, crossstackType.sigma);
            eCrossStack += energyMors * phiFactor;

        }

        if ((dTheta1 >= -M_PI/(cone)) && (dTheta1 <= M_PI/(cone))) {
            if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {

                energyMors = morsAttrEnrgy(invLens[1], alpha, crossstackType.epsi, crossstackType.sigma);
                eCrossStack += phiFactor * energyMors;

            }
            else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || ((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range)))) {
                float cosine2 = cosf(range*dTheta2);
                float sine2 = sinf(range*dTheta2);
                float hbon_cosine = 1.0f - cosine2 * cosine2;

                energyMors = morsAttrEnrgy(invLens[1], alpha, crossstackType.epsi, crossstackType.sigma);
                eCrossStack += phiFactor * hbon_cosine * energyMors;

            }
            else if (((dTheta1 >= M_PI/(cone)) && (dTheta1 <= M_PI/(range))) || ((dTheta1 <= -M_PI/(cone)) && (dTheta1 >= -M_PI/(range)))) {
                float cosine = cosf(range*dTheta1);
                float sine = sinf(range*dTheta1);

                if ((dTheta2 >= -M_PI/(cone)) && (dTheta2 <= M_PI/(cone))) {
                    energyMors = morsAttrEnrgy(invLens[1], alpha, crossstackType.epsi, crossstackType.sigma);
                    eCrossStack += phiFactor * energyMors;

                }
                else if (((dTheta2 >= M_PI/(cone)) && (dTheta2 <= M_PI/(range))) || ((dTheta2 <= -M_PI/(cone)) && (dTheta2 >= -M_PI/(range)))) {
                    float cosine2 = cosf(range*dTheta2);
                    float sine2 = sinf(range*dTheta2);
                    float hbon_cosine2 = 1.0f - cosine2 * cosine2;

                    energyMors = morsAttrEnrgy(invLens[1], alpha, crossstackType.epsi, crossstackType.sigma);
                    eCrossStack += phiFactor * hbon_cosine2 * energyMors;


                }
            }
        }
        return eCrossStack;*/
    }
};

#endif

