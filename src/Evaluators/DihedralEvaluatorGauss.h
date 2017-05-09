#pragma once
#ifndef EVALUATOR_DIHEDRAL_GAUSS
#define EVALUATOR_DIHEDRAL_GAUSS

#include "cutils_math.h"
#include "Dihedral.h"
class DihedralEvaluatorGauss {
    public:
        //dihedralType, phi, c, scValues, invLenSqrs, c12Mags, c0,

                //float3 myForce = evaluator.force(dihedralType, phi, c, scValues, invLenSqrs, c12Mags, c0, c, invMagProds, c12Mags, invLens, directors, myIdxInDihedral);
        inline __device__ float dPotential(DihedralGaussType dihedralType, float phi)  {
            float diffPhi = phi-dihedralType.phi0;
            //float diffPhi = 0.0;

            //Make sure difference follows periodic behavior correctly
            if (diffPhi < -M_PI) diffPhi += 2.0f * M_PI;
            else if (diffPhi > M_PI) diffPhi -= 2.0f * M_PI; 
            //printf("Phi is %f\tPhiref is %f\n", phi*180.0f/M_PI, dihedralType.phi0*180.0f/M_PI);
            float invSigma = 1.0f / (2.0f * dihedralType.sigma * dihedralType.sigma);
            
            //printf("Force is %f\n", - dihedralType.k0 * diffPhi * expf( - diffPhi * diffPhi * invSigma) / (dihedralType.sigma * dihedralType.sigma));
            return  (
                      dihedralType.k0 
                    * diffPhi * expf( - diffPhi * diffPhi * invSigma) * 2.0f * invSigma
                    )
                ;
        }

        inline __device__ float potential(DihedralGaussType dihedralType, float phi) {

            return -dihedralType.k0 * (
                               expf(-(phi - dihedralType.phi0) * (phi - dihedralType.phi0)
                               / (2.0f * dihedralType.sigma * dihedralType.sigma))
                              )
                ;

        }
};

#endif


