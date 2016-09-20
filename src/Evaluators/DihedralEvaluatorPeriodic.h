#pragma once
#ifndef EVALUATOR_DIHEDRAL_PERIODIC
#define EVALUATOR_DIHEDRAL_PERIODIC

#include "cutils_math.h"
#include "Dihedral.h"

//This is the dihedral potential from the Golike interactions developed by Clementi et al.
class DihedralEvaluatorPeriodic {
    public:

                //float3 myForce = evaluator.force(dihedralType, phi, c, scValues, invLenSqrs, c12Mags, c0, c, invMagProds, c12Mags, invLens, directors, myIdxInDihedral);
        inline __device__ float dPotential(DihedralPeriodicType dihedralType, float phi) {
            //The difference between system dihedral and reference (dtor in 3SPN.2 Lammps)
            float diffPhi = phi-dihedralType.phiRef;
    
            //Make sure difference follows periodic behavior correctly
            if (diffPhi < -M_PI) diffPhi += 2.0f * M_PI;
            else if (diffPhi > M_PI) diffPhi -= 2.0f * M_PI;

            
            return (
                    dihedralType.coefs[0] * sinf(diffPhi)
                    + 3.0f * dihedralType.coefs[1] * sinf(3.0f*diffPhi)
                    )
                ;

            /*//Evaluate derivatives
            float fg = dot(directors[0], -directors[1]); 
            float hg = dot(directors[2], -directors[1]);

            //Evaluate cross products
            float3 cVector;
            cVector.x = -directors[0].y*directors[1].z + directors[0].z*directors[1].y;
            cVector.y = -directors[0].z*directors[1].x + directors[0].x*directors[1].z;
            cVector.z = -directors[0].x*directors[1].y + directors[0].y*directors[1].x;
            float cVectorLen = length(cVector);
            float invCLenSqr = 1.0/cVectorLen;
            invCLenSqr *= invCLenSqr;
            float fga =  fg * invLens[1] * invCLenSqr;
    
            float3 dVector;
            dVector.x = -directors[2].y*directors[1].z + directors[2].z*directors[1].y;
            dVector.y = -directors[2].z*directors[1].x + directors[2].x*directors[1].z;
            dVector.z = -directors[2].x*directors[1].y + directors[2].y*directors[1].x;
            float dVectorLen = length(dVector);
            float invDLenSqr = 1.0/dVectorLen;
            invDLenSqr *= invDLenSqr;
            float hgb = hg * invLens[1] * invDLenSqr;

            float gaa = -invCLenSqr * invLens[1];
            float gbb = invDLenSqr * invLens[1];

            //Evaluate the forces 
            if (myIdxInDihedral <= 1) {
                myForce = derivOfPotential * gaa * cVector; 
                if (myIdxInDihedral == 1) {
                    float3 s2;
                    s2 = derivOfPotential * ( fga * cVector - hgb * dVector);
                    myForce = s2 - myForce;
                }
            }
            else {
                myForce = derivOfPotential * gbb * dVector; 
                if (myIdxInDihedral == 2) {
                    float3 s2;
                    s2 = derivOfPotential * ( fga * cVector - hgb * dVector);
                    myForce = s2 - myForce;
                }
            } 
            return myForce;*/
        }
 
        //Compute the energy of this potential
        inline __device__ float potential(DihedralPeriodicType dihedralType, float phi) {
            float phiDiff = phi - dihedralType.phiRef;
            return (
                               dihedralType.coefs[0] * (1.0f - cosf(phiDiff))
                               + dihedralType.coefs[1] * (1.0f - cosf(3.0f*phiDiff)) 
                        )
                ;
        }
};

#endif


