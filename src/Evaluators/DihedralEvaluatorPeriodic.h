#pragma once
#ifndef EVALUATOR_DIHEDRAL_PERIODIC
#define EVALUATOR_DIHEDRAL_PERIODIC

#include "cutils_math.h"
#include "Dihedral.h"

//This is the dihedral potential from the Golike interactions developed by Clementi et al.
class DihedralEvaluatorPeriodic {
    public:

                //float3 myForce = evaluator.force(dihedralType, phi, c, scValues, invLenSqrs, c12Mags, c0, c, invMagProds, c12Mags, invLens, directors, myIdxInDihedral);
        inline __device__ float3 force(DihedralPeriodicType dihedralType, float phi, float scValues[3], float invLenSqrs[3], float c12Mangs[3], float c0, float c, float invMagProds[2], float c12Mags[2], float invLens[3], float3 directors[3], int myIdxInDihedral) {
            float3 myForce;
        
            //The difference between system dihedral and reference (dtor in 3SPN.2 Lammps)
            float diffPhi = phi-dihedralType.phiRef;
    
            //Make sure difference follows periodic behavior correctly
            if (diffPhi < -M_PI) diffPhi += 2.0f * M_PI;
            else if (diffPhi > M_PI) diffPhi -= 2.0f * M_PI;

            
            float derivOfPotential = (
                    - dihedralType.coefs[0] * sinf(diffPhi)
                    - 3.0f * dihedralType.coefs[1] * sinf(3.0f*diffPhi)
                    )
                ;

            //Evaluate derivatives
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
            return myForce;
        }
        inline __device__ void forcesAll(DihedralPeriodicType dihedralType, float phi, float scValues[3], float invLenSqrs[3], float c12Mangs[3], float c0, float c, float invMagProds[2], float c12Mags[2], float invLens[3], float3 directors[3], float3 forces[4]) {
            //The difference between system dihedral and reference (dtor in 3SPN.2 Lammps)
            float diffPhi = phi-dihedralType.phiRef;
    
            //Make sure difference follows periodic behavior correctly
            if (diffPhi < -M_PI) diffPhi += 2.0f * M_PI;
            else if (diffPhi > M_PI) diffPhi -= 2.0f * M_PI;

            
            float derivOfPotential = (
                    dihedralType.coefs[0] * sinf(diffPhi)
                    + 3.0f * dihedralType.coefs[1] * sinf(3.0f*diffPhi)
                    )
                ;

            //Evaluate derivatives
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

            //Only applied to middle two atoms            
            float3 sFloat3 = derivOfPotential * ( fga * cVector - hgb * dVector); 

            forces[0].x = derivOfPotential * gaa * cVector.x;
            forces[0].y = derivOfPotential * gaa * cVector.y;
            forces[0].z = derivOfPotential * gaa * cVector.z;
            forces[1] = sFloat3 - forces[0];


            forces[3].x = derivOfPotential * gbb * dVector.x;
            forces[3].y = derivOfPotential * gbb * dVector.y;
            forces[3].z = derivOfPotential * gbb * dVector.z;
            forces[2] = -sFloat3 - forces[3];


        }
    
        //Compute the energy of this potential
        inline __device__ float energy(DihedralPeriodicType dihedralType, float phi, float scValues[3], float invLenSqrs[3], float c12Mangs[3], float c0, float c, float invMagProds[2], float c12Mags[2], float invLens[3], float3 directors[3], int myIdxInDihedral) {
            float phiDiff = phi - dihedralType.phiRef;
            float eng = (
                               dihedralType.coefs[0] * (1.0f - cosf(phiDiff))
                               + dihedralType.coefs[1] * (1.0f - cosf(3.0f*phiDiff)) 
                        )
                ;
            return (float) (0.25f * eng);

        }
};

#endif


