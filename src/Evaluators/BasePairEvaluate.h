template <class BASEPAIRTYPE, class EVALUATOR, bool COMPUTEVIRIALS> //don't need BasePairGPU, are all BasePairGPU.  Worry about later 
__global__ void compute_force_basepair(int nBasePairs, float4 *xs, float4 *fs, int *idToIdxs, BasePairGPU *basepairs, int *startstops, BoundsGPU bounds, BASEPAIRTYPE *parameters_arg, int nParameters, Virial *virials, bool usingSharedMemForParams, EVALUATOR evaluator) {


    int idx = GETIDX();
    extern __shared__ char all_shr[];
    BASEPAIRTYPE *parameters;
    if (usingSharedMemForParams) {
        parameters = (BASEPAIRTYPE *) (all_shr);
        copyToShared<BASEPAIRTYPE>(parameters_arg, parameters, nParameters);
    } else {
        parameters = parameters_arg;
    }
    __syncthreads();
    if (idx < nBasePairs) {
        Virial sumVirials(0, 0, 0, 0, 0, 0);
        int idxs[4];
        BasePairGPU basepair = basepairs[idx];
    
        uint32_t typeFull = basepair.type;
        int type = (typeFull << 3) >> 3;
        BASEPAIRTYPE basepairType = parameters[type];

        float3 positions[4];

        for (int i=0; i<4; i++) {
            int idxOther = idToIdxs[basepair.ids[i]];
            idxs[i] = idxOther;
            positions[i] = make_float3(xs[idxOther]);
        }
        for (int i=1; i<4; i++) {
            positions[i] = positions[0] + bounds.minImage(positions[i]-positions[0]);
        }
        //for (int i=0; i<4; i++) {
        //    printf("position %d %f %f %f\n", i, positions[i].x, positions[i].y, positions[i].z);
       // }
        float3 directors[3]; //vb_xyz in lammps
        float lenSqrs[3]; //bnmag2 in lammps
        float lens[3]; //bnmag in lammps
        float invLenSqrs[3]; //sb in lammps
        float invLens[3];
        directors[0] = positions[1] - positions[0];
        directors[1] = positions[3] - positions[1];
        directors[2] = positions[3] - positions[2];
        for (int i=0; i<3; i++) {
            //printf("directors %d is %f %f %f\n", i, directors[i].x, directors[i].y, directors[i].z);
            lenSqrs[i] = lengthSqr(directors[i]);
            lens[i] = sqrtf(lenSqrs[i]);
            invLenSqrs[i] = 1.0f / lenSqrs[i];
            invLens[i] = 1.0f / lens[i];
         //   printf("inv len sqrs %d is %f\n", i, invLenSqrs[i]);
        }

     //   printf("c0 is %f\n", c0);
        float c12Mags[2];
        float invMagProds[2]; //r12c1, 2 in lammps
        float dotProd = dot(directors[1], directors[0]);
        invMagProds[0] = invLens[0] * invLens[1];
        c12Mags[0] = -dotProd * invMagProds[0]; //3spn2 indexing 

        //2nd angle calc now
        dotProd = dot(directors[2], directors[1]);
        invMagProds[1] = invLens[2] * invLens[1];
        c12Mags[1] = dotProd * invMagProds[1]; //lammps variable names are opaque

        float thetas[2];
        for (int i =0; i<2; i++) {
            if (c12Mags[i] > 1.0f) {
                c12Mags[i] = 1.0f;
            } else if (c12Mags[i] < -1.0f) {
                c12Mags[i] = -1.0f;
            }
            thetas[i] = acosf(c12Mags[i]);
        }

        //float c0 = dot(directors[0], directors[2]) * invLens[0] * invLens[2];

        //printf("theta0 is %f, theta1 is %f\n", thetas[0] * 180.0f / M_PI, thetas[1] * 180.0f / M_PI);


        //IS THIS EVEN NEEDED?????
        float scValues[2]; //???, is s1, s2, s12 in lammps
        for (int i=0; i<2; i++) {
            float x = max(1.0f - c12Mags[i]*c12Mags[i], 0.0f);
            float sqrtVal = max(sqrtf(x), EPSILON);
            scValues[i] = 1.0f / sqrtVal;
        }
        /*float scValues[2];
        for(int i=0; i<2; i++) {
            scValues[i] = sqrtf(1.0f / (1.0f - c12Mags[i] * c12Mags[i]));
        }*/

       //printf("sc values %f %f\n", scValues[0], scValues[1]);
     
        float3 cVector; //pad_xyz in 3spn2 lammps
        float3 dVector; //pbc_xyz -> will be sent to compute forces
        float3 eVector; //pac_xyz

        cVector.x = directors[0].y*directors[1].z - directors[0].z*directors[1].y;
        cVector.y = directors[0].z*directors[1].x - directors[0].x*directors[1].z;
        cVector.z = directors[0].x*directors[1].y - directors[0].y*directors[1].x;
        cVector *= invMagProds[0];

        dVector.x = directors[2].z*directors[1].y - directors[2].y*directors[1].z;
        dVector.y = directors[2].x*directors[1].z - directors[2].z*directors[1].x;
        dVector.z = directors[2].y*directors[1].x - directors[2].x*directors[1].y;
        dVector *= invMagProds[1];

        eVector.x = cVector.z*dVector.y - cVector.y*dVector.z;
        eVector.y = cVector.x*dVector.z - cVector.z*dVector.x;
        eVector.z = cVector.y*dVector.x - cVector.x*dVector.y;
        float dx = dot(eVector, directors[1]) * invLens[1];
        float c = -dot(cVector, dVector) * scValues[0] * scValues[1];

       //printf("c is %f\n", c);
        if (c > 1.0f) {
            c = 1.0f;
        } else if (c < -1.0f) {
            c = -1.0f;
        }
        float phi = acosf(c);
        //printf("phi is %f\n", phi);
        if (dx < 0) {
            phi = -phi;
        }

        //isb2 and isd2 in 3spn2 lammps
        for (int i=0; i<2; i++) { 
            scValues[i] *= scValues[i]; 
        }

   // printf("Theta0 is %f, theta1 is %f, and phi is %f\n", thetas[0] * 180.0f/ M_PI, thetas[1] * 180.0f / M_PI, phi * 180.0f / M_PI);

        float3 allForces[4];

        for(int i=0; i<4; i++) {
            allForces[i] = 0.0f * allForces[i];
        }

        evaluator.forces(basepairType, phi, thetas, cVector, dVector, scValues, invMagProds, invLens, invLenSqrs, c12Mags, c, directors, allForces);

        //printf("force 1: %f\t%f\t%f\n", allForces[0].x, allForces[0].y, allForces[0].z);
        //printf("force 2: %f\t%f\t%f\n", allForces[1].x, allForces[1].y, allForces[1].z);
        //printf("force 3: %f\t%f\t%f\n", allForces[2].x, allForces[2].y, allForces[2].z);
        //printf("force 4: %f\t%f\t%f\n", allForces[3].x, allForces[3].y, allForces[3].z);

        for (int i=0; i<4; i++) {
            atomicAdd(&(fs[idxs[i]].x), (allForces[i].x));
            atomicAdd(&(fs[idxs[i]].y), (allForces[i].y));
            atomicAdd(&(fs[idxs[i]].z), (allForces[i].z));
            //printf("f %d is %f %f %f\n", i, forces[i].x, forces[i].y, forces[i].z);
        }

        if (COMPUTEVIRIALS) {
            computeVirial(sumVirials, allForces[0], directors[0]);
            computeVirial(sumVirials, allForces[2], directors[1]);
            computeVirial(sumVirials, allForces[3], directors[1] + directors[2]);
            for (int i=0; i<6; i++) {
                //printf("virial %d %f\n", i, sumVirials[i]);
                atomicAdd(&(virials[idxs[0]][i]), sumVirials[i]);
            }            
        }
    }
}


template <class BASEPAIRTYPE, class EVALUATOR>
__global__ void compute_energy_basepair(int nAtoms, float4 *xs, float *perParticleEng, int *idToIdxs, BasePairGPU *basepairs, int *startstops, BoundsGPU bounds, BASEPAIRTYPE *parameters, int nParameters, EVALUATOR evaluator) {


    /*int idx = GETIDX();
    extern __shared__ char all_shr[];
    int idxBeginCopy = startstops[blockDim.x*blockIdx.x];
    int idxEndCopy = startstops[min(nAtoms, blockDim.x*(blockIdx.x+1))];
    BasePairGPU *basepairs_shr = (BasePairGPU *) all_shr;
    int sizeBasePairs = (idxEndCopy - idxBeginCopy) * sizeof(BasePairGPU);
    BASEPAIRTYPE *parameters_shr = (BASEPAIRTYPE *) (all_shr + sizeBasePairs);
    copyToShared<BasePairGPU>(basepairs + idxBeginCopy, basepairs_shr, idxEndCopy - idxBeginCopy);
    copyToShared<BASEPAIRTYPE>(parameters, parameters_shr, nParameters);
    __syncthreads();
    if (idx < nAtoms) {
  //      printf("going to compute %d\n", idx);
        int startIdx = startstops[idx];
        int endIdx = startstops[idx+1];
        //so start/end is the index within the entire bond list.
        //startIdx - idxBeginCopy gives my index in shared memory
        int shr_idx = startIdx - idxBeginCopy;
        int n = endIdx - startIdx;
        if (n) {
            int myIdxInBasePair = basepairs_shr[shr_idx].type >> 29;
            int idSelf = basepairs_shr[shr_idx].ids[myIdxInBasePair];
            
            int idxSelf = idToIdxs[idSelf]; 
        
            float3 pos = make_float3(xs[idxSelf]);
           // printf("I am idx %d and I am evaluating atom with pos %f %f %f\n", idx, pos.x, pos.y, pos.z);
            float energySum = 0;
            for (int i=0; i<n; i++) {
                BasePairGPU basepair = basepairs_shr[shr_idx + i];
                uint32_t typeFull = basepair.type;
                myIdxInBasePair = typeFull >> 29;
                int type = (typeFull << 3) >> 3;
                BASEPAIRTYPE basepairType = parameters_shr[type];
                //USE SHARED AGAIN ONCE YOU FIGURE OUT BUG

                float3 positions[4];
                positions[myIdxInBasePair] = pos;
                int toGet[3];
                if (myIdxInBasePair==0) {
                    toGet[0] = 1;
                    toGet[1] = 2;
                    toGet[2] = 3;
                } else if (myIdxInBasePair==1) {
                    toGet[0] = 0;
                    toGet[1] = 2;
                    toGet[2] = 3;
                } else if (myIdxInBasePair==2) {
                    toGet[0] = 0;
                    toGet[1] = 1;
                    toGet[2] = 3;
                } else if (myIdxInBasePair==3) {
                    toGet[0] = 0;
                    toGet[1] = 1;
                    toGet[2] = 2;
                }


                for (int i=0; i<3; i++) {
                    int idxOther = idToIdxs[basepair.ids[toGet[i]]];
                    positions[toGet[i]] = make_float3(xs[idxOther]);
                }
                for (int i=1; i<3; i++) {
                    positions[i] = positions[0] + bounds.minImage(positions[i]-positions[0]);
                }
                float3 directors[3]; //vb_xyz in lammps
                float lenSqrs[3]; //bnmag2 in lammps
                float lens[3]; //bnmag in lammps
                float invLenSqrs[3]; //sb in lammps
                float invLens[3];
                directors[0] = positions[0] - positions[1];
                directors[1] = positions[3] - positions[1];
                directors[2] = positions[3] - positions[2];
                for (int i=0; i<3; i++) {
                    //printf("directors %d is %f %f %f\n", i, directors[i].x, directors[i].y, directors[i].z);
                    lenSqrs[i] = lengthSqr(directors[i]);
                    lens[i] = sqrtf(lenSqrs[i]);
                    invLenSqrs[i] = 1.0f / lenSqrs[i];
                    invLens[i] = 1.0f / lens[i];
                 //   printf("inv len sqrs %d is %f\n", i, invLenSqrs[i]);
                }

             //   printf("c0 is %f\n", c0);
                float c12Mags[3];
                float invMagProds[3]; //r12c1, 2 in lammps
                float dotProd = dot(directors[1], directors[0]);
                dotProd *= -1;
                invMagProds[0] = invLens[0] * invLens[1];
                c12Mags[0] = dotProd * invMagProds[0]; //3spn2 indexing 
                for (int i=1; i<3; i++) { // 3spn2 indexing
                    float dotProd = dot(directors[2], directors[i]);
              //      printf("ctmp is %f\n", dotProd);
                    invMagProds[i] = invLens[2] * invLens[i];
                    c12Mags[i] = dotProd * invMagProds[i]; //lammps variable names are opaque
              //      printf("c12 mag %d %f\n", i, c12Mags[i]);
                }
                float thetas[2];
                thetas[0] = -c12Mags[0];
                thetas[1] = c12Mags[1];

                for (int i =0; i<2; i++) {
                    if (thetas[i] > 1.0f) {
                        thetas[i] = 1.0f;
                    } else if (thetas[i] < -1.0f) {
                        thetas[i] = -1.0f;
                    }
                    thetas[i] = acosf(thetas[i]);
                }


                float scValues[2]; //???, is s1, s2, s12 in lammps
                for (int i=0; i<2; i++) {
                    float x = max(1 - c12Mags[i]*c12Mags[i], 0.0f);
                    float sqrtVal = max(sqrtf(x), EPSILON);
                    scValues[i] = 1.0 / sqrtVal;
                }


                for (int i=0; i<3; i++) { //Not sure if I want the squared values
                    scValues[i] *= scValues[i]; 
                }
             //   printf("sc values %f %f %f\n", scValues[0], scValues[1], scValues[2]);
             
                float3 cVector; //pad_xyz in 3spn2 lammps
                float3 dVector; //pbc_xyz -> will be sent to compute forces
                float3 eVector; //pac_xyz
                cVector.x = directors[0].y*directors[1].z - directors[0].z*directors[1].y;
                cVector.y = directors[0].z*directors[1].x - directors[0].x*directors[1].z;
                cVector.z = directors[0].x*directors[1].y - directors[0].y*directors[1].x;
                float cVectorLen = length(cVector);
                dVector.x = directors[2].z*directors[0].y - directors[2].y*directors[0].z;
                dVector.y = directors[2].x*directors[0].z - directors[2].z*directors[0].x;
                dVector.z = directors[2].y*directors[0].x - directors[2].x*directors[0].y;
                float dVectorLen = length(dVector);
                eVector.x = cVector.z*dVector.y - cVector.y*dVector.z;
                eVector.y = cVector.x*dVector.z - cVector.z*dVector.x;
                eVector.z = cVector.y*dVector.x - cVector.x*dVector.y;
                float eVectorLen = length(eVector);
                float dx = dot(eVector, directors[1]) * invLens[1] / eVectorLen;
                float c = dot(cVector,dVector) / cVectorLen / dVectorLen;

            //    printf("c is %f\n", c);
                if (c > 1.0f) {
                    c = 1.0f;
                } else if (c < -1.0f) {
                    c = -1.0f;
                }
                float phi = acosf(c);
                //printf("phi is %f\n", phi);
                if (dx < 0) {
                    phi = -phi;
                }
            
                energySum += evaluator.energy(basepairType, phi, alpha, range, thetas, invLens, invLenSqrs, directors);

                


            }
            perParticleEng[idxSelf] += energySum;
        }
    }*/
}
