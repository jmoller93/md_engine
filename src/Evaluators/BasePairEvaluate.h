template <class BASEPAIRTYPE, class EVALUATOR, bool COMPUTEVIRIALS> //don't need BasePairGPU, are all BasePairGPU.  Worry about later 
__global__ void compute_force_basepair(int nAtoms, float4 *xs, float4 *forces, int *idToIdxs, BasePairGPU *basepairs, int *startstops, BoundsGPU bounds, BASEPAIRTYPE *parameters, int nParameters, Virial *virials, EVALUATOR evaluator) {


    int idx = GETIDX();
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
            Virial sumVirials(0, 0, 0, 0, 0, 0);
            int myIdxInBasePair = basepairs_shr[shr_idx].type >> 29;
            int idSelf = basepairs_shr[shr_idx].ids[myIdxInBasePair];
            int idxSelf = idToIdxs[idSelf]; 
        
            float3 pos = make_float3(xs[idxSelf]);
           // printf("I am idx %d and I am evaluating atom with pos %f %f %f\n", idx, pos.x, pos.y, pos.z);
            float3 forceSum = make_float3(0, 0, 0);
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
                //for (int i=0; i<4; i++) {
                //    printf("position %d %f %f %f\n", i, positions[i].x, positions[i].y, positions[i].z);
               // }
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
                    float dotProd = dot(directors[3], directors[i]);
              //      printf("ctmp is %f\n", dotProd);
                    invMagProds[i] = invLens[3] * invLens[i];
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


                float scValues[3]; //???, is s1, s2, s12 in lammps
                for (int i=0; i<3; i++) {
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
                float dx = dot(eVector, directors[1]) * invLens[1] / eVectorLen;
                c = dot(cVector,dVector) / cVectorLen / dVectorLen;

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

                //printf("no force\n");
                if (COMPUTEVIRIALS) {
                    float3 allForces[4];
                    evaluator.forcesAll(basepairType, phi, thetas, cVector, dVector, scValues, invLens, invLenSqrs, c12Mags, c0, c, invMagProds, c12Mags, invLens, directors, allForces);
                    computeVirial(sumVirials, allForces[0], directors[0]);
                    computeVirial(sumVirials, allForces[2], directors[1]);
                    computeVirial(sumVirials, allForces[3], directors[1] + directors[2]);
                    forceSum += allForces[myIdxInBasePair];
                    
                } else {
                    float3 myForce = evaluator.force(basepairType, phi, thetas, cVector, dVector, scValues, invLens, invLenSqrs, c12Mags, c, invMagProds, invLens, directors, myIdxInBasePair);
                    forceSum += myForce;
                }
                


            }
            forces[idxSelf] += forceSum;
            if (COMPUTEVIRIALS) {
                sumVirials *= 0.25f;
                virials[idx] += sumVirials;
            }
        }
    }
}


template <class BASEPAIRTYPE, class EVALUATOR>
__global__ void compute_energy_basepair(int nAtoms, float4 *xs, float *perParticleEng, int *idToIdxs, BasePairGPU *basepairs, int *startstops, BoundsGPU bounds, BASEPAIRTYPE *parameters, int nParameters, EVALUATOR evaluator) {


    int idx = GETIDX();
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
                    float dotProd = dot(directors[3], directors[i]);
              //      printf("ctmp is %f\n", dotProd);
                    invMagProds[i] = invLens[3] * invLens[i];
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
                float dx = dot(eVector, directors[1]) * invLens[1] / eVectorLen;
                c = dot(cVector,dVector) / cVectorLen / dVectorLen;

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
            
                energySum += evaluator.energy(basepairType, phi, thetas, invLens,  invLenSqrs, directors, myIdxInBasePair);

                


            }
            perParticleEng[idxSelf] += energySum;
        }
    }
}
