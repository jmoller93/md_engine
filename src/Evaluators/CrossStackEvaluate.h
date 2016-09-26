template <class CROSSSTACKTYPE, class EVALUATOR, bool COMPUTEVIRIALS>
__global__ void compute_force_crossstack(int nCrossStacks, float4 *xs, float4 *fs, int *idToIdxs, CrossStackGPU *crossstacks, int *startstops, BoundsGPU bounds, CROSSSTACKTYPE *parameters_arg, int nParameters, Virial *virials, bool usingSharedMemForParams, EVALUATOR evaluator) {


    int idx = GETIDX();
    extern __shared__ char all_shr[];
    CROSSSTACKTYPE *parameters;
    if (usingSharedMemForParams) {
        parameters = (CROSSSTACKTYPE *) (all_shr);
        copyToShared<CROSSSTACKTYPE>(parameters_arg, parameters, nParameters);
    } else {
        parameters = parameters_arg;
    }
    __syncthreads();
    if (idx < nCrossStacks) {
        Virial sumVirials(0, 0, 0, 0, 0, 0);
        int idxs[6];
        CrossStackGPU crossstack = crossstacks[idx];
    
        uint32_t typeFull = crossstack.type;
        int type = (typeFull << 3) >> 3;
        CROSSSTACKTYPE crossstackType = parameters[type];

        float3 positions[6];

        //Getting the positions of all particles
        for (int i=0; i<6; i++) {
            int idxOther = idToIdxs[crossstack.ids[i]];
            idxs[i] = idxOther;
            positions[i] = make_float3(xs[idxOther]);
        }
        for (int i=1; i<6; i++) {
            positions[i] = positions[0] + bounds.minImage(positions[i]-positions[0]);
        }

        // Schematic of layout for all sites in a 3spn2 cross stacking interaction

        /* A drawing of the system
        a====b ---- d====c   --> The Hydrogen Bond
              \    /
           Q1  \  /  Q2
                \/  Cross stacking interactions
                /\
               /  \
              /    \
             f      e
        */

        //As division is a longer operation than multiplication it is faster to store the inverses
        float3 directors[5]; //distance between particles
        float lenSqrs[5]; //magnitude of distance btw particles squared
        float lens[5]; //magnitude of distance btw particles
        float invLenSqrs[5]; //inverse magnitude squared
        float invLens[5]; //inverse magnitude

        //Vectors between specific atoms
        directors[0] = positions[1] - positions[0]; //b-a vector
        directors[1] = positions[3] - positions[1]; //d-b vector (hbond vector)
        directors[2] = positions[3] - positions[2]; //d-c vector 
        directors[3] = positions[1] - positions[4]; //b-e vector (1st cstk vector)
        directors[4] = positions[3] - positions[5]; //d-f vector (2nd cstk vector)

        for (int i=0; i<5; i++) {
            lenSqrs[i] = lengthSqr(directors[i]);
            lens[i] = sqrtf(lenSqrs[i]);
            invLenSqrs[i] = 1.0f / lenSqrs[i];
            invLens[i] = 1.0f / lens[i];
        }

        //Evaluation of the angles
        float c12Mags[3]; //Cosines of all relevant angles
        float invMagProds[3]; //Product of inverse magnitudes

        //Angle 0 (b-a = d-c angle) theta3
        //We are using this as angle 0 because it is the same as the lammps notation
        //In other words, blame Dan Hinckley
        float dotProd = dot(directors[0], directors[2]); 
        invMagProds[0] = invLens[0] * invLens[2];
        c12Mags[0] = dotProd * invMagProds[0]; 

        //Angle 1 (a-b-e cstck 1 angle) theta1
        dotProd = dot(directors[0], directors[3]);
        invMagProds[1] = invLens[0] * invLens[3];
        c12Mags[1] = dotProd * invMagProds[1]; 

        //Angle 2 (c-d-f cstck 2 angle) theta2
        dotProd = dot(directors[2], directors[4]);
        invMagProds[1] = invLens[2] * invLens[4];
        c12Mags[2] = dotProd * invMagProds[2]; 

        //Make sure the angles are periodic and evaluate
        float thetas[3];
        for (int i =0; i<2; i++) {
            if (c12Mags[i] > 1.0f) {
                c12Mags[i] = 1.0f;
            } else if (c12Mags[i] < -1.0f) {
                c12Mags[i] = -1.0f;
            }
            thetas[i] = acosf(c12Mags[i]);
        }

        //The forces and initialization
        float3 allForces[4];

        for(int i=0; i<4; i++) {
            allForces[i] = 0.0f * allForces[i];
        }

        //Run the cross stacking evaluator kernel
        evaluator.forces(crossstackType, thetas, invMagProds, invLens, invLenSqrs, c12Mags, directors, allForces);

        //The atomic adds make bonded interactions run quicker
        for (int i=0; i<4; i++) {
            atomicAdd(&(fs[idxs[i]].x), (allForces[i].x));
            atomicAdd(&(fs[idxs[i]].y), (allForces[i].y));
            atomicAdd(&(fs[idxs[i]].z), (allForces[i].z));
        }

        //Compute the virials, NOTE: may not be working at this moment, but running NPT on a single molecule would be unhelpful
        if (COMPUTEVIRIALS) {
            computeVirial(sumVirials, allForces[0], directors[0]);
            computeVirial(sumVirials, allForces[2], directors[1]);
            computeVirial(sumVirials, allForces[3], directors[1] + directors[2]);
            for (int i=0; i<6; i++) {
                atomicAdd(&(virials[idxs[0]][i]), sumVirials[i]);
            }            
        }
    }
}


template <class CROSSSTACKTYPE, class EVALUATOR>
__global__ void compute_energy_crossstack(int nAtoms, float4 *xs, float *perParticleEng, int *idToIdxs, CrossStackGPU *crossstacks, int *startstops, BoundsGPU bounds, CROSSSTACKTYPE *parameters, int nParameters, EVALUATOR evaluator) {


    /*int idx = GETIDX();
    extern __shared__ char all_shr[];
    int idxBeginCopy = startstops[blockDim.x*blockIdx.x];
    int idxEndCopy = startstops[min(nAtoms, blockDim.x*(blockIdx.x+1))];
    CrossStackGPU *crossstacks_shr = (CrossStackGPU *) all_shr;
    int sizeCrossStacks = (idxEndCopy - idxBeginCopy) * sizeof(CrossStackGPU);
    CROSSSTACKTYPE *parameters_shr = (CROSSSTACKTYPE *) (all_shr + sizeCrossStacks);
    copyToShared<CrossStackGPU>(crossstacks + idxBeginCopy, crossstacks_shr, idxEndCopy - idxBeginCopy);
    copyToShared<CROSSSTACKTYPE>(parameters, parameters_shr, nParameters);
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
            int myIdxInCrossStack = crossstacks_shr[shr_idx].type >> 29;
            int idSelf = crossstacks_shr[shr_idx].ids[myIdxInCrossStack];
            
            int idxSelf = idToIdxs[idSelf]; 
        
            float3 pos = make_float3(xs[idxSelf]);
           // printf("I am idx %d and I am evaluating atom with pos %f %f %f\n", idx, pos.x, pos.y, pos.z);
            float energySum = 0;
            for (int i=0; i<n; i++) {
                CrossStackGPU crossstack = crossstacks_shr[shr_idx + i];
                uint32_t typeFull = crossstack.type;
                myIdxInCrossStack = typeFull >> 29;
                int type = (typeFull << 3) >> 3;
                CROSSSTACKTYPE crossstackType = parameters_shr[type];
                //USE SHARED AGAIN ONCE YOU FIGURE OUT BUG

                float3 positions[4];
                positions[myIdxInCrossStack] = pos;
                int toGet[3];
                if (myIdxInCrossStack==0) {
                    toGet[0] = 1;
                    toGet[1] = 2;
                    toGet[2] = 3;
                } else if (myIdxInCrossStack==1) {
                    toGet[0] = 0;
                    toGet[1] = 2;
                    toGet[2] = 3;
                } else if (myIdxInCrossStack==2) {
                    toGet[0] = 0;
                    toGet[1] = 1;
                    toGet[2] = 3;
                } else if (myIdxInCrossStack==3) {
                    toGet[0] = 0;
                    toGet[1] = 1;
                    toGet[2] = 2;
                }


                for (int i=0; i<3; i++) {
                    int idxOther = idToIdxs[crossstack.ids[toGet[i]]];
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
            
                energySum += evaluator.energy(crossstackType, phi, alpha, range, thetas, invLens, invLenSqrs, directors);

                


            }
            perParticleEng[idxSelf] += energySum;
        }
    }*/
}
