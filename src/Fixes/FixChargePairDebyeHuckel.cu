#include "FixChargePairDebyeHuckel.h"
#include "BoundsGPU.h"
#include "GPUData.h"
#include "GridGPU.h"
#include "State.h"

#include "boost_for_export.h"
#include "cutils_func.h"
// #include <cmath>

namespace py=boost::python;
using namespace std;

const std::string chargePairDHType = "ChargePairDH";

//Pairwise Debye Huckel 
//force calculation:
//  F=q_i*q_j/(4*PI*eps_0*eps*r_ij)*exp(-r_ij/lambda_D)*[1/r_ij*(1/r_ij + 1/lambda_D)]



//    compute_cu<<<NBLOCK(nAtoms), PERBLOCK>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), neighborCounts, grid.neighborlist.data(), grid.perBlockArray.d_data.data(), gpd.qs(activeIdx), alpha, r_cut, A, shift, state->boundsGPU, state->devManager.prop.warpSize, 0.5);// state->devManager.prop.warpSize, sigmas.getDevData(), epsilons.getDevData(), numTypes, state->rCut, state->boundsGPU, oneFourStrength);
__global__ void compute_charge_pair_DH_cu(int nAtoms, float4 *xs, float4 *fs, uint16_t *neighborCounts, uint *neighborlist, uint32_t *cumulSumMaxPerBlock, float *qs, float lambdai, float epsi, float qqr_to_eng,  BoundsGPU bounds, int warpSize, float onetwoStr, float onethreeStr, float onefourStr) {

    float multipliers[4] = {1, onetwoStr, onethreeStr, onefourStr};
    int idx = GETIDX();
    if (idx < nAtoms) {
        float4 posWhole = xs[idx];
        float3 pos = make_float3(posWhole);

        float3 forceSum = make_float3(0, 0, 0);
        float qi = qs[idx];//tex2D<float>(qs, XIDX(idx, sizeof(float)), YIDX(idx, sizeof(float)));

        //printf("start, end %d %d\n", start, end);
        int baseIdx = baseNeighlistIdx(cumulSumMaxPerBlock, warpSize);
        int numNeigh = neighborCounts[idx];
        for (int i=0; i<numNeigh; i++) {
            int nlistIdx = baseIdx + warpSize * i;
            uint otherIdxRaw = neighborlist[nlistIdx];
            uint neighDist = otherIdxRaw >> 30;
            uint otherIdx = otherIdxRaw & EXCL_MASK;
            float3 otherPos = make_float3(xs[otherIdx]);
            //then wrap and compute forces!
            float3 dr = bounds.minImage(pos - otherPos);
            float lenSqr = lengthSqr(dr);
            //   printf("dist is %f %f %f\n", dr.x, dr.y, dr.z);
            float multiplier = multipliers[neighDist];
            float len=sqrtf(lenSqr);
            float qj = qs[otherIdx];

            float rinv = 1.0f/len;
            float forceScalar = qi*qj*epsi*rinv*expf(-len*lambdai)*(rinv*(rinv+lambdai)) * multiplier;
    
            float3 forceVec = dr * forceScalar;
            forceSum += forceVec;

        }   
        fs[idx] += forceSum; //operator for float4 + float3

    }

}
FixChargePairDH::FixChargePairDH(SHARED(State) state_, string handle_, string groupHandle_) : FixCharge(state_, handle_, groupHandle_, chargePairDHType, true) {
   setParameters(temp, ionic);
   canOffloadChargePairCalc = true;
};

void FixChargePairDH::setParameters(float temp_, float ionic_)
{
    temp = temp_;
    ionic = ionic_;
    
    //Charge of an electron
    double ec = 1.60217653E-19;
    double temp_eps0 = 8.8541878176E-22;
    
    //Calculate the epsilon parameter as a function of temperature and ionic strength
    double temp_eps = ((249.4 - 0.788 * temp + 7.20E-4 * temp * temp) * 
                     (1.000 - 0.2551 * ionic + 5.151E-2 * ionic * ionic -
                      6.889E-3 * ionic * ionic * ionic)
                     )
                    ;

    //Calculate the debye length from the given temperature and ionic strength calculations
    double lambda = sqrt(temp_eps * temp_eps0 * 1.3806505E-23 * temp * 1.0E30 / (2.0f * 6.0221415E23 * ec * ec * ionic)); 
    //We store the inverse as that is the only form that is used by the Evaluator
    lambdai = 1.0 / lambda;
    //Convert epsilon to kcal/mol from kJ/mol
    epsi = 1.0 / (temp_eps);  
    printf("Epsi is %f, debye length is %f\n", 1.0/epsi, 1.0/lambdai);
}

void FixChargePairDH::compute(bool computeVirials) {
    int nAtoms = state->atoms.size();
    GPUData &gpd = state->gpd;
    GridGPU &grid = state->gridGPU;
    int activeIdx = gpd.activeIdx();
    uint16_t *neighborCounts = grid.perAtomArray.d_data.data();
    float *neighborCoefs = state->specialNeighborCoefs;
    compute_charge_pair_DH_cu<<<NBLOCK(nAtoms), PERBLOCK>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), neighborCounts, grid.neighborlist.data(), grid.perBlockArray.d_data.data(), gpd.qs(activeIdx), lambdai, epsi, state->units.qqr_to_eng, state->boundsGPU, state->devManager.prop.warpSize, neighborCoefs[0], neighborCoefs[1], neighborCoefs[2]);// state->devManager.prop.warpSize, sigmas.getDevData(), epsilons.getDevData(), numTypes, state->rCut, state->boundsGPU, oneFourStrength);
  //  compute_charge_pair_DH_cu<<<NBLOCK(nAtoms), PERBLOCK>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), neighborIdxs, grid.neighborlist.tex, gpd.qs(activeIdx), alpha,r_cut, A,shift, state->boundsGPU, 0.5);


}


std::vector<float> FixChargePairDH::getRCuts() {
    std::vector<float> vals(1, state->rCut);
    return vals;
}

void export_FixChargePairDH() {
    py::class_<FixChargePairDH, SHARED(FixChargePairDH), boost::python::bases<FixCharge> > (
        "FixChargePairDH",
        py::init<SHARED(State), string, string> (
            py::args("state", "handle", "groupHandle"))
    )
    .def("setParameters", &FixChargePairDH::setParameters,
            (py::arg("temp"), py::arg("ionic"))
        )
    ;
}
