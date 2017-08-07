#include "FixChargePairDebyeHuckel.h"
#include "BoundsGPU.h"
#include "GPUData.h"
#include "GridGPU.h"
#include "State.h"

#include "boost_for_export.h"
#include "cutils_func.h"
#include "EvaluatorWrapper.h"
#include "PairEvaluatorNone.h"
// #include <cmath>

namespace py=boost::python;
using namespace std;

const std::string chargePairDHType = "ChargePairDH";

//Pairwise Debye Huckel 
//force calculation:
//  F=q_i*q_j/(4*PI*eps_0*eps*r_ij)*exp(-r_ij/lambda_D)*[1/r_ij*(1/r_ij + 1/lambda_D)]



//    compute_cu<<<NBLOCK(nAtoms), PERBLOCK>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), neighborCounts, grid.neighborlist.data(), grid.perBlockArray.d_data.data(), gpd.qs(activeIdx), alpha, r_cut, A, shift, state->boundsGPU, state->devManager.prop.warpSize, 0.5);// state->devManager.prop.warpSize, sigmas.getDevData(), epsilons.getDevData(), numTypes, state->rCut, state->boundsGPU, oneFourStrength);

FixChargePairDH::FixChargePairDH(SHARED(State) state_, string handle_, string groupHandle_) : FixCharge(state_, handle_, groupHandle_, chargePairDHType, true) {
   setParameters(temp, ionic, r_cut);
   canOffloadChargePairCalc = true;
   setEvalWrapper();
};

void FixChargePairDH::setParameters(float temp_, float ionic_, float r_cut_)
{
    temp = temp_;
    ionic = ionic_;
    r_cut = r_cut_;
    
    //Charge of an electron
    double ec = 1.60217653E-19;
    double temp_eps0 = 8.8541878176E-22;
    double na = 6.0221415E23;
    double kb = 1.3806505E-23;
    
    //Calculate the epsilon parameter as a function of temperature and ionic strength
    double temp_eps = (249.4 - 0.788 * temp + 7.20E-4 * temp * temp);
    temp_eps *= (1.000 - 0.2551 * ionic + 5.151E-2 * ionic * ionic - 6.889E-3 * ionic * ionic * ionic);

    //Calculate the debye length from the given temperature and ionic strength calculations
    double lambda = sqrt(temp_eps * temp_eps0 * kb * temp * 1.0E27 / (2.0f * na * ec * ec * ionic)); 
    //We store the inverse as that is the only form that is used by the Evaluator
    lambdai = 1.0 / lambda;
    //Convert epsilon to kcal/mol from kJ/mol
    epsi = 1.0 / (temp_eps);  
    //printf("Epsi is %f, debye length is %f\n", 1.0/epsi, 1.0/lambdai);
}

void FixChargePairDH::compute(int virialMode) {
    int nAtoms = state->atoms.size();
    int nPerRingPoly = state->nPerRingPoly;
    GPUData &gpd = state->gpd;
    GridGPU &grid = state->gridGPU;
    int activeIdx = gpd.activeIdx();
    uint16_t *neighborCounts = grid.perAtomArray.d_data.data();
    float *neighborCoefs = state->specialNeighborCoefs;
    evalWrap->compute(nAtoms, nPerRingPoly, gpd.xs(activeIdx), gpd.fs(activeIdx),
                  neighborCounts, grid.neighborlist.data(), grid.perBlockArray.d_data.data(),
                  state->devManager.prop.warpSize, nullptr, 0, state->boundsGPU,
                  neighborCoefs[0], neighborCoefs[1], neighborCoefs[2], gpd.virials.d_data.data(), gpd.qs(activeIdx), r_cut, virialMode,  nThreadPerBlock(), nThreadPerAtom());
  //  compute_charge_pair_DH_cu<<<NBLOCK(nAtoms), PERBLOCK>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), neighborIdxs, grid.neighborlist.tex, gpd.qs(activeIdx), alpha,r_cut, A,shift, state->boundsGPU, 0.5);
}

std::vector<float> FixChargePairDH::getRCuts() {
    std::vector<float> vals(1, state->rCut);
    return vals;
}

void FixChargePairDH::setEvalWrapper() {
    if (evalWrapperMode == "offload") {
        if (hasOffloadedChargePairCalc) {
            evalWrap = pickEvaluator<EvaluatorNone, 1, false>(EvaluatorNone(), nullptr); //nParams arg is 1 rather than zero b/c can't have zero sized argument on device
        } else {
            evalWrap = pickEvaluator<EvaluatorNone, 1, false>(EvaluatorNone(), this);
        }
    } else if (evalWrapperMode == "self") {
        evalWrap = pickEvaluator<EvaluatorNone, 1, false>(EvaluatorNone(), this);
    }

}

ChargeEvaluatorDH FixChargePairDH::generateEvaluator() {
    return ChargeEvaluatorDH(lambdai, epsi, state->units.qqr_to_eng);
}

void export_FixChargePairDH() {
    py::class_<FixChargePairDH, SHARED(FixChargePairDH), boost::python::bases<FixCharge> > (
        "FixChargePairDH",
        py::init<SHARED(State), string, string> (
            py::args("state", "handle", "groupHandle"))
    )
    .def("setParameters", &FixChargePairDH::setParameters,
            (py::arg("temp"), py::arg("ionic"), py::arg("r_cut"))
        )
    ;
}
