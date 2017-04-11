#pragma once
#ifndef FIX_CHARGEPAIR_DH_H
#define FIX_CHARGEPAIR_DH_H

//#include "AtomParams.h"
#include "GPUArrayTex.h"
#include "FixCharge.h"
#include "ChargeEvaluatorDH.h"

class State;

void export_FixChargePairDH();

extern const std::string chargePairDHType;
class FixChargePairDH : public FixCharge {

private:
    float epsi;  // temps for compute
    float lambdai;

protected:
    float temp;
    float ionic;
    float r_cut;

public:
    FixChargePairDH(boost::shared_ptr<State> state_,
                     std::string handle_, std::string groupHandle_);

    //Right now temp has to be in K and ionic has to be in M for units of params to work out
    void setParameters(float temp_, float ionic_, float r_cut_);
    void compute(bool);
    ChargeEvaluatorDH generateEvaluator();
    void setEvalWrapper();
    std::vector<float> getRCuts();

};

#endif
