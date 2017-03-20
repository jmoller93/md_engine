#pragma once
#ifndef FIXBASEPAIR_H
#define FIXBASEPAIR_H

#include <boost/python.hpp>

#include "FixPotentialMultiAtom.h"
#include "BasePair.h"
#include "BasePairEvaluator.h"

void export_FixBasePair3SPN2();

class FixBasePair3SPN2 : public FixPotentialMultiAtom<BasePairVariant, BasePair3SPN2, BasePair, BasePairGPU, BasePair3SPN2Type, 4> {

private:
    BasePairEvaluator3SPN2 evaluator;

protected:
    float alpha;
    float range;
    
public:
    //DataSet *eng;
    //DataSet *press;
    std::vector<BondVariant> bonds;

    FixBasePair3SPN2(boost::shared_ptr<State> state_, std::string handle);

    void setParameters(float, float);
    void compute(bool);
    void singlePointEng(float *);
    bool prepareForRun();

    void createBasePair(Atom *, Atom *, Atom *, Atom *, double, double, double, double, double,  int);
    void setBasePairTypeCoefs(int, double, double, double, double, double);

    bool readFromRestart();

    virtual std::vector<BondVariant> *getBonds() {
        return &bonds;
    }

};

#endif