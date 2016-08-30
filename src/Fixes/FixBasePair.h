#pragma once
#ifndef FIXBASEPAIR_H
#define FIXBASEPAIR_H

#include "FixPotentialMultiAtom.h"
#include "BasePair.h"
#include "BasePairEvaluator.h"

void export_FixBasePair();

class FixBasePair : public FixPotentialMultiAtom<BasePairVariant, BasePair3SPN2, BasePair, BasePairGPU, BasePair3SPN2Type, 4> {

private:
    BasePairEvaluator evaluator;
public:
    //DataSet *eng;
    //DataSet *press;

    FixBasePair(boost::shared_ptr<State> state_, std::string handle);

    void compute(bool);
    void singlePointEng(float *);

    void createBasePair(Atom *, Atom *, Atom *, Atom *, double, double, double, double, double, double, double,  int);
    void setBasePairTypeCoefs(int, double, double, double, double, double, double, double);

    bool readFromRestart(pugi::xml_node restData);

};

#endif
