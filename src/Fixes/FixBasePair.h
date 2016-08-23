#pragma once
#ifndef FIXBASEPAIR_H
#define FIXBASEPAIR_H

#include <boost/python.hpp>

#include "FixPotentialMultiAtom.h"
#include "BasePair.h"
#include "BasePairEvaluate.h"

void export_FixBasePair();

class FixBasePair : public FixPotentialMultiAtom<BasePairVariant, BasePair, BasePair, BasePairGPU, BasePairType, 4> {

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

    //std::vector<pair<int, std::vector<int> > > neighborlistExclusions();
    bool readFromRestart(pugi::xml_node restData);

};

#endif
