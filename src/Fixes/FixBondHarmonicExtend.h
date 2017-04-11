#pragma once
#ifndef FIXBONDHARMONICEXTEND_H
#define FIXBONDHARMONICEXTEND_H

#include "Bond.h"
#include "FixBond.h"
#include "BondEvaluatorHarmonicExtend.h"
void export_FixBondHarmonicExtend();

class FixBondHarmonicExtend : public FixBond<BondHarmonicExtend, BondGPU, BondHarmonicExtendType> {

public:
    //int maxBondsPerBlock;
    //DataSet *eng;
    //DataSet *press;

    FixBondHarmonicExtend(boost::shared_ptr<State> state_, std::string handle);

    ~FixBondHarmonicExtend(){};

    void compute(bool);
    void singlePointEng(float *);
    std::string restartChunk(std::string format);
    bool readFromRestart();
    BondEvaluatorHarmonicExtend evaluator;

    // HEY - NEED TO IMPLEMENT REFRESHATOMS
    // consider that if you do so, max bonds per block could change
    //bool refreshAtoms();

    void createBond(Atom *, Atom *, double, double, double, int);  // by ids
    void setBondTypeCoefs(int, double, double, double);

    const BondHarmonicExtend getBond(size_t i) {
        return boost::get<BondHarmonicExtend>(bonds[i]);
    }
    virtual std::vector<BondVariant> *getBonds() {
        return &bonds;
    }

    //std::vector<pair<int, std::vector<int> > > neighborlistExclusions();

};

#endif
