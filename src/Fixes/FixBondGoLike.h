#pragma once
#ifndef FIXGOLIKE_H
#define FIXGOLIKE_H

#include "Bond.h"
#include "FixBond.h"
#include "BondEvaluatorGoLike.h"
void export_FixBondGoLike();

class FixBondGoLike : public FixBond<BondGoLike, BondGPU, BondGoLikeType> {

public:
    //int maxBondsPerBlock;
    //DataSet *eng;
    //DataSet *press;

    FixBondGoLike(boost::shared_ptr<State> state_, std::string handle);

    ~FixBondGoLike(){};

    void compute(bool);
    void singlePointEng(float *);
    std::string restartChunk(std::string format);
    bool readFromRestart(pugi::xml_node restData);
    BondEvaluatorGoLike evaluator;

    // HEY - NEED TO IMPLEMENT REFRESHATOMS
    // consider that if you do so, max bonds per block could change
    //bool refreshAtoms();

    void createBond(Atom *, Atom *, double, double, int);  // by ids
    void setBondTypeCoefs(int, double, double);

    const BondGoLike getBond(size_t i) {
        return boost::get<BondGoLike>(bonds[i]);
    }
    virtual std::vector<BondVariant> *getBonds() {
        return &bonds;
    }

    //std::vector<pair<int, std::vector<int> > > neighborlistExclusions();

};

#endif
