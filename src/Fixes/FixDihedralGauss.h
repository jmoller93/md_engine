#pragma once
#ifndef FIXDIHEDRAL_GAUSS_H
#define FIXDIHEDRAL_GAUSS_H

#include "FixPotentialMultiAtom.h"
#include "Dihedral.h"
#include "DihedralEvaluatorGauss.h"

void export_FixDihedralGauss();

class FixDihedralGauss : public FixPotentialMultiAtom<DihedralVariant, DihedralGauss, Dihedral, DihedralGPU, DihedralGaussType, 4> {

private:
    DihedralEvaluatorGauss evaluator;
public:
    //DataSet *eng;
    //DataSet *press;

    FixDihedralGauss(boost::shared_ptr<State> state_, std::string handle);

    void compute(bool);
    void singlePointEng(float *);

    void createDihedral(Atom *, Atom *, Atom *, Atom *, double, double, double, int type_);
    void setDihedralTypeCoefs(int, double, double, double);

    //std::vector<pair<int, std::vector<int> > > neighborlistExclusions();
    bool readFromRestart();

};

#endif
