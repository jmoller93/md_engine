#pragma once
#ifndef FIXDIHEDRALPERIODIC_H
#define FIXDIHEDRALPERIODIC_H

#include <boost/python.hpp>

#include "FixPotentialMultiAtom.h"
#include "Dihedral.h"
#include "DihedralEvaluatorPeriodic.h"

void export_FixDihedralPeriodic();

class FixDihedralPeriodic : public FixPotentialMultiAtom<DihedralVariant, DihedralPeriodic, Dihedral, DihedralGPU, DihedralPeriodicType, 4> {

private:
    DihedralEvaluatorPeriodic evaluator;
public:
    //DataSet *eng;
    //DataSet *press;
    std::vector<BondVariant> bonds;

    FixDihedralPeriodic(boost::shared_ptr<State> state_, std::string handle);

    void compute(int);
    void singlePointEng(float *);

    void createDihedral(Atom *, Atom *, Atom *, Atom *, double, double, double, int);
    void createDihedralPy(Atom *, Atom *, Atom *, Atom *, boost::python::list, double, int);
    void setDihedralTypeCoefs(int, boost::python::list, double);

    //std::vector<pair<int, std::vector<int> > > neighborlistExclusions();
    bool readFromRestart();

    virtual std::vector<BondVariant> *getBonds() {
        return &bonds;
    }

};

#endif
