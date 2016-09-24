#pragma once
#ifndef FIXCROSSSTACK_H
#define FIXCROSSSTACK_H

#include "FixPotentialMultiAtom.h"
#include "CrossStack.h"
#include "CrossStackEvaluator.h"

void export_FixCrossStack3SPN2();

class FixCrossStack3SPN2 : public FixPotentialMultiAtom<CrossStackVariant, CrossStack3SPN2, CrossStack, CrossStackGPU, CrossStack3SPN2Type, 6> {

private:
    CrossStackEvaluator3SPN2 evaluator;

protected:
    float alpha;
    float range;
    
public:
    //DataSet *eng;
    //DataSet *press;

    FixCrossStack3SPN2(boost::shared_ptr<State> state_, std::string handle);

    void setParameters(float, float);
    void compute(bool);
    void singlePointEng(float *);
    bool prepareForRun();

    void createCrossStack(Atom *, Atom *, Atom *, Atom *, Atom *, Atom *, double, double, double, double, double, double,  int);
    void setCrossStackTypeCoefs(int, double, double, double, double, double, double);

    bool readFromRestart();

};

#endif
