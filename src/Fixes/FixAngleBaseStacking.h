#pragma once
#ifndef FIXANGLEBASESTACKING_H
#define FIXANGLEBASESTACKING_H

#include "FixPotentialMultiAtom.h"
#include "Angle.h"
#include "AngleEvaluatorBaseStacking.h"

void export_FixAngleBaseStacking();

class FixAngleBaseStacking : public FixPotentialMultiAtom<AngleVariant, AngleBaseStacking, Angle, AngleGPU, AngleBaseStackingType, 3> {

private:
    AngleEvaluatorBaseStacking evaluator; 
public:
    //DataSet *eng;
    //DataSet *press;

    FixAngleBaseStacking(boost::shared_ptr<State> state_, std::string handle);

    void compute(bool);
    void singlePointEng(float *);

    void createAngle(Atom *, Atom *, Atom *, double, double, double, double, double, int type_);
    void setAngleTypeCoefs(int, double, double, double, double, double);

    bool readFromRestart(pugi::xml_node restData);

};

#endif
