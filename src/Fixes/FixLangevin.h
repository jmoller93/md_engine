#pragma once

#ifndef FIXLANGEVIN_H
#define FIXLANGEVIN_H

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#include <string>

#include "FixThermostatBase.h"
#include "Fix.h"
#include "GPUArrayDeviceGlobal.h"

#include <curand_kernel.h>

void export_FixLangevin();

class State;

class FixLangevin : public FixThermostatBase, public Fix {
private:
    int seed;
    float gamma;
    void setDefaults();
    GPUArrayDeviceGlobal<curandState_t> randStates;
public:

    FixLangevin(boost::shared_ptr<State> state_, std::string handle_, std::string groupHandle_, double temp_);
    FixLangevin(boost::shared_ptr<State> state_, std::string handle_, std::string groupHandle_, boost::python::list, boost::python::list);
    FixLangevin(boost::shared_ptr<State> state_, std::string handle_, std::string groupHandle_, boost::python::object);
    bool prepareForRun();
    void compute(bool);
    bool postRun();
    void setParams(double seed, double gamma);
};



#endif
