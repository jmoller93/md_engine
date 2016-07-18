#pragma once
#ifndef DATACOMPUTERTEMPERATURE_H
#define DATACOMPUTERTEMPERATURE_H

#include "DataComputer.h"
#include "GPUArrayGlobal.h"
#include "Virial.h"
namespace MD_ENGINE {
    class State;

    class DataComputerTemperature : public DataComputer {
        public:
            void appendScalar(boost::python::list);
            void appendTensor(boost::python::list);

            void computeScalar_GPU(bool, uint32_t);
            void computeTensor_GPU(bool, uint32_t);

            void computeScalar_CPU();
            void computeTensor_CPU();

            DataComputerTemperature(State *, bool, bool);
            void prepareForRun();
            double tempScalar;
            Virial tempTensor;
            //so these are just length 2 arrays.  First value is used for the result of the sum.  Second value is bit-cast to an int and used to cound how many values are present.
            GPUArrayGlobal<float> tempGPUScalar;
            GPUArrayGlobal<Virial> tempGPUTensor;

            void appendScalar(boost::python::list &);
            void appendTensor(boost::python::list &);

            double getScalar();
            Virial getTensor();

    };
}

#endif