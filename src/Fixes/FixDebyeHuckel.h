#pragma once
#ifndef FIXDEBYEHUCKEL_H
#define FIXDEBYEHUCKEL_H

#include "FixPair.h"
#include "PairEvaluatorDH.h"
#include "xml_func.h"
void export_FixDH();

//! Fix for Debye-Huckel pairwise interactions

class FixDH : public FixPair {
    public:
        //! Constructor
        FixDH(SHARED(State), std::string handle);

        //! Compute forces
        void compute(bool);

        //! Compute single point energy
        void singlePointEng(float *);

        //! Prepare Fix
        /*!
         * \returns Always returns True
         *
         * This function needs to be called before simulation run.
         */
        bool prepareForRun();

        //! Run after simulation
        /*!
         * This function needs to be called after simulation run.
         */
        bool postRun();

        //! Create restart string
        /*!
         * \param format Format of the pair parameters.
         *
         * \returns restart chunk string.
         */
        std::string restartChunk(std::string format);

        //! Read parameters from restart file
        /*!
         * \return Always True
         *
         * \param restData XML node containing the restart data.
         */
        bool readFromRestart(pugi::xml_node restData);

        //! Add new type of atoms
        /*!
         * \param handle Not used
         *
         * This function adds a new particle type to the fix.
         */
        void addSpecies(std::string handle);

        //! Return list of cutoff values
        std::vector<float> getRCuts();

    public:
        const std::string 


