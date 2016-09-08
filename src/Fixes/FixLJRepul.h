#pragma once
#ifndef FIXLJREPUL_H
#define FIXLJREPUL_H

#include "FixPair.h"
#include "PairEvaluatorLJRepul.h"
#include "xml_func.h"
void export_FixLJRepul();

//! Fix for repulsive only Lennard-Jones
/*!
 * Fix to calculate only the repulsive part of the Lennard-Jones interactions of particles. The LJ potential
 * is defined as
 * \f[
 * V(r_{ij}) = \varepsilon \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} - 2
 *                               \left(\frac{\sigma}{r_{ij}}\right)^{6}\right] + \varepsilon,
 * \f]
 * where \f$ r \f$ is the distance between two particles and \f$ \varepsilon \f$
 * and \f$ \sigma \f$ are the two relevant parameters. The LJ pair interaction
 * is only calculated for particles closer than \$r_{\text{cut}}\$.
 *
 * From the potential, the force can be derived as
 * \f[
 * F(r_{ij}) = 12 \varepsilon \frac{1}{r} \left[ 
 *                            \left(\frac{\sigma}{r_{ij}}\right)^{12} -
 *                            \left(\frac{\sigma}{r_{ij}}\right)^{6}
 *                          \right].
 * \f]
 */
class FixLJRepul : public FixPair {
    public:
        //! Constructor
        FixLJRepul(SHARED(State), std::string handle);

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
        const std::string epsHandle; //!< Handle for parameter epsilon
        const std::string sigHandle; //!< Handle for parameter sigma
        const std::string rCutHandle; //!< Handle for parameter rCut
        std::vector<float> epsilons; //!< vector storing epsilon values
        std::vector<float> sigmas; //!< vector storing sigma values
        std::vector<float> rCuts; //!< vector storing cutoff distance values

        EvaluatorLJRepul evaluator; //!< Evaluator for generic pair interactions
};

#endif
