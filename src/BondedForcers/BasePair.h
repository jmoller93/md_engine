#pragma once
#ifndef BASEPAIR_H
#define BASEPAIR_H
#include "globalDefs.h"
#include "Atom.h"

#include "cutils_math.h"
#include <boost/variant.hpp>
#include <boost/functional/hash.hpp>
#include <array>
class BasePair3SPN2;
void export_BasePairs();

class BasePair{
    public:
        //going to try storing by id instead.  Makes preparing for a run less intensive
        std::array<int, 4> ids;
        int type;
        void takeIds(BasePair *);
};

class BasePairGPU {
    public:
        int ids[4];
        uint32_t type;
        void takeIds(BasePair *);


};

class BasePair3SPN2Type {
    public:
        float phi0;
        float sigma;
        float epsi;
        float theta1;
        float theta2;
        BasePair3SPN2Type(BasePair3SPN2 *);
        BasePair3SPN2Type(){};
        bool operator==(const BasePair3SPN2Type &) const;
	std::string getInfoString();
};

class BasePair3SPN2 : public BasePair, public BasePair3SPN2Type {
    public: 
        BasePair3SPN2(Atom *a, Atom *b, Atom *c, Atom *d, double phi0_, double sigma_, double epsi_, double theta1_, double theta2_, int type_);
        BasePair3SPN2(double phi0_, double sigma_, double epsi_, double theta1_, double theta2_, int type_);
        BasePair3SPN2(){};
	std::string getInfoString();
};


namespace std {
    template<> struct hash<BasePair3SPN2Type> {
        size_t operator() (BasePair3SPN2Type const& bp) const {
            size_t seed = 0;
            boost::hash_combine(seed, bp.phi0);
            boost::hash_combine(seed, bp.sigma);
            boost::hash_combine(seed, bp.epsi);
            boost::hash_combine(seed, bp.theta1); 
            boost::hash_combine(seed, bp.theta2); 
            return seed;
        }
    };


}

typedef boost::variant<
    BasePair3SPN2,
    BasePair	
> BasePairVariant;
#endif
