#pragma once
#ifndef BASEPAIR_H
#define BASEPAIR_H
#include "globalDefs.h"
#include "Atom.h"

#include "cutils_math.h"
#include <boost/variant.hpp>
#include <boost/functional/hash.hpp>
#include <array>
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

class BasePairType {
    public:
        float phi0;
        float sigma;
        float k;
        float epsi;
        float alpha;
        float theta1;
        float theta2;
        BasePairType(BasePair *);
        BasePairType(){};
        bool operator==(const BasePairType &) const;
	std::string getInfoString();
};

class BasePair : public BasePair, public BasePairType {
    public: 
        BasePair(Atom *a, Atom *b, Atom *c, Atom *d, double phi0_, double sigma_, double k_, double epsi_, double alpha_, double theta1_, double theta2_, int type_);
        BasePair(double phi0_, double sigma_, double k_, double epsi_, double alpha_, double theta1_, double theta2_, int type_);
        BasePair(){};
	std::string getInfoString();
};

namespace std {
    template<> struct hash<BasePairType> {
        size_t operator() (BasePairType const& bp) const {
            size_t seed = 0;
            boost::hash_combine(seed, bp.phi0);
            boost::hash_combine(seed, bp.sigma);
            boost::hash_combine(seed, bp.k);
            boost::hash_combine(seed, bp.epsi);
            boost::hash_combine(seed, bp.alpha); 
            boost::hash_combine(seed, bp.theta1); 
            boost::hash_combine(seed, bp.theta2); 
            return seed;
        }
    };


}

typedef boost::variant<
    BasePair	
> BasePairVariant;
#endif
