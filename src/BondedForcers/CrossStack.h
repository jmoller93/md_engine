#pragma once
#ifndef CROSSSTACK_H
#define CROSSSTACK_H
#include "globalDefs.h"
#include "Atom.h"

#include "cutils_math.h"
#include <boost/variant.hpp>
#include <boost/functional/hash.hpp>
#include <array>
class CrossStack3SPN2;
void export_CrossStacks();
class CrossStack{
    public:
        //going to try storing by id instead.  Makes preparing for a run less intensive
        std::array<int, 6> ids;
        int type;
        void takeIds(CrossStack *);
};

class CrossStack3SPN2Type {
    public:
        //Sigma for the first cross stack interaction
        float sigma1;
        //Sigma for the second one
        float sigma2;
        //Energy scale
        float epsi;
        //Theta for 1st cross stack interaction
        float theta1;
        //Theta for 2nd cross stack interaction
        float theta2;
        //Just known as theta 3
        float theta3;
        CrossStack3SPN2Type(CrossStack3SPN2 *);
        CrossStack3SPN2Type(){};
        bool operator==(const CrossStack3SPN2Type &) const;
	std::string getInfoString();
};

class CrossStack3SPN2 : public CrossStack, public CrossStack3SPN2Type {
    public: 
        CrossStack3SPN2(Atom *a, Atom *b, Atom *c, Atom *d, Atom *e, Atom *f, double sigma1_, double sigma2_, double epsi_, double theta1_, double theta2_, double theta3_, int type_);
        CrossStack3SPN2(double sigma1_, double sigma2_, double epsi_, double theta1_, double theta2_, double theta3_, int type_);
        CrossStack3SPN2(){};
	std::string getInfoString();
};

class CrossStackGPU {
    public:
        int ids[6];
        uint32_t type;
        void takeIds(CrossStack *);


};

namespace std {
    template<> struct hash<CrossStack3SPN2Type> {
        size_t operator() (CrossStack3SPN2Type const& bp) const {
            size_t seed = 0;
            boost::hash_combine(seed, bp.sigma1);
            boost::hash_combine(seed, bp.sigma2);
            boost::hash_combine(seed, bp.epsi);
            boost::hash_combine(seed, bp.theta1); 
            boost::hash_combine(seed, bp.theta2); 
            boost::hash_combine(seed, bp.theta3); 
            return seed;
        }
    };


}

typedef boost::variant<
    CrossStack3SPN2,
    CrossStack	
> CrossStackVariant;
#endif
