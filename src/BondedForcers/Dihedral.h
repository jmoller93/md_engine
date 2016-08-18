#pragma once
#ifndef DIHEDRAL_H
#define DIHEDRAL_H
#include "globalDefs.h"
#include "Atom.h"

#include "cutils_math.h"
#include <boost/variant.hpp>
#include <boost/functional/hash.hpp>
#include <array>
class DihedralOPLS;
class DihedralGauss;
void export_Dihedrals();
class Dihedral{
    public:
        //going to try storing by id instead.  Makes preparing for a run less intensive
        std::array<int, 4> ids;
        int type;
        void takeIds(Dihedral *);
};

class DihedralOPLSType {
    public:
        float coefs[4];
        DihedralOPLSType(DihedralOPLS *);
        DihedralOPLSType(){};
        bool operator==(const DihedralOPLSType &) const;
	std::string getInfoString();
};



class DihedralOPLS : public Dihedral, public DihedralOPLSType {
    public:
        DihedralOPLS(Atom *a, Atom *b, Atom *c, Atom *d, double coefs_[4], int type_);
        DihedralOPLS(double coefs_[4], int type_);
        DihedralOPLS(){};
	std::string getInfoString();
};

class DihedralGPU {
    public:
        int ids[4];
        uint32_t type;
        void takeIds(Dihedral *);


};

//for forcer maps
namespace std {
    template<> struct hash<DihedralOPLSType> {
        size_t operator() (DihedralOPLSType const& dih) const {
            size_t seed = 0;
            for (int i=0; i<4; i++) {
                boost::hash_combine(seed, dih.coefs[i]);
            }
            return seed;
        }
    };


}

//Start 3SPN.2 Gaussian dihedrals
class DihedralGaussType {
    public:
        float phi0;
        float sigma;
        float k0;
        DihedralGaussType(DihedralGauss *);
        DihedralGaussType(){};
        bool operator==(const DihedralGaussType &) const;
	std::string getInfoString();
};

class DihedralGauss : public Dihedral, public DihedralGaussType {
    public: 
        DihedralGauss(Atom *a, Atom *b, Atom *c, Atom *d, double phi0_, double sigma_, double k0_, int type_);
        DihedralGauss(double phi0_, double sigma_, double k0_, int type_);
        DihedralGauss(){};
	std::string getInfoString();
};

namespace std {
    template<> struct hash<DihedralGaussType> {
        size_t operator() (DihedralGaussType const& dih) const {
            size_t seed = 0;
            boost::hash_combine(seed, dih.phi0);
            boost::hash_combine(seed, dih.sigma);
            boost::hash_combine(seed, dih.k0);
            return seed;
        }
    };


}

typedef boost::variant<
	DihedralOPLS,
    DihedralGauss,
    Dihedral	
> DihedralVariant;
#endif
