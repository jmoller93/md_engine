#include "BasePair.h"
#include "boost_for_export.h"
#include "array_indexing_suite.hpp"
namespace py = boost::python;

//3SPN.2 Base Pairs as a four atom bonded potential
BasePair::BasePair(Atom *a, Atom *b, Atom *c, Atom *d, double phi0_, double sigma_, double k_, double epsi = epsi_, double alpha = alpha_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    ids[2] = c->id;
    ids[3] = d->id;

    phi0 = phi0_;
    sigma = sigma_;
    k = k_;
    epsi = epsi_;
    alpha = alpha_;
    type = type_;
}

BasePair::BasePair(double phi0_, double sigma_, double k_, double epsi = epsi_, double alpha = alpha_, int type_) {
    for (int i=0; i<4; i++) {
        ids[i] = -1;
    }
    phi0 = phi0_;
    sigma = sigma_;
    k = k_;
    epsi = epsi_;
    alpha = alpha_;
    type = type_;
}

void BasePair::takeIds(BasePair *other) {
    for (int i=0; i<4; i++) {
        ids[i] = other->ids[i];
    }
}


void BasePairGPU::takeIds(BasePair *other) {
    for (int i=0; i<4; i++) {
        ids[i] = other->ids[i];
    }
}

BasePairType::BasePairType(BasePair *dihedral) {
    phi0  = dihedral->phi0;
    sigma = dihedral->sigma;
    k     = dihedral->k;
    epsi  = dihedral->epsi;
    alpha = dihedral->alpha;

}

std::string BasePair::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' atomID_a='" << ids[0] << "' atomID_b='" << ids[1] << "' atomID_c='" << ids[2] << "' atomID_d='" << ids[3] << "' phi0='" << phi0<< "'sigma='" << sigma << "' k='" << k << "' epsi='" << epsi << "' alpah='" << alpha << "'/>\n";
  return ss.str();
}

std::string BasePairType::getInfoString() {
  std::stringstream ss;
  ss << " phi0='" << phi0<< "' sigma='" << sigma << "' k='" << k << "' epsi='" << epsi << "' alpha='" << alpha;
  return ss.str();
}

bool BasePairType::operator==(const BasePairType &other) const {
    if (phi0 != other.phi0) {
        return false;
    }
    else if (sigma != other.sigma) {
        return false;
    }
    else if (k != other.k) {
        return false;
    }
    else if (epsi != other.epsi) {
        return false;
    }
    else if (alpha != other.alpha) {
        return false;
    }
    return true;
}

void export_BasePairs() {
    py::class_<BasePair, SHARED(BasePair)> ( "SimBasePair", py::init<>())
        .def_readwrite("type", &BasePair::type)
        .def_readonly("phi0", &BasePair::phi0)
        .def_readonly("sigma", &BasePair::sigma)
        .def_readonly("k", &BasePair::k)
        .def_readonly("epsi", &BasePair::epsi)
        .def_readonly("alpha", &BasePair::alpha)
        .def_readonly("ids", &BasePair::ids)

    ;

}
