#include "BasePair.h"
#include "boost_for_export.h"
#include "array_indexing_suite.hpp"
namespace py = boost::python;

//3SPN.2 Base Pairs as a four atom bonded potential
BasePair::BasePair(Atom *a, Atom *b, Atom *c, Atom *d, double phi0_, double sigma_, double k_, double epsi_, double alpha_, double theta1_, double theta2_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    ids[2] = c->id;
    ids[3] = d->id;

    phi0 = phi0_;
    sigma = sigma_;
    k = k_;
    epsi = epsi_;
    alpha = alpha_;
    theta1 = theta1_;
    theta2 = theta2_;
    type = type_;
}

BasePair::BasePair(double phi0_, double sigma_, double k_, double epsi_, double alpha_, double theta1_, double theta2_, int type_) {
    for (int i=0; i<4; i++) {
        ids[i] = -1;
    }
    phi0 = phi0_;
    sigma = sigma_;
    k = k_;
    epsi = epsi_;
    alpha = alpha_;
    theta1 = theta1_;
    theta2 = theta2_;
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

BasePairType::BasePairType(BasePair *basepair) {
    phi0  = basepair->phi0;
    sigma = basepair->sigma;
    k     = basepair->k;
    epsi  = basepair->epsi;
    alpha = basepair->alpha;
    theta1 = basepair->theta1;
    theta2 = basepair->theta2;

}

std::string BasePair::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' atomID_a='" << ids[0] << "' atomID_b='" << ids[1] << "' atomID_c='" << ids[2] << "' atomID_d='" << ids[3] << "' phi0='" << phi0<< "'sigma='" << sigma << "' k='" << k << "' epsi='" << epsi << "' alpha='" << alpha << "' theta1='" << theta1 << "' theta2='" << theta2 <<  "'/>\n";
  return ss.str();
}

std::string BasePairType::getInfoString() {
  std::stringstream ss;
  ss << " phi0='" << phi0<< "' sigma='" << sigma << "' k='" << k << "' epsi='" << epsi << "' alpha='" << alpha << "' theta1='" << theta1 << "' theta2='" << theta2;
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
    else if (theta1 != other.theta1) {
        return false;
    }
    else if (theta2 != other.theta2) {
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
        .def_readonly("theta1", &BasePair::theta1)
        .def_readonly("theta2", &BasePair::theta2)
        .def_readonly("ids", &BasePair::ids)

    ;

}
