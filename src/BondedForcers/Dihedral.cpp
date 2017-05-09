#include "Dihedral.h"
#include <boost/python.hpp>
#include "boost_for_export.h"
#include "array_indexing_suite.hpp"
namespace py = boost::python;

void Dihedral::takeIds(Dihedral *other) {
    for (int i=0; i<4; i++) {
        ids[i] = other->ids[i];
    }
}


void DihedralGPU::takeIds(Dihedral *other) {
    for (int i=0; i<4; i++) {
        ids[i] = other->ids[i];
    }
}


DihedralOPLS::DihedralOPLS(Atom *a, Atom *b, Atom *c, Atom *d, double coefs_[4], int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    ids[2] = c->id;
    ids[3] = d->id;
    for (int i=0; i<4; i++) {
        coefs[i] = coefs_[i];
    }
    type = type_;
}

DihedralOPLS::DihedralOPLS(double coefs_[4], int type_) {
    for (int i=0; i<4; i++) {
        ids[i] = -1;
    }
    for (int i=0; i<4; i++) {
        coefs[i] = coefs_[i];
    }
    type = type_;
}

DihedralOPLSType::DihedralOPLSType(DihedralOPLS *dihedral) {
    for (int i=0; i<4; i++) {
        coefs[i] = dihedral->coefs[i];
    }
}

DihedralPeriodic::DihedralPeriodic(Atom *a, Atom *b, Atom *c, Atom *d, double coefs_[2], double phiRef_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    ids[2] = c->id;
    ids[3] = d->id;
    for (int i=0; i<2; i++) {
        coefs[i] = coefs_[i];
    }
    phiRef = phiRef_;
    type = type_;
}

DihedralPeriodic::DihedralPeriodic(double coefs_[2], double phiRef_, int type_) {
    for (int i=0; i<2; i++) {
        ids[i] = -1;
    }
    for (int i=0; i<2; i++) {
        coefs[i] = coefs_[i];
    }
    phiRef = phiRef_;
    type = type_;
}
DihedralCHARMM::DihedralCHARMM(Atom *atomA, Atom *atomB, Atom *atomC, Atom *atomD, double k_, int n_, double d_,  int type_) {
    ids[0] = atomA->id;
    ids[1] = atomB->id;
    ids[2] = atomC->id;
    ids[3] = atomD->id;
    k = k_;
    n = n_;
    d = d_;
    type = type_;
}

DihedralCHARMM::DihedralCHARMM(double k_, int n_, double d_, int type_) {
    for (int i=0; i<4; i++) {
        ids[i] = -1;
    }
    k = k_;
    n = n_;
    d = d_;
    type = type_;
}


//3SPN.2 Gaussian dihedrals
DihedralGauss::DihedralGauss(Atom *a, Atom *b, Atom *c, Atom *d, double phi0_, double sigma_, double k0_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    ids[2] = c->id;
    ids[3] = d->id;

    phi0 = phi0_;
    sigma = sigma_;
    k0 = k0_;
    type = type_;
}

DihedralGauss::DihedralGauss(double phi0_, double sigma_, double k0_, int type_) {
    for (int i=0; i<4; i++) {
        ids[i] = -1;
    }
    phi0 = phi0_;
    sigma = sigma_;
    k0 = k0_;
    type = type_;
}


DihedralPeriodicType::DihedralPeriodicType(DihedralPeriodic *dihedral) {
    for (int i=0; i<2; i++) {
        coefs[i] = dihedral->coefs[i];
    }
    phiRef = dihedral->phiRef;
}
DihedralCHARMMType::DihedralCHARMMType(DihedralCHARMM *dihedral) {
    k = dihedral->k;
    n = dihedral->n;
    d = dihedral->d;
}


DihedralGaussType::DihedralGaussType(DihedralGauss *dihedral) {
    phi0  = dihedral->phi0;
    sigma = dihedral->sigma;
    k0    = dihedral->k0;
}

std::string DihedralOPLS::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' atomID_a='" << ids[0] << "' atomID_b='" << ids[1] << "' atomID_c='" << ids[2] << "' atomID_d='" << ids[3] << "' coef_a='" << coefs[0]<< "' coef_b='" << coefs[1] << "' coef_c='" << coefs[2] << "' coef_d='" << coefs[3] << "'/>\n";
  return ss.str();
}

std::string DihedralOPLSType::getInfoString() {
    std::stringstream ss;
    ss << " coef_a='" << coefs[0]<< "' coef_b='" << coefs[1] << "' coef_c='" << coefs[2] << "' coef_d='" << coefs[3];
    return ss.str();
}

std::string DihedralPeriodic::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' atomID_a='" << ids[0] << "' atomID_b='" << ids[1] << "' atomID_c='" << ids[2] << "' atomID_d='" << ids[3] << "' coef_a='" << coefs[0]<< "' coef_b='" << coefs[1] << "' phi_ref='"<< phiRef << "'/>\n";
  return ss.str();
}

std::string DihedralPeriodicType::getInfoString() {
    std::stringstream ss;
    ss << " coef_a='" << coefs[0]<< "' coef_b='" << coefs[1] << "' phi_ref='" << phiRef;
    return ss.str();
}

std::string DihedralCHARMM::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' atomID_a='" << ids[0] << "' atomID_b='" << ids[1] << "' atomID_c='" << ids[2] << "' atomID_d='" << ids[3] << " k='" << k << "' n='" << n << "' d='" << d << "'/>\n";
  return ss.str();
}

std::string DihedralGauss::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' atomID_a='" << ids[0] << "' atomID_b='" << ids[1] << "' atomID_c='" << ids[2] << "' atomID_d='" << ids[3] << "' phi0='" << phi0<< "'sigma='" << sigma << "' k0='" << k0 << "'/>\n";
  return ss.str();
}

std::string DihedralGaussType::getInfoString() {
  std::stringstream ss;
  ss << " phi0='" << phi0<< "' sigma='" << sigma << "' k0='" << k0;
  return ss.str();
}

std::string DihedralCHARMMType::getInfoString() {
  std::stringstream ss;
  ss << " k='" << k << "' n='" << n << "' d='" << d;
  return ss.str();
}

bool DihedralOPLSType::operator==(const DihedralOPLSType &other) const {
    for (int i=0; i<4; i++) {
        if (coefs[i] != other.coefs[i]) {
            return false;
        }
    }
    return true;
}

bool DihedralPeriodicType::operator==(const DihedralPeriodicType &other) const {
    for (int i=0; i<2; i++) {
        if (coefs[i] != other.coefs[i]) {
            return false;
        }
    }
    if (phiRef != other.phiRef) {
        return false;
    }
    return true;
}

bool DihedralGaussType::operator==(const DihedralGaussType &other) const {
    if (phi0 != other.phi0) {
        return false;
    }
    else if (sigma != other.sigma) {
        return false;
    }
    else if (k0 != other.k0) {
        return false;
    }
    return true;
}

bool DihedralCHARMMType::operator==(const DihedralCHARMMType &other) const {
    return other.k == k and other.d == d and other.n == n;
}

void export_Dihedrals() {
    py::class_<DihedralOPLS, SHARED(DihedralOPLS)> ( "SimDihedralOPLS", py::init<>())
        .def_readwrite("type", &DihedralOPLS::type)
        .def_readonly("coefs", &DihedralOPLS::coefs)
        .def_readonly("ids", &DihedralOPLS::ids)

    ;
    py::class_<DihedralPeriodic, SHARED(DihedralPeriodic)> ( "SimDihedralPeriodic", py::init<>())
        .def_readwrite("type", &DihedralPeriodic::type)
        .def_readonly("coefs", &DihedralPeriodic::coefs)
        .def_readonly("phiRef", &DihedralPeriodic::phiRef)
        .def_readonly("ids", &DihedralPeriodic::ids)

    ;
    py::class_<DihedralCHARMM, SHARED(DihedralCHARMM)> ( "SimDihedralCHARMM", py::init<>())
        .def_readwrite("type", &DihedralCHARMM::type)
        .def_readwrite("k", &DihedralCHARMM::k)
        .def_readwrite("n", &DihedralCHARMM::n)
        .def_readwrite("d", &DihedralCHARMM::d)
        .def_readonly("ids", &DihedralCHARMM::ids)

    ;

    py::class_<DihedralGauss, SHARED(DihedralGauss)> ( "SimDihedralGauss", py::init<>())
        .def_readwrite("type", &DihedralGauss::type)
        .def_readwrite("phi0", &DihedralGauss::phi0)
        .def_readwrite("sigma", &DihedralGauss::sigma)
        .def_readwrite("k0", &DihedralGauss::k0)
        .def_readonly("ids", &DihedralGauss::ids)

    ;
}


