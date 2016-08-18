#include "Dihedral.h"
#include "boost_for_export.h"
#include "array_indexing_suite.hpp"
namespace py = boost::python;
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

DihedralOPLSType::DihedralOPLSType(DihedralOPLS *dihedral) {
    for (int i=0; i<4; i++) {
        coefs[i] = dihedral->coefs[i];
    }
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

bool DihedralOPLSType::operator==(const DihedralOPLSType &other) const {
    for (int i=0; i<4; i++) {
        if (coefs[i] != other.coefs[i]) {
            return false;
        }
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

void export_Dihedrals() {
    py::class_<DihedralOPLS, SHARED(DihedralOPLS)> ( "SimDihedralOPLS", py::init<>())
        .def_readwrite("type", &DihedralOPLS::type)
        .def_readonly("coefs", &DihedralOPLS::coefs)
        .def_readonly("ids", &DihedralOPLS::ids)

    ;

    py::class_<DihedralGauss, SHARED(DihedralGauss)> ( "SimDihedralGauss", py::init<>())
        .def_readwrite("type", &DihedralGauss::type)
        .def_readonly("phi0", &DihedralGauss::phi0)
        .def_readonly("sigma", &DihedralGauss::sigma)
        .def_readonly("k0", &DihedralGauss::k0)
        .def_readonly("ids", &DihedralGauss::ids)

    ;

}
