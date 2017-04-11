#include "Bond.h"
#include "boost_for_export.h"

namespace py = boost::python;

//BondHarmonicType::BondHarmonicType(BondHarmonic *bond) {
//    k = bond->k;
//    r0 = bond->r0;
//}

bool BondHarmonicType::operator==(const BondHarmonicType &other) const {
    return k == other.k and r0 == other.r0;
}





BondHarmonic::BondHarmonic(Atom *a, Atom *b, double k_, double r0_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    k = k_;
    r0 = r0_;
    type = type_;
}
BondHarmonic::BondHarmonic(double k_, double r0_, int type_) {
    k = k_;
    r0 = r0_;
    type = type_;
}

void BondGPU::takeIds(Bond *b) { 
    myId = b->ids[0];
    otherId = b->ids[1];
}

std::string BondHarmonicType::getInfoString() {
  std::stringstream ss;
  ss << " k='" << k << "' r0='" << r0;
  return ss.str();
}

std::string BondHarmonic::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' k='" << k << "' r0='" << r0 << "' atomID_a='" << ids[0] <<  "' atomID_b='" << ids[1] << "'/>\n";
  return ss.str();
}

void export_BondHarmonic() {
  
    boost::python::class_<BondHarmonic,SHARED(BondHarmonic)> ( "BondHarmonic", boost::python::init<>())
//         .def(boost::python::init<int, int ,double, double,int>())
        .def_readonly("ids", &BondHarmonic::ids)
        .def_readwrite("k", &BondHarmonic::k)
        .def_readwrite("r0", &BondHarmonic::r0)
    ;
}




//extended harmonic bond (3SPN.2 version)
bool BondHarmonicExtendType::operator==(const BondHarmonicExtendType &other) const {
    return k1 == other.k1 and k2 == other.k2 and r0 == other.r0;
}





BondHarmonicExtend::BondHarmonicExtend(Atom *a, Atom *b, double k1_, double k2_, double r0_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    k1 = k1_;
    k2 = k2_;
    r0 = r0_;
    type = type_;
}
BondHarmonicExtend::BondHarmonicExtend(double k1_, double k2_, double r0_, int type_) {
    k1 = k1_;
    k2 = k2_;
    r0 = r0_;
    type = type_;
}


std::string BondHarmonicExtendType::getInfoString() {
  std::stringstream ss;
  ss << "' k1='" << k1 << "' k2='" << k2 << "' r0='" << r0;
  return ss.str();
}

std::string BondHarmonicExtend::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' k1='" <<  k1 << "' k2='" << k2 << "' r0='" << r0 << "' atomID_a='" << ids[0] <<  "' atomID_b='" << ids[1] << "'/>\n";
  return ss.str();
}

void export_BondHarmonicExtend() {
  
    boost::python::class_<BondHarmonicExtend,SHARED(BondHarmonicExtend)> ( "BondHarmonicExtend", boost::python::init<>())
//         .def(boost::python::init<int, int ,double, double,int>())
        .def_readonly("ids", &BondHarmonicExtend::ids)
        .def_readwrite("k1", &BondHarmonicExtend::k1)
        .def_readwrite("k2", &BondHarmonicExtend::k2)
        .def_readwrite("r0", &BondHarmonicExtend::r0)
    ;
}




//bond FENE
bool BondFENEType::operator==(const BondFENEType &other) const {
    return k == other.k and r0 == other.r0 and eps == other.eps and sig == other.sig;
}





BondFENE::BondFENE(Atom *a, Atom *b, double k_, double r0_, double eps_, double sig_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    k = k_;
    r0 = r0_;
    eps = eps_;
    sig = sig_;
    type = type_;
}
BondFENE::BondFENE(double k_, double r0_, double eps_, double sig_, int type_) {
    k = k_;
    r0 = r0_;
    eps = eps_;
    sig = sig_;
    type = type_;
}

std::string BondFENEType::getInfoString() {
  std::stringstream ss;
  ss << " k='" << k << "' r0='" << r0 << "' eps='" << eps << "' sig='" << sig;
  return ss.str();
}

std::string BondFENE::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' k='" << k << "' r0='" << r0 << "' eps='" << eps << "' sig='" << sig << "' atomID_a='" << ids[0] <<  "' atomID_b='" << ids[1] << "'/>\n";
  return ss.str();
}

void export_BondFENE() {
  
    py::class_<BondFENE,SHARED(BondFENE)> ( "BondFENE", py::init<>())
//         .def(py::init<int, int ,double, double,int>())
        .def_readonly("ids", &BondFENE::ids)
        .def_readwrite("k", &BondFENE::k)
        .def_readwrite("r0", &BondFENE::r0)
        .def_readwrite("eps", &BondFENE::eps)
        .def_readwrite("sig", &BondFENE::sig)
    ;
}

//bond GoLike
bool BondGoLikeType::operator==(const BondGoLikeType &other) const {
    return eps == other.eps and sig == other.sig;
}





BondGoLike::BondGoLike(Atom *a, Atom *b, double eps_, double sig_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    eps = eps_;
    sig = sig_;
    type = type_;
}
BondGoLike::BondGoLike(double eps_, double sig_, int type_) {
    eps = eps_;
    sig = sig_;
    type = type_;
}

std::string BondGoLikeType::getInfoString() {
  std::stringstream ss;
  ss << "' eps='" << eps << "' sig='" << sig;
  return ss.str();
}

std::string BondGoLike::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' eps='" << eps << "' sig='" << sig << "' atomID_a='" << ids[0] <<  "' atomID_b='" << ids[1] << "'/>\n";
  return ss.str();
}

void export_BondGoLike() {
  
    py::class_<BondGoLike,SHARED(BondGoLike)> ( "BondGoLike", py::init<>())
//         .def(py::init<int, int ,double, double,int>())
        .def_readonly("ids", &BondGoLike::ids)
        .def_readwrite("eps", &BondGoLike::eps)
        .def_readwrite("sig", &BondGoLike::sig)
    ;
}
