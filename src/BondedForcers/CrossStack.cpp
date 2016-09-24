#include "CrossStack.h"
#include "boost_for_export.h"
#include "array_indexing_suite.hpp"
namespace py = boost::python;

//3SPN.2 Base Pairs as a four atom bonded potential
CrossStack3SPN2::CrossStack3SPN2(Atom *a, Atom *b, Atom *c, Atom *d, Atom *e, Atom *f, double sigma1_, double sigma2_, double epsi_, double theta1_, double theta2_, double theta3_, int type_) {
    ids[0] = a->id;
    ids[1] = b->id;
    ids[2] = c->id;
    ids[3] = d->id;
    ids[4] = e->id;
    ids[5] = f->id;

    sigma = sigma1_;
    sigma = sigma2_;
    epsi = epsi_;
    theta1 = theta1_;
    theta2 = theta2_;
    theta3 = theta3_;
    type = type_;
}

CrossStack3SPN2::CrossStack3SPN2(double sigma1_, double sigma2_,  double epsi_,  double theta1_, double theta2_, double theta3_, int type_) {
    for (int i=0; i<6; i++) {
        ids[i] = -1;
    }
    sigma1 = sigma1_;
    sigma2 = sigma2_;
    epsi = epsi_;
    theta1 = theta1_;
    theta2 = theta2_;
    theta3 = theta3_;
    type = type_;
}

void CrossStack::takeIds(CrossStack *other) {
    for (int i=0; i<6; i++) {
        ids[i] = other->ids[i];
    }
}


void CrossStackGPU::takeIds(CrossStack *other) {
    for (int i=0; i<6; i++) {
        ids[i] = other->ids[i];
    }
}

CrossStack3SPN2Type::CrossStack3SPN2Type(CrossStack3SPN2 *crossstack) {
    sigma1 = crossstack->sigma1;
    sigma2 = crossstack->sigma2;
    epsi  = crossstack->epsi;
    theta1 = crossstack->theta1;
    theta2 = crossstack->theta2;
    theta3 = crossstack->theta3;

}

std::string CrossStack3SPN2::getInfoString() {
  std::stringstream ss;
  ss << "<member type='" << type << "' atomID_a='" << ids[0] << "' atomID_b='" << ids[1] << "' atomID_c='" << ids[2] << "' atomID_d='" << ids[3] << "' atomID_e='" << ids[4] << "' atomID_f='" << ids[5] << "'sigma1='" << sigma1 << "'sigma2='" << sigma2 <<"' epsi='" << epsi << "' theta1='" << theta1 << "' theta2='" << theta2 << "' theta3='" << theta3 << "'/>\n";
  return ss.str();
}

std::string CrossStack3SPN2Type::getInfoString() {
  std::stringstream ss;
  ss << "' sigma1='" << sigma1 << "'sigma2='" << sigma2 <<"' epsi='" << epsi << "' theta1='" << theta1 << "' theta2='" << theta2 << "' theta3='" << theta3;
  return ss.str();
}

bool CrossStack3SPN2Type::operator==(const CrossStack3SPN2Type &other) const {
    if (sigma1 != other.sigma1) {
        return false;
    }
    else if (sigma2 != other.sigma2) {
        return false;
    }
    else if (epsi != other.epsi) {
        return false;
    }
    else if (theta1 != other.theta1) {
        return false;
    }
    else if (theta2 != other.theta2) {
        return false;
    }
    else if (theta3 != other.theta3) {
        return false;
    }
    return true;
}

void export_CrossStacks() {
    py::class_<CrossStack3SPN2, SHARED(CrossStack3SPN2)> ( "SimCrossStack", py::init<>())
        .def_readwrite("type", &CrossStack3SPN2::type)
        .def_readonly("sigma1", &CrossStack3SPN2::sigma1)
        .def_readonly("sigma2", &CrossStack3SPN2::sigma2)
        .def_readonly("epsi", &CrossStack3SPN2::epsi)
        .def_readonly("theta1", &CrossStack3SPN2::theta1)
        .def_readonly("theta2", &CrossStack3SPN2::theta2)
        .def_readonly("theta3", &CrossStack3SPN2::theta2)
        .def_readonly("ids", &CrossStack3SPN2::ids)

    ;

}
