#include "helpers.h"
#include "FixCrossStack.h"
#include "FixHelpers.h"
#include "cutils_func.h"
#include "CrossStackEvaluate.h"
namespace py = boost::python;
using namespace std;

const std::string crossstackType = "CrossStack3SPN2";


FixCrossStack3SPN2::FixCrossStack3SPN2(SHARED(State) state_, string handle) : FixPotentialMultiAtom (state_, handle, crossstackType, true){
      readFromRestart();
}

void FixCrossStack3SPN2::setParameters(float alpha_,float range_)
{
  alpha=alpha_;
  range=range_;
}

bool FixCrossStack3SPN2::prepareForRun() {
    evaluator = CrossStackEvaluator3SPN2(alpha, range);
    FixPotentialMultiAtom::prepareForRun();
    return true;
}


void FixCrossStack3SPN2::compute(bool computeVirials) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();
   
    GPUData &gpd = state->gpd;
    if (forcersGPU.size()) {
        if (computeVirials) {
            compute_force_crossstack<CrossStack3SPN2Type, CrossStackEvaluator3SPN2, true><<<NBLOCK(nAtoms), PERBLOCK, sharedMemSizeForParams>>>(forcersGPU.size(), gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), false, evaluator);
        } else {
            compute_force_crossstack<CrossStack3SPN2Type, CrossStackEvaluator3SPN2, false><<<NBLOCK(nAtoms), PERBLOCK, sharedMemSizeForParams >>>(forcersGPU.size(), gpd.xs(activeIdx), gpd.fs(activeIdx),gpd.idToIdxs.d_data.data(), forcersGPU.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), false, evaluator);
        }
    }
}

//Define specific parameters for all base pairing interactions
void FixCrossStack3SPN2::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();

    GPUData &gpd = state->gpd;
    compute_energy_crossstack<<<forcersGPU.size(), PERBLOCK, sharedMemSizeForParams>>>(nAtoms, gpd.xs(activeIdx), perParticleEng, gpd.idToIdxs.d_data.data(), forcersGPU.data(), state->boundsGPU, parameters.data(), parameters.size(), evaluator);
    

}



void FixCrossStack3SPN2::createCrossStack(Atom *a, Atom *b, Atom *c, Atom *d, Atom *e, Atom*f, double sigma1, double sigma2, double epsi, double theta1, double theta2, double theta3, int type) {
    if (type==-1) {
            assert(sigma1!=COEF_DEFAULT and sigma2!= COEF_DEFAULT and epsi!=COEF_DEFAULT and theta1!=COEF_DEFAULT and theta2!=COEF_DEFAULT and theta3!=COEF_DEFAULT);
    }
    forcers.push_back(CrossStack3SPN2(a, b, c, d, e, f, sigma1, sigma2, epsi, theta1, theta2, theta3, type));
    bonds.push_back(BondHarmonic(a,e, 0.0, 0.0, type));
    //bonds.push_back(BondHarmonic(a,b, 0.0, 0.0, type));
    bonds.push_back(BondHarmonic(b,e, 0.0, 0.0, type));
    bonds.push_back(BondHarmonic(d,f, 0.0, 0.0, type));
    //bonds.push_back(BondHarmonic(d,c, 0.0, 0.0, type));
    bonds.push_back(BondHarmonic(f,c, 0.0, 0.0, type));
    pyListInterface.updateAppendedMember();
}

void FixCrossStack3SPN2::setCrossStackTypeCoefs(int type, double sigma1, double sigma2, double epsi, double theta1, double theta2, double theta3) {
    assert(sigma1>0);
    assert(sigma2>0);
    CrossStack3SPN2 dummy(sigma1, sigma2, epsi, theta1, theta2, theta3, type);
    setForcerType(type, dummy);
}

bool FixCrossStack3SPN2::readFromRestart() {
  auto restData = getRestartNode();
  if(restData) {
      auto curr_node = restData.first_child();
      while (curr_node) {
        string tag = curr_node.name();
        if (tag == "types") {
          for (auto type_node = curr_node.first_child(); type_node; type_node = type_node.next_sibling()) {
            int type;
            double sigma1;
            double sigma2;
            double epsi;
            double theta1;
            double theta2;
            double theta3;
            std::string type_ = type_node.attribute("id").value();
            type = atoi(type_.c_str());
            std::string sigma1_ = type_node.attribute("sigma1").value();
            std::string sigma2_ = type_node.attribute("sigma2").value();
            std::string epsi_ = type_node.attribute("epsi").value();
            std::string theta1_ = type_node.attribute("theta1").value();
            std::string theta2_ = type_node.attribute("theta2").value();
            std::string theta3_ = type_node.attribute("theta3").value();
            sigma1 = atof(sigma1_.c_str());
            sigma2 = atof(sigma2_.c_str());
            epsi  = atof(epsi_.c_str());
            theta1= atof(theta1_.c_str());
            theta2= atof(theta2_.c_str());
            theta3= atof(theta3_.c_str());
            CrossStack3SPN2 dummy(sigma1, sigma2, epsi, theta1, theta2, theta3, type);
            setForcerType(type, dummy);
          }
        } else if (tag == "members") {
          for (auto member_node = curr_node.first_child(); member_node; member_node = member_node.next_sibling()) {
        int type;
        double sigma1;
        double sigma2;
        double epsi;
        double theta1;
        double theta2;
        double theta3;
        int ids[6];
        std::string type_ = member_node.attribute("type").value();
        std::string atom_a = member_node.attribute("atomID_a").value();
        std::string atom_b = member_node.attribute("atomID_b").value();
        std::string atom_c = member_node.attribute("atomID_c").value();
        std::string atom_d = member_node.attribute("atomID_d").value();
        std::string atom_e = member_node.attribute("atomID_e").value();
        std::string atom_f = member_node.attribute("atomID_f").value();
        std::string sigma1_ = member_node.attribute("sigma1").value();
        std::string sigma2_ = member_node.attribute("sigma2").value();
        std::string epsi_  = member_node.attribute("epsi").value();
        std::string theta1_ = member_node.attribute("theta1").value();
        std::string theta2_ = member_node.attribute("theta2").value();
        std::string theta3_ = member_node.attribute("theta3").value();
        type = atoi(type_.c_str());
        ids[0] = atoi(atom_a.c_str());
        ids[1] = atoi(atom_b.c_str());
        ids[2] = atoi(atom_c.c_str());
        ids[3] = atoi(atom_d.c_str());
        ids[4] = atoi(atom_e.c_str());
        ids[5] = atoi(atom_f.c_str());
        sigma1 = atof(sigma1_.c_str());
        sigma2 = atof(sigma2_.c_str());
        epsi  = atof(epsi_.c_str());
        theta1= atof(theta1_.c_str());
        theta2= atof(theta2_.c_str());
        theta3= atof(theta3_.c_str());
        Atom * a = &state->idToAtom(ids[0]);
        Atom * b = &state->idToAtom(ids[1]);
        Atom * c = &state->idToAtom(ids[2]);
        Atom * d = &state->idToAtom(ids[3]);
        Atom * e = &state->idToAtom(ids[4]);
        Atom * f = &state->idToAtom(ids[3]);
        if (a == NULL) {cout << "The first atom does not exist" <<endl; return false;};
        if (b == NULL) {cout << "The second atom does not exist" <<endl; return false;};
        if (c == NULL) {cout << "The third atom does not exist" <<endl; return false;};
        if (d == NULL) {cout << "The fourth atom does not exist" <<endl; return false;};
        if (e == NULL) {cout << "The fifth atom does not exist" <<endl; return false;};
        if (f == NULL) {cout << "The sixth atom does not exist" <<endl; return false;};
        createCrossStack(a, b, c, d, e, f, sigma1, sigma2, epsi, theta1, theta2, theta3, type);
          }
        }
        curr_node = curr_node.next_sibling();
    }
  }
  return true;
}


void export_FixCrossStack3SPN2() {
    py::class_<FixCrossStack3SPN2,
                          SHARED(FixCrossStack3SPN2),
                          py::bases<Fix, TypedItemHolder> > (
        "FixCrossStack3SPN2",
        py::init<SHARED(State), string> (
            py::args("state", "handle")
        )
    )
    .def("createCrossStack", &FixCrossStack3SPN2::createCrossStack,
            (py::arg("sigma1")=-1,
             py::arg("sigma2")=-1,
             py::arg("epsi")=-1,
             py::arg("theta1")=-1,
             py::arg("theta2")=-1,
             py::arg("theta3")=-1,
             py::arg("type")=-1)
        )

    .def("setCrossStackTypeCoefs", &FixCrossStack3SPN2::setCrossStackTypeCoefs, 
            (py::arg("type"), 
             py::arg("sigma1"),
             py::arg("sigma2"),
             py::arg("epsi"),
             py::arg("theta1"),
             py::arg("theta2"),
             py::arg("theta3"))
        )

    .def("setParameters", &FixCrossStack3SPN2::setParameters,
            (py::arg("alpha"), py::arg("range"))
        )
    
    .def_readonly("crossstacks", &FixCrossStack3SPN2::pyForcers)

    ;

}

