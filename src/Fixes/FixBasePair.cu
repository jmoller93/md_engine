#include "helpers.h"
#include "FixBasePair.h"
#include "FixHelpers.h"
#include "cutils_func.h"
#include "BasePairEvaluate.h"
namespace py = boost::python;
using namespace std;

const std::string basepairType = "BasePair3SPN2";


FixBasePair3SPN2::FixBasePair3SPN2(SHARED(State) state_, string handle) : FixPotentialMultiAtom (state_, handle, basepairType, true){
  if (state->readConfig->fileOpen) {
    auto restData = state->readConfig->readFix(type, handle);
    if (restData) {
      std::cout << "Reading restart data for fix " << handle << std::endl;
      readFromRestart(restData);
    }
  }
}


void FixBasePair3SPN2::compute(bool computeVirials) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();


    GPUData &gpd = state->gpd;
    if (computeVirials) {
        compute_force_basepair<BasePair3SPN2Type, BasePairEvaluator3SPN2, true><<<NBLOCK(nAtoms), PERBLOCK, sizeof(BasePairGPU) * maxForcersPerBlock + sizeof(BasePair3SPN2Type) * parameters.size() >>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), evaluator);
    } else {
        compute_force_basepair<BasePair3SPN2Type, BasePairEvaluator3SPN2, false><<<NBLOCK(nAtoms), PERBLOCK, sizeof(BasePairGPU) * maxForcersPerBlock + sizeof(BasePair3SPN2Type) * parameters.size() >>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), evaluator);
    }

}

void FixBasePair3SPN2::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();

    GPUData &gpd = state->gpd;
    compute_energy_basepair<<<NBLOCK(nAtoms), PERBLOCK, sizeof(BasePairGPU) * maxForcersPerBlock + sizeof(BasePair3SPN2Type) * parameters.size() >>>(nAtoms, gpd.xs(activeIdx), perParticleEng, gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), evaluator);
    

}



void FixBasePair3SPN2::createBasePair(Atom *a, Atom *b, Atom *c, Atom *d, double phi0, double sigma, double k, double epsi, double alpha, double theta1, double theta2, int type) {
    if (type==-1) {
            assert(phi0!=COEF_DEFAULT and sigma!=COEF_DEFAULT and k!=COEF_DEFAULT and epsi!=COEF_DEFAULT and alpha!=COEF_DEFAULT and theta1!=COEF_DEFAULT and theta2!=COEF_DEFAULT);
    }
    forcers.push_back(BasePair3SPN2(a, b, c, d, phi0, sigma, k, epsi, alpha, type, theta1, theta2));
    pyListInterface.updateAppendedMember();
}

void FixBasePair3SPN2::setBasePairTypeCoefs(int type, double phi0, double sigma, double k, double epsi, double alpha, double theta1, double theta2) {
    assert(sigma>0);
    BasePair3SPN2 dummy(phi0, sigma, k, epsi, alpha, theta1, theta2, type);
    setForcerType(type, dummy);
}

bool FixBasePair3SPN2::readFromRestart(pugi::xml_node restData) {
  auto curr_node = restData.first_child();
  while (curr_node) {
    string tag = curr_node.name();
    if (tag == "types") {
      for (auto type_node = curr_node.first_child(); type_node; type_node = type_node.next_sibling()) {
        int type;
        double phi0;
        double sigma;
        double k;
        double epsi;
        double alpha;
        double theta1;
        double theta2;
        std::string type_ = type_node.attribute("id").value();
        type = atoi(type_.c_str());
        std::string phi0_ = type_node.attribute("phi0").value();
        std::string sigma_ = type_node.attribute("sigma").value();
        std::string k_ = type_node.attribute("k").value();
        std::string epsi_ = type_node.attribute("epsi").value();
        std::string alpha_ = type_node.attribute("alpha").value();
        std::string theta1_ = type_node.attribute("theta1").value();
        std::string theta2_ = type_node.attribute("theta2").value();
        phi0 = atof(phi0_.c_str());
        sigma = atof(sigma_.c_str());
        k     = atof(k_.c_str());
        epsi  = atof(epsi_.c_str());
        alpha = atof(alpha_.c_str());
        alpha = atof(alpha_.c_str());
        theta1= atof(theta1_.c_str());
        theta2= atof(theta2_.c_str());
        BasePair3SPN2 dummy(phi0, sigma, k, epsi, alpha, theta1, theta2, type);
        setForcerType(type, dummy);
      }
    } else if (tag == "members") {
      for (auto member_node = curr_node.first_child(); member_node; member_node = member_node.next_sibling()) {
	int type;
	double phi0;
    double sigma;
    double k;
    double epsi;
    double alpha;
    double theta1;
    double theta2;
	int ids[4];
	std::string type_ = member_node.attribute("type").value();
	std::string atom_a = member_node.attribute("atomID_a").value();
	std::string atom_b = member_node.attribute("atomID_b").value();
	std::string atom_c = member_node.attribute("atomID_c").value();
	std::string atom_d = member_node.attribute("atomID_d").value();
	std::string phi0_  = member_node.attribute("phi0").value();
	std::string sigma_ = member_node.attribute("sigma").value();
	std::string k_     = member_node.attribute("k").value();
	std::string epsi_  = member_node.attribute("epsi").value();
	std::string alpha_ = member_node.attribute("alpha").value();
	std::string theta1_ = member_node.attribute("theta1").value();
	std::string theta2_ = member_node.attribute("theta2").value();
	type = atoi(type_.c_str());
	ids[0] = atoi(atom_a.c_str());
	ids[1] = atoi(atom_b.c_str());
	ids[2] = atoi(atom_c.c_str());
	ids[3] = atoi(atom_d.c_str());
	phi0  = atof(phi0_.c_str());
	sigma = atof(sigma_.c_str());
	k     = atof(k_.c_str());
	epsi  = atof(epsi_.c_str());
	alpha = atof(alpha_.c_str());
    theta1= atof(theta1_.c_str());
    theta2= atof(theta2_.c_str());
	Atom * a = &state->idToAtom(ids[0]);
	Atom * b = &state->idToAtom(ids[1]);
	Atom * c = &state->idToAtom(ids[2]);
	Atom * d = &state->idToAtom(ids[3]);
	if (a == NULL) {cout << "The first atom does not exist" <<endl; return false;};
	if (b == NULL) {cout << "The second atom does not exist" <<endl; return false;};
	if (c == NULL) {cout << "The third atom does not exist" <<endl; return false;};
	if (d == NULL) {cout << "The fourth atom does not exist" <<endl; return false;};
	createBasePair(a, b, c, d, phi0, sigma, k, epsi, alpha, theta1, theta2, type);
      }
    }
    curr_node = curr_node.next_sibling();
  }
  return true;
}


void export_FixBasePair3SPN2() {
    py::class_<FixBasePair3SPN2,
                          SHARED(FixBasePair3SPN2),
                          py::bases<Fix, TypedItemHolder> > (
        "FixBasePair3SPN2",
        py::init<SHARED(State), string> (
            py::args("state", "handle")
        )
    )
    .def("createBasePair", &FixBasePair3SPN2::createBasePair,
            (py::arg("phi0")=-1,
             py::arg("sigma")=-1,
             py::arg("k")=-1,
             py::arg("epsi")=-1,
             py::arg("alpha")=-1,
             py::arg("theta1")=-1,
             py::arg("theta2")=-1,
             py::arg("type")=-1)
        )

    .def("setBasePairTypeCoefs", &FixBasePair3SPN2::setBasePairTypeCoefs, 
            (py::arg("type"), 
             py::arg("phi0"),
             py::arg("sigma"),
             py::arg("k"),
             py::arg("epsi"),
             py::arg("alpha"),
             py::arg("theta1"),
             py::arg("theta2"))
        )

    .def_readonly("basepairs", &FixBasePair3SPN2::pyForcers)

    ;

}

