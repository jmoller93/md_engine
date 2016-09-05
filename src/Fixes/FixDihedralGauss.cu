#include "helpers.h"
#include "FixDihedralGauss.h"
#include "FixHelpers.h"
#include "cutils_func.h"
#include "DihedralEvaluate.h"
namespace py = boost::python;
using namespace std;

const std::string dihedralGaussType = "DihedralGauss";


FixDihedralGauss::FixDihedralGauss(SHARED(State) state_, string handle) : FixPotentialMultiAtom (state_, handle, dihedralGaussType, true){
  if (state->readConfig->fileOpen) {
    auto restData = state->readConfig->readFix(type, handle);
    if (restData) {
      std::cout << "Reading restart data for fix " << handle << std::endl;
      readFromRestart(restData);
    }
  }
}


void FixDihedralGauss::compute(bool computeVirials) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();


    GPUData &gpd = state->gpd;
    if (computeVirials) {
        compute_force_dihedral<DihedralGaussType, DihedralEvaluatorGauss, true><<<NBLOCK(nAtoms), PERBLOCK, sizeof(DihedralGPU) * maxForcersPerBlock + sizeof(DihedralGaussType) * parameters.size() >>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), evaluator);
    } else {
        compute_force_dihedral<DihedralGaussType, DihedralEvaluatorGauss, false><<<NBLOCK(nAtoms), PERBLOCK, sizeof(DihedralGPU) * maxForcersPerBlock + sizeof(DihedralGaussType) * parameters.size() >>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), evaluator);
    }

}

void FixDihedralGauss::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();

    GPUData &gpd = state->gpd;
    compute_energy_dihedral<<<NBLOCK(nAtoms), PERBLOCK, sizeof(DihedralGPU) * maxForcersPerBlock + sizeof(DihedralGaussType) * parameters.size() >>>(nAtoms, gpd.xs(activeIdx), perParticleEng, gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), evaluator);
    

}



void FixDihedralGauss::createDihedral(Atom *a, Atom *b, Atom *c, Atom *d, double phi0, double sigma, double k0, int type) {
    vector<Atom *> atoms = {a, b, c, d};
    if (type==-1) {
            assert(phi0!=COEF_DEFAULT and sigma!=COEF_DEFAULT and k0!=COEF_DEFAULT);
    }
    forcers.push_back(DihedralGauss(a, b, c, d, phi0, sigma, k0, type));
    pyListInterface.updateAppendedMember();
}

void FixDihedralGauss::setDihedralTypeCoefs(int type, double phi0, double sigma, double k0) {
    assert(sigma>0);
    DihedralGauss dummy(phi0, sigma, k0, type);
    setForcerType(type, dummy);
}

bool FixDihedralGauss::readFromRestart(pugi::xml_node restData) {
  auto curr_node = restData.first_child();
  while (curr_node) {
    string tag = curr_node.name();
    if (tag == "types") {
      for (auto type_node = curr_node.first_child(); type_node; type_node = type_node.next_sibling()) {
        int type;
        double phi0;
        double sigma;
        double k0;
        std::string type_ = type_node.attribute("id").value();
        type = atoi(type_.c_str());
        std::string phi0_ = type_node.attribute("phi0").value();
        std::string sigma_ = type_node.attribute("sigma").value();
        std::string k0_ = type_node.attribute("k0").value();
        phi0 = atof(phi0_.c_str());
        sigma = atof(sigma_.c_str());
        k0    = atof(k0_.c_str());

        setDihedralTypeCoefs(type, phi0, sigma, k0);
      }
    } else if (tag == "members") {
      for (auto member_node = curr_node.first_child(); member_node; member_node = member_node.next_sibling()) {
	int type;
	double phi0;
    double sigma;
    double k0;
	int ids[4];
	std::string type_ = member_node.attribute("type").value();
	std::string atom_a = member_node.attribute("atomID_a").value();
	std::string atom_b = member_node.attribute("atomID_b").value();
	std::string atom_c = member_node.attribute("atomID_c").value();
	std::string atom_d = member_node.attribute("atomID_d").value();
	std::string phi0_  = member_node.attribute("phi0").value();
	std::string sigma_ = member_node.attribute("sigma").value();
	std::string k0_    = member_node.attribute("k0").value();
	type = atoi(type_.c_str());
	ids[0] = atoi(atom_a.c_str());
	ids[1] = atoi(atom_b.c_str());
	ids[2] = atoi(atom_c.c_str());
	ids[3] = atoi(atom_d.c_str());
	phi0  = atof(phi0_.c_str());
	sigma = atof(sigma_.c_str());
	k0    = atof(k0_.c_str());
	Atom * a = &state->idToAtom(ids[0]);
	Atom * b = &state->idToAtom(ids[1]);
	Atom * c = &state->idToAtom(ids[2]);
	Atom * d = &state->idToAtom(ids[3]);
	if (a == NULL) {cout << "The first atom does not exist" <<endl; return false;};
	if (b == NULL) {cout << "The second atom does not exist" <<endl; return false;};
	if (c == NULL) {cout << "The third atom does not exist" <<endl; return false;};
	if (d == NULL) {cout << "The fourth atom does not exist" <<endl; return false;};
	createDihedral(a, b, c, d, phi0, sigma, k0, type);
      }
    }
    curr_node = curr_node.next_sibling();
  }
  return true;
}


void export_FixDihedralGauss() {
    py::class_<FixDihedralGauss,
                          SHARED(FixDihedralGauss),
                          py::bases<Fix, TypedItemHolder> > (
        "FixDihedralGauss",
        py::init<SHARED(State), string> (
            py::args("state", "handle")
        )
    )
    .def("createDihedral", &FixDihedralGauss::createDihedral,
             (py::arg("phi0")=COEF_DEFAULT,
             py::arg("sigma")=COEF_DEFAULT,
             py::arg("k0")=COEF_DEFAULT,
             py::arg("type")=-1)
        )

    .def("setDihedralTypeCoefs", &FixDihedralGauss::setDihedralTypeCoefs, 
            (py::arg("type")=-1, 
             py::arg("phi0")=COEF_DEFAULT,
             py::arg("sigma")=COEF_DEFAULT,
             py::arg("k0")=COEF_DEFAULT
            )
        )
    .def_readonly("dihedrals", &FixDihedralGauss::pyForcers)
    ;

}

