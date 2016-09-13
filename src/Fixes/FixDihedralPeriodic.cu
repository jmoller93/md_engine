#include "helpers.h"
#include "FixDihedralPeriodic.h"
#include "FixHelpers.h"
#include "cutils_func.h"
#include "DihedralEvaluate.h"
namespace py = boost::python;
using namespace std;

const std::string dihedralPeriodicType = "DihedralPeriodic";


FixDihedralPeriodic::FixDihedralPeriodic(SHARED(State) state_, string handle) : FixPotentialMultiAtom (state_, handle, dihedralPeriodicType, true){
    readFromRestart();
}


void FixDihedralPeriodic::compute(bool computeVirials) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();


    GPUData &gpd = state->gpd;
    if (computeVirials) {
        compute_force_dihedral<DihedralPeriodicType, DihedralEvaluatorPeriodic, true><<<NBLOCK(forcersGPU.size()), PERBLOCK, sharedMemSizeForParams>>>(forcersGPU.size(), gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), usingSharedMemForParams, evaluator);
    } else {
        compute_force_dihedral<DihedralPeriodicType, DihedralEvaluatorPeriodic, false><<<NBLOCK(forcersGPU.size()), PERBLOCK, sharedMemSizeForParams>>>(forcersGPU.size(), gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), usingSharedMemForParams, evaluator);
    }

}

void FixDihedralPeriodic::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();

    GPUData &gpd = state->gpd;
    //compute_energy_dihedral<<<NBLOCK(nAtoms), PERBLOCK, sizeof(DihedralGPU) * maxForcersPerBlock + sizeof(DihedralPeriodicType) * parameters.size() >>>(nAtoms, gpd.xs(activeIdx), perParticleEng, gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), evaluator);
    compute_energy_dihedral<<<NBLOCK(nAtoms), PERBLOCK, sharedMemSizeForParams>>>(forcersGPU.size(), gpd.xs(activeIdx), perParticleEng, gpd.idToIdxs.d_data.data(), forcersGPU.data(), state->boundsGPU, parameters.data(), parameters.size(), usingSharedMemForParams, evaluator);
    

}



void FixDihedralPeriodic::createDihedral(Atom *a, Atom *b, Atom *c, Atom *d, double v1, double v2, double phiRef, int type) {
    double vs[2] = {v1, v2};
    if (type==-1) {
        for (int i=0; i<2; i++) {
            assert(vs[i] != COEF_DEFAULT);
        }
        assert(phiRef != COEF_DEFAULT);
    }
    forcers.push_back(DihedralPeriodic(a, b, c, d, vs, phiRef, type));
    pyListInterface.updateAppendedMember();
}

void FixDihedralPeriodic::createDihedralPy(Atom *a, Atom *b, Atom *c, Atom *d, py::list coefs, double phiRef, int type) {
    double coefs_c[2];
    if (type!=-1) {
        createDihedral(a, b, c, d, COEF_DEFAULT, COEF_DEFAULT, COEF_DEFAULT, type);
    } else {
        assert(len(coefs) == 2);
        for (int i=0; i<2; i++) {
            py::extract<double> coef(coefs[i]);
            assert(coef.check());
            coefs_c[i] = coef;
        }
        createDihedral(a, b, c, d, coefs_c[0], coefs_c[1], phiRef, type);

    }
}
void FixDihedralPeriodic::setDihedralTypeCoefs(int type, py::list coefs, double phiRef) {
    assert(len(coefs)==2);
    double coefs_c[2];
    for (int i=0; i<2; i++) {
        py::extract<double> coef(coefs[i]);
        assert(coef.check());
        coefs_c[i] = coef;
    }

    DihedralPeriodic dummy(coefs_c, phiRef, type);
    setForcerType(type, dummy);
}

bool FixDihedralPeriodic::readFromRestart() {
    auto restData = getRestartNode();
    if (restData) {
        auto curr_node = restData.first_child();
        while (curr_node) {
            string tag = curr_node.name();
            if (tag == "types") {
                for (auto type_node = curr_node.first_child(); type_node; type_node = type_node.next_sibling()) {
                    int type;
                    double coefs[4];
                    double phiRef;
                    std::string type_ = type_node.attribute("id").value();
                    type = atoi(type_.c_str());
                    std::string coef_a = type_node.attribute("coef_a").value();
                    std::string coef_b = type_node.attribute("coef_b").value();
                    std::string phi_ref = type_node.attribute("phi_ref").value();
                    coefs[0] = atof(coef_a.c_str());
                    coefs[1] = atof(coef_b.c_str());
                    phiRef = atof(phi_ref.c_str());
                    DihedralPeriodic dummy(coefs, phiRef, type);
                    setForcerType(type, dummy);
                }
            } else if (tag == "members") {
                for (auto member_node = curr_node.first_child(); member_node; member_node = member_node.next_sibling()) {
                    int type;
                    double coefs[4];
                    double phiRef;
                    int ids[4];
                    std::string type_ = member_node.attribute("type").value();
                    std::string atom_a = member_node.attribute("atomID_a").value();
                    std::string atom_b = member_node.attribute("atomID_b").value();
                    std::string atom_c = member_node.attribute("atomID_c").value();
                    std::string atom_d = member_node.attribute("atomID_d").value();
                    std::string coef_a = member_node.attribute("coef_a").value();
                    std::string coef_b = member_node.attribute("coef_b").value();
                    std::string phi_ref = member_node.attribute("phi_ref").value();
                    type = atoi(type_.c_str());
                    ids[0] = atoi(atom_a.c_str());
                    ids[1] = atoi(atom_b.c_str());
                    ids[2] = atoi(atom_c.c_str());
                    ids[3] = atoi(atom_d.c_str());
                    coefs[0] = atof(coef_a.c_str());
                    coefs[1] = atof(coef_b.c_str());
                    phiRef = atof(phi_ref.c_str());
                    Atom * a = &state->idToAtom(ids[0]);
                    Atom * b = &state->idToAtom(ids[1]);
                    Atom * c = &state->idToAtom(ids[2]);
                    Atom * d = &state->idToAtom(ids[3]);
                    if (a == NULL) {cout << "The first atom does not exist" <<endl; return false;};
                    if (b == NULL) {cout << "The second atom does not exist" <<endl; return false;};
                    if (c == NULL) {cout << "The third atom does not exist" <<endl; return false;};
                    if (d == NULL) {cout << "The fourth atom does not exist" <<endl; return false;};
                    createDihedral(a, b, c, d, coefs[0], coefs[1],  phiRef, type);
                }
            }
            curr_node = curr_node.next_sibling();
        }
    }
    return true;
}


void export_FixDihedralPeriodic() {
    py::class_<FixDihedralPeriodic,
                          SHARED(FixDihedralPeriodic),
                          py::bases<Fix, TypedItemHolder> > (
        "FixDihedralPeriodic",
        py::init<SHARED(State), string> (
            py::args("state", "handle")
        )
    )
    .def("createDihedral", &FixDihedralPeriodic::createDihedralPy,
            (py::arg("coefs")=py::list(),
             py::arg("phiRef")=-1,
             py::arg("type")=-1)
        )

    .def("setDihedralTypeCoefs", &FixDihedralPeriodic::setDihedralTypeCoefs, 
            (py::arg("type"), 
             py::arg("coefs"),
             py::arg("phiRef")
            )
        )
    .def_readonly("dihedrals", &FixDihedralPeriodic::pyForcers)

    ;

}

