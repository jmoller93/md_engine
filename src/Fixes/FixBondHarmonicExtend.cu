#include "helpers.h"
#include "FixBondHarmonicExtend.h"
#include "cutils_func.h"
#include "FixHelpers.h"
#include "BondEvaluate.h"
#include "ReadConfig.h"
namespace py = boost::python;
using namespace std;

const std::string bondHarmonicExtendType = "BondHarmonicExtend";

FixBondHarmonicExtend::FixBondHarmonicExtend(SHARED(State) state_, string handle)
    : FixBond(state_, handle, string("None"), bondHarmonicExtendType, true, 1) {
        readFromRestart();
    }



void FixBondHarmonicExtend::createBond(Atom *a, Atom *b, double k1, double k2, double r0, int type) {
    vector<Atom *> atoms = {a, b};
    validAtoms(atoms);
    if (type == -1) {
        assert(k1!=-1 and k2!=-1 and r0!=-1);
    }
    bonds.push_back(BondHarmonicExtend(a, b, k1, k2, r0, type));
    pyListInterface.updateAppendedMember();
    
}

void FixBondHarmonicExtend::setBondTypeCoefs(int type, double k1, double k2, double r0) {
    assert(r0>=0);
    BondHarmonicExtend dummy(k1, k2, r0, type);
    setBondType(type, dummy);
}

void FixBondHarmonicExtend::compute(int virialMode) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();
    GPUData &gpd = state->gpd;
    //cout << "Max bonds per block is " << maxBondsPerBlock << endl;
    if (virialMode) {
        compute_force_bond<BondHarmonicExtendType, BondEvaluatorHarmonicExtend, true> <<<NBLOCK(nAtoms), PERBLOCK, sizeof(BondGPU) * maxBondsPerBlock + sharedMemSizeForParams>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), bondsGPU.data(), bondIdxs.data(), parameters.data(), parameters.size(), state->boundsGPU, gpd.virials.d_data.data(), usingSharedMemForParams, evaluator);
    } else {
        compute_force_bond<BondHarmonicExtendType, BondEvaluatorHarmonicExtend, false> <<<NBLOCK(nAtoms), PERBLOCK, sizeof(BondGPU) * maxBondsPerBlock + sharedMemSizeForParams>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), bondsGPU.data(), bondIdxs.data(), parameters.data(), parameters.size(), state->boundsGPU, gpd.virials.d_data.data(), usingSharedMemForParams, evaluator);
    }
}

void FixBondHarmonicExtend::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();
    //cout << "Max bonds per block is " << maxBondsPerBlock << endl;
    compute_energy_bond<<<NBLOCK(nAtoms), PERBLOCK, sizeof(BondGPU) * maxBondsPerBlock + sharedMemSizeForParams>>>(nAtoms, state->gpd.xs(activeIdx), perParticleEng, state->gpd.idToIdxs.d_data.data(), bondsGPU.data(), bondIdxs.data(), parameters.data(), parameters.size(), state->boundsGPU, usingSharedMemForParams, evaluator);
}

string FixBondHarmonicExtend::restartChunk(string format) {
    stringstream ss;
    ss << "<types>\n";
    for (auto it = bondTypes.begin(); it != bondTypes.end(); it++) {
        ss << "<" << "type id='" << it->first << "'";
        ss << bondTypes[it->first].getInfoString() << "'/>\n";
    }
    ss << "</types>\n";
    ss << "<members>\n";
    for (BondVariant &forcerVar : bonds) {
        BondHarmonicExtend &forcer= boost::get<BondHarmonicExtend>(forcerVar);
        ss << forcer.getInfoString();
    }
    ss << "</members>\n";
    return ss.str();
}

bool FixBondHarmonicExtend::readFromRestart() {
    auto restData = getRestartNode();
    if (restData) {
        auto curr_node = restData.first_child();
        while (curr_node) {
            std::string tag = curr_node.name();
            if (tag == "types") {
                for (auto type_node = curr_node.first_child(); type_node; type_node = type_node.next_sibling()) {
                    int type;
                    double k1;
                    double k2;
                    double r0;
                    std::string type_ = type_node.attribute("id").value();
                    type = atoi(type_.c_str());
                    std::string k1_ = type_node.attribute("k1").value();
                    std::string k2_ = type_node.attribute("k2").value();
                    std::string r0_ = type_node.attribute("r0").value();
                    k1 = atof(k1_.c_str());
                    k2 = atof(k2_.c_str());
                    r0 = atof(r0_.c_str());

                    setBondTypeCoefs(type, k1, k2, r0);
                }
            } else if (tag == "members") {
                for (auto member_node = curr_node.first_child(); member_node; member_node = member_node.next_sibling()) {
                    int type;
                    double k1;
                    double k2;
                    double r0;
                    int ids[2];
                    std::string type_ = member_node.attribute("type").value();
                    std::string atom_a = member_node.attribute("atom_a").value();
                    std::string atom_b = member_node.attribute("atom_b").value();
                    std::string k1_ = member_node.attribute("k1").value();
                    std::string k2_ = member_node.attribute("k2").value();
                    std::string r0_ = member_node.attribute("r0").value();
                    type = atoi(type_.c_str());
                    ids[0] = atoi(atom_a.c_str());
                    ids[1] = atoi(atom_b.c_str());
                    Atom * a = &state->idToAtom(ids[0]);
                    Atom * b = &state->idToAtom(ids[1]);
                    k1 = atof(k1_.c_str());
                    k2 = atof(k2_.c_str());
                    r0 = atof(r0_.c_str());

                    createBond(a, b, k1, k2, r0, type);
                }
            }
            curr_node = curr_node.next_sibling();
        }
    }
    return true;
}

void export_FixBondHarmonicExtend() {
  

  
    py::class_<FixBondHarmonicExtend, SHARED(FixBondHarmonicExtend), py::bases<Fix, TypedItemHolder> >
    (
        "FixBondHarmonicExtend", py::init<SHARED(State), string> (py::args("state", "handle"))
    )
    .def("createBond", &FixBondHarmonicExtend::createBond,
            (py::arg("k1")=-1,
             py::arg("k2")=-1,
             py::arg("r0")=-1,
             py::arg("type")=-1)
        )
    .def("setBondTypeCoefs", &FixBondHarmonicExtend::setBondTypeCoefs,
            (py::arg("type"),
             py::arg("k1"),
             py::arg("k2"),
             py::arg("r0"))
        )
    .def_readonly("bonds", &FixBondHarmonicExtend::pyBonds)    
    ;

}
