#include "helpers.h"
#include "FixBondGoLike.h"
#include "cutils_func.h"
#include "FixHelpers.h"
#include "BondEvaluate.h"
#include "ReadConfig.h"
namespace py = boost::python;
using namespace std;

const std::string bondGoLikeType = "BondGoLike";

FixBondGoLike::FixBondGoLike(SHARED(State) state_, string handle)
    : FixBond(state_, handle, string("None"), bondGoLikeType, true, 1) {
        readFromRestart();
    }



void FixBondGoLike::createBond(Atom *a, Atom *b, double eps, double sig, int type) {
    vector<Atom *> atoms = {a, b};
    validAtoms(atoms);
    if (type == -1) {
        assert(eps!=-1 and sig!=-1);
    }
    bonds.push_back(BondGoLike(a, b, eps, sig, type));
    pyListInterface.updateAppendedMember();
    
}

void FixBondGoLike::setBondTypeCoefs(int type, double eps, double sig) {
    assert(sig>=0);
    BondGoLike dummy(eps, sig, type);
    setBondType(type, dummy);
}

void FixBondGoLike::compute(bool computeVirials) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();
    GPUData &gpd = state->gpd;
    //cout << "Max bonds per block is " << maxBondsPerBlock << endl;
    if (computeVirials) {
        compute_force_bond<BondGoLikeType, BondEvaluatorGoLike, true> <<<NBLOCK(nAtoms), PERBLOCK, sizeof(BondGPU) * maxBondsPerBlock + sharedMemSizeForParams>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), bondsGPU.data(), bondIdxs.data(), parameters.data(), parameters.size(), state->boundsGPU, gpd.virials.d_data.data(), usingSharedMemForParams, evaluator);
    } else {
        compute_force_bond<BondGoLikeType, BondEvaluatorGoLike, false> <<<NBLOCK(nAtoms), PERBLOCK, sizeof(BondGPU) * maxBondsPerBlock + sharedMemSizeForParams>>>(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), bondsGPU.data(), bondIdxs.data(), parameters.data(), parameters.size(), state->boundsGPU, gpd.virials.d_data.data(), usingSharedMemForParams, evaluator);
    }
}

void FixBondGoLike::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();
    //cout << "Max bonds per block is " << maxBondsPerBlock << endl;
    compute_energy_bond<<<NBLOCK(nAtoms), PERBLOCK, sizeof(BondGPU) * maxBondsPerBlock + sharedMemSizeForParams>>>(nAtoms, state->gpd.xs(activeIdx), perParticleEng, state->gpd.idToIdxs.d_data.data(), bondsGPU.data(), bondIdxs.data(), parameters.data(), parameters.size(), state->boundsGPU, usingSharedMemForParams, evaluator);
}

string FixBondGoLike::restartChunk(string format) {
    stringstream ss;
    ss << "<types>\n";
    for (auto it = bondTypes.begin(); it != bondTypes.end(); it++) {
        ss << "<" << "type id='" << it->first << "'";
        ss << bondTypes[it->first].getInfoString() << "'/>\n";
    }
    ss << "</types>\n";
    ss << "<members>\n";
    for (BondVariant &forcerVar : bonds) {
        BondGoLike &forcer= boost::get<BondGoLike>(forcerVar);
        ss << forcer.getInfoString();
    }
    ss << "</members>\n";
    return ss.str();
}

bool FixBondGoLike::readFromRestart() {
    auto restData = getRestartNode();
    if (restData) {
        auto curr_node = restData.first_child();
        while (curr_node) {
            std::string tag = curr_node.name();
            if (tag == "types") {
                for (auto type_node = curr_node.first_child(); type_node; type_node = type_node.next_sibling()) {
                    int type;
                    double eps;
                    double sig;
                    std::string type_ = type_node.attribute("id").value();
                    type = atoi(type_.c_str());
                    std::string eps_ = type_node.attribute("eps").value();
                    std::string sig_ = type_node.attribute("sig").value();
                    eps = atof(eps_.c_str());
                    sig = atof(sig_.c_str());

                    setBondTypeCoefs(type, eps, sig);
                }
            } else if (tag == "members") {
                for (auto member_node = curr_node.first_child(); member_node; member_node = member_node.next_sibling()) {
                    int type;
                    double eps;
                    double sig;
                    int ids[2];
                    std::string type_ = member_node.attribute("type").value();
                    std::string atom_a = member_node.attribute("atom_a").value();
                    std::string atom_b = member_node.attribute("atom_b").value();
                    std::string eps_ = member_node.attribute("eps").value();
                    std::string sig_ = member_node.attribute("sig").value();
                    type = atoi(type_.c_str());
                    ids[0] = atoi(atom_a.c_str());
                    ids[1] = atoi(atom_b.c_str());
                    Atom * a = &state->idToAtom(ids[0]);
                    Atom * b = &state->idToAtom(ids[1]);
                    eps = atof(eps_.c_str());
                    sig = atof(sig_.c_str());

                    createBond(a, b, eps, sig, type);
                }
            }
            curr_node = curr_node.next_sibling();
        }
    }
    return true;
}

void export_FixBondGoLike() {
  

  
    py::class_<FixBondGoLike, SHARED(FixBondGoLike), py::bases<Fix, TypedItemHolder> >
    (
        "FixBondGoLike", py::init<SHARED(State), string> (py::args("state", "handle"))
    )
    .def("createBond", &FixBondGoLike::createBond,
            (py::arg("eps")=-1,
             py::arg("sig")=-1,
             py::arg("type")=-1)
        )
    .def("setBondTypeCoefs", &FixBondGoLike::setBondTypeCoefs,
            (py::arg("type"),
             py::arg("eps"),
             py::arg("sig"))
        )
    .def_readonly("bonds", &FixBondGoLike::pyBonds)    
    ;

}