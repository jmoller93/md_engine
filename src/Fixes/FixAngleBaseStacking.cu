
#include "FixHelpers.h"
#include "helpers.h"
#include "FixAngleBaseStacking.h"
#include "cutils_func.h"
#include "AngleEvaluate.h"
using namespace std;
const string angleBaseStackingType = "AngleBaseStacking";
FixAngleBaseStacking::FixAngleBaseStacking(boost::shared_ptr<State> state_, string handle)
  : FixPotentialMultiAtom(state_, handle, angleBaseStackingType, true)
{
  if (state->readConfig->fileOpen) {
    auto restData = state->readConfig->readFix(type, handle);
    if (restData) {
      std::cout << "Reading restart data for fix " << handle << std::endl;
      readFromRestart(restData);
    }
  }
}

namespace py = boost::python;

void FixAngleBaseStacking::compute(bool computeVirials) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();
    GPUData &gpd = state->gpd;
    if (computeVirials) {
        compute_force_angle<AngleBaseStackingType, AngleEvaluatorBaseStacking, true> <<<NBLOCK(nAtoms), PERBLOCK, sizeof(AngleGPU) * maxForcersPerBlock + parameters.size() * sizeof(AngleBaseStackingType)>>>(nAtoms, state->gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), evaluator);
    } else {
        compute_force_angle<AngleBaseStackingType, AngleEvaluatorBaseStacking, false> <<<NBLOCK(nAtoms), PERBLOCK, sizeof(AngleGPU) * maxForcersPerBlock + parameters.size() * sizeof(AngleBaseStackingType)>>>(nAtoms, state->gpd.xs(activeIdx), gpd.fs(activeIdx), gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), gpd.virials.d_data.data(), evaluator);

    }

}

void FixAngleBaseStacking::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int activeIdx = state->gpd.activeIdx();
    compute_energy_angle<<<NBLOCK(nAtoms), PERBLOCK, sizeof(AngleGPU) * maxForcersPerBlock + parameters.size() * sizeof(AngleBaseStackingType)>>>(nAtoms, state->gpd.xs(activeIdx), perParticleEng, state->gpd.idToIdxs.d_data.data(), forcersGPU.data(), forcerIdxs.data(), state->boundsGPU, parameters.data(), parameters.size(), evaluator);
}
//void cumulativeSum(int *data, int n);
// okay, so the net result of this function is that two arrays (items, idxs of
// items) are on the gpu and we know how many bonds are in bondiest block

void FixAngleBaseStacking::createAngle(Atom *a, Atom *b, Atom *c, double k, double theta0, double epsi, double sigma, double alpha, int type) {
    vector<Atom *> atoms = {a, b, c};
    validAtoms(atoms);
    if (type == -1) {
        assert(k!=COEF_DEFAULT and theta0!=COEF_DEFAULT and epsi!=COEF_DEFAULT and sigma!=COEF_DEFAULT and alpha!=COEF_DEFAULT);
    }
    forcers.push_back(AngleBaseStacking(a, b, c, k, theta0, epsi, sigma, alpha, type));
    pyListInterface.updateAppendedMember();
}

void FixAngleBaseStacking::setAngleTypeCoefs(int type, double k, double theta0, double epsi, double sigma, double alpha) {
    //cout << type << " " << k << " " << theta0 << endl;
    mdAssert(theta0>=0 and theta0 <= M_PI, "Angle theta must be between zero and pi");
    AngleBaseStacking dummy(k, theta0, epsi, sigma, alpha);
    setForcerType(type, dummy);
}


bool FixAngleBaseStacking::readFromRestart(pugi::xml_node restData) {
    auto curr_node = restData.first_child();
    while (curr_node) {
        std::string tag = curr_node.name();
        if (tag == "types") {
            for (auto type_node = curr_node.first_child(); type_node; type_node = type_node.next_sibling()) {
                int type;
                double k;
                double theta0;
                double epsi;
                double sigma;
                double alpha;
                std::string type_ = type_node.attribute("id").value();
                type = atoi(type_.c_str());
                std::string k_ = type_node.attribute("k").value();
                std::string theta0_ = type_node.attribute("theta0").value();
                std::string epsi_ = type_node.attribute("epsi").value();
                std::string sigma_ = type_node.attribute("sigma").value();
                std::string alpha_ = type_node.attribute("alpha").value();
                k = atof(k_.c_str());
                theta0 = atof(theta0_.c_str());
                epsi = atof(epsi_.c_str());
                sigma = atof(sigma_.c_str());
                alpha = atof(alpha_.c_str());
                setAngleTypeCoefs(type, k, theta0, epsi, sigma, alpha);
            }
        } else if (tag == "members") {
            for (auto member_node = curr_node.first_child(); member_node; member_node = member_node.next_sibling()) {
                int type;
                double k;
                double theta0;
                double epsi;
                double sigma;
                double alpha;
                int ids[3];
                std::string type_ = member_node.attribute("type").value();
                std::string atom_a = member_node.attribute("atom_a").value();
                std::string atom_b = member_node.attribute("atom_b").value();
                std::string atom_c = member_node.attribute("atom_c").value();
                std::string k_ = member_node.attribute("k").value();
                std::string theta0_ = member_node.attribute("theta0").value();
                std::string epsi_ = member_node.attribute("epsi").value();
                std::string sigma_ = member_node.attribute("sigma").value();
                std::string alpha_ = member_node.attribute("alpha").value();
                type = atoi(type_.c_str());
                ids[0] = atoi(atom_a.c_str());
                ids[1] = atoi(atom_b.c_str());
                ids[2] = atoi(atom_c.c_str());
                Atom * a = &state->idToAtom(ids[0]);
                Atom * b = &state->idToAtom(ids[1]);
                Atom * c = &state->idToAtom(ids[2]);
                k = atof(k_.c_str());
                theta0 = atof(theta0_.c_str());
                epsi = atof(epsi_.c_str());
                sigma = atof(sigma_.c_str());
                alpha = atof(alpha_.c_str());

                createAngle(a, b, c, k, theta0, epsi, sigma, alpha, type);
            }
        }
        curr_node = curr_node.next_sibling();
    }
    return true;
}

void export_FixAngleBaseStacking() {
    boost::python::class_<FixAngleBaseStacking,
                          boost::shared_ptr<FixAngleBaseStacking>,
                          boost::python::bases<Fix, TypedItemHolder> >(
        "FixAngleBaseStacking",
        boost::python::init<boost::shared_ptr<State>, string>(
                                boost::python::args("state", "handle"))
    )
    .def("createAngle", &FixAngleBaseStacking::createAngle,
            (boost::python::arg("k")=COEF_DEFAULT,
             boost::python::arg("theta0")=COEF_DEFAULT,
             boost::python::arg("epsi")=COEF_DEFAULT,
             boost::python::arg("sigma")=COEF_DEFAULT,
             boost::python::arg("alpha")=COEF_DEFAULT,
             boost::python::arg("type")=-1)
        )
    .def("setAngleTypeCoefs", &FixAngleBaseStacking::setAngleTypeCoefs,
            (boost::python::arg("type")=-1,
             boost::python::arg("k")=COEF_DEFAULT,
             boost::python::arg("theta0")=COEF_DEFAULT,
             boost::python::arg("epsi")=COEF_DEFAULT,
             boost::python::arg("sigma")=COEF_DEFAULT,
             boost::python::arg("alpha")=COEF_DEFAULT)
        )
    .def_readonly("angles", &FixAngleBaseStacking::pyForcers)
    ;
}

