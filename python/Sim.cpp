#include "Python.h"
#include <functional>
#include <map>
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/args.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/operators.hpp>

using namespace boost::python;
#include <tuple>
using namespace std;
#include "State.h"
#include "Atom.h"
#include "Vector.h" 
#include "InitializeAtoms.h"
#include "Bounds.h"
#include "includeFixes.h"
#include "IntegratorVerlet.h"
#include "IntegratorRelax.h"
#include "FixLangevin.h"
#include "boost_stls.h"
#include "PythonOperation.h"
//#include "DataManager.h"
#include "DataSetUser.h"
#include "ReadConfig.h"
#include "TypedItemHolder.h"
#include "Angle.h"
#include "Dihedral.h"
#include "Improper.h"
#include "Molecule.h"
#include "BasePair.h"
//#include "DataTools.h"
BOOST_PYTHON_MODULE(Sim) {
    export_stls();	

    export_Vector();
    export_VectorInt();	
    export_Atom();
    export_Molecule();
    export_Bounds();
    export_Integrator();
    export_IntegratorVerlet();
    export_IntegratorRelax();
    export_TypedItemHolder();
    export_Fix();
    export_FixBondHarmonic();
    export_FixBondHarmonicExtend();
    export_FixBondFENE();
    export_FixBondGoLike();
    export_BondHarmonic();
    export_BondHarmonicExtend();
    export_BondFENE();
    export_BondGoLike();
    
    export_FixAngleBaseStacking();
    export_FixAngleHarmonic();
    export_FixAngleCosineDelta();
    export_AngleHarmonic();
    export_AngleCosineDelta();
    export_AngleBaseStacking();
    export_FixBasePair3SPN2();
    export_FixCrossStack3SPN2();
    export_BasePairs();
    export_CrossStacks();

    export_FixImproperHarmonic();
    export_Impropers();
    export_FixDihedralOPLS();
    export_FixDihedralCHARMM();
    export_FixDihedralPeriodic();
    export_FixDihedralGauss();
    export_Dihedrals();
    export_FixWall();
    export_FixWallHarmonic();
    export_FixWallLJ126();
    export_FixSpringStatic();
    export_Fix2d();
    export_FixLinearMomentum();

    export_FixPair();
    export_FixLJCut(); //make there be a pair base class in boost!
    export_FixLJCutFS();
    export_FixLJRepul();
    export_FixTICG();
    export_FixWCA();
    
    export_FixCharge();
    export_FixChargePairDSF();
    export_FixChargePairDH();
    export_FixChargeEwald();

    export_FixNoseHoover();
    export_FixNVTRescale();
    export_FixLangevin();

    export_FixPressureBerendsen();

    export_AtomParams();
    export_DataManager();
    export_ReadConfig();
    export_PythonOperation();

    export_WriteConfig();
    export_InitializeAtoms();
        
    export_Units();

    export_State(); 	
    export_DeviceManager();
    export_DataSetUser();

}
