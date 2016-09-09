import sys
import matplotlib.pyplot as plt
sys.path = sys.path + ['../build/python/build/lib.linux-x86_64-2.7']
#from Sim import *
from Sim import *
state = State()
state.deviceManager.setDevice(0)
state.bounds = Bounds(state, lo = Vector(-571.389, -571.389, -571.389), hi = Vector(571.389, 571.389, 571.389))
state.rCut = 2.0
state.padding = 0.6
state.periodicInterval = 7
state.shoutEvery = 1000
#I believe this is correct for ps scale
state.dt = 1.26492

#Working directory name (really it's the input directory name, but I'm consistent with Gordo's scripts)
wdir = 'input_3spn2'

#This is a nice tool that allows for debugging/ printing without too much slowdown

def operate(turn):
    fnme = open('force_test.dat', 'a+')
    fnme.write("%10.5f\t\t" % state.atoms[1000].force[0])
    fnme.write("%10.5f\t\t" % state.atoms[1000].force[1])
    fnme.write("%10.5f\n" % state.atoms[1000].force[2])
    fnme.close()

#Activate the python function as a "Fix"
pyOp = PythonOperation(operateEvery=1, handle='hello', operation=operate)
state.activatePythonOperation(pyOp)

#Read species info
f = open('%s/in00_spcs.xml' % wdir).readlines()
for line in f:
    spcsInfo = line.split()
    state.atomParams.addSpecies(handle=str(spcsInfo[1]), mass=float(spcsInfo[3]), atomicNum=float(spcsInfo[0]))

#Read the input configuration
f = open('%s/in00_conf.xml' % wdir).readlines()
for line in f:
    bits = line.split()
    state.addAtom(str(bits[0]), Vector(float(bits[1]), float(bits[2]), float(bits[3])))

The repulsive force may be a challenge
nonbond = FixLJRepul(state, 'excluded')
nonbond.setParameter('sig', 'spc1', 'spc1', 1)
nonbond.setParameter('eps', 'spc1', 'spc1', 1)
#state.activateFix(nonbond)

#Read the bond information
#The 2.5 is the energy scale
f = open('%s/in00_bond.xml' % wdir).readlines()
for line in f:
    bondInfo = line.split()
    bond = FixBondHarmonic(state, 'bond')
    bond.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5])/2.5, float(bondInfo[4]))

#state.activateFix(bond)
#Read and activate the bonded angles
f = open('%s/in00_bend.xml' % wdir).readlines()
for line in f:
    angleInfo = line.split()
    angle = FixAngleHarmonic(state, 'angle')
    angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], float(angleInfo[6])/2.5, float(angleInfo[5]))

#state.activateFix(angle)
#Read and activate the bonded dihedrals (Periodic and Gauss)
f = open('%s/in00_tors.xml' % wdir).readlines()
for line in f:
    diheInfo = line.split()
    if int(diheInfo[5]) == 2:
        diheGauss = FixDihedralGauss(state, 'gauss')
        diheGauss.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],float(diheInfo[6]) * math.pi / 180.0 + math.pi, float(diheInfo[8]), float(diheInfo[7])/2.5)
        dihePeri = FixDihedralPeriodic(state, 'periodic')
        coef = [2.0,0.0,0.0,0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef)
    else:
        dihePeri = FixDihedralPeriodic(state, 'OPLS')
        coef = [float(diheInfo[6])*2.0, 0.0,float(diheInfo[8])*2.0, 0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef)

if diheGauss.exists():
    state.activateFix(diheGauss)
else:
    print("Gaussian dihedrals not found. This run is borked")
state.activateFix(dihePeri)

natv = FixGoLike()
f = open('%s/in00_natv.xml' % wdir).readlines()
for line in f:
    natvInfo = line.split()
    natv = FixBondGoLike(state, 'natv')
    natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])])

InitializeAtoms.initTemp(state, 'all', 1.0)

fixNVT = FixLangevin(state, 'temp', 'all', 1.0)

state.activateFix(fixNVT)

integVerlet = IntegratorVerlet(state)

tempData = state.dataManager.recordTemperature('all', 50)
#pressureData = state.dataManager.recordPressure('all', 100)
#boundsData = state.dataManager.recordBounds(100)
engData = state.dataManager.recordEnergy('all', 50)

writeconfig = WriteConfig(state, fn='test_out', writeEvery=1000, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integVerlet.run(100)
