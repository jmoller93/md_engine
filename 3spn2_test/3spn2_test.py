import sys
import matplotlib.pyplot as plt
import math
sys.path = sys.path + ['../build/python/build/lib.linux-x86_64-2.7']
from Sim import *
state = State()
state.deviceManager.setDevice(0)
state.bounds = Bounds(state, lo = Vector(-571.389, -571.389, -571.389), hi = Vector(571.389, 571.389, 571.389))
state.rCut = 30.0
state.padding = 0.5
state.periodicInterval = 100
state.shoutEvery = 1000
#I believe this is correct for ps scale
state.dt = 0.200

#Working directory name (really it's the input directory name, but I'm consistent with Gordo's scripts)
wdir = 'input_conf_test'

#This is a nice tool that allows for debugging/ printing without too much slowdown

def operate(turn):
    fnme = open('force_test.dat', 'a+')
    fnme.write("%10.5f\t\t" % state.atoms[1000].force[0])
    fnme.write("%10.5f\t\t" % state.atoms[1000].force[1])
    fnme.write("%10.5f\n" % state.atoms[1000].force[2])
    fnme.close()

#Activate the python function as a "Fix"
pyOp = PythonOperation(operateEvery=1000, handle='hello', operation=operate)
#state.activatePythonOperation(pyOp)

#Read species info
f = open('%s/in00_spcs.xml' % wdir).readlines()
for line in f:
    spcsInfo = line.split()
    state.atomParams.addSpecies(handle=str(spcsInfo[1]), mass=float(spcsInfo[3]), atomicNum=int(spcsInfo[0]))

#Read the input configuration
f = open('%s/in00_conf.xml' % wdir).readlines()
for line in f:
    bits = line.split()
    state.addAtom(str(bits[0]), Vector(float(bits[1]), float(bits[2]), float(bits[3])))

#Read the bond information
#The 2.5 is the energy scale
f = open('%s/in00_bond.xml' % wdir).readlines()
bond = FixBondHarmonic(state, 'bond')
for line in f:
    bondInfo = line.split()
    bond.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5]), float(bondInfo[4]))
#state.activateFix(bond)

#Read and activate the bonded angles
f = open('%s/in00_bend.xml' % wdir).readlines()
angle = FixAngleHarmonic(state, 'angle')
for line in f:
    angleInfo = line.split()
    angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], float(angleInfo[6]), float(angleInfo[5]))
state.activateFix(angle)

#Read and activate the bonded dihedrals (Periodic and Gauss)
f = open('%s/in00_tors.xml' % wdir).readlines()
diheGauss = FixDihedralGauss(state, 'gauss')
dihePeri = FixDihedralPeriodic(state, 'periodic')
for line in f:
    diheInfo = line.split()
    if int(diheInfo[5]) == 2:
        diheGauss.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],float(diheInfo[6]) * math.pi / 180.0 + math.pi, float(diheInfo[8]), float(diheInfo[7]))
        coef = [2.0,0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef)
    else:
        coef = [float(diheInfo[6]), float(diheInfo[8])]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef)

#state.activateFix(diheGauss)
#state.activateFix(dihePeri)

#Read and generate the GoLike interactions between the protein sites
f = open('%s/in00_natv.xml' % wdir).readlines()
natv = FixBondGoLike(state, 'natv')
for line in f:
    natvInfo = line.split()
    natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])], float(natvInfo[2]), float(natvInfo[3]))

#state.activateFix(natv)

##Read and generate the BasePair interactions
#f = open('%s/in00_hbon.xml' % wdir).readlines()
#for line in f:
#    bpInfo = line.split()
#    basepair = FixBasePair3SPN2(state, 'basepair')
#    basepair.createBasePair(state.atoms[int(bpInfo[0])], state.atoms[int(bpInfo[1])], float(bpInfo[2]), float(bpInfo[3]))
#
#state.activateFix(basepair)
#
##Read and generate the angle stacking interactions
#f = open('%s/in00_base_stacking.xml' % wdir).readlines()
#for line in f:
#    stackInfo = line.split()
#    stack = FixAngleBaseStacking(state, 'stacking')
#    stack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[3])], float(bpInfo[2]), float(bpInfo[3])/2.5)
#
#state.activateFix(stack)

##Read and generate the cross stacking interactions
#f = open('%s/in00_base_stacking.xml' % wdir).readlines()
#for line in f:
#    stackInfo = line.split()
#    stack = FixAngleBaseStacking(state, 'stacking')
#    stack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[3])], float(bpInfo[2]), float(bpInfo[3])/2.5)
#
#state.activateFix(stack)


#The excluded volume nonbonded potential
#nonbond = FixLJRepul(state, 'excluded')
#nonbond.setParameter('sig', 'spc1', 'spc1', 1)
#nonbond.setParameter('eps', 'spc1', 'spc1', 1)
#state.activateFix(nonbond)


InitializeAtoms.initTemp(state, 'all', 2.479)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
fixNVT = FixLangevin(state, 'temp', 'all', 2.479)

#state.activateFix(fixNVT)


#Other things that we may want to output to make sure the simulation is running correctly
tempData = state.dataManager.recordTemperature('all', 50)
#pressureData = state.dataManager.recordPressure('all', 100)
#boundsData = state.dataManager.recordBounds(100)
engData = state.dataManager.recordEnergy('all', 50)

integRelax = IntegratorRelax(state)
integRelax.run(100000, 0.01)
integVerlet = IntegratorVerlet(state)
#writeconfig = WriteConfig(state, fn='test_out', writeEvery=100, format='xyz', handle='writer')
#state.activateWriteConfig(writeconfig)
integVerlet.run(1000)


