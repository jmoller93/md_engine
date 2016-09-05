import sys
import matplotlib.pyplot as plt
sys.path = sys.path + ['../build/python/build/lib.linux-x86_64-2.7']
#from Sim import *
from Sim import *
state = State()
state.deviceManager.setDevice(0)
state.bounds = Bounds(state, lo = Vector(0, 0, 0), hi = Vector(1142.778000,1142.778000,1142.778000))
state.rCut = 3.0
state.padding = 0.6
state.periodicInterval = 7
state.shoutEvery = 100

wdir = 'input_dihedral'

f = open('%s/in00_spcs.xml' % wdir).readlines()
for line in f:
    spcsInfo = line.split()
    state.atomParams.addSpecies(handle=str(spcsInfo[1]), mass=float(spcsInfo[3]), atomicNum=float(spcsInfo[0]))

f = open('%s/in00_conf.xml' % wdir).readlines()
for line in f:
    bits = line.split()
    state.addAtom(str(bits[0]), Vector(float(bits[1]), float(bits[2]), float(bits[3])))

#The repulsive force may be a challenge
#nonbond = FixLJCut(state, 'cut')
#nonbond.setParameter('sig', 'spc1', 'spc1', 1)
#nonbond.setParameter('eps', 'spc1', 'spc1', 1)
#state.activateFix(nonbond)

f = open('%s/in00_bond.xml' % wdir).readlines()
for line in f:
    bondInfo = line.split()
    bond = FixBondHarmonic(state, 'bond')
    bond.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5]), float(bondInfo[4]))

state.activateFix(bond)

f = open('%s/in00_bend.xml' % wdir).readlines()
for line in f:
    angleInfo = line.split()
    angle = FixAngleHarmonic(state, 'angle')
    angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], float(angleInfo[5]), float(angleInfo[4]))

state.activateFix(angle)

f = open('%s/in00_tors.xml' % wdir).readlines()
for line in f:
    diheInfo = line.split()
    print(line)
    dihe = FixDihedralGauss(state, 'gauss')
    dihe.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],float(diheInfo[6]), float(diheInfo[8]), float(diheInfo[7]))

#state.activateFix(dihe)

#Not yet implemented
#natv = FixGoLike()
#f = open('%s/in00_natv.xml' % wdir).readlines()
#for line in f:
#    natvInfo = line.split()
#    natv = FixBondGoLike(state, 'natv')
#    natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])])

InitializeAtoms.initTemp(state, 'all', 1.2)

fixNVT = FixLangevin(state, 'temp', 'all', 1.2)

state.activateFix(fixNVT)

integVerlet = IntegratorVerlet(state)

tempData = state.dataManager.recordTemperature('all', 100)
#pressureData = state.dataManager.recordPressure('all', 100)
#boundsData = state.dataManager.recordBounds(100)
engData = state.dataManager.recordEnergy('all', 100)

writeconfig = WriteConfig(state, fn='test_out', writeEvery=1, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integVerlet.run(100)
