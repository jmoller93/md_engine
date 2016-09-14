import sys
import matplotlib.pyplot as plt
import math
import numpy as np
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
wdir = 'input_3spn2'

#Define arrays for radii and species names
sigma = []
spcs = []

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
    sigma.append(float(spcsInfo[2]))
    spcs.append(str(spcsInfo[1]))
    state.atomParams.addSpecies(handle=str(spcsInfo[1]), mass=float(spcsInfo[3]))

#The excluded volume nonbonded potential
nonbond = FixLJRepul(state, 'excluded')

for i in range(len(spcs)):
    for j in range(len(spcs)):
        sig = (sigma[i] + sigma[j])*0.5
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)
        nonbond.setParameter('eps', spcs[i], spcs[j], 1)

state.activateFix(nonbond)

#Remember the base pair identity of the particle if it is a base pair for
#later interactions
siteId = []

#Read the input configuration
f = open('%s/in00_conf.xml' % wdir).readlines()
for line in f:
    bits = line.split()
    siteId.append(str(bits[0]))
    state.addAtom(str(bits[0]), Vector(float(bits[1]), float(bits[2]), float(bits[3])))

#Read the bond information
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
        coef = [float(diheInfo[7]), float(diheInfo[8])]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef)

#state.activateFix(diheGauss)
state.activateFix(dihePeri)

#Read and generate the GoLike interactions between the protein sites
f = open('%s/in00_natv.xml' % wdir).readlines()
natv = FixBondGoLike(state, 'natv')
for line in f:
    natvInfo = line.split()
    natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])], float(natvInfo[2]), float(natvInfo[3]))

#state.activateFix(natv)

#All base stacking interactions use functions to help them parse correct inputs

#First up is intrastrand base stacking
#This generates the epsilon values from the intrastrand interactions
def base_stack(atom1, atom2, atom3, siteId):
    if siteId[atom2] == 'A':
        if siteId[atom3] == 'A':
            return 13.810
        elif siteId[atom3] == 'C':
            return 15.820
        elif siteId[atom3] == 'G':
            return 13.320
        elif siteId[atom3] == 'T':
    elif siteId[atom2] == 'C':
        if siteId[atom3] == 'A':
            return 9.250
        elif siteId[atom3] == 'C':
            return 14.010
        elif siteId[atom3] == 'G':
            return 8.830
        elif siteId[atom3] == 'T':
            return 12.420
    elif siteId[atom2] == 'T':
        if siteId[atom3] == 'A':
            return 9.150
        elif siteId[atom3] == 'C':
            return 13.110
        elif siteId[atom3] == 'G':
            return 9.580
        elif siteId[atom3] == 'T':
            return 12.440
    elif siteId[atom2] == 'G':
        if siteId[atom3] == 'A':
            return 13.760
        elif siteId[atom3] == 'C':
            return 15.170
        elif siteId[atom3] == 'G':
            return 14.770
        elif siteId[atom3] == 'T':
            return 14.590

#Read and generate the BasePair interactions
f = open('%s/in00_hbon.xml' % wdir).readlines()

#phi0, sigma, k, epsi, alpha, type, theta1, theta2
basepair = FixBasePair3SPN2(state, 'basepair')
rad = math.pi / 180.0
for line in f:
    bpInfo = line.split()
    if int(bpInfo[0]) == 0:
        epsi = 16.372218 * 0.88
        phi0 = -38.35 * rad
        theta1 = 156.54 * rad
        theta2 = 135.78 * rad
    elif(int(bpInfo[0]) == 1:
        epsi = 16.372218 * 0.88
        phi0 = -38.35 * math.pi / 180.0
        theta1 = 135.78 * rad
        theta2 = 156.54 * rad
    elif int(bpInfo[0]) == 2:
        epsi = 20.727228 * 0.88
        phi0 = -45.81 * math.pi / 180.0
        theta1 = 154.62 * rad
        theta2 = 152.74 * rad
    elif(int(bpInfo[0]) == 3:
        epsi = 16.372218 * 0.88
        phi0 = -45.81 * math.pi / 180.0
        theta1 = 152.74 * rad
        theta2 = 154.62 * rad

    basepair.createBasePair(state.atoms[int(bpInfo[1])], state.atoms[int(bpInfo[2])], (state.atoms[int(bpInfo[3])], state.atoms[int(bpInfo[4])], phi0, epsi, 2.000, 12.000, theta1, theta2)

state.activateFix(basepair)

#Read and generate the angle stacking interactions
f = open('%s/in00_base_stacking.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    stack = FixAngleBaseStacking(state, 'stacking')
    baseStack = base_stack(int(stackInfo[0]), int(stackInfo[1]), int(stackInfo[2]), siteId)
    stack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], baseStack, 3.000, 6.000)

state.activateFix(stack)

##Read and generate the cross stacking interactions
#f = open('%s/in00_base_stacking.xml' % wdir).readlines()
#for line in f:
#    stackInfo = line.split()
#    stack = FixAngleBaseStacking(state, 'stacking')
#    stack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[3])], float(bpInfo[2]), float(bpInfo[3])/2.5)
#
#state.activateFix(stack)

InitializeAtoms.initTemp(state, 'all', 2.479)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
fixNVT = FixLangevin(state, 'temp', 'all', 2.479)

state.activateFix(fixNVT)


#Other things that we may want to output to make sure the simulation is running correctly
tempData = state.dataManager.recordTemperature('all', 50)
#pressureData = state.dataManager.recordPressure('all', 100)
#boundsData = state.dataManager.recordBounds(100)
engData = state.dataManager.recordEnergy('all', 50)

integRelax = IntegratorRelax(state)
integRelax.run(1000000, 10)
integVerlet = IntegratorVerlet(state)
writeconfig = WriteConfig(state, fn='test_out', writeEvery=10, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integVerlet.run(100)


