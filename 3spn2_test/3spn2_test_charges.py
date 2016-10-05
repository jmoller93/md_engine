import sys
import matplotlib.pyplot as plt
import math
import numpy as np
sys.path = sys.path + ['../build/python/build/lib.linux-x86_64-2.7']
from Sim import *

#We use real units here now
state = State()
state.deviceManager.setDevice(1)
state.units.setReal()
state.bounds = Bounds(state, lo = Vector(-571.389, -571.389, -571.389), hi = Vector(571.389, 571.389, 571.389))
state.rCut = 50.0
state.padding = 0.5
state.periodicInterval = 100
state.shoutEvery = 1000
#I believe this is correct for ps scale
state.dt = 20.0

#Working directory name (really it's the input directory name, but I'm consistent with Gordo's scripts)
wdir = 'input_conf_0'

#kJ to kcal converter
kcal = 4.18

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
charge = {}
for line in f:
    spcsInfo = line.split()
    sigma.append(float(spcsInfo[2]))
    spcs.append(str(spcsInfo[1]))
    charge[str(spcsInfo[1])] = float(spcsInfo[4])
    state.atomParams.addSpecies(handle=str(spcsInfo[1]), mass=float(spcsInfo[3]))

#The excluded volume nonbonded potential
nonbond = FixLJRepul(state, 'excluded')

for i in range(len(spcs)):
    for j in range(len(spcs)):
        sig = (sigma[i] + sigma[j])*0.5
        nonbond.setParameter('rCut', spcs[i], spcs[j], sig)
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)
        nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)

state.activateFix(nonbond)

#The fake electric charges for now
electric = FixChargePairDH(state, 'debyeHuckel', 'all')
electric.setParameters(300,0.150)
state.activateFix(electric)

#Remember the base pair identity of the particle if it is a base pair for
#later interactions
siteId = []

#Read the input configuration
f = open('%s/in00_conf_eq.xml' % wdir).readlines()
for line in f:
    bits = line.split()
    qspcs = charge[str(bits[0])]
    siteId.append(str(bits[0]))
    if qspcs != 0:
        state.addAtom(handle=str(bits[0]), pos=Vector(float(bits[1])*10.0, float(bits[2])*10.0, float(bits[3])*10.0),q=qspcs)
    else:
        state.addAtom(handle=str(bits[0]), pos=Vector(float(bits[1])*10.0, float(bits[2])*10.0, float(bits[3])*10.0))

#Read the bond information
f = open('%s/in00_bond.xml' % wdir).readlines()
bondDNA  = FixBondHarmonicExtend(state, 'bondDNA')
bondProt = FixBondHarmonic(state, 'bondProt')
for line in f:
    bondInfo = line.split()
    if float(bondInfo[6]) == 0:
        bondProt.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], 2.0*float(bondInfo[5])/kcal, float(bondInfo[4]))
    else:
        bondDNA.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5])/kcal, float(bondInfo[4]))

state.activateFix(bondDNA)
state.activateFix(bondProt)

#Read and activate the bonded angles
f = open('%s/in00_bend.xml' % wdir).readlines()
angle = FixAngleHarmonic(state, 'angle')
for line in f:
    angleInfo = line.split()
    angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], float(angleInfo[6])/kcal, float(angleInfo[5])*math.pi/180.0)
state.activateFix(angle)

#Read and activate the bonded dihedrals (Periodic and Gauss)
f = open('%s/in00_tors.xml' % wdir).readlines()
diheGauss = FixDihedralGauss(state, 'gauss')
dihePeri = FixDihedralPeriodic(state, 'periodic')
for line in f:
    diheInfo = line.split()
    if int(diheInfo[5]) == 2:
        #Theres only a couple forms of dihedral, the rest equal 0
        if float(diheInfo[7]) != 0:
            diheGauss.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],float(diheInfo[6]) * math.pi / 180.0, float(diheInfo[8]), float(diheInfo[7])/kcal)

        coef = [2.0/kcal,0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,float(diheInfo[6]) * math.pi /180.0)
    else:
        coef = [float(diheInfo[7])/kcal, float(diheInfo[8])/kcal]
#        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,float(diheInfo[6]) * math.pi /180.0)

state.activateFix(diheGauss)
state.activateFix(dihePeri)

#Read and generate the GoLike interactions between the protein sites
f = open('%s/in00_natv.xml' % wdir).readlines()
natv = FixBondGoLike(state, 'natv')
for line in f:
    natvInfo = line.split()
    natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])], float(natvInfo[3])/kcal, float(natvInfo[2]))


state.activateFix(natv)

#All base stacking interactions use functions to help them parse correct inputs

#First up is intrastrand base stacking
#This generates the epsilon values from the intrastrand interactions


################################################################
#I think I want the with curvature variant....
def base_stack(atom1, atom2, atom3, siteId):
    if siteId[atom2] == 'A':
        if siteId[atom3] == 'A':
            return [13.810,3.58,100.13]
        elif siteId[atom3] == 'T':
            return [15.050,3.56,90.48]
        elif siteId[atom3] == 'G':
            return [13.320,3.85,104.39]
        elif siteId[atom3] == 'C':
            return [15.820,3.45,93.23]
    elif siteId[atom2] == 'T':
        if siteId[atom3] == 'A':
            return [9.250,4.15,102.59]
        elif siteId[atom3] == 'T':
            return [12.44,3.93,93.32]
        elif siteId[atom3] == 'G':
            return [9.58,4.32,103.70]
        elif siteId[atom3] == 'C':
            return [13.11,3.87,94.55]
    elif siteId[atom2] == 'G':
        if siteId[atom3] == 'A':
            return [13.76,3.51,95.45]
        elif siteId[atom3] == 'T':
            return [14.59,3.47,87.63]
        elif siteId[atom3] == 'G':
            return [14.77,3.67,106.36]
        elif siteId[atom3] == 'C':
            return [15.17,3.42,83.12]
    elif siteId[atom2] == 'C':
        if siteId[atom3] == 'A':
            return [9.25,4.15,102.69]
        elif siteId[atom3] == 'T':
            return [12.42,3.99,96.05]
        elif siteId[atom3] == 'G':
            return [8.83,4.34,100.46]
        elif siteId[atom3] == 'C':
            return [14.01,3.84,100.68]


#Parameters for the a-b-e first cross stacking interaction
def cross_stack_1(atom2, atom5, siteId):
    if siteId[atom2] == 'A':
        theta3 = 110.92
        if siteId[atom5] == 'A':
            theta1 = 154.04
            sigma = 6.420
            eps = 2.139
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'T':
            theta1 = 158.77
            sigma = 6.770
            eps = 2.714
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'G':
            theta1 = 153.88
            sigma = 6.270
            eps = 2.772
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'C':
            theta1 = 157.69
            sigma = 6.840
            eps = 1.909
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
    elif siteId[atom2] == 'T':
        theta3 = 110.92
        if siteId[atom5] == 'A':
            theta1 = 148.62
            sigma = 6.770
            eps = 2.714
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'T':
            theta1 = 155.05
            sigma = 7.210
            eps = 2.139
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'G':
            theta1 = 147.54
            sigma = 6.530
            eps = 2.485
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'C':
            theta1 = 153.61
            sigma = 7.080
            eps = 2.916
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
    elif siteId[atom2] == 'G':
        theta3 = 120.45
        if siteId[atom5] == 'A':
            theta1 = 153.91
            sigma = 6.270
            eps = 2.772
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'T':
            theta1 = 155.72
            sigma = 6.530
            eps = 2.485
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'G':
            theta1 = 151.84
            sigma = 5.740
            eps = 3.693
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'C':
            theta1 = 157.80
            sigma = 6.860
            eps = 1.104
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
    elif siteId[atom2] == 'C':
        theta3 = 120.45
        if siteId[atom5] == 'A':
            theta1 = 152.04
            sigma = 6.840
            eps = 1.909
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'T':
            theta1 = 157.72
            sigma = 7.080
            eps = 2.916
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'G':
            theta1 = 151.65
            sigma = 6.860
            eps = 1.104
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom5] == 'C':
            theta1 = 154.49
            sigma = 6.790
            eps = 4.699
            return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]


#Second cross stacking interaction parameters
def cross_stack_2(atom4, atom6, siteId):
    if siteId[atom4] == 'A':
        if siteId[atom6] == 'A':
            theta2 = 116.34
            sigma = 5.580
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'T':
            theta2 = 119.61
            sigma = 6.140
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'G':
            theta2 = 115.19
            sigma = 5.630
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'C':
            theta2 = 120.92
            sigma = 6.180
            return [theta2 * math.pi /180.0, sigma]
    elif siteId[atom4] == 'T':
        if siteId[atom6] == 'A':
            theta2 = 207.40
            sigma = 6.140
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'T':
            theta2 = 110.76
            sigma = 6.800
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'G':
            theta2 = 106.33
            sigma = 6.070
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'C':
            theta2 = 111.57
            sigma = 6.640
            return [theta2 * math.pi /180.0, sigma]
    elif siteId[atom4] == 'G':
        if siteId[atom6] == 'A':
            theta2 = 121.61
            sigma = 5.630
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'T':
            theta2 = 124.92
            sigma = 6.070
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'G':
            theta2 = 120.52
            sigma = 5.870
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'C':
            theta2 = 124.88
            sigma = 5.660
            return [theta2 * math.pi /180.0, sigma]
    elif siteId[atom4] == 'C':
        if siteId[atom6] == 'A':
            theta2 = 112.45
            sigma = 6.180
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'T':
            theta2 = 115.43
            sigma = 6.640
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'G':
            theta2 = 110.51
            sigma = 5.660
            return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom6] == 'C':
            theta2 = 115.80
            sigma = 6.800
            return [theta2 * math.pi /180.0, sigma]

################################################################3

#Read and generate the BasePair interactions
f = open('%s/in00_hbon.xml' % wdir).readlines()

#phi0, sigma, k, epsi, alpha, type, theta1, theta2
basepair = FixBasePair3SPN2(state, 'basepair')
basepair.setParameters(2.000,12.000)
rad = math.pi / 180.0
for line in f:
    bpInfo = line.split()
    #Order is 0 = A, 1 = T, 2 = G, 3 = C
    if int(bpInfo[0]) == 0: #A
        epsi = 16.372218 * 0.861
        phi0 = -38.18 * rad
        theta1 = 153.17 * rad
        theta2 = 133.51 * rad
        sigma = 5.82
    elif int(bpInfo[0]) == 1: #T
        epsi = 16.372218 * 0.861
        phi0 = -38.18 * math.pi / 180.0
        theta1 = 133.51 * rad
        theta2 = 153.17 * rad
        sigma = 5.82
    elif int(bpInfo[0]) == 2: #G
        epsi = 20.727228 * 0.861
        phi0 = -35.75 * math.pi / 180.0
        theta1 = 159.50 * rad
        theta2 = 138.08 * rad
        sigma = 5.552
    elif int(bpInfo[0]) == 3: #C
        epsi = 20.727228 * 0.861
        phi0 = -35.75 * math.pi / 180.0
        theta1 = 138.08 * rad
        theta2 = 159.50 * rad
        sigma = 5.552
    #print(siteId[int(bpInfo[1])])
    #print(siteId[int(bpInfo[2])])
    #print(sigma)
    basepair.createBasePair(state.atoms[int(bpInfo[1])], state.atoms[int(bpInfo[2])], state.atoms[int(bpInfo[3])], state.atoms[int(bpInfo[4])], phi0, sigma, epsi/kcal, theta1, theta2)

state.activateFix(basepair)

bstack = FixAngleBaseStacking(state, 'intrastacking')
bstack.setParameters(3.000,6.000)
#Read and generate the angle stacking interactions
f = open('%s/in00_base_stack.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    baseStack = base_stack(int(stackInfo[0]), int(stackInfo[1]), int(stackInfo[2]), siteId)
    bstack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], baseStack[2] * math.pi/180.0, baseStack[0]/kcal , baseStack[1])

state.activateFix(bstack)

##Read and generate the cross stacking interactions
cstack = FixCrossStack3SPN2(state, 'cross-stack')
cstack.setParameters(4.000,8.000)
f = open('%s/in00_cross_stack.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    crossStack1 = cross_stack_1(int(stackInfo[1]), int(stackInfo[4]), siteId)
    crossStack2 = cross_stack_2(int(stackInfo[3]), int(stackInfo[5]), siteId)
    cstack.createCrossStack(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], state.atoms[int(stackInfo[3])], state.atoms[int(stackInfo[4])], state.atoms[int(stackInfo[5])], crossStack1[2], crossStack2[1], crossStack1[3]/kcal, crossStack1[1], crossStack2[0], crossStack1[0])

state.activateFix(cstack)

#Run the relaxation integrator
#integRelax = IntegratorGradientDescent(state)
#integRelax.run(10000,0.001)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
InitializeAtoms.initTemp(state, 'all', 300)
fixNVT = FixLangevin(state, 'temp', 'all', 300)
fixNVT.setParameters(0,0.01)
state.activateFix(fixNVT)

#Other things that we may want to output to make sure the simulation is running correctly
tempData = state.dataManager.recordTemperature('all', 50)
#pressureData = state.dataManager.recordPressure('all', 100)
#boundsData = state.dataManager.recordBounds(100)
#engData = state.dataManager.recordEnergy('all', 50)

#Run the actual system
state.dt = 20.0
integVerlet = IntegratorVerlet(state)
writeconfig = WriteConfig(state, fn='test_out', writeEvery=1000, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integRelax = IntegratorGradientDescent(state)
integRelax.run(100000,0.1)
integVerlet.run(100000)



