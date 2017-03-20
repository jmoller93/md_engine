#!/usr/bin/python

import sys,os,glob
import matplotlib.pyplot as plt
import math
import numpy as np
sys.path = sys.path + ['../build/python/build/lib.linux-x86_64-2.7']
from Sim import *
import Params

#Create parameters instance
params = Params.Params()

#We use real units here now (fs, ang, kcal)
state = State()
state.deviceManager.setDevice(1)
state.units.setReal()
state.bounds = Bounds(state, lo = Vector(-571.389, -571.389, -571.389), hi = Vector(571.389, 571.389, 571.389))
state.rCut = 50.0
state.padding = 0.5
state.periodicInterval = 100
state.shoutEvery = 1000
state.dt = 20.0
state.setSpecialNeighborCoefs(0, 0, 0, 0, 0)

#kJ to kcal converter
kcal = 4.18
#degrees to radians
rad = math.pi / 180.0
natoms = 11


#Working directory name (really it's the input directory name, but I'm consistent with Gordo's scripts)
wdir = 'input_base_pair'

#Define arrays for radii and species names
sigma = []
spcs = []

#This is a nice tool that allows for debugging/ printing without too much slowdown
def operate(turn):
    fnme = open('force_2.dat','w');
    i = 0
    for atom in state.atoms:
        fnme.write("Atom %d\n Force is: %f\t%f\t%f\n" % (i,atom.force[0],atom.force[1],atom.force[2]))
        i = i+1
    fnme.close()

#Activate the python function as a "Fix"

pyOp = PythonOperation(operateEvery=1, handle='force_test', operation=operate)
state.activatePythonOperation(pyOp)

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
nonbond = FixWCA(state, 'excluded')

for i in range(len(spcs)):
    for j in range(len(spcs)):
        sig = (sigma[i] + sigma[j])
        ##nonbond.setParameter('rCut', spcs[i], spcs[j], sig)
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)
        nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)

state.activateFix(nonbond)

#The Debye Huckel electrostatics
electric = FixChargePairDH(state, 'debyeHuckel', 'all')
electric.setParameters(100,0.150)
#state.activateFix(electric)

#Remember the base pair identity of the particle if it is a base pair for
#later interactions
siteId = []

#Read the input configuration
f = open('%s/in00_conf.xml' % wdir).readlines()
for line in f:
    bits = line.split()
    qspcs = charge[str(bits[0])]
    siteId.append(str(bits[0]))
    if qspcs != 0:
        state.addAtom(handle=str(bits[0]), pos=Vector(float(bits[1]), float(bits[2]), float(bits[3])),q=qspcs)
    else:
        state.addAtom(handle=str(bits[0]), pos=Vector(float(bits[1]), float(bits[2]), float(bits[3])))

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

#Activate both kinds of bonded information
state.activateFix(bondDNA)
#state.activateFix(bondProt)

#Read and activate the bonded angles
f = open('%s/in00_bend.xml' % wdir).readlines()
angle = FixAngleHarmonic(state, 'angle')
for line in f:
    angleInfo = line.split()
    angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], 2.0*float(angleInfo[6])/kcal, float(angleInfo[5])*rad)

#Activate the angles
state.activateFix(angle)

#Read the bonded dihedrals (Periodic and Gauss)
f = open('%s/in00_tors.xml' % wdir).readlines()
diheGauss = FixDihedralGauss(state, 'gauss')
dihePeri = FixDihedralPeriodic(state, 'periodic')
for line in f:
    diheInfo = line.split()
    if int(diheInfo[5]) == 2:
        #Theres only a couple forms of dihedral, the rest equal 0
        if float(diheInfo[7]) != 0:
            diheGauss.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])], -float(diheInfo[6])*rad, float(diheInfo[8]), float(diheInfo[7])/kcal)

        coef = [2.0/kcal,0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,-float(diheInfo[6])*rad)
    else:
        coef = [float(diheInfo[7])/kcal, float(diheInfo[8])/kcal]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,-float(diheInfo[6])*rad)

#Active the dihedrals
state.activateFix(diheGauss)
state.activateFix(dihePeri)

#Read and generate the BasePair interactions
f = open('%s/in00_hbon.xml' % wdir).readlines()

#phi0, sigma, k, epsi, alpha, type, theta1, theta2
basepair = FixBasePair3SPN2(state, 'basepair')
basepair.setParameters(2.000,12.000)

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
        phi0 = -38.18 * rad
        theta1 = 133.51 * rad
        theta2 = 153.17 * rad
        sigma = 5.82
    elif int(bpInfo[0]) == 2: #G
        epsi = 20.727228 * 0.861
        phi0 = -35.75 * rad
        theta1 = 159.50 * rad
        theta2 = 138.08 * rad
        sigma = 5.552
    elif int(bpInfo[0]) == 3: #C
        epsi = 20.727228 * 0.861
        phi0 = -35.75 * rad
        theta1 = 138.08 * rad
        theta2 = 159.50 * rad
        sigma = 5.552
    basepair.createBasePair(state.atoms[int(bpInfo[1])], state.atoms[int(bpInfo[2])], state.atoms[int(bpInfo[3])], state.atoms[int(bpInfo[4])], phi0, sigma, epsi/kcal, theta1, theta2)
    #basepair.createBasePair(state.atoms[int(bpInfo[1])], state.atoms[int(bpInfo[2])], state.atoms[int(bpInfo[3])], state.atoms[int(bpInfo[4])], phi0, sigma, 100.0/kcal, theta1, theta2)

#Activate the base pair interactions
state.activateFix(basepair)

bstack = FixAngleBaseStacking(state, 'intrastacking')
bstack.setParameters(3.000,6.000)
#Read and generate the angle stacking interactions
f = open('%s/in00_base_stack.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    baseStack = params.base_stack(int(stackInfo[0]), int(stackInfo[1]), int(stackInfo[2]), siteId)
    bstack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], baseStack[2]*rad, baseStack[0]/kcal , baseStack[1])
    #bstack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], baseStack[2]*rad, 50.0/kcal , baseStack[1])

#Active the base stacking interactions
state.activateFix(bstack)

##Read and generate the cross stacking interactions
cstack = FixCrossStack3SPN2(state, 'cross-stack')
cstack.setParameters(4.000,8.000)
f = open('%s/in00_cross_stack.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    crossStack1 = params.cross_stack_1(int(stackInfo[1]), int(stackInfo[4]), siteId)
    crossStack2 = params.cross_stack_2(int(stackInfo[3]), int(stackInfo[5]), siteId)
    cstack.createCrossStack(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], state.atoms[int(stackInfo[3])], state.atoms[int(stackInfo[4])], state.atoms[int(stackInfo[5])], crossStack1[2], crossStack2[1], crossStack1[3]/kcal, crossStack1[1], crossStack2[0], crossStack1[0])
    #cstack.createCrossStack(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], state.atoms[int(stackInfo[3])], state.atoms[int(stackInfo[4])], state.atoms[int(stackInfo[5])], crossStack1[2], crossStack2[1], 50.0/kcal, crossStack1[1], crossStack2[0], crossStack1[0])

#Active the cross stacking
state.activateFix(cstack)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
InitializeAtoms.initTemp(state, 'all', 100.0)
fixNVT = FixLangevin(state, 'temp', 'all', 100.0)

#500 is the damping coefficient which is ~1/gamma (gamma is 500 for lammps version of 3spn2)
fixNVT.setParameters(0,0.002)
state.activateFix(fixNVT)

#Other things that we may want to output to make sure the simulation is running correctly
tempData = state.dataManager.recordTemperature('all', 50)

#Run the actual system
integVerlet = IntegratorVerlet(state)
integRelax = IntegratorGradientDescent(state)
integRelax2 = IntegratorRelax(state)
#integRelax.run(10000,0.001)
#integRelax2.run(10000,1e-4)
writeconfig = WriteConfig(state, fn='test_out', writeEvery=100, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integVerlet.run(1)
exit()

