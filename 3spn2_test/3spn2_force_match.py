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
#state.deviceManager.setDevice(1)
state.units.setReal()
state.bounds = Bounds(state, lo = Vector(-571.389, -571.389, -571.389), hi = Vector(571.389, 571.389, 571.389))
state.rCut = 30.0
state.padding = 1.0
state.periodicInterval = 100
state.shoutEvery = 5000
state.setSpecialNeighborCoefs(0, 0, 0)

#kJ to kcal converter
kcal = 4.18
#degrees to radians
rad = math.pi / 180.0
#Temperature and salt
temp = 300
salt = 0.150
#Cutoff factor
#natoms = 2310

#Working directory name (really it's the input directory name, but I'm consistent with Gordo's scripts)
wdir = 'input_conf_2'

#Define arrays for radii and species names
sigma = []
spcs = []

#Create groups for both DNA and Protein
#state.createGroup('Protein')
#state.createGroup('DNA')

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
        sig = (sigma[i] + sigma[j])*0.89
        nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)

state.activateFix(nonbond)

#The Debye Huckel electrostatics
electric = FixChargePairDH(state, 'debyeHuckel', 'all')
electric.setParameters(temp,salt,30)
state.activateFix(electric)

#Remember the base pair identity of the particle if it is a base pair for
#later interactions
siteId = []

#Read the input configuration
f = open('%s/in00_conf.dat' % wdir).readlines()
for line in f:
    bits = line.split()
    qspcs = charge[str(bits[0])]
    siteId.append(str(bits[0]))
    if qspcs != 0:
        state.addAtom(handle=str(bits[0]), pos=Vector(float(bits[1]), float(bits[2]), float(bits[3])),q=qspcs,mol_id=int(bits[4]))
    else:
        state.addAtom(handle=str(bits[0]), pos=Vector(float(bits[1]), float(bits[2]), float(bits[3])),mol_id=int(bits[4]))

#Read the bond information
f = open('%s/in00_bond.dat' % wdir).readlines()
bondDNA  = FixBondHarmonicExtend(state, 'bondDNA')
bondProt = FixBondHarmonicExtend(state, 'bondProt')
for line in f:
    bondInfo = line.split()
    if float(bondInfo[6]) == 0:
        bondProt.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5])/kcal, 0.0, float(bondInfo[4]))
    else:
        bondDNA.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5])/kcal, float(bondInfo[5])/kcal, float(bondInfo[4]))

#Activate both kinds of bonded information
state.activateFix(bondDNA)
state.activateFix(bondProt)

#Read and activate the bonded angles
f = open('%s/in00_bend.dat' % wdir).readlines()
angle = FixAngleHarmonic(state, 'angle')
for line in f:
    angleInfo = line.split()
    angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], 2.0*float(angleInfo[6])/kcal, float(angleInfo[5])*rad)

#Activate the angles
state.activateFix(angle)

#Read the bonded dihedrals (Periodic and Gauss)
f = open('%s/in00_tors.dat' % wdir).readlines()
diheGauss = FixDihedralGauss(state, 'gauss')
dihePeri = FixDihedralPeriodic(state, 'periodic')
for line in f:
    diheInfo = line.split()
    if int(diheInfo[5]) == 2:
        #Theres only a couple forms of dihedral, the rest equal 0
        if float(diheInfo[7]) != 0:
            diheGauss.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],-float(diheInfo[6])*rad, float(diheInfo[8]), float(diheInfo[7])/kcal)
        coef = [2.0/kcal,0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,-float(diheInfo[6])*rad)
    else:
        coef = [float(diheInfo[7])/kcal, float(diheInfo[8])/kcal]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,-float(diheInfo[6])*rad)
#Active the dihedrals
state.activateFix(diheGauss)
state.activateFix(dihePeri)

#Read and generate the GoLike interactions between the protein sites
f = open('%s/in00_natv.dat' % wdir).readlines()
natv = FixBondGoLike(state, 'natv')
for line in f:
    natvInfo = line.split()
    natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])], float(natvInfo[3])/kcal, float(natvInfo[2]))

#Activate the native contacts
state.activateFix(natv)

#Read and generate the BasePair interactions
f = open('%s/in00_hbon.xml' % wdir).readlines()

#phi0, sigma, k, epsi, alpha, type, theta1, theta2
basepair = FixBasePair3SPN2(state, 'basepair')
basepair.setParameters(2.000,12.000)

for line in f:
    bpInfo = line.split()
    #Order is 0 = A, 1 = T, 2 = G, 3 = C
    if int(bpInfo[0]) == 0: #A
        epsi = 16.372218 * 0.88
        phi0 = -38.18 * rad
        theta1 = 153.17 * rad
        theta2 = 133.51 * rad
        sigma = 5.82
    elif int(bpInfo[0]) == 1: #T
        epsi = 16.372218 * 0.88
        phi0 = -38.18 * rad
        theta1 = 133.51 * rad
        theta2 = 153.17 * rad
        sigma = 5.82
    elif int(bpInfo[0]) == 2: #G
        epsi = 20.727228 * 0.88
        phi0 = -35.75 * rad
        theta1 = 159.50 * rad
        theta2 = 138.08 * rad
        sigma = 5.552
    elif int(bpInfo[0]) == 3: #C
        epsi = 20.727228 * 0.88
        phi0 = -35.75 * rad
        theta1 = 138.08 * rad
        theta2 = 159.50 * rad
        sigma = 5.552
    basepair.createBasePair(state.atoms[int(bpInfo[1])], state.atoms[int(bpInfo[2])], state.atoms[int(bpInfo[3])], state.atoms[int(bpInfo[4])], phi0, sigma, epsi/kcal, theta1, theta2)

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

#Active the base stacking interactions
state.activateFix(bstack)

##Read and generate the cross stacking interactions
crossstack = FixCrossStack3SPN2(state, 'cross-stack')
crossstack.setParameters(4.000,8.000)
f = open('%s/in00_cross_stack.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    crossStack1 = params.cross_stack_1(int(stackInfo[1]), int(stackInfo[4]), siteId)
    crossStack2 = params.cross_stack_2(int(stackInfo[3]), int(stackInfo[5]), siteId)
    crossstack.createCrossStack(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], state.atoms[int(stackInfo[3])], state.atoms[int(stackInfo[4])], state.atoms[int(stackInfo[5])], crossStack1[2], crossStack2[1],crossStack1[3]/kcal, crossStack1[1], crossStack2[0], crossStack1[0])

#Active the cross stacking
state.activateFix(crossstack)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
InitializeAtoms.initTemp(state, 'all', temp)
fixNVT = FixLangevin(state, 'temp', 'all', temp)

#500 is the damping coefficient which is ~1/gamma (gamma is 500 for lammps version of 3spn2)
fixNVT.setParameters(0,0.002)
state.activateFix(fixNVT)

#Run the actual system
#integRelax = IntegratorGradientDescent(state)
integRelax2 = IntegratorRelax(state)
integVerlet = IntegratorVerlet(state)
writeconfig = WriteConfig(state, fn='traj_2', writeEvery=10000, format='xyz', handle='writer')
writeconfig.write()
state.activateWriteConfig(writeconfig)
#integRelax.run(10000000,0.00001)
state.dt = 0.2
#integRelax2.run(200000,0.0001)
state.dt = 20.0
#integRelax2.run(100000,0.00001)
state.dt = 15.0
integVerlet.run(10000000)

