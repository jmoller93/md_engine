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
state.dt = 20.0
state.setSpecialNeighborCoefs(0, 0, 0)

#kJ to kcal converter
kcal = 4.18
#degrees to radians
rad = math.pi / 180.0
#Temperature and salt
temp = 300
salt = 0.150
#Cutoff factor
natoms = 2310

#Working directory name (really it's the input directory name, but I'm consistent with Gordo's scripts)
wdir = 'input_conf_0'

#Define arrays for radii and species names
sigma = []
spcs = []

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
        #if len(spcs[i]) == 3 and len(spcs[j]) == 3:
        #    nonbond.setParameter('eps', spcs[i], spcs[j], 0.0)
        #else:
        #    nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)
        #nonbond.setParameter('rCut', spcs[i], spcs[j], sig)
        nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)

#The Debye Huckel electrostatics
electric = FixChargePairDH(state, 'debyeHuckel', 'all')
electric.setParameters(temp,salt,30)

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
bondProt = FixBondHarmonicExtend(state, 'bondProt')
for line in f:
    bondInfo = line.split()
    if float(bondInfo[6]) == 0:
        #bondProt.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5])/kcal, 0.0, float(bondInfo[4]))
        bondProt.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], 0.0, 0.0, float(bondInfo[4]))
    else:
        #bondDNA.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], float(bondInfo[5])/kcal, float(bondInfo[5])/kcal, float(bondInfo[4]))
        bondDNA.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], 0.0, 0.0, float(bondInfo[4]))

#Read and activate the bonded angles
f = open('%s/in00_bend.xml' % wdir).readlines()
angle = FixAngleHarmonic(state, 'angle')
for line in f:
    angleInfo = line.split()
    #angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], 2.0*float(angleInfo[6])/kcal, float(angleInfo[5])*rad)
    angle.createAngle(state.atoms[int(angleInfo[1])], state.atoms[int(angleInfo[2])], state.atoms[int(angleInfo[3])], 0.0, float(angleInfo[5])*rad)

#Read the bonded dihedrals (Periodic and Gauss)
f = open('%s/in00_tors.xml' % wdir).readlines()
diheGauss = FixDihedralGauss(state, 'gauss')
dihePeri = FixDihedralPeriodic(state, 'periodic')
for line in f:
    diheInfo = line.split()
    if int(diheInfo[5]) == 2:
        #Theres only a couple forms of dihedral, the rest equal 0
        if float(diheInfo[7]) != 0:
            diheGauss.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],-float(diheInfo[6])*rad, float(diheInfo[8]), float(diheInfo[7])/kcal)
        #coef = [2.0/kcal,0.0]
        coef = [0.0,0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,-float(diheInfo[6])*rad)
    else:
        #coef = [float(diheInfo[7])/kcal, float(diheInfo[8])/kcal]
        coef = [0.0, 0.0]
        dihePeri.createDihedral(state.atoms[int(diheInfo[1])],state.atoms[int(diheInfo[2])],state.atoms[int(diheInfo[3])],state.atoms[int(diheInfo[4])],coef,-float(diheInfo[6])*rad)

#Read and generate the GoLike interactions between the protein sites
f = open('%s/in00_natv.xml' % wdir).readlines()
natv = FixBondGoLike(state, 'natv')
for line in f:
    natvInfo = line.split()
    #natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])], float(natvInfo[3])/kcal, float(natvInfo[2]))
    natv.createBond(state.atoms[int(natvInfo[0])], state.atoms[int(natvInfo[1])], 0.0, float(natvInfo[2]))

#Read and generate the BasePair interactions
f = open('%s/in00_hbon_2.xml' % wdir).readlines()

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
   # basepair.createBasePair(state.atoms[int(bpInfo[1])], state.atoms[int(bpInfo[2])], state.atoms[int(bpInfo[3])], state.atoms[int(bpInfo[4])], phi0, sigma, epsi/kcal, theta1, theta2)
    basepair.createBasePair(state.atoms[int(bpInfo[1])], state.atoms[int(bpInfo[2])], state.atoms[int(bpInfo[3])], state.atoms[int(bpInfo[4])], phi0, sigma, 0.0, theta1, theta2)

bstack = FixAngleBaseStacking(state, 'intrastacking')
bstack.setParameters(3.000,6.000)
#Read and generate the angle stacking interactions
f = open('%s/in00_base_stack.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    baseStack = params.base_stack(int(stackInfo[0]), int(stackInfo[1]), int(stackInfo[2]), siteId)
    #bstack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], baseStack[2]*rad, baseStack[0]/kcal , baseStack[1])
    bstack.createAngle(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], baseStack[2]*rad, 0.0, baseStack[1])

##Read and generate the cross stacking interactions
crossstack = FixCrossStack3SPN2(state, 'cross-stack')
crossstack.setParameters(4.000,8.000)
f = open('%s/in00_cross_stack_2.xml' % wdir).readlines()
for line in f:
    stackInfo = line.split()
    crossStack1 = params.cross_stack_1(int(stackInfo[1]), int(stackInfo[4]), siteId)
    crossStack2 = params.cross_stack_2(int(stackInfo[3]), int(stackInfo[5]), siteId)
    #crossstack.createCrossStack(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], state.atoms[int(stackInfo[3])], state.atoms[int(stackInfo[4])], state.atoms[int(stackInfo[5])], crossStack1[2], crossStack2[1],crossStack1[3]/kcal, crossStack1[1], crossStack2[0], crossStack1[0])
    crossstack.createCrossStack(state.atoms[int(stackInfo[0])], state.atoms[int(stackInfo[1])], state.atoms[int(stackInfo[2])], state.atoms[int(stackInfo[3])], state.atoms[int(stackInfo[4])], state.atoms[int(stackInfo[5])], crossStack1[2], crossStack2[1],0.0, crossStack1[1], crossStack2[0], crossStack1[0])

InitializeAtoms.initTemp(state, 'all', temp)

#Activate all the forces
state.activateFix(electric)
#state.activateFix(nonbond)
state.activateFix(bondDNA)
state.activateFix(bondProt)
state.activateFix(angle)
state.activateFix(diheGauss)
state.activateFix(dihePeri)
state.activateFix(natv)
state.activateFix(bstack)
state.activateFix(basepair)
state.activateFix(crossstack)

#Run the actual system
writeconfig = WriteConfig(state, fn='force_match', writeEvery=1000, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integVerlet = IntegratorVerlet(state)
integVerlet.run(0)

i = 0
fnme = "force.dat"
f = open(fnme,'w')
for atom in state.atoms:
    f.write("Atom %d\n %f %f %f\n" % (i,atom.force[0],atom.force[1],atom.force[2]))
    i = i+1
