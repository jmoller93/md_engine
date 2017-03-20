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
state.bounds = Bounds(state, lo = Vector(-200.389, -200.389, -200.389), hi = Vector(200.389, 200.389, 200.389))
state.rCut = 50.0
state.padding = 0.5
state.periodicInterval = 100
state.shoutEvery = 1000
state.dt = 20.0
#Note!!! Timestep has to be reduced because it seems to crash at a dt of 20 with this code
#NOTE!!! Dihedrals for the GPU code are defined the opposite of that in the 3SPN code.
#All values should be opposite, but do in fact work!

#kJ to kcal converter
kcal = 4.18
#degrees to radians
rad = math.pi / 180.0


#Working directory name (really it's the input directory name, but I'm consistent with Gordo's scripts)
wdir = 'input_base_pair'

#Define arrays for radii and species names
sigma = []
spcs = []

#This is a nice tool that allows for debugging/ printing without too much slowdown
def force(turn):
    fnme = open('force_test.dat', 'a+')
    sum_force = 0.0
    for i in range(0,3):
        for j in range(0,4):
            sum_force = sum_force +state.atoms[j].force[i]
    fnme.write("%10.10f\n" % sum_force)
    fnme.close()

pyOp = PythonOperation(operateEvery=100, handle='force_test', operation=force)
state.activatePythonOperation(pyOp)

state.setSpecialNeighborCoefs(0, 0, 0)

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
        sig = (sigma[i] + sigma[j])*0.89
        nonbond.setParameter('rCut', spcs[i], spcs[j], sig)
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)
        nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)

#state.activateFix(nonbond)

#The Debye Huckel electrostatics
electric = FixChargePairDH(state, 'debyeHuckel', 'all')
electric.setParameters(300,0.150)
#state.activateFix(electric)

#Remember the base pair identity of the particle if it is a base pair for
#later interactions
siteId = []

#Create atoms
state.addAtom(handle='P', pos=Vector(0,0,0),q=-0.6)
state.addAtom(handle='S', pos=Vector(0,0,5))
state.addAtom(handle='P', pos=Vector(0,0,11),q=-0.6)
state.addAtom(handle='S', pos=Vector(0,0,16))


#Read the bond information
bondDNA  = FixBondHarmonicExtend(state, 'bondDNA')
#bondProt = FixBondHarmonic(state, 'bondProt')

bondDNA.createBond(state.atoms[0], state.atoms[1], 1.2/kcal, 4.8635)
#bondDNA.createBond(state.atoms[1], state.atoms[2], 0.6/kcal, 4.8635)
bondDNA.createBond(state.atoms[2], state.atoms[3], 1.2/kcal, 4.8635)
#Activate both kinds of bonded information
state.activateFix(bondDNA)
#state.activateFix(bondProt)

angleDNA = FixAngleHarmonic(state, 'angleDNA')
angleDNA.createAngle(state.atoms[0],state.atoms[1],state.atoms[2],2.0*355.0/kcal,94.784734*rad)
angleDNA.createAngle(state.atoms[1],state.atoms[2],state.atoms[3],2.0*355.0/kcal,108.833652*rad)
state.activateFix(angleDNA)

#Lets test all 3 dihedral potentials
diheDNA = FixDihedralGauss(state, 'diheDNA')
dihePeri = FixDihedralPeriodic(state, 'diheProt')
diheCharm = FixDihedralCHARMM(state, 'diheCharm')
diheBase = FixBasePair3SPN2(state, 'basepair')
#I trust the Gaussian for NVT, but NVE is a little wonky. Energy is conserved, but sometimes the simulation breaks
diheDNA.createDihedral(state.atoms[0],state.atoms[1],state.atoms[2],state.atoms[3], -35.80*rad, 0.3, 7.0/kcal)

#I do not trust the Periodic dihedrals
#dihePeri.createDihedral(state.atoms[0],state.atoms[1],state.atoms[2],state.atoms[3], [10.0/kcal,5.0/kcal], 155.393480*rad)
#diheCharm.createDihedral(state.atoms[0],state.atoms[1],state.atoms[2],state.atoms[3], 7.0/kcal,1, -155.393480*rad)


#I 100% trust the base pair dihedral potential, conserves energy and does not explode
diheBase.createBasePair(state.atoms[0],state.atoms[1],state.atoms[2],state.atoms[3], -38.18*rad,5.82,15/kcal,153.17*rad,133.51*rad)
diheBase.setParameters(2.000,12.000)
#state.activateFix(diheDNA)
state.activateFix(diheBase)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
#InitializeAtoms.initTemp(state, 'all', 1.0)
InitializeAtoms.initTemp(state, 'all', 300.0)
fixNVT = FixLangevin(state, 'temp', 'all', 300.0)

#500 is the damping coefficient which is ~1/gamma (gamma is 500 for lammps version of 3spn2)
fixNVT.setParameters(0,0.002)
state.activateFix(fixNVT)

#Other things that we may want to output to make sure the simulation is running correctly
tempData = state.dataManager.recordTemperature('all', 100)

#state.dataManager.recordEnergy(interval=10, mode='scalar', fixes=[angleDNA, dihePeri])

#Run the actual system
integVerlet = IntegratorVerlet(state)
integRelax  = IntegratorRelax(state)
writeconfig = WriteConfig(state, fn='test_out', writeEvery=100, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
#integRelax.run(100,1e-4)
integVerlet.run(10000)
print tempData.vals
exit()
