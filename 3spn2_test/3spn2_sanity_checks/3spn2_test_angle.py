#!/usr/bin/python

import sys,os,glob
#import matplotlib.pyplot as plt
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
#Note!!! Timestep has to be reduced because it seems to crash at a dt of 20 with this code

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
    #fnme.write("Atom 1\n")
    #fnme.write("%10.5f\t\t" % state.atoms[0].force[0])
    #fnme.write("%10.5f\t\t" % state.atoms[0].force[1])
    #fnme.write("%10.5f\n" % state.atoms[0].force[2])
    #fnme.write("Atom 2\n")
    #fnme.write("%10.5f\t\t" % state.atoms[1].force[0])
    #fnme.write("%10.5f\t\t" % state.atoms[1].force[1])
    #fnme.write("%10.5f\n" % state.atoms[1].force[2])
    #fnme.write("Atom 3\n")
    #fnme.write("%10.5f\t\t" % state.atoms[2].force[0])
    #fnme.write("%10.5f\t\t" % state.atoms[2].force[1])
    #fnme.write("%10.5f\n" % state.atoms[2].force[2])
    sum_force = 0.0
    for i in range(0,3):
        for j in range(0,3):
            sum_force = sum_force +state.atoms[j].force[i]
    fnme.write("%10.5f\n" % sum_force)
    fnme.close()

#Activate the python function as a "Fix"

pyOp = PythonOperation(operateEvery=1, handle='force_test', operation=force)
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
nonbond = FixLJRepul(state, 'excluded')

for i in range(len(spcs)):
    for j in range(len(spcs)):
        sig = (sigma[i] + sigma[j])*0.5
        #nonbond.setParameter('rCut', spcs[i], spcs[j], sig)
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)
        nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)

state.activateFix(nonbond)

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
state.addAtom(handle='P', pos=Vector(0,0,10),q=-0.6)


#Read the bond information
bondDNA  = FixBondHarmonicExtend(state, 'bondDNA')
bondProt = FixBondHarmonic(state, 'bondProt')
#bondProt.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], 2.0*float(bondInfo[5])/kcal, float(bondInfo[4]))
bondDNA.createBond(state.atoms[0], state.atoms[1], 0.6/kcal, 4.8635)
bondDNA.createBond(state.atoms[1], state.atoms[2], 0.6/kcal, 4.8635)

#Activate both kinds of bonded information
state.activateFix(bondDNA)
#state.activateFix(bondProt)

#Both angle potentials work, its just that the base stack is weak and outside the cone, does not bias!
angleDNA = FixAngleHarmonic(state, 'angleDNA')
angleBStack = FixAngleBaseStacking(state, 'intrastacking')
#angleDNA.createAngle(state.atoms[0],state.atoms[1],state.atoms[2],355.0/kcal,94.784734*rad)
angleBStack.createAngle(state.atoms[0],state.atoms[1],state.atoms[2],100.13*rad,13.8/kcal,3.58)
angleBStack.setParameters(3.000,6.000)
state.activateFix(angleBStack)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
#InitializeAtoms.initTemp(state, 'all', 300.0)
fixNVT = FixLangevin(state, 'temp', 'all', 300.0)

#500 is the damping coefficient which is ~1/gamma (gamma is 500 for lammps version of 3spn2)
fixNVT.setParameters(0,0.002)
#state.activateFix(fixNVT)

#Other things that we may want to output to make sure the simulation is running correctly
tempData = state.dataManager.recordTemperature('all', 50)

#Run the actual system
integVerlet = IntegratorVerlet(state)
writeconfig = WriteConfig(state, fn='test_out', writeEvery=10, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integVerlet.run(2000)

#Analyze the output file
fnme = "test_out.xyz"
gnme = "hist.dat"
f = open(fnme, "r")
g = open(gnme, "w")
lines = f.readlines()
f.close()

flag1 = 0
for line in lines:
    if len(line.split()) > 3:
        q = line.split()
        if int(q[0]) == 0 and flag1 == 0:
            r = np.array([float(q[1]),float(q[2]),float(q[3])])
            flag1 = 1;
        elif int(q[0]) == 1:
            p = np.array([float(q[1]),float(q[2]),float(q[3])])
        elif int(q[0]) == 0 and flag1 == 1:
            q = np.array([float(q[1]),float(q[2]),float(q[3])])
            a = np.linalg.norm(q-p)
            b = np.linalg.norm(r-q)
            cost = np.dot((q-p)/a,(p-r)/b)
            theta = (math.pi - math.acos(cost)) * 180.0/math.pi
            g.write("%s\n" % theta)
            flag1 = 0

g.close()
exit()
