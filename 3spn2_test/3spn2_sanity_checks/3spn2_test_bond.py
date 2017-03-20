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
state.bounds = Bounds(state, lo = Vector(-100, -100, -100), hi = Vector(100, 100, 100))
state.rCut = 18.0
state.padding = 0.5
state.periodicInterval = 100
state.shoutEvery = 1000
state.dt = 10.0
state.setSpecialNeighborCoefs(1, 1, 1)
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
        nonbond.setParameter('rCut', spcs[i], spcs[j], sig/0.89)
        nonbond.setParameter('sig', spcs[i], spcs[j], sig)
        nonbond.setParameter('eps', spcs[i], spcs[j], 1.0/kcal)

state.activateFix(nonbond)

#The Debye Huckel electrostatics
#electric = FixChargePairDH(state, 'debyeHuckel', 'all')
#electric.setParameters(300,0.150)
#state.activateFix(electric)

#Remember the base pair identity of the particle if it is a base pair for
#later interactions
siteId = []

#Create atoms
state.addAtom(handle='P', pos=Vector(0,0,0),q=-0.6)
state.addAtom(handle='S', pos=Vector(0,0,10))

#Read the bond information
bondDNA  = FixBondHarmonicExtend(state, 'bondDNA')
bondProt = FixBondHarmonic(state, 'bondProt')
#bondProt.createBond(state.atoms[int(bondInfo[1])], state.atoms[int(bondInfo[2])], 2.0*float(bondInfo[5])/kcal, float(bondInfo[4]))
bondDNA.createBond(state.atoms[0], state.atoms[1], 0.6/kcal, 4.8635)
bondProt.createBond(state.atoms[0], state.atoms[1], 0.6/kcal, 2.0)

#Activate both kinds of bonded information
#state.activateFix(bondDNA)
state.activateFix(bondProt)

#We use a Langevin thermostat (Bussi-Parrinello in original model)
#InitializeAtoms.initTemp(state, 'all', 100.0)
#fixNVT = FixLangevin(state, 'temp', 'all', 300.0)

#500 is the damping coefficient which is ~1/gamma (gamma is 500 for lammps version of 3spn2)
#fixNVT.setParameters(12345,0.02)
#state.activateFix(fixNVT)

#Other things that we may want to output to make sure the simulation is running correctly
tempData = state.dataManager.recordTemperature('all', 50)
eng = state.dataManager.recordEnergy('all',100)
#Run the actual system
integVerlet = IntegratorVerlet(state)
writeconfig = WriteConfig(state, fn='test_out', writeEvery=100, format='xyz', handle='writer')
state.activateWriteConfig(writeconfig)
integVerlet.run(10000)
#for vals in eng.vals:
#    print(vals)

#Analyze the output file
fnme = "test_out.xyz"
gnme = "hist.dat"
f = open(fnme, "r")
g = open(gnme, "w")
lines = f.readlines()
f.close()

for line in lines:
    if len(line.split()) > 3:
        q = line.split()
        if int(q[0]) == 0:
            r = [float(q[1]),float(q[2]),float(q[3])]
        else:
            p = [float(q[1]),float(q[2]),float(q[3])]
            g.write("%s\n" % math.sqrt((r[0]-p[0])**2+(r[1]-p[1])**2+(r[2]-p[2])**2))

g.close()
exit()
