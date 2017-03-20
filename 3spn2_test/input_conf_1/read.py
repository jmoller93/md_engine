import numpy as np
import scipy as sp
import math
import pdb

f = open("in00_conf.xyz", mode="r")
g = open("in00_cvmd.psf", mode='r')

lines = f.readlines()
lines2 = g.readlines()

f.close()
g.close()

x1=0
y1=0
z1=0
i=0

for line in lines:
    if len(line.split())>3:
        p = line.split()
        x1=x1+float(p[1])
        y1=y1+float(p[2])
        z1=z1+float(p[3])
        i=i+1

x1=x1/i; y1=y1/i; z1=z1/i;

print x1, y1, z1
