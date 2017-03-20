import sys
import matplotlib.pyplot as plt
import math
import numpy as np

#Generate the forcefield parameters for 3SPN.2 in three seperate files
fnme = "in00_bstck_params.xml"
gnme = "in00_cstck_params.xml"
hnme = "in00_bspar_params.xml"

f = open(fnme,'w')
f.write("Base1\tBase2\tTheta\tEpsilon\tSigma\n")
f.write("A\tA\t100.13\t13.810\t3.58\n")
f.write("A\tT\t90.48\t15.050\t3.56\n")
f.write("A\tG\t104.39\t13.320\t3.85\n")
f.write("A\tC\t93.23\t15.820\t3.45\n")
f.write("T\tA\t102.59\t9.250\t4.15\n")
f.write("T\tT\t93.32\t12.44\t3.92\n")
f.write("T\tG\t103.70\t9.58\t4.32\n")
f.write("T\tC\t94.55\t13.11\t3.87\n")
f.write("G\tA\t95.45\t13.76\t3.51\n")
f.write("G\tT\t87.63\t14.59\t3.47\n")
f.write("G\tG\t106.36\t14.77\t3.67\n")
f.write("G\tC\t83.12\t15.17\t3.42\n")
f.write("C\tA\t102.69\t9.25\t4.15\n")
f.write("C\tT\t96.05\t12.42\t3.99\n")
f.write("C\tG\t100.46\t8.83\t4.34\n")
f.write("C\tC\t100.68\t14.01\t3.84\n")
f.close()

g = open(gnme, 'w')
g.write("Cross-stack Interactions 1\n")
g.write("A\tA\t
