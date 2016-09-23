import numpy as np
import scipy as sp
import math
import pdb
import copy

#########################################################
#                                                       #
#                                                       #
#         Function that generates in00_bond.dat         #
#                                                       #
#                                                       #
#########################################################

def write_bonds(fnme,fnme2):


    atom = np.zeros(3)
    dont = 0
    length = 0.0
    charge = 0.0
    nope = 0.0
    max = 0

    f = open(fnme,"r")
    g = open(fnme2, "w")

    lines = f.readlines()

    for i in range(2):
        for line in lines:
            if len(line.split())>3:
                p = line.split()
                atom[0] = int(p[0])
                atom[1] = int(p[1])
                atom[2] = int(p[2])
                dont = int(p[3])
                length = float(p[4])
                charge = float(p[5])
                nope = float(p[6])

                if( max < atom[2]):
                    max = atom[2]

                if (i == 0):
                    g.write("%i  %i  %i  %i  %.6f  %.6f  %.6f\n" % (atom[0], atom[1], atom[2], dont, length, charge, nope))

                #make sure that the two nucleosomes are separated
                else:
                    g.write("%i  %i  %i  %i  %.6f  %.6f  %.6f\n" % (atom[0]+num, atom[1]+max+1, atom[2]+max+1, dont, length, charge, nope))
            else:
                if (i == 0):
                    num = int(line.split()[0])
                    num2 = num*2
                    g.write("%i\n" % (num2))
                else:
                    print("Clearly your data parsing skills need work")
    f.close()
    g.close()

#########################################################
#                                                       #
#                                                       #
#         Function that generates in00_bend.dat         #
#                                                       #
#                                                       #
#########################################################


def write_bend(fnme, fnme2):

    atom = np.zeros(4)
    angle = np.zeros(2)
    dont = 0
    max = 0

    f = open(fnme,"r")
    g = open(fnme2, "w")

    lines = f.readlines()

    for i in range(2):
        for line in lines:
            if len(line.split())>3:
                p = line.split()
                atom[0] = int(p[0])
                atom[1] = int(p[1])
                atom[2] = int(p[2])
                atom[3] = int(p[3])
                dont = int(p[4])
                angle[0] = float(p[5])
                angle[1] = float(p[6])

                if( max < atom[3]):
                    max = atom[3]

                if (i == 0):
                    g.write("%i  %i  %i  %i  %i  %.6f  %.6f\n" % (atom[0], atom[1], atom[2], atom[3], dont, angle[0], angle[1]))

                #make sure that the two nucleosomes are separated
                else:
                    g.write("%i  %i  %i  %i  %i  %.6f  %.6f\n" % (atom[0]+num, atom[1]+max+1, atom[2]+max+1, atom[3], dont, angle[0], angle[1]))

            else:
                if (i == 0):
                    num = int(line.split()[0])
                    num2 = num*2
                    g.write("%i\n" % (num2))
                else:
                    print("Clearly your data parsing skills need work")
    f.close()
    g.close()

#########################################################
#                                                       #
#                                                       #
#         Function that generates in00_tors.dat         #
#                                                       #
#                                                       #
#########################################################


def write_tors(fnme, fnme2):

    atom = np.zeros(6)
    data = np.zeros(3)
    max = 0

    f = open(fnme,"r")
    g = open(fnme2, "w")

    lines = f.readlines()

    for i in range(2):
        for line in lines:
            if len(line.split())>3:
                p = line.split()
                atom[0] = int(p[0])
                atom[1] = int(p[1])
                atom[2] = int(p[2])
                atom[3] = int(p[3])
                atom[4] = int(p[4])
                atom[5] = int(p[5])
                data[0] = float(p[6])
                data[1] = float(p[7])
                data[2] = float(p[8])

                if( max < atom[4]):
                    max = atom[4]

                if (i == 0):
                    g.write("%i  %i  %i  %i  %i  %i  %.6f  %.6f  %.6f\n" % (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], data[0], data[1], data[2]))

                else:
                    g.write("%i  %i  %i  %i  %i  %i  %.6f  %.6f  %.6f\n" % (atom[0]+num, atom[1]+max+1, atom[2]+max+1, atom[3]+max+1, atom[4]+max+1, atom[5], data[0], data[1], data[2]))
            else:
                if (i == 0):
                    num = int(line.split()[0])
                    num2 = num*2
                    g.write("%i\n" % (num2))
                else:
                    print("Clearly your data parsing skills need work")
    f.close()
    g.close()


#########################################################
#                                                       #
#                                                       #
#         Function that generates in00_natv.dat         #
#                                                       #
#                                                       #
#########################################################


def write_natv(fnme, fnme2):

    atom = np.zeros(2)
    data = np.zeros(2)
    max = 0

    f = open(fnme,"r")
    g = open(fnme2, "w")

    lines = f.readlines()

    for i in range(2):
        for line in lines:
            if len(line.split())>3:
                p = line.split()
                atom[0] = int(p[0])
                atom[1] = int(p[1])
                data[0] = float(p[2])
                data[1] = float(p[3])

                if( max < atom[1]):
                    max = atom[1]

                if (i == 0):
                    g.write("%i  %i  %.6f  %.6f\n" % (atom[0], atom[1], data[0], data[1]))

                else:
                    g.write("%i  %i  %.6f  %.6f\n" % (atom[0]+max+1, atom[1]+max+1, data[0], data[1]))
            else:
                if (i == 0):
                    num = int(line.split()[0])
                    num2 = num*2
                    g.write("%i\n" % (num2))
                else:
                    print("Clearly your data parsing skills need work")
    f.close()
    g.close()

#########################################################
#                                                       #
#                                                       #
#         Function that generates in00_conf.xyz         #
#                                                       #
#                                                       #
#########################################################


def write_conf(fnme, fnme2):

    cords = np.zeros(3)
    max = 0
    atom = ""

    f = open(fnme,"r")
    g = open(fnme2, "w")

    lines = f.readlines()

    num = int(lines[0])
    text = copy.copy(lines[1])

    g.write("%i\n%s\n" % (num*2,text))

    for i in range(2):
        for line in lines:
            if len(line.split())>3:
                p = line.split()
                atom = copy.copy(p[0])
                cords[0] = float(p[1])
                cords[1] = float(p[2])
                cords[2] = float(p[3])

                if (i == 0):
                    g.write("%s     %.6f      %.6f      %.6f\n" % (atom, cords[0], cords[1], cords[2]))

                else:
                    g.write("%s     %.6f      %.6f      %.6f\n" % (atom,  cords[0], cords[1], cords[2]+30.0))

    f.close()
    g.close()


def main():

    bond  = "in00_bond.dat"
    bond2 = "in00_bond_2.dat"
    bend = "in00_bend.dat"
    bend2 = "in00_bend_2.dat"
    tors = "in00_tors.dat"
    tors2 = "in00_tors_2.dat"
    natv = "in00_natv.dat"
    natv2 = "in00_natv_2.dat"
    conf = "in00_conf.xyz"
    conf2 = "in00_conf_2.xyz"

    write_bonds(bond, bond2)
    write_bend(bend, bend2)
    write_tors(tors, tors2)
    write_natv(natv, natv2)
    write_conf(conf, conf2)

    print("Your dirty work is done, signore")

if __name__=="__main__":
    main()
