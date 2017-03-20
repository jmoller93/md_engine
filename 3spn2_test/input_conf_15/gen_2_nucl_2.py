import numpy as np
import scipy as sp
import math
import pdb
import copy
import sys
import argparse

class Nucleosome:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
   # def write_psf(self,fnme):


class Atom:
    def __init__(self, type, index, pos):
        self.residue_info = []
        self.index = index
        self.type = copy.copy(type)
        self.pos = pos
    #def write_xyz(self,fnme):

class Bond:
    def __init__(self, index, atom1, atom2, length ):
        self.index = index
        self.atom1 = atom1
        self.atom2 = atom2
        self.length = length
    #def write_bonds(self,fnme):


class Angle:
    def __init__(self,atom1,atom2,atom3):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
    #def write_bends(self,fnme):

class Dihedral:
    def __init__(self,atom1,atom2,atom3,atom4):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
    #def write_tors(self,fnme)

class Residue:
    def __init__(self,seg_name, resID, res_name, charge, mass):
        self.seg_name = copy.copy(seg_name)
        self.resID = resID
        self.res_name = copy.copy(res_name)
        self.charge = charge
        self.mass = mass


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('outputdir', type=str, help = "Which directory are we putting this rubbish into?")
    parser.add_argument('inputdir', type=str, help = "Where is everything coming from?")
    parser.add_argument('nucl_num', default = 2, type = int, help = "How many nucleosomes would you liek?")

    args = parser.parse_args()

    p_out_dir = args.outputdir
    p_in_dir = args.inputdir
    p_num_nucl = arg.nucl_num

    #Now it is time to read the xyz and natv files to input the atom info into the main Nucleosome class

    #Maybe just want to make this easier on myself and not use a dictionary
    d = {}
    for i in p_num_nucl:
        d["nucleosome" + str(i)] = Nucleosome()


    #Read in the xyz file first
    f = open("%s/in00_conf.xyz" % p_in_dir)
    lines = f.readlines()
    f.close()

    total_num = int(lines[0])
    index = 0
    pos = np.zeros(3)
    dna_num = lines.find('ALA')
    dna_num = dna_num-2
    prot_num = total_num-dna_num


    for i in range (p_num_nucl):
        for line in lines:
            if len(line.split())>3:
                p = line.split()
                type = copy.copy(p[0])
                pos[0] = float(p[1])
                pos[1] = float(p[2])
                pos[2] = float(p[3])
                d["nucleosome" + str(i)].atoms.append(Atom(type, index, pos))
                index = index+1


if __name__=="__main__":
    main()


