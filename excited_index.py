#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Mar 26 11:36:12 2024

@author: adufou05
"""

#-----------------------------------------------------------------------------#
import sys
import numpy as np
from sys import exit, argv
import argparse, copy

__version__ = 2.0
bohr = 0.52917
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Classe pour les fichiers Cube Electron Density (Code 1 - CP2K)
class CubeCP2K():

    def __init__(self, fname=None):
        if fname != None:
            try:
                self.read_cube(fname)
            except IOError as e:
                print( "File used as input: %s" % fname )
                print( "File error ({0}): {1}".format(e.errno, e.strerror))
                self.terminate_code()
        else:
            self.default_values()
        return None

    def terminate_code(self):
        print("Code terminating now")
        exit()
        return None

    def default_values(self):
        self.natoms = 0
        self.comment1 = 0
        self.comment2 = 0
        self.origin = np.array([0, 0, 0])
        self.NX = 0
        self.NY = 0
        self.NZ = 0
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.atoms = ['0']
        self.atomsXYZ = [0, 0, 0]
        self.data = [0]
        return None

    def read_cube(self, fname):
        """
        Method to read cube file. Just needs the filename
        """
        with open(fname, 'r') as fin:
            self.filename = fname
            self.comment1 = fin.readline()  # Save 1st comment
            self.comment2 = fin.readline()  # Save 2nd comment
            nOrigin = fin.readline().split()  # Number of Atoms and Origin
            self.natoms = int(nOrigin[0])  # Number of Atoms
            self.origin = np.array([float(nOrigin[1]), float(nOrigin[2]), float(nOrigin[3])])  # Position of Origin
            nVoxel = fin.readline().split()  # Number of Voxels
            self.NX = int(nVoxel[0])
            self.X = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])
            nVoxel = fin.readline().split()  #
            self.NY = int(nVoxel[0])
            self.Y = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])
            nVoxel = fin.readline().split()  #
            self.NZ = int(nVoxel[0])
            self.Z = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])
            self.atoms = []
            self.atomsXYZ = []
            for atom in range(self.natoms):
                line = fin.readline().split()
                self.atoms.append(line[0])
                self.atomsXYZ.append(list(map(float, [line[2], line[3], line[4]])))
            self.data = np.zeros((self.NX, self.NY, self.NZ))
            i = int(0)
            for s in fin:
                for v in s.split():
                    self.data[int(i / (self.NY * self.NZ)), int((i / self.NZ) % self.NY), int(i % self.NZ)] = float(v)
                    i += 1
        return None

    def write_cube(self, fname, comment='Cube file written by CubeToolz_M\nCubeToolz_M %3.1f' % __version__):
        """
        Write out a Gaussian Cube file
        """
        try:
            with open(fname, 'w') as fout:
                if len(comment.split('\n')) != 2:
                    print('Comment line NEEDS to be two lines!')
                    self.terminate_code()
                fout.write('%s\n' % comment)
                fout.write("%5d  %0.6f  %.6f  %.6f\n" % (self.natoms, self.origin[0], self.origin[1], self.origin[2]))
                fout.write("%5d   %.6f   %.6f   %.6f\n" % (self.NX, self.X[0], self.X[1], self.X[2]))
                fout.write("%5d   %.6f   %.6f   %.6f\n" % (self.NY, self.Y[0], self.Y[1], self.Y[2]))
                fout.write("%5d   %.6f   %.6f   %.6f\n" % (self.NZ, self.Z[0], self.Z[1], self.Z[2]))
                for atom, xyz in zip(self.atoms, self.atomsXYZ):
                    fout.write("%5s   %.6f   %.6f   %.6f %.6f\n" % (atom, int(atom), xyz[0], xyz[1], xyz[2]))
                for ix in range(self.NX):
                    for iy in range(self.NY):
                        for iz in range(self.NZ):
                            fout.write("%.5e " % self.data[ix, iy, iz]),
                            if (iz % 6 == 5): fout.write('\n')
                    fout.write("\n")
        except IOError as e:
            print("File used as output does not work: %s" % fname)
            print("File error ({0}): {1}".format(e.errno, e.strerror))
            self.terminate_code()
        return None

    def square_cube(self, power=2):
        """
        Function to raise cube data to a power. Squares cube data by default.
        """
        self.data = self.data ** power
        return None

    def cube_int(self):
        """
        Integrate the entire cube data.
        """
        vol = np.linalg.det(np.array([self.X, self.Y, self.Z]))
        edensity = np.sum(self.data)
        nelectron = vol * edensity
        return nelectron

#-----------------------------------------------------------------------------#
# Classe pour les fichiers OM Cube (Code 2 - Orca)
class CubeOrca():
    def __init__(self, fname=None):
        if fname != None:
            try:
                self.read_cube_OM(fname)
            except IOError as e:
                print("File used as input: %s" % fname)
                print("File error ({0}): {1}".format(e.errno, e.strerror))
                self.terminate_code_OM()
        else:
            self.default_values_OM()
        return None

    def terminate_code_OM(self):
        print("Code terminating now")
        exit()
        return None

    def default_values_OM(self):
        self.natoms = 0
        self.comment1 = 0
        self.comment2 = 0
        self.origin = np.array([0, 0, 0])
        self.NX = 0
        self.NY = 0
        self.NZ = 0
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.atoms = ['0']
        self.atomsXYZ = [0, 0, 0]
        self.data = [0]
        return None

    def read_cube_OM(self,fname):
        """
        Method to read cube file. Just needs the filename
        """

        with open(fname, 'r') as fin:
            self.filename = fname
            self.comment1 = fin.readline() #Save 1st comment
            self.comment2 = fin.readline() #Save 2nd comment
            nOrigin = fin.readline().split() # Number of Atoms and Origin
            self.natoms = - int(nOrigin[0]) #Number of Atoms
            self.origin = np.array([float(nOrigin[1]),float(nOrigin[2]),float(nOrigin[3])]) #Position of Origin
            nVoxel = fin.readline().split() #Number of Voxels
            self.NX = int(nVoxel[0])
            self.X = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
            nVoxel = fin.readline().split() #
            self.NY = int(nVoxel[0])
            self.Y = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
            nVoxel = fin.readline().split() #
            self.NZ = int(nVoxel[0])
            self.Z = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
            self.atoms = []
            self.atomsXYZ = []
            for atom in range(self.natoms):
                line= fin.readline().split()
                self.atoms.append(line[0])
                # self.atomicmass.append(line[1])
                self.atomsXYZ.append(list(map(float,[line[2], line[3], line[4]])))
            self.data = np.zeros((self.NX,self.NY,self.NZ))
            no_orb = fin.readline().split()
            self.Orbital = int(no_orb[1])
            i= int(0)
            for s in fin:
                for v in s.split() :
                    if int(i/(self.NY*self.NZ)) < int(nVoxel[0]) :
                        self.data[int(i/(self.NY*self.NZ)), int((i/self.NZ)%self.NY), int(i%self.NZ)] = float(v)
                        i+=1
        return None      

    def read_cube(self, fname):
        """
        Method to read cube file. Just needs the filename
        """
        with open(fname, 'r') as fin:
            self.filename = fname
            self.comment1 = fin.readline()  # Save 1st comment
            self.comment2 = fin.readline()  # Save 2nd comment
            nOrigin = fin.readline().split()  # Number of Atoms and Origin
            self.natoms = int(nOrigin[0])  # Number of Atoms
            self.origin = np.array([float(nOrigin[1]), float(nOrigin[2]), float(nOrigin[3])])  # Position of Origin
            nVoxel = fin.readline().split()  # Number of Voxels
            self.NX = int(nVoxel[0])
            self.X = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])
            nVoxel = fin.readline().split()  #
            self.NY = int(nVoxel[0])
            self.Y = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])
            nVoxel = fin.readline().split()  #
            self.NZ = int(nVoxel[0])
            self.Z = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])
            self.atoms = []
            self.atomsXYZ = []
            for atom in range(self.natoms):
                line = fin.readline().split()
                self.atoms.append(line[0])
                self.atomsXYZ.append(list(map(float, [line[2], line[3], line[4]])))
            self.data = np.zeros((self.NX, self.NY, self.NZ))
            i = int(0)
            for s in fin:
                for v in s.split():
                    self.data[int(i / (self.NY * self.NZ)), int((i / self.NZ) % self.NY), int(i % self.NZ)] = float(v)
                    i += 1
        return None

    def write_cube_OM(self,fname,comment='Cube file written by CubeToolz_M\nCubeToolz_M %3.1f' % __version__):
        '''
        Write out a Gaussian Cube file
        '''
        try:
            with open(fname,'w') as fout:
                if len(comment.split('\n')) != 2:
                    print( 'Comment line NEEDS to be two lines!')
                    self.terminate_code()
                fout.write('%s\n' % comment)
                fout.write("%5d  %0.6f  %.6f  %.6f\n" % (- self.natoms, self.origin[0], self.origin[1], self.origin[2]))
                fout.write("%5d   %.6f   %.6f   %.6f\n" % (self.NX, self.X[0], self.X[1], self.X[2]))
                fout.write("%5d   %.6f   %.6f   %.6f\n" % (self.NY, self.Y[0], self.Y[1], self.Y[2]))
                fout.write("%5d   %.6f   %.6f   %.6f\n" % (self.NZ, self.Z[0], self.Z[1], self.Z[2]))
                for atom,xyz in zip(self.atoms,self.atomsXYZ):
                    fout.write("%5s   %0.6f   %6.3f   %6.3f %6.3f\n" % (atom, int(atom), xyz[0], xyz[1], xyz[2]))
                fout.write("%5s %4s\n" % (1, self.Orbital))
                for ix in range(self.NX):
                   for iy in range(self.NY):
                       for iz in range(self.NZ):
                           fout.write("%.5e " % self.data[ix,iy,iz]),
                           if (iz % 6 == 5): fout.write('\n')
                       fout.write("\n")
        except IOError as e:
            print( "File used as output does not work: %s" % fname)
            print( "File error ({0}): {1}".format(e.errno, e.strerror))
            self.terminate_code()
        return None

    def write_cube(self, fname, comment='Cube file written by CubeToolz\nCubeToolz %3.1f' % __version__):
        """
        Write out a Gaussian Cube file
        """
        try:
            with open(fname, 'w') as fout:
                if len(comment.split('\n')) != 2:
                    print('Comment line NEEDS to be two lines!')
                    self.terminate_code()
                fout.write('%s\n' % comment)
                fout.write("%4d %.6f %.6f %.6f\n" % (self.natoms, self.origin[0], self.origin[1], self.origin[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NX, self.X[0], self.X[1], self.X[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NY, self.Y[0], self.Y[1], self.Y[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NZ, self.Z[0], self.Z[1], self.Z[2]))
                for atom, xyz in zip(self.atoms, self.atomsXYZ):
                    fout.write("%s %d %6.3f %6.3f %6.3f\n" % (atom, 0, xyz[0], xyz[1], xyz[2]))
                for ix in range(self.NX):
                    for iy in range(self.NY):
                        for iz in range(self.NZ):
                            fout.write("%.5e " % self.data[ix, iy, iz]),
                            if (iz % 6 == 5): fout.write('\n')
                    fout.write("\n")
        except IOError as e:
            print("File used as output does not work: %s" % fname)
            print("File error ({0}): {1}".format(e.errno, e.strerror))
            self.terminate_code()
        return None

    def square_cube_OM(self, power=2):
        """
        Function to raise cube data to a power. Squares cube data by default.
        """
        self.data = self.data ** power
        return None

    def cube_OM_int(self):
        """
        Integrate the entire cube data.
        """
        vol = np.linalg.det(np.array([self.X, self.Y, self.Z]))
        edensity = np.sum(self.data)
        nelectron = vol * edensity
        return nelectron

#-----------------------------------------------------------------------------#
# Fonctions communes fusionnées : CP2K et Orca
def diff_cubes(files, name):
    cubes = [CubeCP2K(fin) for fin in files]
    print("====== Subtracting cube files ======")
    cube_out = copy.deepcopy(cubes[0])
    for ctmp in cubes[1:]:
        cube_out.data -= ctmp.data
    print("====== Writing output cube as diff.cube ======")
    cube_out.write_cube(name)
    return cube_out

def add_cubes(files, name):
    cubes = [CubeCP2K(fin) for fin in files]
    print("====== Adding cube files ======")
    cube_out = copy.deepcopy(cubes[0])
    for ctmp in cubes[1:]:
        cube_out.data += ctmp.data
    print("====== Writing output cube as diff.cube ======")
    cube_out.write_cube(name)
    return cube_out

def square_cubes(files, power):
    cubes = [CubeCP2K(fin) for fin in files]
    print("====== Squaring cube files ======")
    [ctmp.square_cube(power) for ctmp in cubes]
    print("====== Writing output cubes as squareN.cube ======")
    if len(cubes) == 1:
        cubes[0].write_cube('square.cube{}'.format(files[0]))
    else:
        for ind, cout in enumerate(cubes):
            cout.write_cube('square%d.cube' % ind)
    return cubes[0]

def positive_cubes(files, name):
    cubes = [CubeCP2K(fin) for fin in files]
    print("====== Adding cube files ======")
    cube_out = copy.deepcopy(cubes[0])
    for ctmp in cubes[1:]:
        cube_out.data += ctmp.data
        cube_out.data = cube_out.data * 0.5
    print("====== Writing output cube as diff.cube ======")
    cube_out.write_cube(name)
    return cube_out

def mult_cubes(files):
    cubes = [CubeCP2K(fin) for fin in files]
    print("====== Multiplying cube files ======")
    cube_out = copy.deepcopy(cubes[0])
    for ctmp in cubes[1:]:
        cube_out.data *= ctmp.data
    print("====== Writing output cube as mult.cube ======")
    cube_out.write_cube('mult.cube')
    return cube_out

def barycentre(file):
    cube_input = CubeCP2K(file)
    temp = cube_input.data
    d = 0
    b = [0, 0, 0]
    for i in range(len(temp[:, 0, 0])):
        for j in range(len(temp[0, :, 0])):
            for k in range(len(temp[0, 0, :])):
                d += temp[i, j, k]
                b[0] += temp[i, j, k] * i
                b[1] += temp[i, j, k] * j
                b[2] += temp[i, j, k] * k
    b[0] = (b[0] / d)
    b[1] = (b[1] / d)
    b[2] = (b[2] / d)
    print(f'Barycentre position: {b} Angstrom')
    return b

def rcubes(files):
    cubes = [CubeCP2K(fin) for fin in files]
    for i in range(len(cubes)):
        cube = cubes[i]
        temp = cube.data
        grid_size = [cube.NX, cube.NY, cube.NZ]
        X, Y, Z = np.indices(grid_size)
        center = np.array(cube.origin)
        distances = np.sqrt((X - center[0])**2 + (Y - center[1])**2 + (Z - center[2])**2)
        cube.data = distances
        cube.write_cube(f'radial_{i}.cube')

def cube_integrate(file):
    cube = CubeCP2K(file)
    return cube.cube_int()

def cube_integrate_OM(file):
    cube = CubeOrca(file)
    return cube.cube_OM_int()

def main_dct(input_file, initial_state, final_state, bohr):
    if initial_state != 0 and initial_state < 10 and final_state < 10:
        diff_cubes(["{}.cisdp0{}.cube".format(input_file, final_state), "{}.cisdp0{}.cube".format(input_file, initial_state)], 'diff.cube')
    elif initial_state != 0 and initial_state < 10 and final_state >= 10:
        diff_cubes(["{}.cisdp0{}.cube".format(input_file, final_state), "{}.cisdp{}.cube".format(input_file, initial_state)], 'diff.cube')
    elif initial_state != 0 and initial_state >= 10:
        diff_cubes(["{}.cisdp{}.cube".format(input_file, final_state), "{}.cisdp{}.cube".format(input_file, initial_state)], 'diff.cube') 
    square_cubes(["diff.cube"], 2)
    square_cubes(["square.cubediff.cube"], 0.5)
    positive_cubes(["diff.cube", "square.cubesquare.cubediff.cube"], 'positive.cube')
    diff_cubes(["diff.cube", "positive.cube"], 'negative.cube')
    bary_plus = barycentre('positive.cube')
    bary_moins = barycentre('negative.cube')
    Dct = distance(bary_plus, bary_moins) * bohr
    qct = cube_integrate("positive.cube")
    muct = Dct * qct
    print("")
    print("========================================")
    print("")
    print(" Le Bahers Index : {} Aangstroem".format(Dct))
    print("")
    print("========================================")
    return [Dct , qct , muct, bary_plus, bary_moins]
   
def main_double_overlap(input_file, initial_state, final_state, bohr):
    print("Work in Progress")
    if initial_state == 0:
        add_cubes(["ground_state.cube", "final_state.cube"], "Final_Excited_State.eldens.cube")
    else : 
        add_cubes(["ground_state.cube", "initial_state.cube"], "Initial_Excited_State.eldens.cube")
        add_cubes(["ground_state.cube", "final_state.cube"], "Finale_Excited_State.eldens.cube")
    mult_cubes(["Initial_State_Excited.eldens.cube", "Final_Excited_State.eldens.cube"])
    double_int = cube_integrate("mult.cube")
    print("========================================")
    print("")
    print(" Double Overlaping between State {} and {} : {}".format(initial_state, final_state, double_int))
    print("")
    print("========================================")

def main_double_overlap_cube(initial_state, final_state):
    mult_cubes([initial_state, final_state])
    double_int = cube_integrate("mult.cube")
    int_ini = cube_integrate(initial_state)
    int_fin = cube_integrate(final_state)
    P_if = (2*double_int)/(int_ini**2 + int_fin**2)
    print("========================================")
    print("")
    print(" Double Overlaping between the two State : ", P_if)
    print("")
    print("========================================")
    output = open("Double_Integral.txt", "a")
    output.write("Double Integral of Electronic Density : ")
    output.write(str(P_if))
    output.close()

#-----------------------------------------------------------------------------#
# Cube OM Class Function : Tozer Index and Omega Index

def mult_cubes_OM(files, name):
    cubes = [CubeOrca(fin) for fin in files]
    print( "====== Multiplying cube files ======")
    cube_out = copy.deepcopy(cubes[0])
    for ctmp in cubes[1:]:
        cube_out.data *= ctmp.data
    print( "====== Writing output cube as mult.cube ======")
    cube_out.write_cube_OM(name)
    return cube_out

def square_cubes_OM(files,power):
    cubes = [CubeOrca(fin) for fin in files]
    print( "====== Squaring cube files ======")
    [ctmp.square_cube_OM(power) for ctmp in cubes]
    print( "====== Writing output cubes as squareN.cube ======")
    if len(cubes) == 1:
        cubes[0].write_cube_OM('square.cube{}'.format(files[0]))
    else:
        for ind,cout in enumerate(cubes):
            cout.write_cube_OM('square%d.cube' % ind)
    return cubes[0]

def mult_cubes_CP2K(files, name):
    cubes = [CubeCP2K(fin) for fin in files]
    print( "====== Multiplying cube files ======")
    cube_out = copy.deepcopy(cubes[0])
    for ctmp in cubes[1:]:
        cube_out.data *= ctmp.data
    print( "====== Writing output cube as mult.cube ======")
    cube_out.write_cube(name)
    return cube_out

def square_cubes_CP2K(files,power):
    cubes = [CubeCP2K(fin) for fin in files]
    print( "====== Squaring cube files ======")
    [ctmp.square_cube(power) for ctmp in cubes]
    print( "====== Writing output cubes as squareN.cube ======")
    if len(cubes) == 1:
        cubes[0].write_cube('square.cube{}'.format(files[0]))
    else:
        for ind,cout in enumerate(cubes):
            cout.write_cube('square%d.cube' % ind)
    return cubes[0]

def cube_OM_integrate(file):
    cube_in = CubeOrca(file)
    cube_int_total = cube_in.cube_OM_int()
    return cube_int_total


def cube_CP2K_integrate(file):
    cube_in = CubeCP2K(file)
    cube_int_total = cube_in.cube_int()
    return cube_int_total

def rephase_cube_OM(file, signe, output_name):
    cube_in = CubeOrca(file)
    cube_out = copy.deepcopy(cube_in)
    cube_out.data *= signe
    cube_out.write_cube_OM(output_name)

def rephase_OM(orbitals, transi, signe, input_file):
    coeff_k = transi**2
    n = coeff_k.argmax()
    si = sign(transi[n])
    if si==signe :
        return transi
    else :
        for i in range(len(transi)):
            om = int(orbitals[i])
            rephase_cube_OM("{}.mo{}a.cube".format(input_file, om), -1, "{}.mo{}a.cube".format(input_file, om))
        return -transi

def main_tozer_CP2K(input_file, cutoff):
    MO_occ = molecular_decomposition_CP2K('occupied', False, "transition.txt")
    MO_virt = molecular_decomposition_CP2K('virtual', False, "transition.txt")
    cut = float(cutoff)
    overlapping = np.zeros(len(MO_occ))
    coeff_k = molecular_decomposition_CP2K("Partition", False, "transition.txt")**2
    coeff_k_test = molecular_decomposition_CP2K("Partition", False, "transition.txt")
    for transi_index in range(len(MO_occ)):
        if abs(coeff_k_test[transi_index])>cut :
            occ = int(MO_occ[transi_index])
            virt = int(MO_virt[transi_index])
            print("")
            print("-----------------")
            print("    Transition {} -> {}".format(occ, virt))
            print("")
            square_cubes_CP2K(["{}-WFN_{:05d}_1-1_0.cube".format(input_file, occ)], 2)
            square_cubes_CP2K(["square.cube{}-WFN_{:05d}_1-1_0.cube".format(input_file, occ)], 0.5)
            square_cubes_CP2K(["{}-WFN_{:05d}_1-1_0.cube".format(input_file, virt)], 2)
            square_cubes_CP2K(["square.cube{}-WFN_{:05d}_1-1_0.cube".format(input_file, virt)], 0.5)
            mult_cubes_CP2K(["square.cubesquare.cube{}-WFN_{:05d}_1-1_0.cube".format(input_file, occ), "square.cubesquare.cube{}-WFN_{:05d}_1-1_0.cube".format(input_file, virt)], "mult_{}.cube".format(transi_index))
            overlapping[transi_index] = cube_CP2K_integrate("mult_{}.cube".format(transi_index))
                 
            
    tozer_tmp = 0
    normalization = 0
    for i in range(len(coeff_k)):
        
        if abs(coeff_k_test[i]) > cut :
            tozer_tmp += coeff_k[i] * overlapping[i]
            normalization += coeff_k[i] 
    tozer = tozer_tmp/normalization
    print("")
    print("========================================")
    print("")
    print(" Tozer Index : ", tozer)
    print("")
    print("========================================")
    output_print_tozer_CP2K(tozer)
    

def main_tozer_Orca(input_file, cutoff):
    MO_occ = molecular_decomposition('occupied', False, "transition.txt")
    MO_virt = molecular_decomposition('virtual', False, "transition.txt")
    cut = float(cutoff)
    overlapping = np.zeros(len(MO_occ))
    coeff_k = molecular_decomposition("Partition", False, "transition.txt")
    for transi_index in range(len(MO_occ)):
        if coeff_k[transi_index]>cut :
            occ = int(MO_occ[transi_index])
            virt = int(MO_virt[transi_index])
            print("")
            print("-----------------")
            print("    Transition {} -> {}".format(occ, virt))
            print("")
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, occ)], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, occ)], 0.5)
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, virt)], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, virt)], 0.5)
            mult_cubes_OM(["square.cubesquare.cube{}.mo{}a.cube".format(input_file, occ), "square.cubesquare.cube{}.mo{}a.cube".format(input_file, virt)], "mult_{}.cube".format(transi_index))
            overlapping[transi_index] = cube_OM_integrate("mult_{}.cube".format(transi_index))
    tozer_tmp = 0
    normalization = 0
    for i in range(len(coeff_k)):
        if coeff_k[i] > cut :
            tozer_tmp += coeff_k[i] * overlapping[i]
            normalization += coeff_k[i]
    tozer = tozer_tmp/normalization
    print("")
    print("========================================")
    print("")
    print(" Tozer Index : ", tozer)
    print("")
    print("========================================")
    output_print_tozer(tozer)

def main_omega_soc(input_file, cutoff, initial_state, final_state):
    transi_ini = "transition_S{}.txt".format(initial_state)
    transi_fin = "transition_T{}.txt".format(final_state)
    MO_occ_s = molecular_decomposition('occupied', False, transi_ini)
    MO_virt_s = molecular_decomposition('virtual', False, transi_ini)
    MO_occ_t = molecular_decomposition('occupied', False, transi_fin)
    MO_virt_t = molecular_decomposition('virtual', False, transi_fin)
    coeff_k_s_raw = molecular_decomposition("Partition", False, transi_ini)
    coeff_k_t_raw = molecular_decomposition("Partition", False, transi_fin)
    coeff_k_s = coeff_k_s_raw**2
    indice_transi_maj_s = coeff_k_s.argmax()
    phase = sign(coeff_k_s_raw[indice_transi_maj_s])
    coeff_k_t_phase = rephase_OM(MO_occ_t, coeff_k_t_raw, phase, input_file)
    rephase_OM(MO_virt_t, coeff_k_t_raw, phase, input_file)
    Occ_temp = list(MO_occ_s) + list(MO_occ_t)
    Virt_temp = list(MO_virt_s) + list(MO_virt_t)
    temp_list = []
    Occ = []
    Virt = []
    for i in range(len(Occ_temp)):
        temp = str(Occ_temp[i])+'-'+str(Virt_temp[i])
        if temp not in temp_list:
            temp_list.append(temp)
    for transi in temp_list:
        temp = transi.split('-')
        Occ.append(float(temp[0]))
        Virt.append(float(temp[1]))
    MO_occ = np.array(Occ)
    MO_virt = np.array(Virt)
    cut = float(cutoff)
    overlapping = np.zeros(len(MO_occ))
    coeff_k = np.zeros(len(MO_occ))
    S_temp = []
    T_temp = []
    for i in range(len(MO_occ_s)):
        S_temp.append(str(MO_occ_s[i])+'-'+str(MO_virt_s[i]))
    for i in range(len(MO_occ_t)):
        T_temp.append(str(MO_occ_t[i])+'-'+str(MO_virt_t[i]))
    for i in range(len(temp_list)):
        for j in range(len(S_temp)):
            if temp_list[i]==S_temp[j]:
                coeff_k[i]+=coeff_k_s_raw[j]
        for j in range(len(T_temp)):
            if temp_list[i]==T_temp[j]:
                coeff_k[i]-=coeff_k_t_phase[j]
    coeff_k = np.array(coeff_k)**2
    for transi_index in range(len(MO_occ)):
        if coeff_k[transi_index]>cut :
            occ = int(MO_occ[transi_index])
            virt = int(MO_virt[transi_index])
            print("")
            print("-----------------")
            print("    Transition {} -> {}".format(occ, virt))
            print("")
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, str(occ).replace('.0', ""))], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, str(occ).replace('.0', ""))], 0.5)
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, str(virt).replace('.0', ""))], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, str(virt).replace('.0', ""))], 0.5)
            mult_cubes_OM(["square.cubesquare.cube{}.mo{}a.cube".format(input_file, str(occ).replace('.0', "")), "square.cubesquare.cube{}.mo{}a.cube".format(input_file, str(virt).replace('.0', ""))], "mult_{}.cube".format(transi_index))
            overlapping[transi_index] = cube_OM_integrate("mult_{}.cube".format(transi_index))

    overlapping_s = np.zeros(len(MO_occ_s))
    for transi_index in range(len(MO_occ_s)):
        if coeff_k_s[transi_index]>cut:
            occ = int(MO_occ_s[transi_index])
            virt = int(MO_virt_s[transi_index])
            print("")
            print("-----------------")
            print("    Transition {} -> {}".format(occ, virt))
            print("")
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, occ)], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, occ)], 0.5)
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, virt)], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, virt)], 0.5)
            mult_cubes_OM(["square.cubesquare.cube{}.mo{}a.cube".format(input_file, occ), "square.cubesquare.cube{}.mo{}a.cube".format(input_file, virt)], "mult.cube")
            overlapping_s[transi_index] = cube_OM_integrate("mult.cube")
    
    overlapping_t = np.zeros(len(MO_occ_t))
    coeff_k_t = coeff_k_t_raw**2
    for transi_index in range(len(MO_occ_t)):
        if coeff_k_t[transi_index]>cut:
            occ = int(MO_occ_t[transi_index])
            virt = int(MO_virt_t[transi_index])
            print("")
            print("-----------------")
            print("    Transition {} -> {}".format(occ, virt))
            print("")
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, occ)], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, occ)], 0.5)
            square_cubes_OM(["{}.mo{}a.cube".format(input_file, virt)], 2)
            square_cubes_OM(["square.cube{}.mo{}a.cube".format(input_file, virt)], 0.5)
            mult_cubes_OM(["square.cubesquare.cube{}.mo{}a.cube".format(input_file, occ), "square.cubesquare.cube{}.mo{}a.cube".format(input_file, virt)], "mult.cube")
            overlapping_t[transi_index] = cube_OM_integrate("mult.cube")

    omega_tmp = 0
    normalization = 0
    for i in range(len(coeff_k)):
        if coeff_k[i] > cut :
            omega_tmp += coeff_k[i] * overlapping[i]
            normalization += coeff_k[i]
    omega = omega_tmp/normalization

    sigma_a = 0
    for i in range(len(coeff_k_s)):
        sigma_a += coeff_k_s[i] * overlapping_s[i]
    sigma_b = 0
    for i in range(len(coeff_k_t)):
        sigma_b += coeff_k_t[i] * overlapping_t[i]
    sigma = (sigma_a + sigma_b)/normalization
    gamma = (sigma_a - sigma_b)/normalization
    
    print("")
    print("     ========================================")
    print("")
    print(" Omega Index for ISC between state S{} and state T{}: ".format(initial_state, final_state), omega)
    print("")
    print("     ========================================")
    omega_output = open("Omega.txt","a")
    omega_output.write("Omega Index between state S{} and T{}: ".format(initial_state, final_state))
    omega_output.write(str(omega))
    omega_output.write("\n")
    omega_output.write("Sigma : ")
    omega_output.write(str(sigma))
    omega_output.write("\n")
    omega_output.write("Gamma : ")
    omega_output.write(str(gamma))
    omega_output.write("\n")
    omega_output.write("-- Virtual Transition Decomposition -- ")
    omega_output.write("\n")
    for i in range(len(MO_occ)):
        omega_output.write(str(int(MO_occ[i])))
        omega_output.write(" -> ")
        omega_output.write(str(int(MO_virt[i])))
        omega_output.write(" : ")
        omega_output.write(str(coeff_k[i]))
        omega_output.write(" ( {} )".format(overlapping[i]))
        omega_output.write("\n")
    
#-----------------------------------------------------------------------------#
# General Function

def find_nb_state_calc(name_file): 
    input_file = open("{}".format(name_file),"r").readlines()
    for i in range(len(input_file)):
        line = input_file[i].split()
        input_file[i] = line   
    i = 1
    while 'nroots' not in input_file[i]:
        i+=1
    nb_states = int(input_file[i][3])
    j = 1
    while 'DOSOC' not in input_file[j] and j<(len(input_file)-1):
        j+=1 
    if j!=(len(input_file)-1):
        triplets = input_file[j]
    else:
        triplets = 'NOOO'
    if triplets[3]=='TRUE': 
        print(nb_states, True)
    else:
        print(nb_states, False)
    
def find_nb_state_calc_CP2K(name_file): 
    input_file = open("{}".format(name_file),"r").readlines()
    for i in range(len(input_file)):
        line = input_file[i].split()
        input_file[i] = line   
    i = 1
    while 'NSTATES' not in input_file[i]:
        i+=1
    nb_states = int(input_file[i][1])
    print(nb_states)
    
def molecular_decomposition(type_OM, plot, name_file):
    input_file = open("{}".format(name_file),"r").readlines()
    for i in range(len(input_file)):
        line = input_file[i].split()
        input_file[i] = line
    OM = np.zeros((len(input_file), 3))
    for i in range(len(input_file)):
        OM[i][0] = int(input_file[i][0].replace("a",""))
        OM[i][1] = int(input_file[i][2].replace("a",""))
        OM[i][2] = float(input_file[i][4])
    if type_OM == "occupied" and plot==True:
        print(OM[:,0])
    if type_OM == "occupied" and plot==False:
        return OM[:,0]
    if type_OM == 'virtual' and plot==True:
        print(OM[:,1])
    if type_OM == 'virtual' and plot==False:
        return OM[:,1]
    if type_OM == "Partition":
        return OM[:,2]

def molecular_decomposition_CP2K(type_OM, plot, name_file):
    input_file = open("{}".format(name_file),"r").readlines()
    for i in range(len(input_file)):
        line = input_file[i].split()
        input_file[i] = line
    OM = np.zeros((len(input_file), 3))
    for i in range(len(input_file)):
        OM[i][0] = input_file[i][0]
        OM[i][1] = input_file[i][1]
        OM[i][2] = input_file[i][2]
    if type_OM == "occupied" and plot==True:
        print(OM[:,0])
    if type_OM == "occupied" and plot==False:
        return OM[:,0]
    if type_OM == 'virtual' and plot==True:
        print(OM[:,1])
    if type_OM == 'virtual' and plot==False:
        return OM[:,1]
    if type_OM == "Partition":
        return OM[:,2]

    
def distance(point_A, point_B):
    pas_x = 1
    pas_y = 1
    pas_z = 1
    x = ((point_A[0] - point_B[0])*pas_x)**2
    y = ((point_A[1] - point_B[1])*pas_y)**2
    z = ((point_A[2] - point_B[2])*pas_z)**2
    return np.sqrt(x + y + z)    

def sign(num):
    return -1 if num < 0 else 1

def output_print_tozer(tozer_index):
    tozer_output = open("Tozer.txt","a")
    tozer_output.write("Tozer Index : ")
    tozer_output.write(str(tozer_index))
    tozer_output.write("\n")
    MO_occ = molecular_decomposition('occupied', False, "transition.txt")
    MO_virt = molecular_decomposition('virtual', False, "transition.txt")
    coeff_k = molecular_decomposition("Partition", False, "transition.txt") ** 2
    tozer_output.write("-- Transition Decomposition -- ")
    tozer_output.write("\n")
    for i in range(len(MO_occ)):
        tozer_output.write(str(int(MO_occ[i])))
        tozer_output.write(" -> ")
        tozer_output.write(str(int(MO_virt[i])))
        tozer_output.write(" : ")
        tozer_output.write(str(coeff_k[i]))
        tozer_output.write("\n")

def output_print_tozer_CP2K(tozer_index):
    tozer_output = open("Tozer.txt","a")
    tozer_output.write("Tozer Index : ")
    tozer_output.write(str(tozer_index))
    tozer_output.write("\n")
    MO_occ = molecular_decomposition_CP2K('occupied', False, "transition.txt")
    MO_virt = molecular_decomposition_CP2K('virtual', False, "transition.txt")
    coeff_k = molecular_decomposition_CP2K("Partition", False, "transition.txt") ** 2
    tozer_output.write("-- Transition Decomposition -- ")
    tozer_output.write("\n")
    for i in range(len(MO_occ)):
        tozer_output.write(str(int(MO_occ[i])))
        tozer_output.write(" -> ")
        tozer_output.write(str(int(MO_virt[i])))
        tozer_output.write(" : ")
        tozer_output.write(str(coeff_k[i]))
        tozer_output.write("\n")

def output_print_Dct(dct_list):
    dct_output = open("Dct.txt","a")
    dct_output.write("Le Bahers Index : ")
    dct_output.write(str(dct_list[0]))
    dct_output.write(" Aangstroem \n")
    dct_output.write("-----------")
    dct_output.write("\n")
    dct_output.write("   Transfered Charge : ")
    dct_output.write(str(dct_list[1]))
    dct_output.write(" |e-| \n")
    dct_output.write("   Dipole Moment : ")
    dct_output.write(str(dct_list[2]))
    dct_output.write(" Debye \n")
    dct_output.write("-----------")
    dct_output.write("\n")
    dct_output.write("   Positive Barycenter Position : ")
    dct_output.write(str(dct_list[3][0]))
    dct_output.write("  ")
    dct_output.write(str(dct_list[3][1]))
    dct_output.write("  ")
    dct_output.write(str(dct_list[3][2]))
    dct_output.write("\n")
    dct_output.write("   Negative Barycenter Position : ")
    dct_output.write(str(dct_list[4][0]))
    dct_output.write("  ")
    dct_output.write(str(dct_list[4][1]))
    dct_output.write("  ")
    dct_output.write(str(dct_list[4][2]))
    dct_output.write("\n")
    dct_output.write("-----------")
    dct_output.close()

#-----------------------------------------------------------------------------#

if __name__ == '__main__':
    index = sys.argv[1]
    if index not in ['tozer', 'dct', 'omega_soc', 'double_overlap', 'double_overlap_cube']:
        sys.exit()
    elif index=='tozer':
        name = sys.argv[2]
        cutoff = sys.argv[3]
        code = sys.argv[4]
        if code == 'CP2K':     
            main_tozer_CP2K(name, cutoff)
        else:
            main_tozer_Orca(name, cutoff)
            
    elif index=='dct':
        input_file = sys.argv[2]
        initial_state = int(sys.argv[3])
        final_state = int(sys.argv[4])
        dct_list = main_dct(input_file, initial_state, final_state, bohr)
        output_print_Dct(dct_list)
    elif index=='omega_soc':
        input_file = sys.argv[2]
        cutoff = sys.argv[3]
        initial_state = int(sys.argv[4])
        final_state = int(sys.argv[5])
        omega = main_omega_soc(input_file, cutoff, initial_state, final_state)
    elif index=="double_overlap":
        input_file = sys.argv[2]
        initial_state = int(sys.argv[3])
        final_state = int(sys.argv[4])
        main_double_overlap(input_file, initial_state, final_state, bohr)
    elif index=="double_overlap_cube":
        initial_state = sys.argv[2]
        final_state = sys.argv[3]
        main_double_overlap_cube(initial_state, final_state)

#-----------------------------------------------------------------------------#
