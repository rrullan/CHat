#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script contains the functions used to compute and plot spectral data.
The script was written by Rémi Grincourt (M2 Student) based on a script of
Prof. Tangui Le Bahers
"""

import numpy as np
import codecs
import sys

convert_au_to_cm = 219474.63
convert_au_to_ev = 27.211386246
convert_J_to_ev = 6.24150907e18
convert_J_to_au = 2.29371228e17
convert_ev_to_cm = 8065.54
h = 6.626068e-34
c = 299792458
e = 1.609217733e-19


def extract_absorption_spectra_cp2k(file):
    """
    extract_absorption_pepctra_cp2k(file)

    Extract the Absorption spectra from a CP2K output
    """


    flag_states = 0
    flag_dipole = 0
    state_transition = []
    transition_energy = []
    abs_velo = []


    for line in codecs.open(file, 'r',encoding="utf-8"):
        if "Transition dipoles" in line:
            flag_states = 4

        elif flag_states:
            flag_states -= 1
            if not flag_states:
                flag_dipole = 1

        elif flag_dipole:
            if len(line) <= 4:
                break
            lsplit = line.split()
            state_transition.append(int(lsplit[1]))
            transition_energy.append(float(lsplit[2]) * convert_ev_to_cm)
            abs_velo.append([float(lsplit[3]),float(lsplit[4]),float(lsplit[5])])

    transition_energy = np.array(transition_energy)

    return state_transition, transition_energy, abs_velo



def extract_absorption_spectra_orca(file):
    """
    extract_absorption_spectra_orca(file)

    Extract the Absorption spectra from an ORCA output
    """

    flag_states = 0
    flag_completing_state = 0
    flag_completed_state = 0

    flag_spec = 0
    flag_spec_in = 0

    abs_elec,abs_velo,cd_elec, cd_velo = [],[],[],[]

    for line in codecs.open(file, 'r',encoding="utf-8"):
        if "Program Version" in line:
            version = int(line.split()[2][0])

        elif "Center of mass = " in line:
            lsplit = line.split()
            num_1,num_2,num_3 = lsplit[-3:]
            num_1, num_2, num_3 = num_1.replace("(",""), num_2.replace("(",""), num_3.replace("(","")
            num_1, num_2, num_3 = num_1.replace(")",""), num_2.replace(")",""), num_3.replace(")","")
            num_1, num_2, num_3 = num_1.replace(",",""), num_2.replace(",",""), num_3.replace(",","")
            num_1,num_2,num_3 = float(num_1),float(num_2),float(num_3)

            center_of_mass = [num_1,num_2,num_3]

        elif "SPECTRUM VIA TRANSITION" in line:
            flag_spec = 4
            if "ABSORPTION" in line and "ELECTRIC DIPOLE" in line:
                abs_ed,abs_vd,cd_ed,cd_vd = 1,0,0,0
                state_transition = []
                transition_energy = []

            elif "ABSORPTION" in line and "VELOCITY DIPOLE" in line:
                abs_ed,abs_vd,cd_ed,cd_vd = 0,1,0,0

            elif "CD" in line and "ELECTRIC DIPOLE" in line:
                abs_ed,abs_vd,cd_ed,cd_vd = 0,0,1,0
            elif "CD" in line and "VELOCITY DIPOLE" in line:
                abs_ed,abs_vd,cd_ed,cd_vd = 0,0,0,1
            elif "CD" in line and version <= 5:
                abs_ed,abs_vd,cd_ed,cd_vd = 0,0,1,0


        elif flag_spec:
            flag_spec -=1
            if flag_spec == 0: flag_spec_in = 1

        elif flag_spec_in:
            lsplit = line.split()
            if len(line)<3: flag_spec_in = 0
            else:
                if abs_ed:
                    if version >= 6:
                        state_transition.append(lsplit[0]+" "+lsplit[1]+" "+lsplit[2])
                        transition_energy.append(float(lsplit[3]))
                        abs_elec.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])
                    else:
                        state_transition.append(lsplit[0])
                        transition_energy.append(float(lsplit[1]))
                        abs_elec.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])
                elif abs_vd:
                    abs_velo.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])
                elif cd_ed:
                    cd_elec.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])
                elif cd_vd:
                    cd_velo.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])

    transition_energy = np.array(transition_energy)
    return center_of_mass, state_transition, transition_energy, abs_elec, abs_velo, cd_elec, cd_velo




def extract_absorption_spectra_gaussian(file):
    """
    extract_absorption_spectra_gaussian(file)

    Extract the Absorption spectra from a gaussian output
    """

    flag_edm = 0
    flag_vdm = 0
    flag_mdm = 0
    flag_excs = 0
    transition_energy = []
    state_transition = []
    for line in codecs.open(file, 'r',encoding="utf-8"):
        if "Ground to excited state" in line:
            flag_edm = 0
            flag_vdm = 0
            flag_mdm = 0

        if "electric dipole moment" in line:
            dipole = []
            flag_edm = 2

        elif "velocity dipole moment" in line:
            dipole = []
            flag_vdm = 2

        elif "magnetic dipole moment" in line:
            dipole = []
            flag_mdm = 2

        elif flag_edm:
            if flag_edm == 2:
                flag_edm = 1
                dm_edm = []
                state_transition = []
            else:
                lsplit = line.split()
                dm_edm.append([float(lsplit[1]),float(lsplit[2]),float(lsplit[3])])
                state_transition.append(int(lsplit[0]))

        elif flag_vdm:
            if flag_vdm == 2:
                flag_vdm = 1
                dm_vdm = []
            else:
                lsplit = line.split()
                dm_vdm.append([float(lsplit[1]),float(lsplit[2]),float(lsplit[3])])

        elif flag_mdm:
            if flag_mdm == 2:
                flag_mdm = 1
                dm_mdm = []
            else:
                lsplit = line.split()
                dm_mdm.append([float(lsplit[1]),float(lsplit[2]),float(lsplit[3])])

        if "Excited State  " in line: transition_energy.append(float(line.split()[4]))

    transition_energy = np.array(transition_energy[-len(state_transition):])
    state_transition = np.array(state_transition)
    dm_edm = np.array(dm_edm)
    dm_vdm = np.array(dm_vdm)
    dm_mdm = np.array(dm_mdm)

    return transition_energy, state_transition, dm_edm, dm_vdm, dm_mdm





def extract_absorption_spectra(file):
    """
    Extract data from an orca or a cp2k file
    """

    orca_file = 0
    gaussian_file = 0
    count = 0
    for line in open(file):
        count+=1
        if "Gaussian" in line:
            gaussian_file = 1
            break

        if count==3:
            if  "O   R   C   A" in line: orca_file = 1
            break

    if orca_file: center_of_mass, state_transition, transition_energy, abs_elec, abs_velo, cd_elec, cd_velo = extract_absorption_spectra_orca(file)

    elif gaussian_file:
        transition_energy, state_transition, abs_elec, abs_velo, abs_mag = extract_absorption_spectra_gaussian(file)
        center_of_mass = np.array([0,0,0])
        cd_elec = None
        cd_velo = None

    else:
        state_transition, transition_energy, abs_velo = extract_absorption_spectra_cp2k(file)
        center_of_mass = np.array([0,0,0])
        abs_elec = None
        cd_elec = None
        cd_velo = None

    return center_of_mass, state_transition, transition_energy, abs_elec, abs_velo, cd_elec, cd_velo




def oscillator_force(transition_energy, moments):
    """
    Computes the oscillator force from the transition energy and the moments
    """

    N_trans = len(moments)
    fosc = np.zeros((N_trans,4)) #Store X*X, Y*Y, Z*Z and sum

    for energy, mom, fosc_i in zip(transition_energy, moments, range(N_trans)):
        fosc_i_x = 2 * energy / convert_au_to_ev * mom[0]**2 /3
        fosc_i_y = 2 * energy / convert_au_to_ev * mom[1]**2 /3
        fosc_i_z = 2 * energy / convert_au_to_ev * mom[2]**2 /3

        fosc[fosc_i] = [fosc_i_x, fosc_i_y, fosc_i_z, fosc_i_x + fosc_i_y + fosc_i_z]

    return fosc



def compute_spectra(transition_energy, moments, DO=1, lambda_min=300, lambda_max=900, n_points=500, gauss=0.3, save=True, plot=False, show=False, file="",spectra="all",fmt="png"):
    """
    Computes the spectra
    The gauss factor provided is the enlargement of the peaks using gaussian fitting.
    If a list is provided, it will apply each term to each transition
    Lambda are in nm

    If plot is True and show is false, save in the same folder as file
                        show is true, show every spectra in one graph
    spectra can be : "all", "x", "xy" ... or a list of them
    """
    fosc = oscillator_force(transition_energy, moments)

    FWHM = 2*(2*np.log(2))**(1/2)
    rac_pi = (2*np.pi)**(1/2)
    if type(gauss) is list:
        gauss = np.array(gauss)
    gauss_corr = gauss / FWHM / convert_au_to_ev

    lambda_list = np.linspace(lambda_min,lambda_max,n_points,endpoint=True, dtype=int)
    lambda_energy = 1239.8 /(lambda_list) 
    distance = lambda_energy.reshape((len(lambda_energy),1)) -transition_energy


    if type(gauss) is np.ndarray:
        gauss_x = np.einsum("ij,j->ij",-distance**2 / 2, 1/gauss_corr**2)

    else:
        gauss_x = -distance**2 / 2 / gauss_corr**2

    spectra_x = np.einsum("j,ij->i",fosc[:,0]/gauss_corr,np.exp(gauss_x)) / rac_pi /3
    spectra_y = np.einsum("j,ij->i",fosc[:,1]/gauss_corr,np.exp(gauss_x)) / rac_pi /3
    spectra_z = np.einsum("j,ij->i",fosc[:,2]/gauss_corr,np.exp(gauss_x)) / rac_pi /3

    spectra_xy = spectra_x + spectra_y
    spectra_xz = spectra_x + spectra_z
    spectra_yz = spectra_y + spectra_z

    spectra_xyz = spectra_x + spectra_y + spectra_z


    norm = np.max(spectra_xyz)
    spectra_x = spectra_x / norm*DO
    spectra_y = spectra_y / norm*DO
    spectra_z = spectra_z / norm*DO
    spectra_xy = spectra_xy / norm*DO
    spectra_xz = spectra_xz / norm*DO
    spectra_yz = spectra_yz / norm*DO
    spectra_xyz = spectra_xyz / norm*DO


    folder = file.split("/")[:-1]
    dirr = "".join(folder)
    if dirr == "":
        dirr = "."
    file_name = (file.split("/")[-1]).split(".")[0]



    #x y z xy xz yz xyz
    Choice_spectra = np.zeros(7)

    if spectra == "all":
        Choice_spectra = np.ones(7)
        print_spectra = "x, y, z, xy, xz, yz, xyz"
    elif type(spectra) is str:
        if spectra == "x":
            Choice_spectra[0] = 1
        if spectra == "y":
            Choice_spectra[1] = 1
        if spectra == "z":
            Choice_spectra[2] = 1
        if spectra == "xy":
            Choice_spectra[3] = 1
        if spectra == "xz":
            Choice_spectra[4] = 1
        if spectra == "yz":
            Choice_spectra[5] = 1
        if spectra == "xyz":
            Choice_spectra[6] = 1
        print_spectra = spectra
    else:
        print_spectra = ""
        if "x" in spectra:
            Choice_spectra[0] = 1
            print_spectra += "x, "
        if "y" in spectra:
            Choice_spectra[1] = 1
            print_spectra += "y, "
        if "z" in spectra:
            Choice_spectra[2] = 1
            print_spectra += "z, "
        if "xy" in spectra:
            Choice_spectra[3] = 1
            print_spectra += "xy, "
        if "xz" in spectra:
            Choice_spectra[4] = 1
            print_spectra += "xz, "
        if "yz" in spectra:
            Choice_spectra[5] = 1
            print_spectra += "yz, "
        if "xyz" in spectra:
            Choice_spectra[6] = 1
            print_spectra += "xyz, "
        print_spectra = print_spectra[:-2]
    if save:
        spectra_all = np.array([spectra_x,spectra_y,spectra_z,spectra_xy,spectra_xz,spectra_yz,spectra_xyz])
        array_to_save = np.append(lambda_list.reshape((1,len(lambda_list))),spectra_all[Choice_spectra.astype("bool")],axis=0)
        # print(np.shape(array_to_save))
        np.savetxt(dirr+ "/" + file_name + "_spectra_data.txt",array_to_save.transpose(),header="Wavelength (nm) Absorption value for {}".format(print_spectra),fmt="%3.4f")

    if plot:
        import matplotlib.pyplot as plt
        if not show:
            fig,ax = plt.subplots(figsize=(8,6),dpi=200)
            plt.rcParams.update({'font.size': 15})
            plt.rcParams['svg.fonttype'] = 'none'
        colors = ["red","blue","green","darkred","cyan","lime","magenta","teal","purple","darkorange"]
        if Choice_spectra[0]:
            plt.plot(lambda_list,spectra_x,'-',label="X",color="cyan")
        if Choice_spectra[1]:
            plt.plot(lambda_list,spectra_y,'-',label="Y",color="magenta")
        if Choice_spectra[2]:
            plt.plot(lambda_list,spectra_z,'-',label="Z",color="yellow")
        if Choice_spectra[3]:
            plt.plot(lambda_list,spectra_xy,'-',label="XY",color="blue")
        if Choice_spectra[4]:
            plt.plot(lambda_list,spectra_xz,'-',label="XZ",color="green")
        if Choice_spectra[5]:
            plt.plot(lambda_list,spectra_yz,'-',label="YZ",color="red")
        if Choice_spectra[6]:
            plt.plot(lambda_list,spectra_xyz,'-',label="XYZ",color="black")
        plt.legend()
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Absorption")

        if show:
            plt.show()
        else:
            plt.savefig(dirr + "/" + file_name + "_spectra." + fmt)

    return lambda_list, spectra_x, spectra_y, spectra_z, spectra_xy, spectra_xz, spectra_yz, spectra_xyz




if __name__ == "__main__":

    # center_of_mass, state_transition, transition_energy, abs_elec, P, M, L = extract_absorption_spectra_orca("Demo/Fcenter_TDA_PBE_TZVP.out")
    # fosc = oscillator_force(transition_energy, D)
    # compute_spectra(transition_energy, abs_elec, plot=True,spectra="all",show=True)
    # compute_spectra(transition_energy, D, plot=True,)

    file = sys.argv[1]
    type_spectra = sys.argv[2]
    chosen_spectra = sys.argv[3]
    OD = float(sys.argv[4])
    gauss = float(sys.argv[5])
    min_lambda = float(sys.argv[6])
    max_lambda = float(sys.argv[7])
    num_points = int(sys.argv[8])
    fmt = sys.argv[9]

    # file = "Demo/Fcenter_TDA_PBE_TZVP.out"
    # type_spectra="abs_elec"
    # chosen_spectra="xyz"
    # OD = 1
    # gauss = 0.3
    # min_lambda=300
    # max_lambda=900
    # num_points= 601
    # fmt = "png"

    center_of_mass, state_transition, transition_energy, abs_elec, abs_velo, cd_elec, cd_velo = extract_absorption_spectra(file)

    if type_spectra == "abs_elec": moments = abs_elec
    if type_spectra == "abs_velo": moments = abs_velo
    if type_spectra == "cd_elec": moments = cd_elec
    if type_spectra == "cd_velo": moments = cd_velo

    if chosen_spectra != "all": chosen_spectra = chosen_spectra.split(",")

    compute_spectra(transition_energy, moments, DO=OD, lambda_min=min_lambda, lambda_max=max_lambda, n_points=num_points, gauss=gauss, plot=True, show=False, file=file,spectra=chosen_spectra,fmt=fmt)

