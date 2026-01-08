#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hello
"""

import numpy as np
import codecs

convert_au_to_cm = 219474.63
convert_au_to_ev = 27.211386246
convert_J_to_ev = 6.24150907e18
convert_J_to_au = 2.29371228e17
h = 6.626068e-34
c = 299_792_458
e = 1.609217733e-19

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

    D,P,M = [],[],[]

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
                        D.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])
                    else:
                        state_transition.append(lsplit[0])
                        transition_energy.append(float(lsplit[1]))
                        D.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])
                elif abs_vd:
                    P.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])
                elif cd_ed:
                    M.append([float(lsplit[-3]),float(lsplit[-2]),float(lsplit[-1])])

    transition_energy = np.array(transition_energy)
    return center_of_mass, state_transition, transition_energy, D, P, M


def oscillator_force(transition_energy, moments):
    """
    Computes the oscillator force from the transition energy and the moments
    """

    N_trans = len(moments)
    fosc = np.zeros((N_trans,4)) #Store X*X, Y*Y, Z*Z and sum

    for energy, mom, fosc_i in zip(transition_energy, moments, range(N_trans)):
        fosc_i_x = 2 * energy / convert_au_to_cm * mom[0]**2 /3
        fosc_i_y = 2 * energy / convert_au_to_cm * mom[1]**2 /3
        fosc_i_z = 2 * energy / convert_au_to_cm * mom[2]**2 /3

        fosc[fosc_i] = [fosc_i_x, fosc_i_y, fosc_i_z, fosc_i_x + fosc_i_y + fosc_i_z]

    return fosc

def compute_spectra(transition_energy, moments, lambda_min=300, lambda_max=900, n_points=500, gauss=0.03, plot=False):
    """
    Computes the spectra
    The gauss factor provided is the enlargement of the peaks using gaussian fitting.
    If a list is provided, it will apply each term to each transition
    Lambda are in nm
    """
    fosc = oscillator_force(transition_energy, moments)


    FWHM = 2*(2*np.log(2))**(1/2)
    rac_pi = (2*np.pi)**(1/2)
    if type(gauss) is list:
        gauss = np.array(gauss)
    gauss_corr = gauss / FWHM

    lambda_list = np.linspace(lambda_min,lambda_max,n_points,endpoint=True)
    lambda_energy = h * c /(lambda_list*1e-9) * convert_J_to_au
    distance = lambda_energy.reshape((len(lambda_energy),1)) -transition_energy/convert_au_to_cm


    if type(gauss) is np.ndarray:
        gauss_x = np.einsum("ij,j->ij",-distance**2 / 2, 1/gauss_corr**2)

    else:
        gauss_x = -distance**2 / 2 / gauss_corr**2


    spectra_x = np.einsum("j,ij->i",fosc[:,0]/gauss_corr,np.exp(gauss_x)) / rac_pi
    spectra_y = np.einsum("j,ij->i",fosc[:,1]/gauss_corr,np.exp(gauss_x)) / rac_pi
    spectra_z = np.einsum("j,ij->i",fosc[:,2]/gauss_corr,np.exp(gauss_x)) / rac_pi

    spectra_xy = spectra_x + spectra_y
    spectra_xz = spectra_x + spectra_z
    spectra_yz = spectra_y + spectra_z

    spectra_xyz = spectra_x + spectra_y + spectra_z




    if plot:
        import matplotlib.pyplot as plt
        plt.plot(lambda_list,spectra_x)
        plt.show()


if __name__ == "__main__":

    center_of_mass, state_transition, transition_energy, D, P, M = extract_absorption_spectra_orca("Demo/Fcenter_TDA_PBE_TZVP.out")
    # fosc = oscillator_force(transition_energy, D)
    compute_spectra(transition_energy, D, plot=True)
