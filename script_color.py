#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script contains the function used to compute the color using an output file
"""

import numpy as np
import codecs
import sys
from script_spectra import *

#Dicts for eye and luminants
eye = {
    "2deg" : "data_color/CIE_1931_2deg.txt",
    "10deg" : "data_color/CIE_1964_10deg.txt"
}


lum = {
    "A" : "A",
    "B" : "data_color/lum_C.txt",
    "D65" : "data_color/lum_D65.txt"
}


def read_txt(file):
    return np.loadtxt(file,unpack=True)


def illuminant_A(lambda_list):
    """computes the CIE illuminant A (https://en.wikipedia.org/wiki/Standard_illuminant)"""
    return 100 * (560/lambda_list)**5 * (np.exp(1.435e7/2848/560)-1) / (np.exp(1.435e7/2848/lambda_list)-1)

def _illuminant_D_deprecated(T):
    """
    Returns the x and y of the illuminant D for a temperature
    The temperature must be between 4000K and 25000K
    https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D

    This function is not actually used as, to my knowledge, the spectra cannot be determined from this
    """

    if T >= 4000 and T <= 7000:
        x = 0.244063 + 0.09911e3 / T + 2.9678e6 /T/T - 4.6070e9 /T/T/T
    elif T >= 7000 and T <= 25000:
        x = 0.237040 + 0.24748e3 /T + 1.9018e6 /T/T - 2.0064e9 /T/T/T
    else: raise ValueError("Temperature not between 4000k and 25000K")
    y = -3*x*x + 2.87*x - 0.275
    return x,y


def illuminant_D65():
    """
    Returns the spectrum of illuminant D65
    """

    lambda_list, spectrum = np.loadtxt("data_color/lum_D65.txt",unpack=True)
    return lambda_list, spectrum

def illuminant_C():
    """
    Returns the spectrum of illuminant C
    """

    lambda_list, spectrum = np.loadtxt("data_color/lum_C.txt",unpack=True)
    return lambda_list, spectrum


def convert_spectrum_to_XYZ(spectrum,lamp,obs,norm=False):
    """
    Convert a spectrum to the XYZ color space
    """
    K = 1/np.sum(lamp * obs[2])
    X = K * np.sum(lamp * obs[1] * spectrum)
    Y = K * np.sum(lamp * obs[2] * spectrum)
    Z = K * np.sum(lamp * obs[3] * spectrum)
    if norm:
        X = X * 100/Y
        Z = Z * 100/Y
        Y = 100
    return X,Y,Z


def cieluv_perf_diff(obs,lamp):
    """
    Computes the luv coefficient for a perfect reflecting diffuser
    This is used to compute the cieluv in other cases
    See cieluv for references
    """

    X,Y,Z = convert_spectrum_to_XYZ(lamp,1,obs,norm=True)

    Yn = 100
    L = 116 * Y/Yn ** (1/3) - 16
    if L < 8: L = (29/3)**3 * Y/Yn


    upn = 4*X / (X+15*Y+3*Z)
    vpn = 9*Y / (X+15*Y+3*Z)
    return upn,vpn



def cieluv(X,Y,Z,obs,lamp):
    """
    Compute in the CIELUV color space
    https://en.wikipedia.org/wiki/CIELUV
    https://doi.org/10.1002/col.22873
    """
    Yn = 100
    L = 116 * (Y/Yn) ** (1/3) - 16
    if L < 8: L = (29/3)**3 * Y/Yn

    up = 4*X / (X+15*Y+3*Z)
    vp = 9*Y / (X+15*Y+3*Z)

    upn,vpn = cieluv_perf_diff(obs,lamp)

    us = 13*L * (up - upn)
    vs = 13*L * (vp - vpn)

    C = (us**2 + vs**2)**(1/2)
    h = np.arctan2(vs,us)

    return L,us,vs,C,h


def cielab(X,Y,Z,obs,lamp):
    """
    compute in the CIELAB color space
    https://en.wikipedia.org/wiki/CIELAB_color_space
    """

    def f(t):
        d = 6/29

        if t > d**3:
            return t**(1/3)
        return 1/3 * t /d**2  + 4/29

    Xn,Yn,Zn = convert_spectrum_to_XYZ(lamp,1,obs,norm=True)

    L = 116 * f(Y/Yn) - 16
    a = 500 * (f(X/Xn) - f(Y/Yn))
    b = 200 * (f(Y/Yn) - f(Z/Zn))
    C = (a**2 + b**2)**(1/2)
    h = np.arctan(b/a)

    return L,a,b,C,h

def srgb(X,Y,Z):
    """
    Compute in the CIE 1999 sRGB
    https://en.wikipedia.org/wiki/SRGB
    """

    def transfer_function(color):
        if color <= 0.0031308: return color*12.92
        return 1.055*(color)**(1/2.4) - 0.055

    X = X/100
    Y = Y/100
    Z = Z/100

    R =  3.2406255 * X - 1.5372080 * Y - 0.4986286 * Z
    G = -0.9689307 * X + 1.8758561 * Y + 0.0415175 * Z
    B =  0.0557101 * X - 0.2040211 * Y + 1.0569959 * Z

    R = transfer_function(R)
    G = transfer_function(G)
    B = transfer_function(B)

    if R>1: R=1
    if R<0: R=0
    if G>1: G=1
    if G<0: G=0
    if B>1: B=1
    if B<0: B=0

    R = int(round(R*255))
    G = int(round(G*255))
    B = int(round(B*255))
    return R,G,B


def compute_all_colors(file,obs,lamp,spectrum,print_data=True,save_data=True,save_color=True):
    """
    Computes and prints the color of a spectrum in all available spaces
    """

    X,Y,Z = convert_spectrum_to_XYZ(spectrum,lamp,obs)

    x = X / (X+Y+Z)
    y = Y / (X+Y+Z)
    z = 1 - x - y

    L_uv,us,vs,C_uv,h_uv = cieluv(X,Y,Z,obs,lamp)
    L_ab,a,b,C_ab,h_ab = cielab(X,Y,Z,obs,lamp)
    R,G,B = srgb(X,Y,Z)

    if print_data:
        print("########################################################")
        print('# Colors in the spaces :                               #')
        print('# XYZ : X = {:<10.4f} Y = {:^10.4f} Z = {:>10.4f}   #'.format(X,Y,Z))
        print('# xyz : x = {:<10.4f} y = {:^10.4f} z = {:>10.4f}   #'.format(x,y,z))
        print('# Luv : L = {:<10.4f} u = {:^10.4f} v = {:>10.4f}   #'.format(L_ab,us,vs))
        print('# Lab : L = {:<10.4f} a = {:^10.4f} b = {:>10.4f}   #'.format(L_ab,a,b))
        print('# RGB : R = {:<10.4f} G = {:^10.4f} B = {:>10.4f}   #'.format(R,G,B))
        print("########################################################")

    folder = file.split("/")[:-1]
    dirr = "".join(folder)
    if dirr == "":
        dirr = "."
    file_name = (file.split("/")[-1]).split(".")[0]

    if save_data:
        with open(dirr + "/" + file_name + "_color.txt","w") as f:

            f.write('Colors in the spaces :\n')
            f.write('XYZ : X = {:<10.4f} Y = {:^10.4f} Z = {:>10.4f}\n'.format(X,Y,Z))
            f.write('xyz : x = {:<10.4f} y = {:^10.4f} z = {:>10.4f}\n'.format(x,y,z))
            f.write('Luv : L = {:<10.4f} u = {:^10.4f} v = {:>10.4f}\n'.format(L_ab,us,vs))
            f.write('Lab : L = {:<10.4f} a = {:^10.4f} b = {:>10.4f}\n'.format(L_ab,a,b))
            f.write('RGB : R = {:<10.4f} G = {:^10.4f} B = {:>10.4f}\n'.format(R,G,B))

    if save_color:
        import matplotlib.image as mpimg
        image = np.zeros((500,500,3)) + np.array([R/255,G/255,B/255])
        mpimg.imsave(dirr + "/" + file_name +  "_color.png",image)



def compute_all_colors_from_file(file,moment="abs_elec",spectrum="xyz",eye_dict="10deg",lamp_dict="D65",DO=1,gauss=0.3):
    """
    Compute all the colors from an input file

    moments :
        abs_elec, absorption via transition electric dipole
        abs_velo, absorption via transition velocity dipole
    spectrum :
        x, y, z, xy, xz, yz, xyz

    eye_dict :
        2deg : the CIE 1931 2 degree observer standard
        10deg : the CIE 1964 10 degree observer standard

    lamp_dict :
        "A" : the CIE A standard illuminant
        "C" : the CIE C standard illuminant
        "D65" : the CIE D65 standard illuminant
    """


    center_of_mass, state_transition, transition_energy, abs_elec, abs_velo, cd_elec, cd_velo = extract_absorption_spectra(file)

    if moment == "abs_elec":
        lambda_list, spectra_x, spectra_y, spectra_z, spectra_xy, spectra_xz, spectra_yz, spectra_xyz = compute_spectra(transition_energy, abs_elec, lambda_min=360, lambda_max=830, n_points=471, gauss=gauss, DO=OD)

    elif moment == "abs_velo":
        lambda_list, spectra_x, spectra_y, spectra_z, spectra_xy, spectra_xz, spectra_yz, spectra_xyz = compute_spectra(transition_energy, abs_velo, lambda_min=360, lambda_max=830, n_points=471, gauss=gauss, DO=OD)

    elif moment == "cd_elec":
        lambda_list, spectra_x, spectra_y, spectra_z, spectra_xy, spectra_xz, spectra_yz, spectra_xyz = compute_spectra(transition_energy, cd_elec, lambda_min=360, lambda_max=830, n_points=471, gauss=gauss, DO=OD)

    elif moment == "cd_velo":
        lambda_list, spectra_x, spectra_y, spectra_z, spectra_xy, spectra_xz, spectra_yz, spectra_xyz = compute_spectra(transition_energy, cd_velo, lambda_min=360, lambda_max=830, n_points=471, gauss=gauss, DO=OD)

    if eye_dict in eye:
        eye_chosen = read_txt(eye[eye_dict])
    else:
        eye_chosen = read_txt(eye)


    if lamp_dict == "A":
        lamp = illuminant_A(eye[0])
    elif lamp_dict == "C":
        lum_lambda_list, lamp = read_txt(lum[lamp_dict])
        lambda_list = lum_lambda_list[12:]
        lamp = lamp[12:]
        eye_chosen = eye_chosen[:,:-50]
        eye_chosen = eye_chosen[:,::5]
    elif lamp_dict == "D65":
        lum_lambda_list, lamp = read_txt(lum[lamp_dict])
        lamp = lamp[60:]
    else:
        lum_lambda_list, lamp = read_txt(lamp_dict)


    if spectrum == "x": spectrum = 100 * 10**(-spectra_x)
    elif spectrum == "y": spectrum = 100 * 10**(-spectra_y)
    elif spectrum == "z": spectrum = 100 * 10**(-spectra_z)
    elif spectrum == "xy": spectrum = 100 * 10**(-spectra_xy)
    elif spectrum == "xz": spectrum = 100 * 10**(-spectra_xz)
    elif spectrum == "yz": spectrum = 100 * 10**(-spectra_yz)
    elif spectrum == "xyz": spectrum = 100 * 10**(-spectra_xyz)

    compute_all_colors(file,eye_chosen,lamp,spectrum,print_data=True,save_data=True,save_color=True)




if __name__ == "__main__":
    # obs = read_txt("data_color/CIE_1964_10deg.txt")
    # obs = read_txt("tangui_eye.txt")
    # obs = read_txt("data_color/CIE_1931_2deg.txt")

    # lambda_list, spectrum = illuminant_D65()
    # lambda_list = lambda_list[80:-50]
    # lamp = spectrum[80:-50]
    # print(lambda_list)
    # lambda_list,lamp = read_txt("tangui_d65.txt")

    # lambda_list = np.linspace(360,830,471)
    # lamp = illuminant_A(lambda_list)
    # lambda_list, lamp = illuminant_C()
    # lambda_list = lambda_list[12:]
    # lamp = lamp[12:]
    # eye = eye[:,:-50]
    # eye = eye[:,::5]
    # upn,vpn = cieluv_perf_diff(eye,lamp)
    # center_of_mass, state_transition, transition_energy, D, P, M = extract_absorption_spectra_orca("Demo/Fcenter_TDA_PBE_TZVP.out")
    # lambda_list, spectra_x, spectra_y, spectra_z, spectra_xy, spectra_xz, spectra_yz, spectra_xyz = compute_spectra(transition_energy, D, lambda_min=360, lambda_max=830, n_points=471)
    # compute_all_colors(eye,lamp,spectra_xyz,print_data=True)

    # transition_energy, MX, MY, MZ, gauss = np.loadtxt("electric_dipole_CAMB3LYP.txt",unpack=True)
    # transition_energy = transition_energy * 8065.54
    # Moments = np.array([MX,MY,MZ]).transpose()
    # #
    # E,A,B = np.loadtxt("electric_dipole_CAMB3LYP_XYZ.txt",unpack=True)
    # #
    # lambda_list, spectra_x, spectra_y, spectra_z, spectra_xy, spectra_xz, spectra_yz, spectra_xyz = compute_spectra(transition_energy, Moments, lambda_min=380, lambda_max=780, n_points=401,gauss=0.3,save=False)
    # spectrum = 100 * 10**(-spectra_xyz)
    # compute_all_colors("electric_dipole_CAMB3LYP.txt",obs,lamp,spectrum,print_data=True)
    #
    # import matplotlib.pyplot as plt
    # plt.plot(lambda_list,spectra_xyz,"o-r",label="CHat")
    # plt.plot(E,B,"o-b",label="spectre.c")
    # plt.legend()
    # plt.show()


    # import matplotlib.pyplot as plt
    # plt.plot(lambda_list,spectrum,"o-r",label="CHat")
    # plt.plot(E,A,"o-b",label="spectre.c")
    # plt.legend()
    # plt.show()

    # compute_all_colors_from_file("Demo/Fcenter_TDA_PBE_TZVP.out",moment="abs_elec",spectrum="xyz",eye_dict="10deg",lamp_dict="D65")

    file = sys.argv[1]
    eye_chosen = sys.argv[2]
    lamp_chosen = sys.argv[3]
    type_spectra = sys.argv[4]
    chosen_spectra = sys.argv[5]
    OD = float(sys.argv[6])
    gauss = float(sys.argv[7])


    compute_all_colors_from_file(file,moment=type_spectra,spectrum=chosen_spectra,eye_dict=eye_chosen,lamp_dict=lamp_chosen,DO=OD,gauss=gauss)
