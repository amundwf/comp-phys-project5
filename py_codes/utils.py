# This file contains useful functions to be used in other .py files.

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import pandas as pd
import os
import re

def energy_EV_analytical_2x2(T, J, kB):
    # Returns the analytical expectation value for the 2x2 case,
    # as a function of temperature and the interaction term J (and
    # the Boltzmann constant kB).
    beta = 1/(kB*T)
    x = 8*J*beta
    print("x: "); print(x)
    print(np.sinh(x)); print(np.cosh(x))
    # The analytical expression:
    return -8*J*(np.sinh(x) / (np.cosh(x) + 3))

'''
def read_matrix_from_CSV(filename, directory):
    directory = "../results/3g_three_body/"

    # Read the planet names from the planet_names.csv file:
    
    #planetNamesData = np.loadtxt(directory + "planet_names.csv", skiprows=0, delimiter=",")
    #planetNamesData = pd.DataFrame(planetNamesData, columns=["planet_names"])
    #planetNames = planetNamesData["planet_names"]
    
    # Read in the file as a dataframe:
    dataFrame = pd.read_csv(directory + filename, header=None)
'''

def plot_spinMatrix(spinMatrix):
    # input: spinMatrix: A numpy array containing values -1 and 1.
    # This function plots a colormap of black (-1) and white (+1) squares.

    tempColorMatrix = spinMatrix
    tempColorMatrix[spinMatrix==-1] = 0
    plt.imshow(tempColorMatrix, cmap = 'gray') # The 'gray' colormap is binary:
    # black for values below 1 and white for values >= 1.
    plt.show()

def plot_4c_ising():
    
    directory = "../results/4c_ising/"
    filename = "2x2_T=1.0_N=100000.csv"
    
    data = np.loadtxt(directory +filename, skiprows=1, delimiter=",")
    data = pd.DataFrame(data, columns=["MC_cycles", "E", "E2", "M", "Mabs", "M2"])
    
    spinState = "ordered"
    T = 1.0
    L = 2.0
    L2 = L*L
    MC = data["MC_cycles"]
    N_MC = len(data['MC_cycles']) -1
    E = data["E"]
    E2 = data["E2"]
    M = data["M"]
    M2 = data["M2"]
    Mabs = data["Mabs"]
    
    f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
    
    plt.subplot(211)
    plt.plot(MC, E/L2, '.')
    plt.xscale("log")
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Energy")
    
    plt.subplot(212)
    plt.plot(MC, Mabs/L2, '.')
    plt.xscale("log")
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Magnetisation")
    
    f.tight_layout(pad=1.0)
    plt.savefig(directory + "{}_cycles_T={}_{}.pdf".format(N_MC, T, spinState))
    
    E_mean = np.mean(E)
    M_mean = np.mean(Mabs)
    E2_mean = np.mean(E2)
    M2_mean = np.mean(M2)
    CV = (E2_mean - E_mean*E_mean)*(1/(T*T))/L2
    chi = (M2_mean - M_mean*M_mean)*(1/(T))/L2
    
    print("N", N_MC)
    print("E mean", E_mean/L2)
    print("M mean", M_mean/L2)
    print("CV", CV)
    print("chi", chi)

def plot_4d_histogram():
    directory = "../results/4d/results/"
 
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            print(filename)
            data = np.loadtxt(directory +filename, skiprows=1, delimiter=",")
            data = pd.DataFrame(data, columns=["MC_cycles", "E", "E2", "M", "M_abs",
                                               "M2", "flipsAccepted"])
            if re.search('random', filename):
                spinState = "random"
            if re.search('ordered', filename):
                spinState = "ordered"
                
            T = re.findall('(?<=T=)[0-9].[0-9]+',filename)
            T = T[0] # take the only value in list.
            L = 20.0
            L2 = L*L
            MC = data["MC_cycles"]
            N_MC = len(data['MC_cycles']) - 1
            E = data["E"]
            E2 = data["E2"]
            M = data["M"]
            M2 = data["M2"]
            Mabs = data["M_abs"]
            
            
            
            f = plt.figure(figsize=(13,13), dpi=80, facecolor='w', edgecolor='k')
            E = E[20000:]
            E2 = E2[20000:]
            E_mean = np.mean(E)
            E_var = np.var(E)
            E2_mean = np.mean(E2)
            sigma2 = E2_mean - E_mean*E_mean
            
            # divide by L2
            E_mean /= L2
            E /= L2
            sigma2 /= L2
            E_var /= L2
            
            weights = np.ones_like(E) / (len(E))
            n, bins, patches = plt.hist(E, bins = 20, weights = weights, label = "T = {}".format(T))
            plt.xlabel("Energy")
            plt.ylabel("P[E]")
            plt.legend(loc="best")
            
            f.suptitle("Probability Density of the Energy\n with mean = {:.3f} and variance = {:.6f}.\n $E[E^2] - E[E]^2 = \u03C3^2_E$ = {:.6f}".format(E_mean, E_var, sigma2))
            f.tight_layout(pad= 6.0)
            plt.savefig(directory + "{}_cycles_T={}_{}.pdf".format(N_MC, T, spinState))
            plt.savefig(directory + "{}_cycles_T={}_{}.png".format(N_MC, T, spinState))
            plt.close()
    return


def plot_4d():
    
    directory = "../results/4d/mean_results/"
    
    ### MAKE PLOTS
    fig, ax = plt.subplots(1, 2, figsize=(20, 10), dpi = 80)
    i = 0
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
    
            data = np.loadtxt(directory +filename, skiprows=1, delimiter=",")
            data = pd.DataFrame(data, columns=["MC_cycles", "E_mean", "E2_mean", "M_mean", "M_abs_mean",
                                               "M2_mean", "CV", "chi"])
            if re.search('random', filename):
                spinState = "random"
            if re.search('ordered', filename):
                spinState = "ordered"
                
            T = re.findall('(?<=T=)[0-9].[0-9]+',filename)
            T = T[0] # take the only value in list.
            T = float(T)
            L = 20.0
            L2 = L*L
            MC = data["MC_cycles"]
            N_MC = len(data['MC_cycles'])
            E_mean = data["E_mean"]

            
            ################ ENERGY 
            ax[i].plot(MC[20000:], E_mean[20000:]/L2, '.')
            #plt.xscale("log")
            ax[i].set_xlabel("Monte Carlo cycles")
            ax[i].set_ylabel("Expectation Energy")
            ax[i].ticklabel_format(axis = 'y', style = 'sci', useOffset = False)
            if T == 1.0:
                ax[i].set_ylim(-1.9972, -1.9965)
            else:
                ax[i].set_ylim(-1.244, -1.235)
                
            axins = ax[i].inset_axes([0.6, 0.6, 0.37, 0.37])
            axins.plot(MC, E_mean/L2, '.')
            axins.set_xscale("log")
            axins.set_yticklabels('')
            axins.set_xticklabels('')
            
            ax[i].set_title("Spin {} configuration at T = {}".format(spinState, T))
            i += 1
    fig.tight_layout(pad=2.0)
    plt.savefig(directory + "{}_cycles_T={:.1f}_Energy.pdf".format(N_MC, T))
    plt.savefig(directory + "{}_cycles_T={:.1f}_Energy.png".format(N_MC, T))
    plt.close()
     
    ### MAKE PLOTS
    fig, ax = plt.subplots(1, 2, figsize=(20, 10), dpi = 80)
    i = 0
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
    
            data = np.loadtxt(directory +filename, skiprows=1, delimiter=",")
            data = pd.DataFrame(data, columns=["MC_cycles", "E_mean", "E2_mean", "M_mean", "M_abs_mean",
                                               "M2_mean", "CV", "chi"])
            if re.search('random', filename):
                spinState = "random"
            if re.search('ordered', filename):
                spinState = "ordered"
                
            T = re.findall('(?<=T=)[0-9].[0-9]+',filename)
            T = T[0] # take the only value in list.
            T = float(T)
            L = 20.0
            L2 = L*L
            MC = data["MC_cycles"]
            N_MC = len(data['MC_cycles'])
            Mabs_mean = data["M_abs_mean"]
            ########## MAGNETISATION
        
            ax[i].plot(MC[20000:], Mabs_mean[20000:]/L2, '.')
            
            ax[i].set_xlabel("Monte Carlo cycles")
            ax[i].set_ylabel("Expectation Magnetisation")
            ax[i].ticklabel_format(axis = 'y', style = 'sci', useOffset = False)
            
            if T == 1.0:
                ax[i].set_ylim(0.99875, 0.99929)
            else:
                ax[i].set_ylim(0.456, 0.472)
                
            axins = ax[i].inset_axes([0.6, 0.01, 0.37, 0.32])
            axins.plot(MC, Mabs_mean/L2, '.')
            axins.set_xscale("log")
            axins.set_yticklabels('')
            axins.set_xticklabels('')
                

            ax[i].set_title("Spin {} configuration at T = {}".format(spinState, T))
            i +=1
    fig.tight_layout(pad=2.0)
            
    plt.savefig(directory + "{}_cycles_T={:.1f}_Magnetisation.pdf".format(N_MC, T))
    plt.savefig(directory + "{}_cycles_T={:.1f}_Magnetisation.png".format(N_MC, T))
    plt.close()

def plot_4d_flips():
    ''' Plot the ration of flips accepted at each monte carlo
    cycle.
    '''
    directory = "../results/4d/results/"
    f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
    
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            data = np.loadtxt(directory + filename, skiprows=1, delimiter=",")
            data = pd.DataFrame(data, columns=["MC_cycles", "E_mean", "E2_mean", "M_mean", "M_abs_mean",
                                       "M2_mean", "flipsAccepted"])
            
            if re.search('random', filename):
                spinState = "random"
            if re.search('ordered', filename):
                spinState = "ordered"
                
            T = re.findall('(?<=T=)[0-9].[0-9]+',filename)
            T = T[0] # take the only value in list.
            L = 20
            L2 = L*L
            MC = data["MC_cycles"]
            N_MC = len(data['MC_cycles'])
            flipsAccepted = data["flipsAccepted"]

            plt.plot(MC, flipsAccepted/L2, '.-', label="T = {}".format(T))
            plt.xscale("log")
            plt.xlabel("Monte Carlo cycles")
            plt.ylabel("Accepted spins flips / Total possible flips")
            plt.legend()
            #f.tight_layout(pad=1.0)
            f.suptitle("The number of accepted spins to be flipped divided by the total possible.\n The 20 x 20 spin matrix was initialised with {} orientations".format(spinState))
    plt.savefig(directory + "accepted_flips_{}.pdf".format(spinState))
    plt.savefig(directory + "accepted_flips_{}.png".format(spinState))
    plt.close()
    return

def plot_4f_python():
        
    directory = "C:/github/Ising/"
    
    f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
    for filename in os.listdir(directory):
            if filename.endswith(".csv"):
    
                data = np.loadtxt(directory + filename, skiprows=1, delimiter=",")
                data = pd.DataFrame(data, columns=["Empty", "T", "E", "M", "CV", "Susc"])
                
                L = re.findall('(?<=N=)[0-9]+',filename)
                L = L[0]
                #L2 = L*L
                T = data["T"]
                E = data["E"]
                M = data["M"]
                CV = data["CV"]
                chi = data["Susc"]
                
                plt.subplot(221)
                plt.plot(T, E, '.-', label="{}".format(L))
                plt.xlabel("Temperature")
                plt.ylabel("Expectation Energy")
                plt.legend()
                
                plt.subplot(222)
                plt.plot(T, M, '.-', label="{}".format(L))
                plt.xlabel("Temperature")
                plt.ylabel("Expectation Magnetisation")
                plt.legend()
                
                plt.subplot(223)
                plt.plot(T, CV/float(L) , '.-', label="{}".format(L))
                plt.xlabel("Temperature")
                plt.ylabel("Specific Heat")
                plt.legend()
                
                plt.subplot(224)
                plt.plot(T, chi/float(L) , '.-', label="{}".format(L))
                plt.xlabel("Temperature")
                plt.ylabel("Magnetic Susceptibility")
                plt.legend()
                f.tight_layout(pad=1.0)
                
    plt.savefig(directory + "python_2.0_2.5_10.pdf")
    
def plot_4f_cpp():
        
    directory = "../results/4f/"
    f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
    i=0
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            print(filename)
            data = np.loadtxt(directory + filename, skiprows=1, delimiter=",")
            data = pd.DataFrame(data, columns=["T", "E_mean", "E2_mean", "M_mean", "M_abs_mean",
                                               "M2_mean", "CV", "chi"])
            
            L = re.findall('(?<=L=)[0-9]+',filename)
            L = L[0]
            T = data["T"]
            E_mean = data["E_mean"]
            M_mean = data["M_abs_mean"]
            CV = data["CV"]
            chi = data["chi"]

            
            plt.subplot(221)
            plt.plot(T, E_mean, '.-', label="{}".format(L))
            plt.xlabel("Temperature")
            plt.ylabel("Expectation Energy")
            plt.legend(loc="best")
            plt.xlim(2.0, 2.4)
            
            plt.subplot(222)
            plt.plot(T, M_mean, '.-', label="{}".format(L))
            plt.xlabel("Temperature")
            plt.ylabel("Expectation Magnetisation")
            plt.legend(loc="best")
            plt.xlim(2.0, 2.4)
            
            plt.subplot(223)
            plt.plot(T, CV, '.-', label="{}".format(L))
            plt.xlabel("Temperature")
            plt.ylabel("Heat Capacity")
            plt.legend(loc="best")
            plt.xlim(2.0, 2.4)
            
            plt.subplot(224)
            plt.plot(T, chi, '.-', label="{}".format(L))
            plt.xlabel("Temperature")
            plt.ylabel("Magnetic Susceptibility")
            plt.legend(loc="best")
            plt.xlim(2.0, 2.4)
            
            #f.suptitle("Expectation values for 100,000 Monte Carlo cycles for a 2x2 spin matrix. An additional 10% were added for equilibrium")
            f.tight_layout(pad=1.0)
            i+= 1
    plt.savefig(directory + "MC=10^5_start=2.0_End=2.3_Steps=16.png")
    plt.savefig(directory + "MC=10^5_start=2.0_End=2.3_Steps=16.pdf")
    return

def calc_4g():
    
    directory = "../results/4f/"
    
    TC_list = np.zeros((4))
    L_list = np.zeros((4))
    i = 0
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            data = np.loadtxt(directory + filename, skiprows=1, delimiter=",")
            data = pd.DataFrame(data, columns=["T", "E_mean", "E2_mean", "M_mean", "M_abs_mean",
                                               "M2_mean", "CV", "chi"])
            
            L = re.findall('(?<=L=)[0-9]+',filename)
            L = L[0]
            T = data["T"]
            E_mean = data["E_mean"]
            M_mean = data["M_abs_mean"]
            CV = data["CV"]
            chi = data["chi"]
            
            
            chiMax_ind = np.argmax(chi[1:])
            print(chiMax_ind)
            TC_list[i] = T[chiMax_ind]
            L_list[i] = L
            i+=1
            
    print("L list is :", L_list)
    print("TC list is : ", TC_list)
    
    ind = np.argsort(L_list)
    L_list = np.sort(L_list)
    TC_list = TC_list[ind]
    
    print("L list is :", L_list)
    print("TC list is : ", TC_list)
    
    x = 1/L_list
    y = TC_list

    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    print ("Gradient = {:.3f} and intercept = {:.3f}".format(m, c))
    
    f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(x, y, '.', label='$T_C  data$', markersize = 20)
    plt.plot(x, m*x + c, 'r', label='Fitted line')
    plt.legend()
    plt.xlabel("$L^{-1}$")
    plt.ylabel("$T_C$")
    plt.title("Gradient = {:.4f} and intercept = {:.4f}".format(m, c))
    
    
    plt.savefig(directory + "TC_in_limit.png")
    plt.savefig(directory + "TC_in_limit.pdf")
    
    return

def Z_analytical_2x2(J, kB, T):
    # The analytical expression of Z for a 2x2 spin lattice.
    beta = 1/float(kB*T); x = float(8*J*beta)
    Z = 4*np.cosh(x) + 12
    return Z
    
def E_mean_analytical_2x2(J, kB, T):
    beta = 1/float(kB*T)
    x = float(8*J*beta)
    Z = float(Z_analytical_2x2(J, kB, T))
    #E_mean = -8*J*np.sinh(x)/(np.cosh(x)+3)
    E_mean = -(1/Z)*32*J*np.sinh(x)
    return E_mean

def E2_mean_analytical_2x2(J, kB, T):
    beta = 1/float(kB*T); x = float(8*J*beta); Z = Z_analytical_2x2(J, kB, T)
    E2_mean = (1/Z)*256*J*J*np.cosh(x)
    return E2_mean

def M_abs_mean_analytical_2x2(J, kB, T):
    beta = 1/float(kB*T); x = float(8*J*beta); Z = Z_analytical_2x2(J, kB, T)
    M_abs_mean = (1/Z)*8*(np.exp(x)+2)
    return M_abs_mean

def M2_mean_analytical_2x2(J, kB, T):
    beta = 1/float(kB*T); x = float(8*J*beta); Z = Z_analytical_2x2(J, kB, T)
    M2_mean = (1/Z)*32*(np.exp(x)+1)
    return M2_mean

def C_V_analytical_2x2(J, kB, T):
    beta = 1/float(kB*T); x = float(8*J*beta); Z = Z_analytical_2x2(J, kB, T)
    E2_mean = E2_mean_analytical_2x2(J, kB, T)
    E_mean = E_mean_analytical_2x2(J, kB, T)
    #C_V = (256*J*J/(kB*T*T*Z)) * (np.cosh(x) - (1/Z)*((np.sinh(x))**2))
    C_V = (1/(kB*T*T))*(E2_mean - E_mean**2)
    return C_V

def chi_analytical_2x2(J, kB, T):
    beta = 1/float(kB*T); x = float(8*J*beta); Z = Z_analytical_2x2(J, kB, T)
    M2_mean = M2_mean_analytical_2x2(J, kB, T)
    M_abs_mean = M_abs_mean_analytical_2x2(J, kB, T)
    chi = beta*(M2_mean - M_abs_mean**2)
    return chi
