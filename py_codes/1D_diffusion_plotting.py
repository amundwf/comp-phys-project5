# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 13:19:01 2020

@author: jamesbc
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
import pandas as pd
import os
import re

directory = "../results/1D_diffusion/"

fig, ax = plt.subplots(2,2, figsize=(9,9), dpi = 80)
ax = ax.flatten()
j=0
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        
        
        
        data = np.loadtxt(directory + filename, delimiter=",")
        data = pd.DataFrame(data)
        print(filename)
        print(data.shape)
       
        
        if re.search('Explicit', filename):
            print("plotting explicit")
            pcm = ax[j].pcolormesh(data)
            ax[j].set_xlabel("N Points")
            ax[j].set_ylabel("Time Points")
            ax[j].set_yscale("linear")
            ax[j].set_title("{}".format(filename))
            
            fig.colorbar(pcm, ax = ax[j])
            j+=1
            
        
        elif re.search("v_xt", filename):
            print("plotting analytic")
            N_x = data.shape[1]
            xList = np.linspace(0,1,N_x)
            u = np.zeros((data.shape))
            
            for i in range(0, data.shape[0]-1):
                u[i,:] = data.iloc[i,:] + xList
            
            pcm = ax[j].pcolormesh(u, vmin=0.0, vmax=1.0)
            ax[j].set_xlabel("N Points")
            ax[j].set_ylabel("Time Points")
            ax[j].set_title("Analytic {} + $x_i/L$".format(filename))
            fig.colorbar(pcm, ax = ax[j])
            j+=1
        
        else:
            pcm = ax[j].pcolormesh(data, vmin=0.0, vmax=1.0)
            ax[j].set_xlabel("N Points")
            ax[j].set_ylabel("Time Points")
            ax[j].set_title("{}".format(filename))
            fig.colorbar(pcm, ax = ax[j])
            j+=1
        