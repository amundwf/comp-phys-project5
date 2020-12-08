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
            name = "Explicit"
        elif re.search("Implicit", filename):
            name = "Implicit"
        elif re.search("Crank", filename):
            name = "Crank Nicolson"
                 
        pcm = ax[j].pcolormesh(data, vmin=0.0, vmax=1.0)
        ax[j].set_xlabel("N Points")
        ax[j].set_ylabel("Time Points")
        ax[j].set_title("{}".format(name))
        fig.colorbar(pcm, ax = ax[j])
        j+=1

'''
fig, ax = plt.subplots(2,2, figsize=(9,9), dpi = 80)
ax = ax.flatten()
j=0
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        
        data = np.loadtxt(directory + filename, delimiter=",")
        data = pd.DataFrame(data)
        print(filename)
        print(data.shape)
        
        ax[j].plot(data.iloc[5000,:])
        ax[j].set_xlabel("N Points")
        ax[j].set_ylabel("$u(x,t)$")
        ax[j].set_title("{}".format(filename))
        
        j+=1
        
'''