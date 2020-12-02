# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 10:35:41 2020

@author: jamesbc
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import pandas as pd
import os
import re

directory = "../results/2D_diffusion/"

fig, ax = plt.subplots(2, 2, figsize=(20, 10), dpi = 80)
    
for filename in os.listdir(directory):
    if filename.endswith(".txt"):
    
        cube = np.loadtxt(directory + filename)
        
        cube = cube.reshape(40,40,-1)
        
        t_dim = cube.shape[2]
        print(t_dim)
        step = int(t_dim/4)
        i = 0
        j=0
        ax = ax.flatten()
        while i < t_dim:
            ax[j].imshow(cube[:,:,i])
            ax[j].set_title("solution at time={}".format(i))
            #ax[j].set_colorbar()
            
            j+=1
            i+=step

