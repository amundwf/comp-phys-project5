# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 10:35:41 2020

@author: jamesbc
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import pandas as pd
import os
import re

directory = "../results/2D_diffusion/"

for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        fig, ax = plt.subplots(4, 5, figsize=(9,9), dpi = 80)
        ax = ax.flatten()
        cube = np.loadtxt(directory + filename)
        
        cube = cube.reshape(-1,100,100)
        
        print("Cube shape is: ",cube.shape)
        
       
        t_dim = cube.shape[0]
        
        step = int(t_dim / len(ax))
        print(step)
        i=0
        j=0
        
        while i < t_dim:
            pcm = ax[j].imshow(cube[i,:,:], cmap="viridis")
            ax[j].set_title("solution \nat time={:.2f}".format(i*0.01))
            #ax[j].set_colorbar()
            
            j+=1
            i+=step
        fig.tight_layout(pad = 2.0)
        fig.colorbar(pcm, ax = ax)
        #fig.set_title("Tpoints = {}".format(t_dim))
        
