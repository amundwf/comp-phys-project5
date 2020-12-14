# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:34:52 2020

@author: jamesbc
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams.update({'font.size': 12})
import pandas as pd
import os
import re

#directory = "../results/2D_diffusion_before_enrichment/"
#directory = "../results/2D_diffusion_after_enrichment/"
   

for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        
        if re.search("Ana", filename):
            name = "Analytic"
        else:
            name = "Numerical"
            
            
        fig, axes = plt.subplots(3,2, figsize=(15, 15), dpi = 80)
        axes = axes.flatten()
        cube = np.loadtxt(directory + filename)
        
        cube = cube.reshape(-1,x_dim,x_dim)
        print("The cube shape is: ",cube.shape)
        
        t_dim = cube.shape[0]
        step = int(t_dim / len(axes))
        
        x = np.arange(x_dim)
        y = np.ones((x_dim))*80
        y1 = np.ones((x_dim))*100
        
        i=0
        for ax in axes:
            
            # For some reason the y axis needs flipping.
            pcm = ax.pcolormesh(cube[i,::-1,:], cmap="viridis")#, vmin = 0.0, vmax=1.0)
            ax.plot(x, y, 'w--')
            ax.plot(x, y1, 'w--')
            
            i+=step
            cbar = fig.colorbar(pcm, ax = ax, format="%.1e")
            cbar.set_label("Temperature /Kelvin")
            
        fig.suptitle("{} solution to two-dimensional diffusion equation".format(name))
        
        plt.savefig(directory + "{}.png".format(filename))
        plt.savefig(directory + "{}.pdf".format(filename))

