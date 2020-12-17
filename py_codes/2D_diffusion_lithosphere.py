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

directory = "../results/2D_diffusion_before_enrichment/"
#directory = "../results/2D_diffusion_after_enrichment/"
   
x_dim = 120
dt = 0.0001
dx = 1
enrich = "before"

for filename in os.listdir(directory):
    if filename.endswith(".txt"):
    
        fig, axes = plt.subplots(3,2, figsize=(9,9), dpi = 80)
        axes = axes.flatten()
        cube = np.loadtxt(directory + filename)
        
        cube = cube.reshape(-1,x_dim,x_dim)
        print("The cube shape is: ",cube.shape)
        
        t_dim = cube.shape[0]
        step = int(t_dim / len(axes))
        
        x = np.arange(x_dim)
        y = np.ones((x_dim))*80
        y1 = np.ones((x_dim))*100
        
        if re.search("Ana", filename):
            name = "Analytic"
        else:
            name = "Numerical"
        
        i=0
        for ax in axes:
            # For some reason the y axis needs flipping.
            pcm = ax.pcolormesh(cube[i,::-1,:], cmap="viridis")
            cbar = fig.colorbar(pcm, ax = ax, format="%.1e")
            cbar.set_label("Temperature / Kelvin")
            ax.set_title("t={:.2e}".format(i*dt))
            
            # Add white dashed lines for crust and mantle edge.
            ax.plot(x, y, 'w--')
            ax.plot(x, y1, 'w--')
            
            # Axes labels
            ax.set_ylabel("Depth / $km$")
            ax.set_xlabel("Width / $km$")
            
            # Next step.
            i+=step
            
            
        fig.suptitle("{} solution to the lithosphere problem {} enrichment".format(name,enrich))
        fig.subplots_adjust(wspace=0.7, hspace=0.7)
        
        plt.savefig(directory + "{}_solution_{}_enrichment.png".format(name,enrich))
        plt.savefig(directory + "{}_solution_{}_enrichment.pdf".format(name, enrich))

