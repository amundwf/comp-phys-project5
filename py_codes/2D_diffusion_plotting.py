# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 10:35:41 2020

@author: jamesbc
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams.update({'font.size': 15})
import pandas as pd
import os
import re

#directory = "../results/2D_diffusion_before_enrichment/"
#directory = "../results/2D_diffusion_after_enrichment/"
directory = "../results/2D_diffusion/"


# PLease enter the x dimension (Npoints)from the file you want
# to load. 
x_dim = 10
        
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
        
        i=0
        for ax in axes:
            print(i)
            pcm = ax.pcolormesh(cube[i,:,:], cmap="viridis")#, vmin = 0.0, vmax=1.0)
            #ax.set_title("time ind={}".format(i))
            #ax.ticklabel_format(style = 'sci', scilimits = (0,0))
            #ax.set_xticklabels([])
            #ax.set_yticklabels([])
            #ax.set_xticks([])
            #ax.set_yticks([])
            
            i+=step
            cbar = fig.colorbar(pcm, ax = ax, format="%.1e")
            cbar.set_label("Temperature /Kelvin")
            
        fig.suptitle("{} Solution to the two dimensional diffusion equation".format(name))
        
        plt.savefig(directory + "{}.png".format(name))
        plt.savefig(directory + "{}.pdf".format(name))

        
        
        
        
        
        