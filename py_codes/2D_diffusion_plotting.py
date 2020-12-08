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
directory = "../results/2D_diffusion/"

for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        if re.search("ana", filename):
            name = "Analytic"
        else:
            name = "Numerical"
        fig, axes = plt.subplots(4,5, figsize=(9,9), dpi = 80)
        axes = axes.flatten()
        cube = np.loadtxt(directory + filename)
        
        # PLease enter the x dimension in the csv file.
        x_dim = 100
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
            
            i+=1
            cbar = fig.colorbar(pcm, ax = axes, format="%.2e")
            cbar.set_label("Temperature / Kelvin")
        fig.suptitle("{} Solution".format(name))
        plt.savefig(directory + "BeforeEnrichment.png")
        #plt.savefig(directory + "BeforeEnrichment.pdf")


'''
for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        fig, axes = plt.subplots(3, 7, figsize=(9,9), dpi = 80)
        axes = axes.flatten()
        cube = np.loadtxt(directory + filename)
        
        # PLease enter the x dimension in the csv file.
        x_dim = 100
        cube = cube.reshape(-1,x_dim,x_dim)
        
        print("Cube shape is: ",cube.shape)
        
       
        t_dim = cube.shape[0]
        
        step = int(t_dim / len(axes))
        print(step)
        i=0
      
        row = int(x_dim/2)
        
        for ax in axes:
            pcm = ax.plot(cube[i,50,:])
            ax.set_title("time ind={}".format(i))
            i+=step
'''    
        
        
        
        
        
        