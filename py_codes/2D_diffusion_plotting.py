# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 10:35:41 2020

@author: jamesbc
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams.update({'font.size': 12})
import pandas as pd
import os
import re

directory = "../results/2D_diffusion/"
             
# PLease enter the x dimension (Npoints) from the file you want
# to load. 
x_dim = 10
dt = 1e-3
dx = 0.1

fig, axes = plt.subplots(3,2, figsize=(9,9), dpi = 80)
axes = axes.flatten()

for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        
        cube = np.loadtxt(directory + filename)
        
        cube = cube.reshape(-1,x_dim,x_dim)
        print("The cube shape is: ",cube.shape)
        
        t_dim = cube.shape[0]
        x = np.arange(x_dim)/x_dim
        step = int(t_dim / len(axes))
        
        if re.search("Analytic", filename):
            name = "Analytic"
        else:
            name = "Numerical"
            
        # Save computation and plot mesh grid too. 
        ################
        fig1, axes1 = plt.subplots(3,2, figsize=(9,9), dpi = 80)
        axes1 = axes1.flatten()

        i=0
        for ax in axes1:
            
            # For some reason the y axis needs flipping.
            pcm = ax.pcolormesh(cube[i,::-1,:], cmap="viridis")#, vmin = 0.0, vmax=1.0)
            #ax.plot(x, y, 'w--')
            #ax.plot(x, y1, 'w--')
            ax.set_title("t={:.2e}".format(i*dt))
            
            i+=step
            cbar = fig.colorbar(pcm, ax = ax, format="%.1e")
            cbar.set_label("u(x,y,t)")
            
        fig1.suptitle("{} solution to two-dimensional diffusion equation".format(name))
        fig1.subplots_adjust(wspace=0.7, hspace=0.7)
        
        plt.savefig(directory + "{}.png".format(filename))
        plt.savefig(directory + "{}.pdf".format(filename))
        
        
        ################
            
        i = 0 
        
        for ax in axes:
            ax.plot(x, cube[i,:,5], label = "{}".format(name))
            i+=step

            ax.set_ylabel("$u(x,y,t)$")
            ax.ticklabel_format(axis='y',style='sci')
            ax.set_xlabel("x")
            ax.set_title("t = {:.2e},\ny = 0.5".format(i*dt))
            ax.legend()
            
fig.subplots_adjust(wspace=0.7, hspace=0.7)
    

plt.savefig(directory + "dx={}_dt={}.png".format(dx, dt))
plt.savefig(directory + "dx={}_dt={}.pdf".format(dx, dt))     
        