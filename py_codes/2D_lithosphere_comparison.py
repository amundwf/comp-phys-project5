# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:57:49 2020

@author: jamesbc
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams.update({'font.size': 12})
import pandas as pd
import os
import re


directory = "../results/Temperature difference before and after/"


# PLease enter the x dimension (Npoints) from the file you want
# to load. 
x_dim = 120
dx = 1 #km

fig, ax = plt.subplots(1,1, figsize=(9,9), dpi = 80)

for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        
        cube = np.loadtxt(directory + filename)
        
        cube = cube.reshape(-1,x_dim,x_dim)
        print("The cube shape is: ",cube.shape)
        
        t_dim = cube.shape[0]
        
        x = np.arange(x_dim)
        x1 = np.ones((x_dim))*20
        x2 = np.ones((x_dim))*40
        
        y = np.linspace(281, 1800, x_dim)
        
        cube = cube[:,:,:]
        
        if re.search("Before", filename):
            name = "without enrichment"
            
            ax.plot(x, cube[1,:,60], label = 'Steady-State\nwithout enrichment')
            
        elif re.search("After", filename):
            name = "with enrichment"
            
            times = [1,99]
            for t in times:
    
                ax.plot(x, cube[t,:,60], label = 'Time = {} $Gy$\nwith enrichment'.format(t*0.01))

        ax.plot(x1, y, 'k--')
        ax.plot(x2, y, 'k--')
        
        
        ax.set_ylabel("Temperature /Kelvin")
        ax.set_xlabel("Depth /km")
        ax.text(5, 1600, "Upper\nCrust")
        ax.text(25, 1600, "Lower\nCrust")
        ax.text(60, 1600, "Mantle")
     
plt.legend()       
            
fig.suptitle("A cross-section of the lithosphere comparing temperature \ndistribution both with and without radioactive enrichment")

plt.savefig(directory + "temp_diff.png")
plt.savefig(directory + "temp_diff.pdf")

'''
fig, ax = plt.subplots(1,1, figsize=(9,9), dpi = 80)

for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        
        if re.search("Before", filename):
            name = "without enrichment"
            colour = "r"
        elif re.search("After", filename):
            name = "with enrichment"
            colour = 'b'
            
        cube = np.loadtxt(directory + filename)
        
        cube = cube.reshape(-1,x_dim,x_dim)
        print("The cube shape is: ",cube.shape)
        
        t_dim = cube.shape[0]
        
        x = np.arange(x_dim)
        x1 = np.ones((x_dim))*20
        x2 = np.ones((x_dim))*40
        
        y = np.linspace(281, 1800, x_dim)
        
        cube = cube[:,:,:]
        
        label_added = False
        
        for t in range(1, t_dim):
            
            if label_added == False:
                # For some reason the y axis needs flipping.
                ax.plot(x, cube[t,:,60], '{}-'.format(colour), label="{}".format(name))
                label_added = True
            else:
                ax.plot(x, cube[t,:,60], '{}-'.format(colour))
                
        ax.plot(x1, y, 'k--')
        ax.plot(x2, y, 'k--')
        
        
        ax.set_ylabel("Temperature /Kelvin")
        ax.set_xlabel("Depth /km")
        ax.text(5, 1600, "Upper\nCrust")
        ax.text(25, 1600, "Lower\nCrust")
        ax.text(60, 1600, "Mantle")
     
plt.legend()       
            
fig.suptitle("A cross-section of the lithosphere comparing temperature \ndistribution both with and without radioactive enrichment")

plt.savefig(directory + "temp_diff.png")
plt.savefig(directory + "temp_diff.pdf")
'''