# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 13:19:01 2020

@author: jamesbc
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import pandas as pd
import os
import re

directory = "../results/1D_diffusion/"

fig, ax = plt.subplots(1, 3, figsize=(9,9), dpi = 80)
ax = ax.flatten()
j=0
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        
        
        
        data = np.loadtxt(directory + filename, delimiter=",")
        data = pd.DataFrame(data)
        print(filename)
        print(data.shape)
       
        
        if re.search('explicit', filename):
            d = 1
            
        pcm = ax[j].pcolor(data)
        ax[j].set_xlabel("distance")
        ax[j].set_ylabel("time")
        ax[j].set_yscale("linear")
        ax[j].set_title("{}".format(filename))
        fig.tight_layout(pad = 2.0)
        fig.colorbar(pcm, ax = ax[j])
        j+=1
        
