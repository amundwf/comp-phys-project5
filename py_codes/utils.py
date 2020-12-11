# This file contains useful functions to be used in other .py files.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def f(x,L): # Defined by v(x,t) = u(x,t) + f(x)
    x = np.array(x) # Make x numpy compatible
    return -x/L

def get_extended_array(arr):
    arr_ext = np.zeros(np.array(arr.shape)+1)
    arr_ext[:-1,:-1] = arr
    return arr_ext

def get_extended_array_and_meshgrid(arr, xList, yList):
    # Takes in an array and increases its dimensions by 1 (adds one column 
    # and one row), and returns corresponding extended meshgrids from xList
    # and yList.
    arr_ext = get_extended_array(arr)
    
    dx = xList[1]-xList[0]
    dy = yList[1]-yList[0]

    xList_ext = np.append(xList, xList[-1]+dx)
    yList_ext = np.append(yList, yList[-1]+dy)

    [X_ext, Y_ext] = np.meshgrid(xList_ext, yList_ext)

    return arr_ext, X_ext, Y_ext


