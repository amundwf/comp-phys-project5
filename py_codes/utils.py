# This file contains useful functions to be used in other .py files.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def f(x,L): # Defined by v(x,t) = u(x,t) + f(x)
    x = np.array(x) # Make x numpy compatible
    return -x/L
