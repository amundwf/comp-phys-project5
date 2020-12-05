import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def v_n(x,t,n,k,k2): # Calculates n-th element in the sum/series.
    return (-1)**(n+1)*np.sin(k*n*x)*np.exp(-k2*(n**2)*t)

def f(x,L): # Defined by v(x,t) = u(x,t) + f(x)
    return -x/L

calculate_function = True
if calculate_function == True:
    N_x = 100
    dt = 1
    N_t = 20

    L = 1
    k = np.pi/L
    k2 = k**2

    xList = np.linspace(0,1,N_x)
    tList = np.arange(0, N_t*dt, dt)
    v_array = np.zeros((len(tList), len(xList)))

    N_sum = 1000 # The number of terms to include in the sum. The analytical 
    # expression for v(x,t) has N_sum = infinity, but it has to be approximated
    # here. The higher N_sum is the better approximation.

    for i_t in range(len(tList)):
        t = tList[i_t]
        print("t: " + str(t))
        for i_x in range(len(xList)):
            x = xList[i_x]
            v_xt = 0 # v(x,t)
            for n in range(1,N_sum+1):
                v_xt +=v_n(x,t,n,k,k2)
            # Put the new calculated vn value
            v_array[i_t, i_x] = v_xt

    # Save the results (.csv):
    np.savetxt("../results/1D_diffusion/v_xt.csv", v_array, delimiter=",")


# Load the results:
#data = np.loadtxt("../results/1D_diffusion/v_xt.csv", skiprows=0, delimiter=",")
# Get the columns of data as lists:
#data = pd.DataFrame(data, rows=["MC_cycle", "E_mean", "E2_mean", "M_mean", \
#    "M_abs_mean", "M2_mean", "C_V", "chi"])
dataFrame = pd.read_csv("../results/1D_diffusion/v_xt.csv", sep=',',header=None)
v_array = dataFrame.values


t_index = 0
v_t_list = v_array[t_index,:] # plots v(x) at some time t.
plot1, = plt.plot(xList, v_t_list)
#plt.legend(plot1, )
plt.legend()
plt.grid()
plt.xlabel(r'x')   # r means 'render'?
plt.ylabel(r'$v(x,t)$') # $$ for latex
#plt.xlim(0,max(xliste)*1.05)  # Sets the limits for the x axis.
#plt.ylim(0,max(yliste)*(1.05)) # Sets the limits for the y axis.
plt.suptitle('')
plt.show()