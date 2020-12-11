import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import re
import os
import utils as ut

runCppCode = 0
# Set to false if you have the results files and just want to plot the results.
if runCppCode == True: 
    # Compile and run the C++ files (this is exactly what is in the makefile):
    os.system("echo compiling C++ codes...")
    #os.system("g++ -O3 -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/ising_solver.cpp -larmadillo")
    os.system("g++ -O3 -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp -larmadillo -fopenmp")
    os.system("echo executing...")
    os.system("./main.out")

resultsDirectory = "../results/1D_diffusion/"

# Get the axes values (time t and space x):
df_xList = pd.read_csv(resultsDirectory + "analytical_1D_xList.csv", sep=',',header=None)
df_tList = pd.read_csv(resultsDirectory + "analytical_1D_tList.csv", sep=',',header=None)
xList = np.array(df_xList.values)
tList = np.array(df_tList.values)

# Read the analytical results:
df_u_xt = pd.read_csv(resultsDirectory + "analytical_1D.csv", sep=',',header=None)
u_xt_array = np.array(df_u_xt.values)
df_v_xt = pd.read_csv(resultsDirectory + "v_xt.csv", sep=',',header=None)
v_xt_array = np.array(df_v_xt.values)

# Read the explicit results:
df_explicit = pd.read_csv(resultsDirectory + "explicit_1D.csv", sep=',',header=None)
u_explicit_array = np.array(df_explicit.values)

# Read the implicit results:
df_implicit = pd.read_csv(resultsDirectory + "implicit_1D.csv", sep=',',header=None)
u_implicit_array = np.array(df_implicit.values)

# Read the Crank Nicolson results:
df_crankNicolson = pd.read_csv(resultsDirectory + "crankNicolson_1D.csv", sep=',',header=None)
u_crankNicolson_array = np.array(df_crankNicolson.values)

#print(u_xt_array)

X, T = np.meshgrid(xList, tList)
#plt.pcolormesh(dataFrame)
#plt.pcolor(xList, tList, u_xt_array)#, cmap=cm)

plotCodeWord = 'colormesh'
#plotCodeWord = 'oneFrame' # Plots u for one specified time.
#plotCodeWord = 'animate' # Animates the time steps using matplotlib.animate.
t_index = 10 # Which timestep should be plotted?
time = tList[t_index]

if plotCodeWord == 'colormesh':
    # Get a figure with 2x2 subplots:
    fig, ax = plt.subplots(2,2)#, figsize=(9,9), dpi = 80)

    # Plot the analytical solution (top left subplot):
    # Get extended arrays for plotting using pcolormesh (pcolormesh is bugged; it
    # doesn't show the right and bottom sides of the plotted array):
    u_ext, X_ext, T_ext = ut.get_extended_array_and_meshgrid(u_xt_array, xList, tList)
    pcm = ax[0,0].pcolormesh(X_ext,T_ext,u_ext)
    ax[0,0].set_xlabel(r"$x$")
    ax[0,0].set_ylabel(r"$t$")
    ax[0,0].set_yscale("linear")
    ax[0,0].set_title("Analytical")
    fig.colorbar(pcm, ax=ax[0,0])
    #ax[0,0].invert_yaxis()
    
    # Plot the results from the explicit scheme:
    u_ext, X_ext, T_ext = ut.get_extended_array_and_meshgrid(u_explicit_array, xList, tList)
    pcm = ax[0,1].pcolormesh(X_ext,T_ext,u_ext)
    ax[0,1].set_xlabel(r"$x$")
    ax[0,1].set_ylabel(r"$t$")
    ax[0,1].set_yscale("linear")
    ax[0,1].set_title("Explicit scheme")
    fig.colorbar(pcm, ax=ax[0,1])

    # Plot the results from the implicit scheme:
    u_ext, X_ext, T_ext = ut.get_extended_array_and_meshgrid(u_implicit_array, xList, tList)
    pcm = ax[1,0].pcolormesh(X_ext,T_ext,u_ext)
    ax[1,0].set_xlabel(r"$x$")
    ax[1,0].set_ylabel(r"$t$")
    ax[1,0].set_yscale("linear")
    ax[1,0].set_title("Implicit scheme")
    fig.colorbar(pcm, ax=ax[1,0])

    # Plot the results from the Crank-Nicolson scheme:
    u_ext, X_ext, T_ext = ut.get_extended_array_and_meshgrid(u_crankNicolson_array, xList, tList)
    pcm = ax[1,1].pcolormesh(X_ext,T_ext,u_ext)
    ax[1,1].set_xlabel(r"$x$")
    ax[1,1].set_ylabel(r"$t$")
    ax[1,1].set_yscale("linear")
    ax[1,1].set_title("Crank-Nicolson scheme")
    fig.colorbar(pcm, ax=ax[1,1])


elif plotCodeWord == 'oneFrame':
    t = tList[t_index]
    u_xt_list = u_xt_array[t_index, :] # u(x,t) at the specific time t.
    v_xt_list = v_xt_array[t_index, :]
    minus_f_x_list = (-ut.f(xList,1))[:,0]

    u_xt_list_simple_sum = v_xt_list + minus_f_x_list

    #plt.plot(xList, u_xt_list, label=r'$u(x,t) = v(x,t)+x/L$')
    plt.plot(xList, v_xt_list, label=r'$v(x,t)$', color='b')
    plt.plot(xList, -ut.f(xList,1), label=r'$x/L = -f(x)$', color='y')
    plt.plot(xList, u_xt_list_simple_sum, label=r'$u(x,t) = v(x,t)+x/L$', color='r')
    
    #plt.ylabel("u(x,t=" + str(t) + ")")

elif plotCodeWord == 'animate':
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    line, = plt.plot([], [], 'ro')

    def init():
        ax.set_xlim(-0.01, 0.99)
        #ax.set_ylim(-1, 40)
        ax.set_ylim(-1, v_xt_array.max()*1.3)
        return line,

    print('tList:'); print(tList)

    def update(frame): # frame = t_index basically
        t = tList[frame]

        u_xt_list = u_xt_array[frame, :] # u(x,t) at the specific time t.

        v_xt_list = v_xt_array[frame, :]
        minus_f_x_list = (-ut.f(xList,1))[:,0]
        u_xt_list_simple_sum = v_xt_list + minus_f_x_list

        xdata.append(t)
        #ydata.append(u_xt_list_simple_sum[frame])
        ydata.append(minus_f_x_list[frame])
        line.set_data(xdata, ydata)
        return line,
    
    def animate(i):
        x = xList

        v_xt_list = v_xt_array[i, :]
        minus_f_x_list = (-ut.f(xList,1))[:,0]
        u_xt_list_simple_sum = v_xt_list + minus_f_x_list
        y = u_xt_list_simple_sum

        line.set_data(x, y)
        return line,
    #ani = FuncAnimation(fig, update, frames=np.array(range(len(tList))),
                        #init_func=init, blit=True)
    ani = FuncAnimation(fig, animate, init_func=init, frames=len(tList), interval=20, blit=True)
    plt.ylabel(r'$u(x,t)$')

plt.xlabel(r'$x$')
plt.suptitle('Analytical 1D solution, diffusion equation')
plt.legend()
plt.show()

 
