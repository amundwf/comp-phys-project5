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

#data = np.loadtxt(directory + filename, delimiter=",")
#data = pd.DataFrame(data)

df_u_xt = pd.read_csv("../results/1D_diffusion/analytical_1D.csv", sep=',',header=None)
u_xt_array = np.array(df_u_xt.values)
df_v_xt = pd.read_csv("../results/1D_diffusion/v_xt.csv", sep=',',header=None)
v_xt_array = np.array(df_v_xt.values)

# Get the axis values (time t and space x):
df_xList = pd.read_csv("../results/1D_diffusion/analytical_1D_xList.csv", sep=',',header=None)
df_tList = pd.read_csv("../results/1D_diffusion/analytical_1D_tList.csv", sep=',',header=None)
xList = np.array(df_xList.values)
tList = np.array(df_tList.values)

X, T = np.meshgrid(xList, tList)
#plt.pcolormesh(dataFrame)
#plt.pcolor(xList, tList, u_xt_array)#, cmap=cm)

#plotCodeWord = 'colormesh'
#plotCodeWord = 'oneFrame' # Plots u for one specified time.
plotCodeWord = 'animate' # Animates the time steps using matplotlib.animate.
t_index = 10 # Which timestep should be plotted?

if plotCodeWord == 'colormesh':
    plt.pcolor(X, T, u_xt_array)#, cmap=cm)
    #plt.pcolormesh(u_xt_array)
    #plot1, = plt.plot(xList, v_t_list)
    #plt.legend(plot1, )
    plt.legend()
    plt.grid()
    plt.ylabel(r'$t$')

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
    plt.ylabel("u(x,t=" + format(float(t), '.7f') + ")")

elif plotCodeWord == 'animate':
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    line, = plt.plot([], [], 'ro')

    def init():
        ax.set_xlim(-0.01, 1.01)
        #ax.set_ylim(-1, 40)
        ax.set_ylim(-1, v_xt_array.max()*1.3)
        return line,

    print('tList:'); print(tList)

    def update(frame): # frame = t_index basically
        t = tList[frame]

        u_xt_list = u_xt_array[frame, :] # u(x,t) at the specific time t.
        #u_xt_list = u_xt_array[frame+1, :]

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
    ani = FuncAnimation(fig, animate, init_func=init, frames=len(tList), interval=50, blit=True)
    plt.ylabel(r'$u(x,t)$')
    

plt.xlabel(r'$x$')
plt.suptitle('Analytical 1D solution, diffusion equation')
plt.legend()
plt.show()

if plotCodeWord=='animate': # Save the animation as a gif:
    ani.save('../results/1D_diffusion/analytical_1D.gif',writer='imagemagick') 