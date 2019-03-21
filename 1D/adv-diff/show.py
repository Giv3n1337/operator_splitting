import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

res = str(sys.argv[1])

data_dir = os.getcwd() + '/tmp/data/convergence_test_states/'
fname_prefix = 'advection_n-' + res + '_t-'
fname_suffix = '.txt'

fnumbers = []
dt = 10.0000 * 0.5000 * 1.0000 / int(res)

for t in np.arange(0.0000+dt, 5.0000, dt):
    number = ('%1.4f' %t)
    fnumbers.append(number)


fnumbers = np.array(fnumbers)
fnumbers.sort()

fig = plt.figure()

def get_data_from_file(i):

    global data_dir
    global fname_prefix
    global fnumbers
    global fname_suffix

    if fnumbers.size < i+1:
        return np.array([])

    fname = data_dir + fname_prefix + fnumbers[i] + fname_suffix
    return np.loadtxt(fname)

def anim(i):
    global res

    data = get_data_from_file(i)

    if data.size == 0:
        return

    plt.clf()
    plt.plot(data)
    axes = plt.gca()
    axes.set_ylim([-1,1])
    

ani = animation.FuncAnimation(fig, anim, interval=1)
plt.show()




