import numpy
from math import *
from matplotlib import pyplot as plt

parameters = {}
with open("config.ini", "r") as config:
    for s in config:
        if '#' in s:
            s1 = s.split('#')
        else:
            s1 = s
        if '=' in s1:
            parameter = s1.split('=')[0].strip()
            value = s1.split('=')[1].strip()
            parameters[parameter] = value
periods = float(parameters["periods"])
dt = float(parameters["delta_t"])
n_steps = floor(periods * 2.0 * pi / dt)
m = int(parameters["m"])
n = int(n_steps / m)
t = numpy.zeros(n)
z = numpy.zeros(n)
z_dot = numpy.zeros(n)
with open("result.dat", "r") as inp:
    k = 0
    for s in inp:
        if '#' in s:
            continue
        s_list = s.split()
        t[k] = float(s_list[0])
        z[k] = float(s_list[1])
        z_dot[k] = float(s_list[2])
        k = k + 1
plt.xlabel("$z$")
plt.ylabel("$\\dot{z}$")
plt.plot(z, z_dot, color="blue", linewidth=0.75)
fig = plt.gcf()
fig.savefig("result.png", dpi=300)
plt.cla()
