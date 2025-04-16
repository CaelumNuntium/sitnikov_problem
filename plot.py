import numpy
from math import *
from matplotlib import pyplot as plt

parameters = {}
with open("config.ini", "r") as config:
    for s in config:
        if '#' in s:
            s1 = s.split('#')[0]
        else:
            s1 = s
        if '=' in s1:
            parameter = s1.split('=')[0].strip()
            value = s1.split('=')[1].strip()
            parameters[parameter] = value
if "poincare_map" not in parameters.keys():
    poincare = False
elif parameters["poincare_map"] == "false":
    poincare = False
elif parameters["poincare_map"] == "true":
    poincare = True
else:
    raise ValueError(f"Invalid value: poincare_map = {parameters['poincare_map']}")
periods = float(parameters["periods"])
if poincare:
    z_lbound = float(parameters["z_limits"].split(':')[0])
    z_rbound = float(parameters["z_limits"].split(':')[1])
    dz_lbound = float(parameters["z_dot_limits"].split(':')[0])
    dz_rbound = float(parameters["z_dot_limits"].split(':')[1])
    n_z = int(parameters["n_z"])
    n_z_dot = int(parameters["n_z_dot"])
    n = (n_z + 1) * (n_z_dot + 1) * floor(periods)
    z = numpy.zeros(n)
    z_dot = numpy.zeros(n)
    with open("poincare_map.dat", "r") as inp:
        k = 0
        for s in inp:
            if '#' in s:
                continue
            s_list = s.split()
            z[k] = float(s_list[0])
            z_dot[k] = float(s_list[1])
            k = k + 1
    plt.xlim(z_lbound, z_rbound)
    plt.ylim(dz_lbound, dz_rbound)
    plt.xlabel("$z$")
    plt.ylabel("$\\dot{z}$")
    plt.scatter(z, z_dot, marker="o", color="black", s=0.05, linewidths=0)
    fig = plt.gcf()
    fig.savefig("poincare_map.png", dpi=1200)
    plt.cla()
else:
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
