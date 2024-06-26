#! /usr/bin/env/python3

# Script pout tracer les sorties timehistory avec matplotlib

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pathlib


parser = argparse.ArgumentParser(description="Plot the TimeHistory file in argument")
parser.add_argument("filepath", type=str, help="Path to the file to plot")
args = parser.parse_args()

title = pathlib.PurePath(args.filepath).name
var_name = title.split(sep="_")[0]
table = np.loadtxt(args.filepath)
plt.figure()
plt.plot(table[:, 0], table[:, 1])
plt.xlabel("Time [s]")
plt.ylabel(var_name)
plt.title(title, weight="bold")
plt.show()
