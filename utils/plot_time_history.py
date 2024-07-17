#! /usr/bin/env/python3

# Copyright 2024 CEA

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Script pout tracer les sorties timehistory avec matplotlib

# Faire un module load python3 avant de lancer le script

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pathlib


parser = argparse.ArgumentParser(description="Plot the TimeHistory file in argument. Launch `module load python3` before launching the script.")
parser.add_argument("filepath", type=str, help="Path to the file to plot")
parser.add_argument("--reference", type=str, help="Path to the reference file to plot")
args = parser.parse_args()

title = pathlib.PurePath(args.filepath).name
var_name = title.split(sep="_")[0] if "Time" not in title else title

table = np.loadtxt(args.filepath)

plt.figure()
plt.plot(table[:, 0], table[:, 1], label="execution")
plt.xlabel("Time [s]")
plt.ylabel(var_name)
plt.title(title, weight="bold")

if args.reference is not None:
  ref_table = np.loadtxt(args.reference)
  plt.plot(ref_table[:, 0], ref_table[:, 1], label="reference")

plt.legend(loc='best')
plt.show()
