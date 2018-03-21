#!/usr/bin/python3

#
#   Copyright 2017-2018 Simon Raschke
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#

import os
import sys
import argparse
import pprint
import numpy as np
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import plot_helper_functions as plthelp

pp = pprint.PrettyPrinter(indent=4, compact=False)

parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--file", type=str, help="path to config_files.json file")
parser.add_argument("--time", type=float, nargs=2, default=[0,1e10], help="path to config_files.json file")
parser.add_argument("--con", type=str, nargs='*', help="constrain to parameters as {...}")
args = parser.parse_args()

fig = plt.figure()


constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})

for t in plthelp.getMatchedValues(args.file, "temperature", constraints):
    newconstraints = constraints
    newconstraints.update({"temperature":str(t)})

    rho = [float(x) for x in plthelp.getMatchedValues(args.file,"density", newconstraints)]
    order = []

    for density in rho:
        order.append(plthelp.getOrder(args.file, {**newconstraints,**{"density":str(density)}}, args.time))

    rho,order = plthelp.removeBadEntries2D(rho,order)

    print("rho    order")
    for r, rf in zip(rho,order):
        try:
            print("  {:.3f}".format(r), "  {:.5f}".format(rf))
        except TypeError:
            print("  {:.3f}".format(r), "  None")
    print()
    label = "T="+str(t)
    plt.plot(rho,order, label=label)

plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\Theta_{s}^{}$')
plt.legend(loc='0')
plt.xlim(0.0, float(max(plthelp.getMatchedValues(args.file,"density",constraints))))
plt.ylim(0.0, 0.025)
fig.tight_layout()
plt.show()