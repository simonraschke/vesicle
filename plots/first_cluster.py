#!/usr/bin/python3.6

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
import matplotlib.pyplot as plt
import plot_helper_functions as plthelp
import paper_style as style

pp = pprint.PrettyPrinter(indent=4, compact=False)

parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--file", type=str, help="path to config_files.json file")
parser.add_argument("--time", type=float, nargs=2, default=[0,1e10], help="path to config_files.json file")
parser.add_argument("--con", type=str, nargs='*', default=[], help="constrain to parameters as {...}")
parser.add_argument("--size", type=int, default=20, help="size of cluster")
parser.add_argument("--out", type=str, default="first_cluster_", help="output file name")
args = parser.parse_args()

style.setStyle()
fig = plt.figure()


constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})

max_time = 0
for t in plthelp.getMatchedValues(args.file, "temperature", constraints):
    newconstraints = constraints
    newconstraints.update({"temperature":str(t)})

    rho = [float(x) for x in plthelp.getMatchedValues(args.file,"density", newconstraints)]
    time = []

    for density in rho:
        time.append(plthelp.getTimeOfFirstCluster(args.file, {**newconstraints,**{"density":str(density)}}, args.time, args.size))

    rho,time = plthelp.removeBadEntries2D(rho,time)
    if max(time) > max_time: max_time = max(time)

    print("rho    time")
    for r, rf in zip(rho,time):
        try:
            print("  {:.3f}".format(r), "  {:.5f}".format(rf))
        except TypeError:
            print("  {:.3f}".format(r), "  None")
    print()
    label = "T="+str(t)
    plt.semilogy(rho,time, label=label)

plt.plot([-1,1],[-1,1], color="black",zorder=1)
plt.legend(loc='best', frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.xlim(0.0, float(max(plthelp.getMatchedValues(args.file,"density",constraints))))
plt.ylim(9000, int(np.ceil(max_time / 10000000.0)) * 10000000)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$t_\textrm{first}$')
fig.tight_layout()
plt.plot(rasterized=False)
plt.savefig(args.out+str(args.size)+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+str(args.size)+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)