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
mpl.use('Agg')
import matplotlib.pyplot as plt
import plot_helper_functions as plthelp
import paper_style as style

pp = pprint.PrettyPrinter(indent=4, compact=False)

parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--file", type=str, help="path to config_files.json file")
parser.add_argument("--time", type=float, nargs=2, default=[0,1e10], help="path to config_files.json file")
parser.add_argument("--con", type=str, nargs='*', default=[], help="constrain to parameters as {...}")
parser.add_argument("--min", type=int, default=1, help="min size of clusters to analyze")
parser.add_argument("--out", type=str, default="N_max", help="output file name")
args = parser.parse_args()

style.setStyle()
fig = plt.figure()


constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})

max_n_max = 0
for t in plthelp.getMatchedValues(args.file, "temperature", constraints):
    newconstraints = constraints
    newconstraints.update({"temperature":str(t)})

    rho = [float(x) for x in plthelp.getMatchedValues(args.file,"density", newconstraints)]
    n_max = []

    for density in rho:
        n_max.append(plthelp.getLargestCluster(args.file, {**newconstraints,**{"density":str(density)}}, args.time, args.min))
    
    for i in range(len(n_max)):
        if n_max[i] < args.min:
            n_max[i] = np.NaN

    rho,n_max = plthelp.removeBadEntries2D(rho,n_max)
    print(max(n_max))
    if max(n_max) > max_n_max: max_n_max = max(n_max)

    print("  rho    n_max")
    for r, n in zip(rho,n_max):
        try:
            print("  {:.3f}".format(r), "  {:.5f}".format(n))
        except TypeError:
            print("  {:.3f}".format(r), "  None")
    print()
    label = "T="+str(t)
    plt.plot(rho,n_max, label=label)

plt.legend(loc='best', frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.xlim(0.0, float(max(plthelp.getMatchedValues(args.file,"density",constraints))))
plt.ylim(0.0, int(np.ceil(max_n_max / 10.0)) * 10)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\mathrm{N}_{\mathrm{max}}$')
fig.tight_layout()
plt.plot(rasterized=False)
plt.savefig(args.out+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)