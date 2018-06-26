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
parser.add_argument("--time", type=float, nargs=2, default=[1e4,1e8], help="path to config_files.json file")
parser.add_argument("--dens", type=float, nargs=2, default=[0,1], help="density range")
parser.add_argument("--con", type=str, nargs='*', default=[], help="constrain to parameters as {...}")
parser.add_argument("--min", type=int, default=20, help="min size of clusters to analyze")
parser.add_argument("--out", type=str, default="gamma", help="output file name")
args = parser.parse_args()

style.setStyle()
fig = plt.figure()


constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})

max_gamma = 0
max_rho_free = 0

print("densities", plthelp.getValuesInRange(args.file, "density", args.dens))
for t in plthelp.getMatchedValues(args.file, "temperature", constraints):
    newconstraints = constraints
    newconstraints.update({"temperature":str(t)})
    
    rho_free,gamma = plthelp.getGamma(args.file, newconstraints, args.time, args.min)
    rho_free,gamma = plthelp.removeBadEntries2D(rho_free,gamma)

    if max(gamma) > max_gamma: max_gamma = max(gamma)
    if max(rho_free) > max_rho_free: max_rho_free = max(rho_free)
    
    print("rho_free   gamma")
    for rf, g in zip(rho_free,gamma):
        try:
            print("  {:.3f}".format(rf), "  {:.5f}".format(g))
        except TypeError:
            print("  {:.3f}".format(rf), "  None")
    print()
    label = "T="+str(t)
    plt.scatter(rho_free,gamma, label=label)

plt.legend(loc='best', frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.ylim(0.0, max_rho_free+0.001)
plt.ylim(0.0, max_gamma+0.001)
plt.xlabel(r'$\rho_\mathrm{free}^{}$')
plt.ylabel(r'$\Gamma_{(\mathrm{c_{first}^{}})}^{}$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fig.tight_layout()
plt.plot(rasterized=False)
plt.savefig(args.out+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)