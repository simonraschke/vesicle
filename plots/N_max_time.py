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
mpl.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import plot_helper_functions as plthelp
import paper_style as style


pp = pprint.PrettyPrinter(indent=4, compact=False)

parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--file", type=str, help="path to config_files.json file")
parser.add_argument("--time", type=float, nargs=2, default=[1e4,1e8], help="path to config_files.json file")
parser.add_argument("--dens", type=float, nargs=2, default=[0,1], help="density range")
parser.add_argument("--opt", type=float, default=100, help="Add optimum line at N")
parser.add_argument("--con", type=str, nargs='*', default=[], help="constrain to parameters as {...}")
parser.add_argument("--out", type=str, default="N_max_time", help="output file name")
args = parser.parse_args()

style.setStyle()
fig = plt.figure()


constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})

max_N_max = 0

print("densities", plthelp.getValuesInRange(args.file, "density", args.dens))
for density in plthelp.getValuesInRange(args.file, "density", args.dens):
    newconstraints = constraints
    newconstraints.update({"density":str(density)})
    
    timepoints,N_max = plthelp.getLargestClusterTimeEvolution(args.file, newconstraints, args.time)
    timepoints,N_max = plthelp.removeBadEntries2D(timepoints,N_max)

    if max(N_max) > max_N_max: max_N_max = max(N_max)
    
    print("timepoints   N_max")
    for t, rf in zip(timepoints,N_max):
        try:
            print("  {:.3f}".format(t), "  {:.5f}".format(rf))
        except TypeError:
            print("  {:.3f}".format(t), "  None")
    print()
    label = r'$\rho=$'
    plt.semilogx(timepoints,N_max, label=label+str(round(density,4)))

plt.plot([min(args.time)+1e4, max(args.time)],[args.opt, args.opt], color="grey",zorder=0, linestyle='--', linewidth=.5)
plt.legend(loc='best', frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.xlim(min(args.time)+1e4, max(args.time))
plt.ylim(0.0, max_N_max+5)
plt.xlabel(r'simulation steps')
plt.ylabel(r'$\mathrm{N}_{\mathrm{max}}$')
fig.tight_layout()
plt.plot(rasterized=False)
plt.savefig(args.out+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)