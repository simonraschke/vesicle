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
parser.add_argument("--out", type=str, default="exchange_time", help="output file name")
args = parser.parse_args()

style.setStyle()
fig = plt.figure()


constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})

max_exchange = 0

print("densities", plthelp.getValuesInRange(args.file, "density", args.dens))
for density in plthelp.getValuesInRange(args.file, "density", args.dens):
    newconstraints = constraints
    newconstraints.update({"density":str(density)})
    
    timepoints,exchange = plthelp.getExchangeTimeEvolution(args.file, newconstraints, args.time)
    timepoints,exchange = plthelp.removeBadEntries2D(timepoints,exchange)
    # exchange = [e*100.0 for e in exchange]
    try:
        if max(exchange) > max_exchange: max_exchange = max(exchange)
    except:
        pass

    print("timepoints   exchange")
    for t, rf in zip(timepoints,exchange):
        try:
            print("  {:.3f}".format(t), "  {:.5f}".format(rf))
        except TypeError:
            print("  {:.3f}".format(t), "  None")
    print()
    label = r'$\rho=$'
    plt.semilogx(timepoints,exchange, label=label+str(round(density,4)))

plt.legend(loc='best', frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.xlim(min(args.time), max(args.time))
plt.ylim(0.0, max_exchange+0.05)
plt.xlabel(r'simulation steps')
plt.ylabel(r'exchange')
plt.gca().yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1,decimals=0))
fig.tight_layout()
plt.plot(rasterized=False)
plt.savefig(args.out+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)