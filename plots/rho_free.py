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
parser.add_argument("--out", type=str, default="rho_free", help="output file name")
parser.add_argument("--read", type=str, default=None, help="input .json to read data from")
parser.add_argument("--extra", action='store_true', help="show the curve for small times")
args = parser.parse_args()

style.setStyle()
fig = plt.figure()

constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})


datadict = {}
if args.read:
    datadict = plthelp.readJson("rho_free.json")

    base_keys = [ str(key) for key in datadict.keys() if not "extra" in key ]
    print(base_keys)
    extra_keys = [ str(key) for key in datadict.keys() if "extra" in key and key.split("extra")[0] in datadict.keys() ]
    print(extra_keys)

    for t in base_keys:
        label = "T="+str(t)
        line = plt.plot(datadict[t]["x"], datadict[t]["y"], label=label, marker="s", zorder=5) 
        if args.extra:
            color = line[-1].get_color()
            line = plt.plot(datadict[t+"extra"]["x"], datadict[t+"extra"]["y"], linestyle='--', color=color, marker="s", markeredgecolor=color, markerfacecolor='white', zorder=4.9)

else: 
    max_rho_free = 0
    for t in plthelp.getMatchedValues(args.file, "temperature", constraints):
        newconstraints = constraints
        newconstraints.update({"temperature":str(t)})

        rho = [float(x) for x in plthelp.getMatchedValues(args.file,"density", newconstraints)]
        rho_free = []

        for density in rho:
            rho_free.append(plthelp.getFreeParticleDensity(args.file, {**newconstraints,**{"density":str(density)}}, args.time))

        rho,rho_free = plthelp.removeBadEntries2D(rho,rho_free)
        if max(rho_free) > max_rho_free: max_rho_free = max(rho_free)

        print("rho    rho_free")
        for r, rf in zip(rho,rho_free):
            try:
                print("  {:.3f}".format(r), "  {:.5f}".format(rf))
            except TypeError:
                print("  {:.3f}".format(r), "  None")
        print()

        label = "T="+str(t)
        line = plt.plot(rho,rho_free, label=label) 

        datadict.update( { str(t) : newconstraints  } )
        datadict.update( { str(t) : {"time":args.time, "x":rho , "y":rho_free}  } )

        if args.extra:
            color = line[-1].get_color()
            rho_free = []

            for density in rho:
                rho_free.append(plthelp.getFreeParticleDensity(args.file, {**newconstraints,**{"density":str(density)}}, [90000,110000]))

            rho,rho_free = plthelp.removeBadEntries2D(rho,rho_free)
            if max(rho_free) > max_rho_free: max_rho_free = max(rho_free)

            print("rho    rho_free")
            for r, rf in zip(rho,rho_free):
                try:
                    print("  {:.3f}".format(r), "  {:.5f}".format(rf))
                except TypeError:
                    print("  {:.3f}".format(r), "  None")
            print()
            line = plt.plot(rho,rho_free, linestyle='--', color=color)

            datadict.update( { str(t)+"extra" : newconstraints  } )
            datadict.update( { str(t)+"extra" : {"time":args.time, "x":rho , "y":rho_free}  } )

        plthelp.writeJson("rho_free.json", datadict)



plt.legend(loc='best', frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
# plt.xlim(0.0, float(max(plthelp.getMatchedValues(args.file,"density",constraints))))
# plt.ylim(0.0, max_rho_free+0.001)
plt.xlim(plt.xlim())
ymin,ymax = plt.ylim()
plt.ylim(ymin, ymax+0.005)
plt.plot([-1,1],[-1,1], color="black",zorder=1)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\rho_\mathrm{free}^{}$')
fig.tight_layout()
plt.plot(rasterized=False)
plt.savefig(args.out+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)