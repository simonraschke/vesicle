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
from scipy.stats import norm
import plot_helper_functions as plthelp
import paper_style as style

pp = pprint.PrettyPrinter(indent=4, compact=False)

parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--file", type=str, help="path to config_files.json file")
parser.add_argument("--time", type=float, nargs='*', default=[1e5,1e8], help="time range of cluster size distributions")
parser.add_argument("--num", type=int, default=10, help="number of time points")
parser.add_argument("--con", type=str, nargs='*', default=[], help="constrain to parameters as {...}")
# parser.add_argument("--opt", type=float, default=100, help="Add optimum line at N")
parser.add_argument("--max", type=float, default=0.05, help="max possibility on y axis")
parser.add_argument("--out", type=str, default="hist", help="output file name")
args = parser.parse_args()

style.setStyle()

args.time = [ round(x,-4) for x in np.logspace(np.log10(min(args.time)), np.log10(max(args.time)), num=args.num, endpoint=True, base=10.0, dtype=float) ]
print("timepoints", args.time)

"""
datadict
  -> temperature
       -> timepoints
       -> means
       -> variances
"""
datadict = {}

constraints = {}
for i in range(len(args.con[::2])):
    constraints.update({args.con[i*2]:args.con[i*2+1]})

for t in plthelp.getMatchedValues(args.file, "temperature", constraints):
    newconstraints = constraints
    newconstraints.update({"temperature":str(t)})

    datadict.update( {t : {}} )
    datadict[t].update( {"timepoints" : []} )
    datadict[t].update( {"means" : []} )
    datadict[t].update( {"stds" : []} )
    datadict[t].update( {"variances" : []} )

    fig = plt.figure(1)
    gs = mpl.gridspec.GridSpec(len(args.time), 1, wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.17, right=0.845) 

    subplot_counter = 0
    max_size = 0
    for timepoint in args.time:
        timepoint = int(timepoint)

        data = plthelp.getClusterSizeDistribution(args.file, {**newconstraints}, timepoint)
        if max(data) > max_size: max_size = max(data)
        bins = np.arange(1,max(data),1,dtype=int)

        datadict[t]["timepoints"].append(timepoint)
        datadict[t]["means"].append(np.mean(list(filter(lambda x: x > 30, data))))
        datadict[t]["stds"].append(np.std(list(filter(lambda x: x > 30, data))))
        datadict[t]["variances"].append(np.var(list(filter(lambda x: x > 30, data))))

        label = "step "+'{0:.2e}'.format(timepoint)
        plt.subplot(gs[subplot_counter,0])

        plt.gca().tick_params(axis='both', which='both', direction='in')
        plt.xlim(0,200)
        plt.ylim(0,args.max)
        plt.gca().locator_params(axis='y', nbins=3)

        data, bins, patches = plt.hist( data, bins=bins, weights=data, label=label, density=True)

        plt.gca().legend(loc='best', frameon=False, fontsize='xx-small', handlelength=0)

        if subplot_counter+1 < len(args.time): plt.gca().set_xticklabels([])
        if subplot_counter+1 != len(args.time): plt.gca().set_yticks(plt.gca().get_yticks()[1:])
        if subplot_counter+1 == int((len(args.time)+1)/2): plt.ylabel(r'\textrm{probability}')
        
        subplot_counter+=1
        # print()

    pp.pprint(datadict[t])

    # plt.plot([0,0.1],[args.opt, args.opt], color="grey",zorder=0, linestyle='--', linewidth=.5)
    plt.xlabel(r'\textrm{N}')
    # plt.tight_layout()
    plt.plot(rasterized=False)
    plt.savefig(args.out+'_T'+str(t)+'.'+'eps'.format(), bbox_inches='tight', format='eps')
    plt.plot(rasterized=True)
    plt.savefig(args.out+'_T'+str(t)+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)

print("\n\nComplete Data")
pp.pprint(datadict)

fig = plt.figure(2)
gs = mpl.gridspec.GridSpec(len(list(datadict.keys())), 1, wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.17, right=0.845) 
subplot_counter = 0

for t,d in datadict.items():
    x = d["timepoints"]
    y = d["means"]
    e = d["stds"]
    v = d["variances"]
    label = "T="+str(t)

    plt.subplot(gs[subplot_counter,0])
    plt.gca().tick_params(axis='both', which='both', direction='in')
    plt.ylim(60,120)
    plt.gca().locator_params(axis='y', nbins=3)
    plt.gca().set_xscale("log", nonposx='clip')
    plt.errorbar(x,y,yerr=e, label=label, elinewidth=.5, ecolor="black", capsize=2, linestyle='None', marker=".", markersize=5, capthick=.5)

    for i in range(len(x)):
        if i == len(x)-1:
            plt.gca().text(s=r'$\bar N={}$'.format(str(round(y[i],1))), x=x[i]*0.98, y=y[i], fontsize=2, stretch='extra-condensed', ha='right', va='bottom')
            plt.gca().text(s=r'$\sigma={}$'.format(str(round(e[i],1))), x=x[i]*0.98, y=y[i], fontsize=2, stretch='extra-condensed', ha='right', va='top')
        else:
            plt.gca().text(s=r'$\bar N={}$'.format(str(round(y[i],1))), x=x[i]*1.02, y=y[i], fontsize=2, stretch='extra-condensed', va='bottom')
            plt.gca().text(s=r'$\sigma={}$'.format(str(round(e[i],1))), x=x[i]*1.02, y=y[i], fontsize=2, stretch='extra-condensed', va='top')

    plt.text(1300000, 118, label, ha='center', va='top', fontsize='xx-small')
    # plt.gca().legend(loc='best', frameon=True, fontsize='xx-small', handlelength=0, markerscale=0)
    # handles, labels = plt.gca().get_legend_handles_labels()
    # handles = [h[0] for h in handles]
    # plt.gca().legend(handles, labels, loc='best', frameon=True, fontsize='xx-small', handlelength=0, markerscale=0)

    if subplot_counter+1 < len(list(datadict.keys())): plt.gca().set_xticklabels([])
    if subplot_counter+1 != len(list(datadict.keys())): plt.gca().set_yticks(plt.gca().get_yticks()[1:])
    if subplot_counter+1 == int((len(list(datadict.keys()))+1)/2): plt.ylabel(r'\textrm{<N>}')
    subplot_counter+=1

plt.xlabel(r'\textrm{time step}')
plt.plot(rasterized=True)
plt.savefig(args.out+'_complete.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)
plt.plot(rasterized=False)
plt.savefig(args.out+'_complete.'+'eps'.format(), bbox_inches='tight', format='eps')