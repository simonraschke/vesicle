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
parser.add_argument("--chi", type=float, nargs='*', default=[0,0.25,0.5,1,2], help="chi values")
parser.add_argument("--out", type=str, default="angular_lj", help="output file name")
args = parser.parse_args()

style.setStyle()
plt.rc('figure', figsize=style.makesize(222,ratio=0.8))
fig = plt.figure()


rr = np.arange(0.8, 3, 0.01)
chi_values = args.chi
  
def lj_with_chi(r,chi):
    #return r+float(chi)
    term2 = (1.0/r)**6
    U = 4.0*( term2*term2 - (1.0-float(chi))*term2 )
    return U

potential = np.vectorize(lj_with_chi)

#for _chi in chi_values:
for i, _chi in list(enumerate(chi_values)):    
    labelstring = "$\\mathrm{\chi=\;}$"+str(_chi)
    plt.plot(rr,potential(rr,_chi), label=labelstring)

plt.legend(frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.plot([-10,10], [0,0], color='grey', linestyle='--', linewidth=.5, zorder=0)
plt.ylim(-1.1, 2.1)
plt.xlim( .8, 2.5)
plt.locator_params(axis='x', nbins=4)
plt.locator_params(axis='y', nbins=4)
plt.xlabel(r"$r / \sigma$")
plt.ylabel(r"$U(r,\chi) / \epsilon$")
plt.plot(rasterized=False)
plt.savefig(args.out+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)