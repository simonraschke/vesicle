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
# parser.add_argument("--file", type=str, help="path to config_files.json file")
# parser.add_argument("--time", type=float, nargs=2, default=[0,1e10], help="path to config_files.json file")
# parser.add_argument("--con", type=str, nargs='*', default=[], help="constrain to parameters as {...}")
# parser.add_argument("--min", type=int, default=1, help="min size of clusters to analyze")
parser.add_argument("--out", type=str, default="tau", help="output file name")
args = parser.parse_args()

style.setStyle()
fig = plt.figure()

label = "T=0.26"
x = [0.001, 0.01, 0.015, 0.02,  0.03, 0.04, 0.05, 0.06, 0.07,   0.08,   0.09,   0.1]
y = [0,     0,    1.1e7, 1.1e6, 8e4,  4e4,  3e4,  2e4,  1.15e4, 1.15e4, 1.15e4, 1.15e4]
plt.semilogy(x[2:], y[2:], label=label)

label = "T=0.27"
x = [0.001, 0.01, 0.015, 0.02,   0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
y = [0,     0,    0,     1.05e7, 6e5,  5e4,  5e4, 3e4, 3e4, 2e4, 1e4, 1e4]
plt.semilogy(x[3:], y[3:], label=label)

label = "T=0.28"
x = [0.001, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
y = [0,     0,    0,     0,    9e5,  9e4,  7e4,  5e4,  4e4,  2e4, 1.5e4, 1.5e4]
plt.semilogy(x[4:], y[4:], label=label)

label = "T=0.29"
x = [0.001, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
y = [0,     0,    0,     0,    0,    7e4,  8e4,  3e4,  2.5e4, 2e4, 1.5e4, 1.5e4]
plt.semilogy(x[5:], y[5:], label=label)


plt.legend(loc='best', frameon=False)
plt.gca().tick_params(axis='both', which='both', direction='in')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\tau$')
plt.plot(rasterized=False)
plt.savefig(args.out+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig(args.out+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)