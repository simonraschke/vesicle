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
args = parser.parse_args()

constraints = {"temperature": "0.26", "mobile": "1000"}
rho = [float(x) for x in plthelp.getMatchedValues(args.file,"density", constraints)]
order = []

for density in rho:
    rho_free.append(plthelp.getFreeParticleDensity(args.file, {**constraints,**{"density":str(density)}}, args.time))

print(rho)
print(rho_free)

plt.style.use('seaborn-paper')
fig = plt.figure()
plt.plot(rho,rho_free)
fig.tight_layout()
plt.show()