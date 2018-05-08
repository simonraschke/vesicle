#!/usr/bin/python3.6

import os
import sys
import argparse
import pprint
import numpy as np
import MDAnalysis as mda
import matplotlib as mpl
import matplotlib.pyplot as plt
import paper_style as style

parser = argparse.ArgumentParser()
parser.add_argument("--top", type=str, help="path to topology file")
parser.add_argument("--gro", type=str, help="path to trajectory file")
args = parser.parse_args()

# a paper_style for plots / can be removed
style.setStyle()
# possible alternative
# plt.style.use('seaborn-paper')

universe = mda.Universe(args.top, args.gro)
print(universe)

# divide universe into groups
solvent = universe.select_atoms("resname OSMOT and name O")
nonsolvent = universe.select_atoms("resname FRAME") + universe.select_atoms("resname MOBIL")

# calc center of mass per residue
centers_of_masses = [ res.atoms.center_of_geometry() for res in nonsolvent.residues ]
# divide all coords by 10, to get the original units back
centers_of_masses = np.divide(centers_of_masses, 10)

# center of micelle from all centers of mass from residues
center_of_micelle =  np.mean(centers_of_masses, axis=0)
print(center_of_micelle)

#  list of distances from residue to center of micelle
dist_hist = [ np.linalg.norm(c-center_of_micelle) for c in centers_of_masses ]

print("average", np.average(dist_hist), "  variance", np.var(dist_hist))

# a new matplotlib figure
fig = plt.figure()

# make a histogram from distances
# bins=auto -> no need to pass bin array
# density=True -> norm to 1
n, bins, patches = plt.hist(dist_hist, bins='auto', density=True)

# label axis
plt.xlabel('r')
plt.ylabel('probability')

#limit x axis to covered area, y will be adjusted automatically
plt.xlim(min(dist_hist) - 0.1, max(dist_hist) + 0.1)

#output
fig.tight_layout()
plt.plot(rasterized=False)
plt.savefig("center_rdf"+'.'+'eps'.format(), bbox_inches='tight', format='eps')
plt.plot(rasterized=True)
plt.savefig("center_rdf"+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)