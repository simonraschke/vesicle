#!/usr/bin/python3.6

import os
import sys
import argparse
import pprint
import numpy as np
import MDAnalysis as mda
import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
import paper_style as style

from sklearn.cluster import DBSCAN
from sklearn import metrics

parser = argparse.ArgumentParser()
parser.add_argument("--top", type=str, help="path to topology file")
parser.add_argument("--traj", type=str, help="path to trajectory file")
args = parser.parse_args()

# a paper_style for plots / can be removed
style.setStyle()
# possible alternative
# plt.style.use('seaborn-paper')

    # a new matplotlib figure

variances = []
averages = []
times = []

universe = mda.Universe(args.top, args.traj)
# divide universe into groups
solvent = universe.select_atoms("resname OSMOT and name O")
nonsolvent = universe.select_atoms("resname FRAME") + universe.select_atoms("resname MOBIL")

for snapshot in universe.trajectory:
    fig = plt.figure()

    '''
    getting all coordinates of residues in the frame guided cluster
    '''
    all_com = [ res.atoms.center_of_geometry() for res in nonsolvent.residues ]
    db = DBSCAN(min_samples=1,eps=15).fit(all_com)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print('Estimated number of clusters: %d' % n_clusters_)
    print(labels)

    frame_guided_cluster = [ ]
    for i in range(len(all_com)):
        if labels[i] == 0:
            frame_guided_cluster.append(all_com[i])
    '''
    getting all coordinates of residues in the frame guided cluster
    '''
    
    # center of micelle from all centers of mass from residues
    center_of_micelle = np.mean(frame_guided_cluster, axis=0)
    print(center_of_micelle)

    #  list of distances from residue to center of micelle
    dist_hist = np.divide([ np.linalg.norm(c-center_of_micelle) for c in frame_guided_cluster ],10)

    print("average", np.average(dist_hist), "  variance", np.var(dist_hist))

    # make a histogram from distances
    # bins=auto -> no need to pass bin array
    # density=True -> norm to 1
    n, bins, patches = plt.hist(dist_hist, bins='auto', density=1)

    # best fit of data
    (mu, sigma) = norm.fit(dist_hist)
    y_fit = mpl.mlab.normpdf(bins, mu, sigma)
    plt.plot(bins, y_fit)

    # label axis
    plt.xlabel('r')
    plt.ylabel('probability')

    #limit x axis to covered area, y will be adjusted automatically
    plt.xlim(8,11)

    #output
    fig.tight_layout()
    plt.plot(rasterized=True)
    plt.savefig("histogram"+str(int(snapshot.time))+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)
    plt.close()

    times.append(snapshot.time)
    averages.append(np.average(dist_hist))
    variances.append(np.var(dist_hist)/2)



fig = plt.figure()
plt.errorbar(times,averages,variances)
plt.savefig("center_rdf"+'.'+'png'.format(), bbox_inches='tight', format='png', dpi=600)
# plt.show()