#!/usr/bin/python3

import sys
import os
import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda 
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.lib.distances import distance_array
import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt
import scipy
import subprocess
import math
import shutil
import h5py
import pprint
import itertools
import time
import multiprocessing
import analysis_helper_functions as HELP

from sklearn.cluster import DBSCAN
from sklearn.preprocessing import normalize

parser = argparse.ArgumentParser()
parser.add_argument("--top", type=str, default="trajectory.gro", help="path to topology file")
parser.add_argument("--traj", type=str, default=None, help="path to trajectory file")
parser.add_argument("--config", type=str, default="config.ini", help="path to config file")
parser.add_argument("--solvent", type=str, nargs='*', default=["resname OSMOT"], help="solvent selection rule")
parser.add_argument("--nonsolvent", type=str, nargs='*', default=["resname MOBIL","resname FRAME"], help="nonsolvent selection rule")
parser.add_argument("--clstr_eps", type=float, default=15, help="max distance for cluster algorithm")
parser.add_argument("--forcenew", action='store_true', help="force new hdf5 file")
args = parser.parse_args()



topology = args.top
trajectory = args.traj



# prepare input files
if os.path.exists(topology):
    print("found topology:", topology)
else:
    raise Exception("unable to find topology "+str(topology))

if trajectory != None:
    if os.path.exists(trajectory):
        print("found trajectory:", trajectory)
    else:
        raise Exception("unable to find trajectory "+str(trajectory))
else:
    print("try to convert topology to .xtc trajectory")
    if shutil.which("gmx") != None:
        print("found gmx")
        trajectory = "trajectory.xtc"
        cmd = "gmx trjconv -f "+topology+" -o "+trajectory
        pprint.pprint(subprocess.getstatusoutput(cmd))
    else:
        raise Exception("gmx not found")



# loading the universe
universe = mda.Universe(topology, trajectory)
print("\n\nloaded universe", universe)
#split into residue groups
solvent = sum(universe.select_atoms(x) for x in args.solvent)
nonsolvent = sum(universe.select_atoms(x) for x in args.nonsolvent)
print("got", len(solvent.residues), "solvent residues and", len(nonsolvent.residues), "nonsolvent residues")



# creating the storage file object
if args.forcenew and os.path.exists("data.h5"):
    os.remove("data.h5")
datafile = pd.HDFStore("data.h5")

for snapshot in universe.trajectory[-5:]:
    t_start = time.perf_counter()
    print("\n",snapshot)
    dimensions = universe.dimensions[:3]

    t_prep = time.perf_counter()
    #get all positions
    df = pd.DataFrame(snapshot.positions, columns=['x','y','z'])
    # calculate center_of_geometries
    coms = df.groupby(np.arange(len(df))//2).mean()
    # calculate orientations (upper atom minus com)
    orientations = pd.DataFrame(normalize(df[0::2].reset_index(drop=True).sub(coms)), columns=['ux','uy','uz'])
    # concatenate both
    particledata = pd.concat([coms,orientations], axis=1).reset_index()

    # scan for clusters
    distances_array = distance_array(coms.values, coms.values, box=dimensions)
    dbscan = DBSCAN(min_samples=2, eps=14, metric="precomputed", n_jobs=-1).fit(distances_array)
    labels = pd.DataFrame(dbscan.labels_, columns=['cluster'])

    # add to data and sort for cluster id
    particledata["cluster"] = labels
    particledata.sort_values('cluster', inplace=True)

    # create cluster dataframe
    clusterdata = particledata.groupby(["cluster"]).size().reset_index(name='particles').drop('cluster', axis=1)
    print(f"prep took     {time.perf_counter()-t_prep:.4f} seconds")

    t_sub = time.perf_counter()
    # subcluster identification
    subcluster_labels = []
    for ID, group in particledata.groupby(["cluster"], as_index=False):
        #skip if noise
        if ID == -1:
            subcluster_labels.extend(np.zeros( len(group.index), dtype=int ) )
            continue
        # arange a DBSCAN without PBC to get subclusters
        coms_subcluster = pd.concat([group['x'], group['y'], group['z']], axis=1)
        distances_array_subcluster = distance_array(coms_subcluster.values, coms_subcluster.values, box=None)
        dbscan_subcluster = DBSCAN(min_samples=1, eps=14, metric="precomputed", n_jobs=-1).fit(distances_array_subcluster)
        subcluster_labels.extend(dbscan_subcluster.labels_)
    # add the subcluster IDs
    particledata["subcluster"] = subcluster_labels
    print(f"subclstr took {time.perf_counter()-t_sub:.4f} seconds")

    t_shift = time.perf_counter()
    particledata["shiftx"] = particledata["x"]
    particledata["shifty"] = particledata["y"]
    particledata["shiftz"] = particledata["z"]
    # shift subclusters towards largest subcluster
    for clusterID, group in particledata.groupby("cluster"):
        # get largest subcluster
        unique, counts = np.unique(group["subcluster"], return_counts=True)
        if len(unique) == 1:
            continue
        max_subclusterID = unique[counts == np.max(counts)][0]
        # calculate shifts per subcluster
        centers = group.groupby("subcluster")['x','y','z'].mean()
        shifts = np.round(( -centers + centers.loc[max_subclusterID] )/dimensions).astype(int)
        shifts *= dimensions
        # calculate new coordinates based on shift
        newx = np.add(group["x"], shifts.loc[group["subcluster"]]["x"])
        newy = np.add(group["y"], shifts.loc[group["subcluster"]]["y"])
        newz = np.add(group["z"], shifts.loc[group["subcluster"]]["z"])
        # assign to main data
        particledata.loc[newx.index, "shiftx"] = newx.values
        particledata.loc[newy.index, "shifty"] = newy.values
        particledata.loc[newz.index, "shiftz"] = newz.values
    print(f"shift took    {time.perf_counter()-t_shift:.4f} seconds")


    t_volume = time.perf_counter()
    particledata["volume"] = 0.0
    for clusterID, group in particledata.groupby(["cluster"]):
        # print(clusterID)
        # print(subclusterID)
        x_vector = np.arange(np.min(group["shiftx"])-14, np.max(group["shiftx"])+14, 3.3, dtype=np.float32)
        y_vector = np.arange(np.min(group["shifty"])-14, np.max(group["shifty"])+14, 3.3, dtype=np.float32)
        z_vector = np.arange(np.min(group["shiftz"])-14, np.max(group["shiftz"])+14, 3.3, dtype=np.float32)
        xx,yy,zz = np.meshgrid(x_vector, y_vector, z_vector)
        # print(len(x_vector))
        # print(len(y_vector))
        # print(len(z_vector))
        meshgrid = np.stack((xx.ravel(), yy.ravel(), zz.ravel()), axis=1)
        flags = np.zeros_like((len(meshgrid)), dtype=bool)
        coms_cluster = pd.concat([group['shiftx'], group['shifty'], group['shiftz']], axis=1)
        # print()
        # print(meshgrid.shape)
        # print(coms_cluster.values.shape)
        # distances_array_volume = distance_array(meshgrid, coms_cluster.values, box=None, backend="OpenMP")
        # isclose = np.where(distances_array_volume < 14, True, False)
        # print(isclose)
        if len(group.index) > 3:
            print(f"volume {scipy.spatial.ConvexHull(coms_cluster.values).volume:.4f}")

        # print()
        # print()
        # print()
        # print()
    print(f"volume took   {time.perf_counter()-t_volume:.4f} seconds")

    t_write = time.perf_counter()
    print(particledata)
    datafile["time"+str(int(snapshot.time))+"/particles"] = particledata
    datafile["time"+str(int(snapshot.time))+"/cluster"] = clusterdata
    print(f"write took    {time.perf_counter()-t_write:.4f} seconds")
    
    t_end = time.perf_counter()
    print(f"complete took {t_end-t_start:.4f} seconds")