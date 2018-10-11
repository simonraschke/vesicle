#!/usr/bin/python3

import sys
import os
import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt
import h5py
import pprint
import shutil
import subprocess
import time
import analysis_helper_functions as helper

from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import normalize
from MDAnalysis.lib.distances import distance_array

parser = argparse.ArgumentParser()
parser.add_argument("--top", type=str, default="trajectory.gro", help="path to topology file")
parser.add_argument("--traj", type=str, default=None, help="path to trajectory file")
parser.add_argument("--config", type=str, default="config_analysis.ini", help="path to config file")
parser.add_argument("--solvent", type=str, nargs='*', default=["resname OSMOT"], help="solvent selection rule")
parser.add_argument("--nonsolvent", type=str, nargs='*', default=["resname MOBIL","resname FRAME"], help="nonsolvent selection rule")
parser.add_argument("--clstr_eps", type=float, default=12, help="max distance for cluster algorithm")
parser.add_argument("--start", type=float, default=-1, help="starting time of analysis")
parser.add_argument("--stop", type=float, default=10e20, help="starting time of analysis")
parser.add_argument("--forcenew", action='store_true', help="force new hdf5 file")
parser.add_argument("--lowmem", action='store_true', help="dont save resname, saves memory BIG TIME")
args = parser.parse_args()


pp = pprint.PrettyPrinter(indent=4, compact=False)
np.set_printoptions(suppress=True)


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
        pp.pprint(subprocess.getstatusoutput(cmd))
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
try:
    # datafile = pd.HDFStore("data.h5", "a")
    datafile = pd.HDFStore('data.h5', 'a', complevel=9, complib='zlib')
except Exception as e:
    print(e)
    sys.exit()



# saving attributes from args.config file
attributes = pd.DataFrame(helper.getAttributeDict(args.config, universe.trajectory[0].dimensions[:3]))
print(attributes)
datafile["attributes"] = attributes

epot = helper.EpotCalculator(attributes)


t_start = time.perf_counter()
for snapshot in universe.trajectory:
    # print("\n",snapshot)
    if snapshot.time < args.start:
        continue
    elif snapshot.time > args.stop:
        continue
    dimensions = universe.dimensions

    t_prep = time.perf_counter()
    #get all positions
    df = pd.DataFrame(snapshot.positions, columns=['x','y','z'])
    # calculate center_of_geometries
    coms = df.groupby(np.arange(len(df))//2).mean()
    # calculate orientations (upper atom minus com)
    orientations = pd.DataFrame(normalize(df[0::2].reset_index(drop=True).sub(coms)), columns=['ux','uy','uz'])
    # concatenate both
    particledata = pd.concat([coms,orientations], axis=1).reset_index()

    # add particle (residue) names
    particledata["resname"] = universe.atoms.residues.resnames
    particledata["resid"] = universe.atoms.residues.resids

    # scan for clusters
    distances_array = distance_array(coms.values, coms.values, box=dimensions)
    dbscan = DBSCAN(min_samples=2, eps=args.clstr_eps, metric="precomputed", n_jobs=-1).fit(distances_array)
    labels = pd.DataFrame(dbscan.labels_, columns=['cluster'])

    # add to data and sort for cluster id
    particledata["cluster"] = labels
    # particledata.sort_values('cluster', inplace=True)
    unique, counts = np.unique(labels, return_counts=True)
    particledata["clustersize"] = particledata["cluster"].apply( lambda x: counts[np.where(unique == x)][0] )
    particledata.loc[particledata["cluster"] == -1, "clustersize"] = 1

    # create cluster dataframe
    # clusterdata = particledata.groupby(["cluster"]).size().reset_index(name='particles').drop('cluster', axis=1)
    # print(f"prep took     {time.perf_counter()-t_prep:.4f} seconds")

    t_sub = time.perf_counter()
    # subcluster identification
    subcluster_labels = []
    for ID, group in particledata.groupby(["cluster"], as_index=False):
        subcluster_labels.extend(helper.getSubclusterLabels(ID, group, args.clstr_eps))
        
    # add the subcluster IDs
    particledata["subcluster"] = subcluster_labels
    # print(f"subclstr took {time.perf_counter()-t_sub:.4f} seconds")

    t_shift = time.perf_counter()
    particledata["shiftx"] = particledata["x"]
    particledata["shifty"] = particledata["y"]
    particledata["shiftz"] = particledata["z"]
    # shift subclusters towards largest subcluster
    for ID, group in particledata.groupby("cluster"):
        newx, newy, newz = helper.getShiftedCoordinates(ID, group, args.clstr_eps, dimensions[:3])
        particledata.loc[newx.index, "shiftx"] = newx.values
        particledata.loc[newy.index, "shifty"] = newy.values
        particledata.loc[newz.index, "shiftz"] = newz.values
    # print(f"shift took    {time.perf_counter()-t_shift:.4f} seconds")

    # t_order = time.perf_counter()
    # get the order of particle in cluster
    particledata["order"] = 0.0
    for ID, group in particledata.groupby("cluster"):
        orders = helper.getOrder(ID, group)
        particledata.loc[group.index, "order"] = orders
    # print(particledata.groupby("clustersize")["order"].mean())
    # print(f"shift took    {time.perf_counter()-t_order:.4f} seconds")

    t_volume = time.perf_counter()
    particledata["volume"] = 0.0
    # get the volume per cluster
    for ID, group in particledata.groupby(["cluster"]):
        volume = helper.getClusterVolume(ID, group, args.clstr_eps, 4)
        particledata.loc[group.index, "volume"] = volume
        if volume / np.cumprod(dimensions[:3])[-1] > 0.5:
            raise Exception(f"volume of cluster {ID} is {volume / np.cumprod(dimensions[:3])[-1]} of box volume")

    # calculate the potential energy per particle
    particledata["epot"] = epot.get(coms, orientations, dimensions, particledata)

    # plt.show()
    # plt.savefig("cluster.png", dpi=600)
    # print(f"volume took   {time.perf_counter()-t_volume:.4f} seconds")

    if args.lowmem:
        particledata = particledata.drop(columns=["x","y","z","ux","uy","uz","shiftx","shifty","shiftz","resname"])
    else:
        particledata = particledata.drop(columns=["resname"])

    t_write = time.perf_counter()
    datafile[f"time{int(snapshot.time)}"] = particledata
    print(datafile[f"time{int(snapshot.time)}"])
    print(particledata.groupby("clustersize")["epot","volume"].sum())
    
    t_end = time.perf_counter()
    print(f"time {snapshot.time} took {t_end-t_start:.4f} seconds")
    t_start = time.perf_counter()

datafile.close()