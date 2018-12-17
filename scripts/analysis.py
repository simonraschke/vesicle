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
import sklearn
import time
import analysis_helper_functions as helper

from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from MDAnalysis.lib.distances import distance_array
from collections import defaultdict

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
parser.add_argument("--hdbscan", action='store_true', help="additionaly perform alternative clustering by HDBSCAN")
parser.add_argument("--timestats", action='store_true', help="show timer statistics")
parser.add_argument("--reanalyze", action='store_true', help="reanalze from data.h5 file instead of trajectory")
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
print()
for column in attributes:
    print(f"{column:<25} {attributes[column].values[0]}")
print()
datafile["attributes"] = attributes

epot_calc = helper.EpotCalculator(attributes)
resname_map = defaultdict(lambda : -1, {
    "MOBIL" : 0,
    "FRAME" : 1,
    "OSMOT" : 2
})

complete_data = pd.DataFrame()

t_start = time.perf_counter()
for snapshot in universe.trajectory:
    if snapshot.time < args.start:
        continue
    elif snapshot.time > args.stop:
        continue
    dimensions = universe.dimensions

    t_prep = time.perf_counter()
    #get all positions
    df = pd.DataFrame(universe.atoms.positions, columns=['x','y','z'])
    # calculate center_of_geometries
    coms = df.groupby(np.arange(len(df))//2).mean()
    # calculate orientations (upper atom minus com)
    orientations = pd.DataFrame(normalize(df[0::2].reset_index(drop=True).sub(coms)), columns=['ux','uy','uz'])
    # concatenate both
    particledata = pd.concat([coms,orientations], axis=1).reset_index()

    """
    add particle (residue) names
    """
    particledata["resname"] = pd.Series(universe.atoms.residues.resnames).map(resname_map).astype(np.int16)
    particledata["resid"] = universe.atoms.residues.resids

    relevant_cond = particledata["resname"] != 2
    coms = particledata[relevant_cond].filter(["x","y","z"])
    orientations = particledata[relevant_cond].filter(["ux","uy","uz"])
    particledata["ux"] = np.where(relevant_cond, particledata["ux"], np.nan)
    particledata["uy"] = np.where(relevant_cond, particledata["uy"], np.nan)
    particledata["uz"] = np.where(relevant_cond, particledata["uz"], np.nan)


    """
    scan for clusters
    """
    distances_array = distance_array(coms.values, coms.values, box=dimensions)
    dbscan = DBSCAN(min_samples=2, eps=args.clstr_eps, metric="precomputed", n_jobs=-1).fit(distances_array)
    labels = pd.DataFrame(np.append(dbscan.labels_, np.full(np.count_nonzero(~relevant_cond), -2)), columns=['cluster'])
    # add to data and sort for cluster id
    particledata["cluster"] = labels

    unique, counts = np.unique(labels, return_counts=True)
    particledata["clustersize"] = particledata["cluster"].apply( lambda x: counts[np.where(unique == x)][0] if x >= 0 else 1 )
    particledata.loc[particledata["cluster"] == -1, "clustersize"] = 1

    if args.timestats: print(f"prep took     {time.perf_counter()-t_prep:.4f} seconds")



    """
    subcluster identification
    """
    t_sub = time.perf_counter()
    particledata["subcluster"] = -1
    for ID, group in particledata[relevant_cond].groupby(["cluster"], as_index=False):
        subclusters = helper.getSubclusterLabels(ID, group, args.clstr_eps)
        particledata.loc[group.index, "subcluster"] = subclusters

    if args.timestats: print(f"subclstr took {time.perf_counter()-t_sub:.4f} seconds")


    """
    shift subclusters towards largest subcluster
    """
    t_shift = time.perf_counter()
    particledata["shiftx"] = np.where(relevant_cond, particledata["x"], np.nan)
    particledata["shifty"] = np.where(relevant_cond, particledata["y"], np.nan)
    particledata["shiftz"] = np.where(relevant_cond, particledata["z"], np.nan)
    for ID, group in particledata[relevant_cond].groupby("cluster"):
        newx, newy, newz = helper.getShiftedCoordinates(ID, group, args.clstr_eps, dimensions[:3])
        particledata.loc[newx.index, "shiftx"] = newx.values
        particledata.loc[newy.index, "shifty"] = newy.values
        particledata.loc[newz.index, "shiftz"] = newz.values
    
    if args.timestats: print(f"shift took    {time.perf_counter()-t_shift:.4f} seconds")



    """
    get the order of particle in cluster
    """
    t_order = time.perf_counter()
    particledata["order"] = np.nan
    for ID, group in particledata[relevant_cond].groupby("cluster"):
        orders = helper.getOrder(ID, group)
        particledata.loc[group.index, "order"] = orders
    
    if args.timestats: print(f"shift took    {time.perf_counter()-t_order:.4f} seconds")



    """
    get the volume per cluster DBSCAN
    """
    t_volume = time.perf_counter()
    particledata["volume"] = np.nan
    for ID, group in particledata[relevant_cond].groupby(["cluster"]):
        volume = helper.getClusterVolume(ID, group, args.clstr_eps, 4)
        particledata.loc[group.index, "volume"] = volume
        if volume / np.cumprod(dimensions[:3])[-1] > 0.9:
            raise Exception(f"volume of cluster {ID} is {volume / np.cumprod(dimensions[:3])[-1]} of box volume")
    
    if args.timestats: print(f"volume took   {time.perf_counter()-t_volume:.4f} seconds")


    
    """
    calculate the potential energy per particle
    """
    t_epot = time.perf_counter()
    epot, chi = epot_calc.get(coms, orientations, dimensions, ret="epot+chi")
    particledata["epot"] = np.append(epot, np.full(np.count_nonzero(~relevant_cond), np.nan))
    particledata["chi"] = np.append(chi, np.full(np.count_nonzero(~relevant_cond), np.nan))
    if args.timestats: print(f"epot and chi took     {time.perf_counter()-t_epot:.4f} seconds")



    """
    calculate the curvature of the structure for every particle
    """
    t_curvature = time.perf_counter()
    particledata["curvature"] = np.append(helper.getCurvature(particledata[relevant_cond], dimensions, cutoff=13), np.full(np.count_nonzero(~relevant_cond), np.nan))
    if args.timestats: print(f"curvature took        {time.perf_counter()-t_curvature:.4f} seconds")



    if args.lowmem:
        particledata = particledata.drop(columns=["shiftx","shifty","shiftz"])

    t_write = time.perf_counter()
    datafile[f"time{int(snapshot.time)}"] = particledata
    if args.timestats: print(f"write took    {time.perf_counter()-t_epot:.4f} seconds")
    
    t_end = time.perf_counter()
    print(f"time {snapshot.time} took {t_end-t_start:.4f} seconds")
    t_start = time.perf_counter()

    # print(datafile[f"time{int(snapshot.time)}"])

datafile.close()