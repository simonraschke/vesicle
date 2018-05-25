#!/usr/bin/python3

import sys
import os
import argparse
import numpy as np
import MDAnalysis as mda 
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt
import scipy
import vtk
import subprocess
import math
import shutil
import h5py
import pprint
import itertools
import time
import multiprocessing
# import tbb

from sklearn.cluster import DBSCAN

parser = argparse.ArgumentParser()
parser.add_argument("--top", type=str, default="trajectory.gro", help="path to topology file")
parser.add_argument("--traj", type=str, default=None, help="path to trajectory file")
parser.add_argument("--solvent", type=str, nargs='*', default=["resname OSMOT"], help="solvent selection rule")
parser.add_argument("--nonsolvent", type=str, nargs='*', default=["resname MOBIL","resname FRAME"], help="nonsolvent selection rule")
parser.add_argument("--clstr_eps", type=float, default=15, help="max distance for cluster algorithm")
parser.add_argument("--eps", type=float, default=1, help="LJ epsilon")
parser.add_argument("--sigma", type=float, default=1, help="LJ sigma")
parser.add_argument("--kappa", type=float, default=1, help="ALJ kappa")
parser.add_argument("--gamma", type=float, help="ALJ gamma")
parser.add_argument("--render", action='store_true', help="enable rendering of pictures")
parser.add_argument("--render_dir", type=str, default="render", help="path to rendered pictures directory")
parser.add_argument("--rdf", action='store_true', help="enable rdf calculation")
args = parser.parse_args()

topology = args.top
trajectory = args.traj
DATAFILE = h5py.File("data.h5",'a')
if "clusters" in DATAFILE:
    del DATAFILE["clusters"]
CLUSTERGROUP = DATAFILE.create_group("clusters")

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

# prepare for rendering
if args.render:
    if os.path.isdir(args.render_dir):
        print("found render directory", args.render_dir, "...removing")
        shutil.rmtree(args.render_dir)
    else:
        print("render directory not found", args.render_dir)
    print("make dir", args.render_dir)
    os.mkdir(args.render_dir)

np.set_printoptions(precision=2)

# make the universe
universe = mda.Universe(topology, trajectory)
print("\n")
print("loaded universe", universe)
#split into residue groups
solvent = sum(universe.select_atoms(x) for x in args.solvent)
nonsolvent = sum(universe.select_atoms(x) for x in args.nonsolvent)
print("got", len(solvent.residues), "solvent residues and", len(nonsolvent.residues), "nonsolvent residues")


'''
The function to calc a distance with periodic boundary conditions
'''
def PBC_distance(a,b):
    return np.sqrt(np.sum(np.square( (b-a) - universe.dimensions[:3] * np.round((b-a)/universe.dimensions[:3]) )))


'''
The function to calc the potential energie
'''
a = args.sigma + args.kappa*np.sin(args.gamma*math.pi/180)
b = args.sigma - args.kappa*np.sin(args.gamma*math.pi/180)
c = np.linalg.norm(np.array([a,0,0]) + np.array([a-args.sigma, args.kappa*np.cos(args.gamma*math.pi/180),0]))
print("a",a)
print("b",b)
print("c",c)


def __detail_Epot(c1, c2, o1, o2, dimensions):
    distance_vec = np.subtract(c2, c1) 
    distance_vec_pbc = np.subtract(distance_vec, np.multiply(dimensions, np.round(distance_vec/dimensions)))
    distance_vec_pbc = np.divide(distance_vec_pbc, 10)
    norm = np.linalg.norm(distance_vec_pbc)
    if norm < 3*args.sigma:
        r2 = args.sigma*args.sigma/norm**2
        r6 = r2**3
        normed = np.divide(distance_vec_pbc,norm)
        res1_u = ( o1 - c1 )
        res1_u = res1_u/np.linalg.norm(res1_u)*(args.kappa/2)
        res2_u = ( o2 - c2 )
        res2_u = res2_u/np.linalg.norm(res2_u)*(args.kappa/2)
        chi = (np.linalg.norm(-res1_u+normed+res2_u) - a)**2 + (np.linalg.norm(res1_u+normed-res2_u) - b)**2 + (np.linalg.norm(-res1_u+normed-res2_u) - c)**2 + (np.linalg.norm(res1_u+normed+res2_u) - c)**2 
        return 4.0*args.eps*(r6*r6-(1.0-chi)*r6)
    else:
        return 0.0


def Epot(residues):
    epot = 0
    dims = universe.dimensions[:3]
    # results = []
    for i in range(len(residues)):
        for j in range(i):
            if i==j: continue
            # if PBC_distance(residues[i].atoms.center_of_geometry(), residues[h].atoms.center_of_geometry()) < 30: continue
            epot += __detail_Epot(residues[i].atoms.center_of_geometry(), residues[j].atoms.center_of_geometry(), residues[i].atoms[0].position, residues[j].atoms[0].position, dims)
            # results.append(pool.apply_async(__detail_Epot, (residues[i].atoms.center_of_geometry(), residues[j].atoms.center_of_geometry(), residues[i].atoms[0].position, residues[j].atoms[0].position, dims,)))
    return epot
    # return sum( [r.get() for r in results] )


'''
The function to calc the order
'''
def Order(residues):
    if len(residues) <= 1: return 0.0
    center = np.mean([ res.atoms.center_of_geometry() for res in residues ], axis=0)
    connections = np.array([ res.atoms.center_of_geometry()-center for res in residues ])
    connections = np.array([ con/np.linalg.norm(con) for con in connections ])
    orientations = np.array([ np.divide(np.subtract(res.atoms[0].position, res.atoms.center_of_geometry()), np.linalg.norm(np.subtract(res.atoms[0].position, res.atoms.center_of_geometry()))) for res in residues ])
    order_values = np.array([np.dot(c,o) for c,o in zip(connections, orientations)])
    return np.mean(order_values)



'''
The function to scan for subclusters an return labels
'''
def SubclusterDBSCAN(cogs):
    if len(cogs) == 1: return np.array([0])
    dbscan_internal = DBSCAN(min_samples=1,eps=args.clstr_eps).fit(cogs)
    core_samples_mask = np.zeros_like(dbscan_internal.labels_, dtype=bool)
    core_samples_mask[dbscan_internal.core_sample_indices_] = True
    return np.array(dbscan_internal.labels_)


'''
The function to calc rdf
'''
def RDF(residues, start, stop, skip=1):
    nbins = 1000
    dmin, dmax = 0.0, 100.0
    rdf, edges = np.histogram([0], bins=nbins, range=(dmin, dmax))
    rdf *= 0
    rdf = rdf.astype(np.float64)
    try:
        n = residues.n_residues
    except:
        residues = mda.core.groups.AtomGroup(sum([res.atoms for res in residues])).residues
        n = residues.n_residues
    dist = np.zeros(int(n * (n - 1) / 2,), dtype=np.float64)
    boxvolume = 0
    for ts in universe.trajectory[start:stop:skip]:
        boxvolume += ts.volume
        coor = np.array([ res.atoms.center_of_geometry() for res in residues.residues ])
        box = ts.dimensions[:3]
        mda.lib.distances.self_distance_array(reference=coor, box=box, result=dist, backend="OpenMP")
        new_rdf, edges = np.histogram(dist, bins=nbins, range=(dmin, dmax))
        rdf += new_rdf
    numframes = (stop-start) / skip
    boxvolume /= numframes # average volume
    # Normalize RDF
    radii = 0.5 * (edges[1:] + edges[:-1])
    vol = (4. / 3.) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    # normalization to the average density n/boxvolume in the simulation
    density = n / boxvolume
    norm = density * (n - 1) / 2 * numframes
    rdf /= norm * vol
    return radii, rdf

max_time = universe.trajectory[-1].time
max_frame = universe.trajectory[-1].frame
universe.trajectory[0]

if args.rdf:
    bins, rdf = RDF(nonsolvent.residues, int(max_frame/2), max_frame)
    bins = np.divide(bins,10)
    fig = plt.figure()
    plt.plot([min(bins), max(bins)], [1,1], color='black')
    plt.plot(bins, rdf)
    plt.show()

pool = multiprocessing.Pool(maxtasksperchild=1000)

for snapshot in universe.trajectory[-5:]:
    print("\n",snapshot)
    t_start = time.time()
    dimensions = universe.dimensions[:3]

    # first general DBSCAN
    all_com = [ res.atoms.center_of_geometry() for res in nonsolvent.residues ]
    dbscan = DBSCAN(min_samples=1,eps=args.clstr_eps, metric=PBC_distance, n_jobs=-1, algorithm='brute').fit(all_com)
    t_dbscan = time.time()-t_start
    print("DBSCAN took", t_dbscan, "seconds")
    core_samples_mask = np.zeros_like(dbscan.labels_, dtype=bool)
    core_samples_mask[dbscan.core_sample_indices_] = True
    labels = dbscan.labels_
    num_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print("found", num_clusters, "clusters")
    clusters = {}
    for ID in np.unique(labels):
        clusters.update({ID : [nonsolvent.residues[i] for i in range(len(nonsolvent.residues)) if labels[i] == ID]})

    volume = 0
    surface_area = 0
    epot = 0

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0,1.0,1.0)

    DSET_clusters = CLUSTERGROUP.create_dataset("time"+str(int(snapshot.time)), (len(clusters.keys()),len(clusters[max(clusters, key= lambda x: len(set(clusters[x])))])+5), dtype=np.float16, fillvalue=-1)
    
    t_basic = time.time()-t_start
    print("basic analysis took", t_basic, "seconds")

    results = {}
    for ID, cluster in clusters.items():
        if len(cluster) > 1:
            internal_com = [ res.atoms.center_of_geometry() for res in cluster ]
            results[ID] = (pool.apply_async(SubclusterDBSCAN,(internal_com,)))
    
    for ID, cluster in clusters.items():
        #internal DBSCAN per cluster
        if len(cluster) > 1:
            labels_internal = results[ID].get()
            num_clusters_internal = len(set(labels_internal)) - (1 if -1 in labels_internal else 0)
            subclusters = {}
            for ID_internal in np.unique(labels_internal):
                subclusters.update({ID_internal : [cluster[i] for i in range(len(cluster)) if labels_internal[i] == ID_internal]})

            largest_subcluster_ID = max(subclusters, key= lambda x: len(set(subclusters[x])))
            cluster_centers = [ np.mean([ res.atoms.center_of_geometry() for res in subcluster], axis=0) for ID_internal,subcluster in subclusters.items() ]
            np.set_printoptions(precision=1)
            center_shifts = [ np.round((cluster_centers[largest_subcluster_ID]-center)/dimensions).astype(np.int) for center in cluster_centers ]

            points = vtk.vtkPoints()
            for ID_internal, subcluster in subclusters.items():
                for residue in subcluster:
                    norm = np.linalg.norm(np.subtract(residue.atoms[0].position, residue.atoms.center_of_geometry()))
                    com = np.add(np.add(residue.atoms.center_of_geometry(), np.multiply(center_shifts[ID_internal], dimensions)), np.divide(np.subtract(residue.atoms[0].position, residue.atoms.center_of_geometry()), norm/10))
                    points.InsertNextPoint(*com)

            polyData = vtk.vtkPolyData()
            polyData.SetPoints(points)

            surface = vtk.vtkDataSetSurfaceFilter()
            properties = vtk.vtkMassProperties()
            if len(cluster) >= 6: 
                delaunay = vtk.vtkDelaunay3D()
                delaunay.SetInputData(polyData)
                delaunay.Update()
                surface.SetInputConnection(delaunay.GetOutputPort())
                properties.SetInputConnection(surface.GetOutputPort())
                volume += properties.GetVolume()
                surface_area += properties.GetSurfaceArea()
            
            if args.render:
                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInputData(surface.GetOutput())
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.GetProperty().SetColor(min([max([0.2,len(cluster)/100]), 1]), min([max([0.2,len(cluster)/100])/2, 0.5]), 0.2)
                renderer.AddActor(actor)

            epot_temp = Epot(cluster)
            # epot_temp = -len(cluster)
            epot += epot_temp


        DSET_clusters[ID,0] = len(cluster)
        DSET_clusters[ID,1] = Order(cluster)
        if len(cluster) >= 6: 
            DSET_clusters[ID,2] = properties.GetVolume()/10**3
            DSET_clusters[ID,3] = properties.GetSurfaceArea()/10**2
        else:
            DSET_clusters[ID,2] = len(cluster) * 4/3*math.pi*(1.1224/2)**3
            DSET_clusters[ID,3] = len(cluster) * 4*math.pi*(1.1224/2)**2
        DSET_clusters[ID,4] = epot
        for i,res in enumerate(cluster):
            DSET_clusters[ID,i+5] = res.ix

        
        # print("single cluster analysis took", time.time()-t_start_internal, "seconds")

    print("volume", volume)
    print("surface area", surface_area)
    print("epot", epot)

    if args.render:
        # a badly programmed box actor
        boxPoints = vtk.vtkPoints()
        boxPoints.InsertNextPoint(0,0,0)
        boxPoints.InsertNextPoint(universe.dimensions[0],0,0)
        boxPoints.InsertNextPoint(universe.dimensions[0],universe.dimensions[1],0)
        boxPoints.InsertNextPoint(0,universe.dimensions[1],0)
        boxPoints.InsertNextPoint(0,0,0)
        boxPoints.InsertNextPoint(0,0,universe.dimensions[2])
        boxPoints.InsertNextPoint(universe.dimensions[0],0,universe.dimensions[2])
        boxPoints.InsertNextPoint(universe.dimensions[0],0,0)
        boxPoints.InsertNextPoint(universe.dimensions[0],0,universe.dimensions[2])
        boxPoints.InsertNextPoint(universe.dimensions[0],universe.dimensions[1],universe.dimensions[2])
        boxPoints.InsertNextPoint(universe.dimensions[0],universe.dimensions[1],0)
        boxPoints.InsertNextPoint(universe.dimensions[0],universe.dimensions[1],universe.dimensions[2])
        boxPoints.InsertNextPoint(0,universe.dimensions[1],universe.dimensions[2])
        boxPoints.InsertNextPoint(0,universe.dimensions[1],0)
        boxPoints.InsertNextPoint(0,universe.dimensions[1],universe.dimensions[2])
        boxPoints.InsertNextPoint(0,0,universe.dimensions[2])
        boxPolyLine = vtk.vtkPolyLine()
        boxPolyLine.GetPointIds().SetNumberOfIds(16)
        for i in range(16):
            boxPolyLine.GetPointIds().SetId(i,i)
        boxCells = vtk.vtkCellArray()
        boxCells.InsertNextCell(boxPolyLine)
        boxPolyData = vtk.vtkPolyData()
        boxPolyData.SetPoints(boxPoints)
        boxPolyData.SetLines(boxCells)
        boxMapper = vtk.vtkPolyDataMapper()
        boxMapper.SetInputData(boxPolyData)
        boxActor = vtk.vtkActor()
        boxActor.SetMapper(boxMapper)
        boxActor.GetProperty().SetColor(0,0,0)
        renderer.AddActor(boxActor)

        camera = vtk.vtkCamera()
        camera.SetPosition(dimensions[0]/2,dimensions[1]/2,dimensions[2] + dimensions[2]*2)
        camera.SetFocalPoint(dimensions[0]/2,dimensions[1]/2,dimensions[2]/2)
        renderer.SetActiveCamera(camera)

        # the render window and actual render
        renderer_window = vtk.vtkRenderWindow()
        renderer_window.SetWindowName("Cluster")
        renderer_window.SetSize(2000,2000)
        renderer_window.AddRenderer(renderer)
        renderer_window.SetOffScreenRendering(1)
        renderer_window.Render()

        win2image = vtk.vtkWindowToImageFilter()
        win2image.SetInput(renderer_window)
        win2image.Update()

        writer = vtk.vtkJPEGWriter()
        writer.SetFileName(os.path.join(args.render_dir, str(int(snapshot.time)).zfill(len(str(int(max_time))))+'.jpg'))
        writer.SetInputConnection(win2image.GetOutputPort())
        writer.Write()
        renderer.Clear()

    print("frame took", time.time()-t_start, "seconds")

DATAFILE.close()