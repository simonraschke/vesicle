#!/usr/bin/python3

import sys
import os
import argparse
import numpy as np
import MDAnalysis as mda
import scipy
import vtk
import subprocess
import math
import shutil
import h5py
import pprint

from sklearn.cluster import DBSCAN

parser = argparse.ArgumentParser()
parser.add_argument("--top", type=str, default="trajectory.gro", help="path to topology file")
parser.add_argument("--traj", type=str, default=None, help="path to trajectory file")
parser.add_argument("--solvent", type=str, nargs='*', default=["resname OSMOT"], help="solvent selection rule")
parser.add_argument("--nonsolvent", type=str, nargs='*', default=["resname MOBIL","resname FRAME"], help="nonsolvent selection rule")
parser.add_argument("--eps", type=float, default=15, help="max distance for cluster algorithm")
parser.add_argument("--render", action='store_true', help="enable rendering of pictures")
parser.add_argument("--render_dir", type=str, default="render", help="path to rendered pictures directory")
args = parser.parse_args()

topology = args.top
trajectory = args.traj
DATAFILE = h5py.File("data.h5",'w')
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

max_time = universe.trajectory[-1].time
universe.trajectory[0]

for snapshot in universe.trajectory[-5:]:
    dimensions = universe.dimensions[:3]

    # first general DBSCAN
    all_com = [ res.atoms.center_of_geometry() for res in nonsolvent.residues ]
    dbscan = DBSCAN(min_samples=1,eps=args.eps, metric=PBC_distance, n_jobs=-1).fit(all_com)
    core_samples_mask = np.zeros_like(dbscan.labels_, dtype=bool)
    core_samples_mask[dbscan.core_sample_indices_] = True
    labels = dbscan.labels_
    num_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print("\n",snapshot)
    print("found", num_clusters, "clusters")
    clusters = {}
    for ID in np.unique(labels):
        clusters.update({ID : [nonsolvent.residues[i] for i in range(len(nonsolvent.residues)) if labels[i] == ID]})
    
    volume = 0
    surface_area = 0

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0,1.0,1.0)

    DSET_clusters = CLUSTERGROUP.create_dataset("time"+str(int(snapshot.time)), (len(clusters.keys()),len(clusters[max(clusters, key= lambda x: len(set(clusters[x])))])+5), dtype=np.float16, fillvalue=-1)

    for ID, cluster in clusters.items():
        #internal DBSCAN per cluster
        internal_com = [ res.atoms.center_of_geometry() for res in cluster ]
        dbscan_internal = DBSCAN(min_samples=1,eps=args.eps).fit(internal_com)
        core_samples_mask = np.zeros_like(dbscan_internal.labels_, dtype=bool)
        core_samples_mask[dbscan_internal.core_sample_indices_] = True
        labels_internal = dbscan_internal.labels_
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
        
        print(ID)
        DSET_clusters[ID,0] = len(cluster)
        DSET_clusters[ID,1] = 1.5
        if len(cluster) >= 6: 
            DSET_clusters[ID,2] = properties.GetVolume()/10**3
            DSET_clusters[ID,3] = properties.GetSurfaceArea()/10**2
        else:
            DSET_clusters[ID,2] = len(cluster) * 4/3*math.pi*(1.1224/2)**3
            DSET_clusters[ID,3] = len(cluster) * 4*math.pi*(1.1224/2)**2
        for i,res in enumerate(cluster):
            DSET_clusters[ID,i+5] = res.ix

    print("volume", volume)
    print("surface area", surface_area)

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
