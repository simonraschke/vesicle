import sys
import numpy as np
import pandas as pd
import scipy
import MDAnalysis as mda 

from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.lib.distances import distance_array
from sklearn.cluster import DBSCAN
from scipy import ndimage
from sklearn.preprocessing import normalize


# iterate a file linewise
# check line for keyword
# if keyword is found return value after seperator
def fileValueFromKeyword(filepath, keyword, seperator='='):
    found_counter = 0
    assert(filepath)
    value = [np.NaN]
    try:
        with open(filepath) as FILE:
            for line in FILE:
                if keyword in line:
                    found_counter += 1
                    value = [float(line.split(seperator,1)[1].rstrip())]
                    break
    except:
        pass
    if found_counter == 0:
        print(Warning("unable to find keyword "+keyword+" in file "+filepath+"\n"))
        # raise Warning("unable to find keyword "+keyword+" in file "+filepath+"\n")
    # if found_counter >= 2:
        # raise Exception("multiple keywords "+keyword+" in file "+filepath+"\n")
    return value



def getAttributeDict(filename, dimensions):
    return {
        'mobile': fileValueFromKeyword(filename, 'mobile', '='),
        'density': fileValueFromKeyword(filename, 'density', '='),
        'boxx': dimensions[0],
        'boxy': dimensions[1],
        'boxz': dimensions[2],
        'temperature': fileValueFromKeyword(filename, 'temperature', '='),
        'time_max': fileValueFromKeyword(filename, 'time_max', '='),
        'kappa': fileValueFromKeyword(filename, 'kappa', '='),
        'gamma': fileValueFromKeyword(filename, 'gamma', '='),
        'guiding_elements_each': fileValueFromKeyword(filename, 'guiding_elements_each', '='),
        'frame_guides_grid_edge': fileValueFromKeyword(filename, 'frame_guides_grid_edge', '='),
        'osmotic_density_inside': fileValueFromKeyword(filename, 'osmotic_density_inside', '='),
        'sw_position_min': fileValueFromKeyword(filename, 'sw_position_min', '='),
        'sw_position_max': fileValueFromKeyword(filename, 'sw_position_max', '='),
        'sw_position_target': fileValueFromKeyword(filename, 'sw_position_target', '='),
        'sw_orientation_min': fileValueFromKeyword(filename, 'sw_orientation_min', '='),
        'sw_orientation_max': fileValueFromKeyword(filename, 'sw_orientation_max', '='),
        'sw_orientation_target': fileValueFromKeyword(filename, 'sw_orientation_target', '='),
        'cell_min_edge': fileValueFromKeyword(filename, 'cell_min_edge', '='),
        'max_cells_dim': fileValueFromKeyword(filename, 'max_cells_dim', '='),
        'skip': fileValueFromKeyword(filename, 'skip', '='),
        'cluster_distance_threshold': fileValueFromKeyword(filename, 'cluster_distance_threshold', '=')
    }



def getSubclusterLabels(ID, group, eps):
    #skip if noise
    if ID == -1:
        return np.zeros( len(group.index), dtype=int )
    # arange a DBSCAN without PBC to get subclusters
    coms_subcluster = pd.concat([group['x'], group['y'], group['z']], axis=1)
    distances_array_subcluster = distance_array(coms_subcluster.values, coms_subcluster.values, box=None)
    dbscan_subcluster = DBSCAN(min_samples=1, eps=eps, metric="precomputed", n_jobs=-1).fit(distances_array_subcluster)
    # subcluster_labels.extend(dbscan_subcluster.labels_)
    return dbscan_subcluster.labels_



def getShiftedCoordinates(ID, group, eps, dimensions):
    # get largest subcluster
    unique, counts = np.unique(group["subcluster"], return_counts=True)
    if len(unique) == 1:
        return group["x"], group["y"], group["z"]
    max_subclusterID = unique[counts == np.max(counts)][0]
    # calculate shifts per subcluster
    centers = group.groupby("subcluster")['x','y','z'].mean()
    shifts = np.round(( -centers + centers.loc[max_subclusterID] )/dimensions[:3]).astype(int)
    shifts *= dimensions[:3]
    # calculate new coordinates based on shift
    newx = np.add(group["x"], shifts.loc[group["subcluster"]]["x"])
    newy = np.add(group["y"], shifts.loc[group["subcluster"]]["y"])
    newz = np.add(group["z"], shifts.loc[group["subcluster"]]["z"])
    # assign to main data
    return newx, newy, newz


def getClusterVolume(ID, group, eps, pps):
    # volume = 0
    if ID == -1:
        # single particles as spheres
        return np.pi * (eps**3) * 4/3
    else:
        # generate meshgrid around cluster particles plus threshold
        x_vector = np.arange(np.min(group["shiftx"])-eps, np.max(group["shiftx"])+eps+10/pps, 10/pps, dtype=np.float32)
        y_vector = np.arange(np.min(group["shifty"])-eps, np.max(group["shifty"])+eps+10/pps, 10/pps, dtype=np.float32)
        z_vector = np.arange(np.min(group["shiftz"])-eps, np.max(group["shiftz"])+eps+10/pps, 10/pps, dtype=np.float32)
        xx,yy,zz = np.meshgrid(x_vector, y_vector, z_vector)
        # stack them together as array of 3D points
        meshgrid = np.stack((xx.ravel(), yy.ravel(), zz.ravel()), axis=1)
        #calculate the distance array with centres of masses of particles
        coms_cluster = pd.concat([group['shiftx'], group['shifty'], group['shiftz']], axis=1)
        distances_array_volume = distance_array(meshgrid, coms_cluster.values, box=None)
        # check if any point in distance array row is close enough, then reshape to meshgrid
        # result is a binary meshgrid with 1 for the cluster shell region
        isclose = np.where(distances_array_volume <= eps, True, False).any(axis=1).reshape(xx.shape[0], yy.shape[1], zz.shape[2])
        # fill hole inside the shell region
        isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
        # calc volum from all points inside cluster
        return np.count_nonzero(isclose)*((10/pps)**3)

        # z,x,y = isclose.nonzero()
        # ax.scatter(x+,y,z, s=2)



def getOrder(ID, group):
    if ID == -1:
        return group["order"].values
    else:
        shifted_coms = pd.concat([group['shiftx'], group['shifty'], group['shiftz']], axis=1)
        normalized_orientations = pd.concat([group['ux'], group['uy'], group['uz']], axis=1)
        return (normalized_orientations.values*normalize(shifted_coms.sub(shifted_coms.mean().values))).sum(axis=1)