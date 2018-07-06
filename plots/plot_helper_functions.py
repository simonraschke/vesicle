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

import sys
import os
import re
import fileinput
import string
import math
import shutil
import subprocess
import time
import json
import pprint
import numpy as np
import h5py
import multiprocessing
import contextlib

pp = pprint.PrettyPrinter(indent=4, compact=True)


def makesize(width,ratio=None, scale=1):
    fig_width_pt = width  # Get this from LaTeX using \the\XXXXwidth
    inches_per_pt = 1.0/72.27                                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0
    if ratio == None:
        ratio = golden_mean
    #golden_mean = 0.54  if twocolumn == True else 0.9
    fig_width = fig_width_pt*inches_per_pt*scale                    # width in inches
    fig_height = fig_width*ratio                              # height in inches
    return [fig_width,fig_height]



class DummyFile(object):
    def write(self, x): pass



@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = DummyFile()
    yield
    sys.stdout = save_stdout



def averageNestedLists(nested_vals):
    """
    Averages a 2-D array and returns a 1-D array of all of the columns
    averaged together, regardless of their dimensions.
    """
    output = []
    maximum = 0
    for lst in nested_vals:
        if len(lst) > maximum:
            maximum = len(lst)
    for index in range(maximum): # Go through each index of longest list
        temp = []
        for lst in nested_vals: # Go through each list
            if index < len(lst): # If not an index error
                temp.append(lst[index])
        if len(temp) > 0:
            try:
                output.append(np.nanmean(temp))
            except:
                output.append(np.NaN)
        else:
            output.append(np.NaN)
    return output



def numbersListFromString(string):
    return re.findall('\d+', string)


# iterate a file linewise
# check line for keyword
# if keyword is found return value after seperator
def fileValueFromKeyword(filepath, keyword, seperator='='):
    found_counter = 0
    assert(filepath)
    value = None
    with open(filepath) as FILE:
        for line in FILE:
            if keyword in line:
                found_counter += 1
                value = line.split(seperator,1)[1].rstrip()
                break
    if found_counter == 0:
        raise Exception("unable to find keyword "+keyword+" in file "+filepath+"\n")
    # if found_counter >= 2:
        # raise Exception("multiple keywords "+keyword+" in file "+filepath+"\n")
    return value



def writeJson(filepath, obj):
    with open(filepath, 'w') as fout:
        json.dump(obj, fout, indent=4, sort_keys=True)



def readJson(filepath):
    assert(os.path.exists(filepath))
    with open(filepath) as json_data:
        d = json.load(json_data)
        return d



def getValues(filepath, key):
    values = []
    for paramdict in readJson(filepath):
        values.append(paramdict[key])
    return sorted(list(set(values)))



def getValuesInRange(overviewfilepath, key, range):
    return sorted(list(set( [float(x) for x in getValues(overviewfilepath, key) if float(min(range)) <= float(x) <= float(max(range))] )))



def getDirs(filepath):
    values = []
    for paramdict in readJson(filepath):
        values.append(paramdict["dir_path"])
    return sorted(list(set(values)))



# return all dicts wich match 'constraints'
# constraints mus be dict
def getMatched(filepath, constraints):
    chosen_dicts = []
    for paramdict in readJson(filepath):
        if all( paramdict[param] == value for param, value in constraints.items() ):
            chosen_dicts.append(paramdict)
    return chosen_dicts



# return all values for 'param'
# in dicts that match all in 'constraints'
# constraints mus be dict
def getMatchedValues(filepath, param, constraints):
    values = []
    for paramdict in getMatched(filepath, constraints):
        values.append(paramdict[param])
    return sorted(list(set(values)))



# return all values for 'dir_path' 
# in dicts that match all in 'constraints'
# constraints mus be dict
def getMatchedDirs(overviewfilepath, constraints):
    return getMatchedValues(overviewfilepath, "dir_path", constraints)



# print an overview of a HDF5 file
def fileOverviewHDF5(filepath):
    ph5 = pprint.PrettyPrinter(indent=4, compact=True, width=1)
    file = h5py.File(filepath, 'r')
    keys = file.keys()
    groups = list(keys)
    print("<HDF5 file", filepath+">")
    for groupname in groups:
        temp = file[groupname]
        print(temp)



# open a HDF5 dataset in @filepath
# return np.array
# if filepath doesnt exist return None
# if Exception return empty list []
def getHDF5Dataset(filepath, datasetname, column="all"):
    if os.path.exists(filepath):
        print("HDF5 file exists:", filepath)

        # try opening file
        try:
            FILE = h5py.File(filepath, 'r')
        except:
            print("cannot open file", filepath)
            return []

        # try getting dataset
        try:
            if column == "all":
                dataset = FILE.get(datasetname).value.transpose()
            else:
                dataset = FILE.get(datasetname).value[:,column].transpose()
        except:
            print("cannot open dataset", datasetname)
            return []

        # check for data in dataset
        try:
            if dataset.ndim > 1:
                return dataset
            else:
                pass
        except:
            return []
    else:
        print("filepath", filepath, "does not exist")
        pass



# open HDF5 dataset @datasetname at @filepath 
# return attribute value for @attributename
def getHDF5DatasetAttributes(filepath, datasetname, attributename):
    if os.path.exists(filepath):
        return h5py.File(filepath, 'r').get(datasetname).attrs[attributename]
    else:
        pass



def save_dict_to_hdf5(dic, filename):
    """
    Save a dictionary whose contents are only strings, np.float64, np.int64,
    np.ndarray, and other dictionaries following this structure
    to an HDF5 file. These are the sorts of dictionaries that are meant
    to be produced by the ReportInterface__to_dict__() method.
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    Take an already open HDF5 file and insert the contents of a dictionary
    at the current path location. Can call itself recursively to fill
    out HDF5 files with the contents of a dictionary.
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes)):
            h5file[path + key] = item
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type'%type(item))

def load_dict_from_hdf5(filename):
    """
    Load a dictionary whose contents are only strings, floats, ints,
    numpy arrays, and other dictionaries following this structure
    from an HDF5 file. These dictionaries can then be used to reconstruct
    ReportInterface subclass instances using the
    ReportInterface.__from_dict__() method.
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    """
    Load contents of an HDF5 group. If further groups are encountered,
    treat them like dicts and continue to load them recursively.
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans



# open @overviewfilepath to look for dir_paths with given @constraints
# open all @datasetname datasets (which must be 2D arrays) and average every field
# return averaged [np.array]s
def getAveragedData2D(overviewfilepath, datasetname, constraints):
    print(getAveragedData2D.__name__, "  get dataset", datasetname, "  with constraints", constraints)
    # get dir_path's
    dirs = getMatchedDirs(overviewfilepath,constraints)
    # prepare datafilenames
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # load everythin into list of np.arrays
    all_data = [getHDF5Dataset(path,datasetname) for path in paths]
    # put all columns of same number in buffers
    x_buffer = np.array([ data[0] for data in all_data ])
    y_buffer = np.array([ data[1] for data in all_data ])
    x, y = [], []
    # iterate range from 0 to length of longest array in buffer
    for i in range(max([ len(array) for array in x_buffer ])):
        x_temp = []
        # if possible add all elements of same number to temp
        for array in x_buffer:
            try:
                x_temp.append(array[i])
            except IndexError:
                pass
        # and add the average to x
        x.append(np.average(x_temp))
    # iterate range from 0 to length of longest array in buffer
    for i in range(max([ len(array) for array in y_buffer ])):
        y_temp = []
        # if possible add all elements of same number to temp
        for array in y_buffer:
            try:
                y_temp.append(array[i])
            except IndexError:
                pass
        # and add the average to y
        y.append(np.average(y_temp))
    return [x,y]



# input 1D array and range of values
# return all indexes to elements in value range
def indexOfValueRange(column,value_range):
    return [ column.index(x) for x in column if x >= min(value_range) and x <= max(value_range) ]



# input 2 1D arrays and get indexes to value range of x
# return 2 arrays with the elements for given indexes
def extractXValueRange(x,y,value_range):
    indexes = indexOfValueRange(x,value_range)
    new_x = [x[i] for i in indexes]
    new_y = [y[i] for i in indexes]
    return [new_x, new_y]



# read HDF5 file @datafilepath
# try and read system.box parameters from potential_energies dataset
# if possible return volume else NaN
def getVolume(datafilepath,elementname="cluster_self_assembled"):
    if not os.path.exists(datafilepath):
        raise AssertionError("{} does not exist".format(datafilepath))
    else:
        try:
            x = float(getHDF5DatasetAttributes(datafilepath,elementname,"system.box.x"))
            y = float(getHDF5DatasetAttributes(datafilepath,elementname,"system.box.y"))
            z = float(getHDF5DatasetAttributes(datafilepath,elementname,"system.box.z"))
        except:
            return np.NaN
        return x*y*z



# unused
def getVolumes(filepath, constraints):
    print(getVolumes.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(filepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    return [getVolume(path) for path in paths]



# input hdf5 file object, datasetname and the searched size
# return the number of occurencies of certain value (here: cluster size)
def getClustersOfSize(datafile, datasetname, size):
    dataset = datafile.get(datasetname).value.transpose()
    unique, counts = np.unique(dataset[0], return_counts=True)
    return dict(zip(unique, counts))[size]



# input hdf5 file object, datasetname and the searched size
# return the number of occurencies greater than certain value (here: cluster size)
def getClustersGreaterThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname).value.transpose()
    unique, counts = np.unique(dataset[0], return_counts=True)
    return sum([c for u,c in zip(unique,counts) if u > size])



# input hdf5 file object, datasetname and the searched size
# return the number of occurencies smaller than certain value (here: cluster size)
def getClustersSmallerThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname).value.transpose()
    unique, counts = np.unique(dataset[0], return_counts=True)
    return sum([c for u,c in zip(unique,counts) if u < size])
    


def getLargestClusterGreaterEqual(datafile, datasetname, min_size=1):
    dataset = datafile.get(datasetname).value.transpose()
    unique, counts = np.unique(dataset[0], return_counts=True)
    return max(unique) if max(unique) >= min_size else 0



def getNumParticlesInClustersGreaterThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname).value.transpose()
    unique, counts = np.unique(dataset[0], return_counts=True)
    return sum([c*u for u,c in zip(unique,counts) if u > size])



def getNumParticlesInClustersSmallerThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname).value.transpose()
    unique, counts = np.unique(dataset[0], return_counts=True)
    return sum([c*u for u,c in zip(unique,counts) if u < size])



def getVolumeOfClustersGreaterThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname).value.transpose()
    return sum([ dataset[2,i] for i in range(len(dataset[0])) if dataset[0,i] > size ])



def getVolumeOfClustersSmallerThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname).value.transpose()
    return sum([ dataset[2,i] for i in range(len(dataset[0])) if dataset[0,i] < size ])



def getOrderOfClustersGreaterThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname)
    if dataset == None:
        return np.NaN
    else:
        dataset = datafile.get(datasetname).value.transpose()
    num_particles = getNumParticlesInClustersGreaterThan(datafile,datasetname,size)
    if np.absolute(num_particles) < 1e-6:
        return np.NaN
    else:
        return sum([ dataset[1,i]*dataset[0,i] for i in range(len(dataset[0])) if dataset[0,i] > size ]) / getNumParticlesInClustersGreaterThan(datafile,datasetname,size)



def getOrderOfClustersSmallerThan(datafile, datasetname, size):
    dataset = datafile.get(datasetname)
    if dataset == None:
        return np.NaN
    else:
        dataset = datafile.get(datasetname).value.transpose()
    num_particles = getNumParticlesInClustersGreaterThan(datafile,datasetname,size)
    if np.absolute(num_particles) < 1e-6:
        return np.NaN
    else:
        return sum([ dataset[1,i]*dataset[0,i] for i in range(len(dataset[0])) if dataset[0,i] < size ]) / getNumParticlesInClustersGreaterThan(datafile,datasetname,size)



def getAverageParticlesInClustersizeGreaterThan(datafile, datasetname, size):
    return getNumParticlesInClustersGreaterThan(datafile,datasetname,size) / getClustersGreaterThan(datafile,datasetname,size)



def getAverageParticlesInClustersizeSmallerThan(datafile, datasetname, size):
    return getNumParticlesInClustersSmallerThan(datafile,datasetname,size) / getClustersSmallerThan(datafile,datasetname,size)



def __detail_getTimeOfFirstCluster_singleSimulationDatafile(datafilepath, time_range, size, time_increment=10000, clustergroup="/cluster_self_assembled"):
    try:
        file = h5py.File(datafilepath, 'r')
    except:
        return []
    # predict dataset names
    # WARNING: this is curcial for performance, time_range input is critical or no data will be found
    datasetnames = [clustergroup+"/time"+str(int(x)) for x in np.arange(int(min(time_range)),int(max(time_range))+1, time_increment)]
    if int(min(time_range)) == int(max(time_range)):
        datasetnames = [clustergroup+"/time"+str(int(time_range[0]))]
    for datasetname in datasetnames:
        try:
            if getClustersGreaterThan(file, datasetname, size) > 0:
                print("got", numbersListFromString(datasetname)[0], "in", datasetname)
                return float(numbersListFromString(datasetname)[0])
        except:
            pass
    # return float(numbersListFromString(datasetnames[-1])[0])
    return None



# read hdf5 file @datafilepath
# count occurencies of FUNCTOR @size
# average for @time_range
# return averaged value
def getTimeAverageNum_FUNCTOR(datafilepath, size, time_range, functor, time_increment=10000, clustergroup="/cluster_self_assembled"):
    try:
        file = h5py.File(datafilepath, 'r')
    except:
        print("cannot open file", datafilepath)
        return []
    values = []
    # predict dataset names
    # WARNING: this is curcial for performance, time_range input is critical or no data will be found
    datasetnames = [clustergroup+"/time"+str(int(x)) for x in np.arange(int(min(time_range)),int(max(time_range))+1, time_increment)]
    if int(min(time_range)) == int(max(time_range)):
        datasetnames = [clustergroup+"/time"+str(int(time_range[0]))]
    for datasetname in datasetnames:
        result = functor(file, datasetname, size)
        if result != None:
            values.append(result)
    # return average if possible, else return zero
    cleaned_data = np.ma.masked_array(values,mask=np.isnan(values))
    if len(values) > 0:
        print("time_averaged",np.ma.average(cleaned_data))
        print(cleaned_data.compressed())
        return np.ma.average(cleaned_data)
    else:
        print("time_averaged", 0)
        return 0



# get the free particle density averaged over @time_range from @datafilepath
def __detail_getFreeParticleDensity_single_simulation_datafile(datafilepath, time_range):
    try:
        volume = getVolume(datafilepath)
    except AssertionError:
        return np.NaN
    inaccessible_volumes = []
    try:
        inaccessible_volumes = getTimeAverageNum_FUNCTOR(datafilepath, 5, time_range, getVolumeOfClustersGreaterThan)
    except:
        return np.NaN
    inaccessible_volume = np.average(inaccessible_volumes)
    free_particles = getTimeAverageNum_FUNCTOR(datafilepath, 5, time_range, getNumParticlesInClustersSmallerThan)
    try:
        print("free particles:", int(free_particles), "  vol:", round(volume), "  vol_in:", round(inaccessible_volume), "  == rho_free:", round(float(free_particles)/(volume - inaccessible_volume),5))
    except:
        print("possible TypeError: free particles:          ", free_particles )
        print("possible TypeError: volume:                  ", volume )
        print("possible TypeError: free inaccessible_volume:", inaccessible_volume )
        return np.NaN
    return float(free_particles)/(volume - inaccessible_volume)



# remove all None and !np.isfinite fields of array
def removeBadEntries1D(x):
    if len(x) == 0:
        return []
    for i in reversed(range(len(x))):
        if np.isfinite(x[i]):
            pass
        else:
            del x[i]
    return x



# remove all None and !np.isfinite fields of array
def removeBadEntries2D(x,y):
    assert(len(x)==len(y))
    for i in reversed(range(len(x))):
        if x[i] == None or y[i] == None:
            # pass
            # FIXME: MAYBE REVERT
            del x[i]
            del y[i]
        elif np.isfinite(x[i]) and np.isfinite(y[i]):
            pass
        else:
            del x[i]
            del y[i]
    return x,y



# get rho free from all files with given @constraints in @time_range
def getFreeParticleDensity(overviewfilepath, constraints, time_range):
    print(getFreeParticleDensity.__name__,"of time_range", time_range, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # get theses values in parallel
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getFreeParticleDensity_single_simulation_datafile,(path, time_range,)))
    rho_free_values = removeBadEntries1D([r.get() for r in results if r.get() != None])
    if len(rho_free_values) > 0:
        print(getFreeParticleDensity.__name__, "  result", "{:.5f}".format(np.average(rho_free_values)))
        print()
        return np.average(rho_free_values)
    else:
        print(getFreeParticleDensity.__name__, "  result", "None")
        print()
        return None



def __detail_getFreeParticleDensityFullEvolution_single_simulation_datafile(datafilepath, size):
    print(__detail_getFreeParticleDensityFullEvolution_single_simulation_datafile.__name__, "in", datafilepath )
    rho_free = []
    FILE = None
    group = None
    volume = None
    try:
        volume = getVolume(datafilepath)
    except AssertionError:
        print("cannot get volume of", datafilepath)
        return [],[]
    try:
        FILE = h5py.File(datafilepath, 'r')
    except:
        print("cannot open file", datafilepath)
        return [],[]
    try:
        group = FILE.get("cluster_self_assembled")
    except:
        print("cannot find group /cluster_self_assembled in", datafilepath)
        return [],[]

    timepoints = []
    inaccessible_volumes = []
    free_particles = []
    datasets = sorted([ d.name for d in group.values() ])
    # only read every n'th dataset
    # datasets = datasets[:100:1]+datasets[100::10]
    for datasetname in datasets:
        timepoints.append(float(numbersListFromString(datasetname)[0]))
        dataset = group.get(datasetname).value.transpose()
        inaccessible_volumes.append(sum([ dataset[2,i] for i in range(len(dataset[0])) if dataset[0,i] > size ]))
        unique, counts = np.unique(dataset[0], return_counts=True)
        free_particles.append(sum([c*u for u,c in zip(unique,counts) if u < size]))
        # if len(free_particles) == len(timepoints) + 1:
    rho_free = [ float(fp)/(volume - vi) for fp,vi in zip(free_particles, inaccessible_volumes) ]
    assert(len(timepoints) == len(rho_free))
    # print("got",len(timepoints),"timepoints and",len(rho_free),"rho_frees")
    return timepoints, rho_free



def getFreeParticleDensityTimeEvolution(overviewfilepath, constraints, time_range):
    print(getFreeParticleDensityTimeEvolution.__name__,"of time_range", time_range, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    timepoints = np.logspace(np.log10(min(time_range)), np.log10(max(time_range)), num=20, endpoint=True, base=10.0, dtype=float)
    timepoints = np.insert(timepoints, 0, 0)
    timepoints = sorted(list(set(np.around(timepoints, decimals=-4).tolist())))
    rho_free_values = []
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getFreeParticleDensityFullEvolution_single_simulation_datafile,(path, 5,)))
    for r in results:
        timepoints_single, rho_free_single = r.get()
        assert(len(timepoints_single) == len(rho_free_single))
        rho_free_average_single = []
        for i,tp in enumerate(timepoints):
            accumulator = [ rho_free_single[j] for j,tps in enumerate(timepoints_single) if timepoints[i-1] < tps <= timepoints[i] ]
            if(len(accumulator) > 0):
                rho_free_average_single.append(np.average(accumulator))
            # except AssertionError:
            #     print("ASSERTION ERROR: between timepoints", timepoints[i-1], "and", timepoints[i] , "accumulator", accumulator)
        if len(rho_free_average_single) > 0:
            rho_free_values.append(rho_free_average_single)
    if len(rho_free_values) > 0:
        rho_free_values = averageNestedLists(rho_free_values)
    timepoints.pop(0)
    assert(len(timepoints) == len(rho_free_values))
    print("got",len(timepoints),"timepoints and",len(rho_free_values),"rho_free_values")
    return timepoints, rho_free_values



# get the order averaged over @time_range from @datafilepath
def __detail_getOrder_single_simulation_datafile(datafilepath, time_range, min_size):
    order_values = []
    # try:
    order_values = getTimeAverageNum_FUNCTOR(datafilepath, min_size-1, time_range, getOrderOfClustersGreaterThan)
    #     print("order_values",order_values)
    # except:
    #     print("order_values EXCEPTION",order_values)
    #     return np.NaN
    return np.average([order_values])



# get order from all files with given @constraints in @time_range
def getOrder(overviewfilepath, constraints, time_range, min_size, part="full"):
    print(getOrder.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # get theses values in parallel
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getOrder_single_simulation_datafile,(path, time_range, min_size,)))
    order_values = removeBadEntries1D([r.get() for r in results if r.get() != None])
    if len(order_values) > 0:
        print(getOrder.__name__, "  result", "{:.5f}".format(np.average(order_values)))
        print()
        return np.average(order_values)
    else:
        print(getOrder.__name__, "  result", "None")
        print()
        return None



# get the order averaged over @time_range from @datafilepath
def __detail_getAverageClustersize_single_simulation_datafile(datafilepath, time_range, min_size):
    average_values = []
    try:
        average_values = getTimeAverageNum_FUNCTOR(datafilepath, min_size-1, time_range, getAverageParticlesInClustersizeGreaterThan)
    except:
        return np.NaN
    return np.average(average_values)


# get order from all files with given @constraints in @time_range
def getAverageClustersize(overviewfilepath, constraints, time_range, min_size, part="full"):
    print(getAverageClustersize.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # get theses values in parallel
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getAverageClustersize_single_simulation_datafile,(path, time_range, min_size)))
    average_values = removeBadEntries1D([r.get() for r in results if r.get() != None])
    if len(average_values) > 0:
        print(getAverageClustersize.__name__, "  result", "{:.5f}".format(np.average(average_values)))
        print()
        return np.average(average_values)
    else:
        print(getAverageClustersize.__name__, "  result", "None")
        print()
        return None



# get the order averaged over @time_range from @datafilepath
def __detail_getLargestCluster_single_simulation_datafile(datafilepath, time_range, min_size):
    average_values = []
    try:
        average_values = getTimeAverageNum_FUNCTOR(datafilepath, min_size, time_range, getLargestClusterGreaterEqual)
    except:
        return np.NaN
    return np.average(average_values)


# get order from all files with given @constraints in @time_range
def getLargestCluster(overviewfilepath, constraints, time_range, min_size, part="full"):
    print(getLargestCluster.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # get theses values in parallel
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getLargestCluster_single_simulation_datafile,(path, time_range, min_size)))
    average_values = removeBadEntries1D([r.get() for r in results if r.get() != None])
    if len(average_values) > 0:
        print(getLargestCluster.__name__, "  result", "{:.5f}".format(np.average(average_values)))
        print()
        return np.average(average_values)
    else:
        print(getLargestCluster.__name__, "  result", "None")
        print()
        return None



def __detail_getlargestClusterFullEvolution_single_simulation_datafile(datafilepath):
    print(__detail_getlargestClusterFullEvolution_single_simulation_datafile.__name__, "in", datafilepath )
    N_max = []
    FILE = None
    group = None
    try:
        FILE = h5py.File(datafilepath, 'r')
    except:
        print("cannot open file", datafilepath)
        return [],[]
    try:
        group = FILE.get("cluster_self_assembled")
    except:
        print("cannot find group /cluster_self_assembled in", datafilepath)
        return [],[]

    timepoints = []
    datasets = sorted([ d.name for d in group.values() ])
    for datasetname in datasets:
        timepoints.append(float(numbersListFromString(datasetname)[0]))
        dataset = group.get(datasetname).value.transpose()
        unique, counts = np.unique(dataset[0], return_counts=True)
        N_max.append(max(unique))
    assert(len(timepoints) == len(N_max))
    return timepoints, N_max



def getLargestClusterTimeEvolution(overviewfilepath, constraints, time_range):
    print(getLargestClusterTimeEvolution.__name__,"of time_range", time_range, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    timepoints = np.logspace(np.log10(min(time_range)), np.log10(max(time_range)), num=20, endpoint=True, base=10.0, dtype=float)
    timepoints = np.insert(timepoints, 0, 0)
    timepoints = sorted(list(set(np.around(timepoints, decimals=-4).tolist())))
    print(timepoints)
    N_max_values = []
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getlargestClusterFullEvolution_single_simulation_datafile,(path,)))
    for r in results:
        timepoints_single, N_max_single = r.get()
        assert(len(timepoints_single) == len(N_max_single))
        N_max_average_single = []
        for i,tp in enumerate(timepoints):
            accumulator = [ N_max_single[j] for j,tps in enumerate(timepoints_single) if timepoints[i-1] < tps <= timepoints[i] ]
            if(len(accumulator) > 0):
                N_max_average_single.append(np.average(accumulator))
        if len(N_max_average_single) > 0:
            N_max_values.append(N_max_average_single)
    if len(N_max_values) > 0:
        N_max_values = averageNestedLists(N_max_values)
    timepoints.pop(0)
    for i,j in zip(timepoints, N_max_values): print(i,j)
    assert(len(timepoints) == len(N_max_values))
    print("got",len(timepoints),"timepoints and",len(N_max_values),"N_max_values")
    return timepoints, N_max_values



# get the order averaged over @time_range from @datafilepath
def __detail_getGamma_single_simulation_datafile(datafilepath, time_range, min_size):
    average_values = []
    try:
        average_values = __detail_getTimeOfFirstCluster_singleSimulationDatafile(datafilepath, min_size-1, time_range, getClustersGreaterThan)
    except:
        return np.NaN, np.NaN
    time_of_first_cluster = np.average(average_values)
    rounded_time = np.around(time_of_first_cluster,-4)
    rho_free = np.NaN
    with nostdout():
        rho_free = __detail_getFreeParticleDensity_single_simulation_datafile(datafilepath, [rounded_time, rounded_time])
    return rho_free, time_of_first_cluster



# get order from all files with given @constraints in @time_range
def getGamma(overviewfilepath, constraints, time_range, min_size, part="full"):
    print(getGamma.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # get theses values in parallel
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getGamma_single_simulation_datafile,(path, time_range, min_size)))
    rho_free,gamma = [],[]
    for r in results:
        rho_free_temp, gamma_temp = r.get()
        if np.isfinite(rho_free_temp) and np.isfinite(gamma_temp) and gamma_temp < (max(time_range)-1e-5): 
            rho_free.append(float(rho_free_temp))
            gamma.append(float(gamma_temp))
    return rho_free, gamma



def getParticleHistory(overviewfilepath, constraints, time_range, particle, part="full"):
    print(getParticleHistory.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    for path in paths:
        result = load_dict_from_hdf5(path)
        pp.pprint(result)
        sys.exit()


def getTimeOfFirstCluster(overviewfilepath, constraints, time_range, size, part="full"):
    print(getTimeOfFirstCluster.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getTimeOfFirstCluster_singleSimulationDatafile,(path, time_range, size, 10000)))
    time = []
    for r in results:
        time_temp = r.get()
        if time_temp != None:
            if np.isfinite(time_temp) and time_temp != None: 
                time.append(float(time_temp))
    if len(time) > 0:
        return np.average(time)



# get the order averaged over @time_range from @datafilepath
def __detail_getPotentialEnergy_single_simulation_datafile(datafilepath, time_range):
    dataset = []
    try:
        dataset = getHDF5Dataset(datafilepath,"potential_energies")
    except:
        return None
    epot_values = []
    if isinstance(dataset, (list, np.generic)):
        if len(dataset) == 2:
            for time, epot in zip(dataset[0],dataset[1]):
                if min(time_range) <= time <= max(time_range):
                    epot_values.append(epot)
            return np.average(epot_values)
    return None



# get order from all files with given @constraints in @time_range
def getPotentialEnergy(overviewfilepath, constraints, time_range):
    print(getPotentialEnergy.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # get theses values in parallel
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getPotentialEnergy_single_simulation_datafile,(path, time_range,)))
    epot_values = removeBadEntries1D([r.get() for r in results if r.get() != None])
    if len(epot_values) > 0:
        print(getPotentialEnergy.__name__, "  result", "{:.5f}".format(np.average(epot_values)))
        print()
        return np.average(epot_values)
    else:
        print(getPotentialEnergy.__name__, "  result", "None")
        print()
        return None



def __detail_getExchangeFullEvolution_single_simulation_datafile(datafilepath, time_range, min_size=50, time_increment=1e4):
    print(__detail_getExchangeFullEvolution_single_simulation_datafile.__name__, "in", datafilepath )
    FILE = None
    group = None
    try:
        FILE = h5py.File(datafilepath, 'r')
    except:
        print("cannot open file", datafilepath)
        return [],[]
    try:
        group = FILE.get("cluster_self_assembled")
    except:
        print("cannot find group /cluster_self_assembled in", datafilepath)
        return [],[]
    data_before = []
    timepoints = [x for x in np.arange(min(time_range),max(time_range)+1,time_increment, dtype=int) ]
    exchanges = [np.NaN]*len(timepoints)
    datasets = ["time"+str(x) for x in timepoints ]
    for i, timepoint in enumerate(timepoints):
        dataset = group.get(datasets[i])
        data_now = []
        if dataset != None:
            data_now = dataset.value
        if len(data_before) == 0 or len(data_now) == 0:
            data_before = data_now
            continue
        clusters = []
        for row in data_now:
            if row[0] >= min_size+4:
                temp_list = [int(x) for x in row[4:] if x >= 0]
                try:
                    assert(len(temp_list) == int(row[0]))
                except AssertionError:
                    assert(False)
                clusters.append(sorted(temp_list))
        if len(clusters) > 0:
            data_now = sorted(clusters)
        else:
            data_now = []
        exchange_collector = []
        if i == 0:
            continue
        for cluster in data_now:
            biggest_overlap = 0
            try:
                biggest_overlap = max( [ len(set(cluster).intersection(set(old_cluster))) for old_cluster in data_before] )
            except:
                pass
            old_cluster = []
            for old in data_before:
                if len(set(cluster).intersection(set(old))) == biggest_overlap:
                    old_cluster = old
                    break
            size_before = len(cluster)
            number_of_old_in_new = biggest_overlap
            exchange = (size_before-number_of_old_in_new)/size_before
            exchange_collector.append(exchange)
        if len(exchange_collector) > 0:
            exchanges[i] = np.average(exchange_collector)
        data_before = data_now
    print("will return", len(timepoints),"timepoints and",len(exchanges),"exchange_values")
    return timepoints, exchanges



def getExchangeTimeEvolution(overviewfilepath, constraints, time_range):
    print(getExchangeTimeEvolution.__name__,"of time_range", time_range, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    min_size = 20
    time_increment = 1e4
    timepoints = np.logspace(np.log10(min(time_range)), np.log10(max(time_range)), num=20, endpoint=True, base=10.0, dtype=float)
    timepoints = np.insert(timepoints, 0, 0)
    timepoints = sorted(list(set(np.around(timepoints, decimals=-4).tolist())))
    print(timepoints)
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getExchangeFullEvolution_single_simulation_datafile,(path, time_range, min_size, time_increment,)))
    timepoint_values = []
    exchange_values = []
    for r in results:
        timepoints_single, exchange_single = r.get()
        exchange_average_single = []
        for i, tp in enumerate(timepoints):
            accumulator = [ exchange_single[j] for j,tps in enumerate(timepoints_single) if timepoints[i-1] < tps <= timepoints[i] ]
            if(len(accumulator) > 0):
                exchange_average_single.append(np.average(accumulator))
        if len(exchange_average_single) > 0:
            exchange_values.append(exchange_average_single)
    if len(exchange_values) > 0:
        exchange_values = averageNestedLists(exchange_values)
    timepoints.pop(0)
    print("got",len(timepoints),"timepoints and",len(exchange_values),"exchange_values")
    assert(len(timepoints) == len(exchange_values))
    return timepoints, exchange_values



def __detail_getClusterSizeDistribution_single_simulation_datafile(datafilepath, timepoint, clustergroup):
    # print(__detail_getClusterSizeDistribution_single_simulation_datafile.__name__, "in", datafilepath )
    datasetname = str(clustergroup+"/time"+str(int(timepoint)))
    print("get size distribution from ", datafilepath+":"+datasetname)
    try:
        data = getHDF5Dataset(datafilepath, datasetname)[0,:]
        return data
    except:
        return []
    # print (data)


def getClusterSizeDistribution(overviewfilepath, constraints, timepoint, clustergroup="/cluster_self_assembled"):
    print(getClusterSizeDistribution.__name__,"of time", timepoint, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getClusterSizeDistribution_single_simulation_datafile,(path, timepoint, clustergroup,)))
    hist_array = []
    for r in results:
        hist_array.extend(list(r.get()))
    # print(np.array(hist_array, dtype=int))
    return np.array(hist_array, dtype=int)



#MUST be at the EOF
# pool = multiprocessing.Pool(1)

import traceback
from multiprocessing.pool import Pool
# Shortcut to multiprocessing's logger
def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)

class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise

        # It was fine, give a normal answer
        return result

class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)

multiprocessing.log_to_stderr()
pool = LoggingPool(processes=10,maxtasksperchild=1)