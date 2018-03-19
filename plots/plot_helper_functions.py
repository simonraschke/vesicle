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

pp = pprint.PrettyPrinter(indent=4, compact=True)



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
def getHDF5Dataset(filepath, datasetname):
    if os.path.exists(filepath):
        print("HDF5 file exists:", filepath)

        # try opening file
        try:
            FILE = h5py.File(filepath, 'r')
        except:
            print("cannot open file", filepath)
            return []
        # if datasetname in h5py.File(filepath, 'r'):

        # try getting dataset
        try:
            dataset = FILE.get(datasetname).value.transpose()
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
        pass



# open HDF5 dataset @datasetname at @filepath 
# return attribute value for @attributename
def getHDF5DatasetAttributes(filepath, datasetname, attributename):
    if os.path.exists(filepath):
        return h5py.File(filepath, 'r').get(datasetname).attrs[attributename]
    else:
        pass



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
def getVolume(datafilepath):
    if not os.path.exists(datafilepath):
        raise AssertionError("{} does not exist".format(datafilepath))
    else:
        try:
            x = float(getHDF5DatasetAttributes(datafilepath,"potential_energies","system.box.x"))
            y = float(getHDF5DatasetAttributes(datafilepath,"potential_energies","system.box.y"))
            z = float(getHDF5DatasetAttributes(datafilepath,"potential_energies","system.box.z"))
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



# read hdf5 file @datafilepath
# count occurencies of cluster @size
# average for @time_range
# return averaged value
def getClustersOfSizeTimeAverage(datafilepath, size, time_range, clustergroup="/cluster_self_assembled"):
    # print(getClustersOfSizeTimeAverage.__name__, "in", datafilepath,"in time range", time_range)
    try:
        file = h5py.File(datafilepath, 'r')
    except:
        return []
    values = []
    # predict dataset names
    # TODO: this is curcial for performance, time_range input is critical or no data will be found
    datasetnames = [clustergroup+"/time"+str(int(x)) for x in np.arange(int(min(time_range)),int(max(time_range))+1, 10000)]
    if int(min(time_range)) == int(max(time_range)):
        datasetnames = [clustergroup+"/time"+str(int(time_range[0]))]
    for datasetname in datasetnames:
        values.append(getClustersOfSize(file,datasetname,size)) 
    # return average if possible, else return zero
    if len(values) > 0:
        return np.average(values)
    else:
        return 0



# get the free particle density averaged over @time_range from @datafilepath
def __detail_getFreeParticleDensity_single_simulation_datafile(datafilepath, time_range):
    # print(getFreeParticleDensitySingleFile.__name__, "in", datafilepath, "in time range", time_range)
    try:
        volume = getVolume(datafilepath)
    except AssertionError:
        return np.NaN
    data = getHDF5Dataset(datafilepath, "cluster_volumes")
    # if len(data) >= 1:
        # print("dataset has length", len(data))
    inaccessible_volumes = extractXValueRange(list(data[0]),list(data[1]),time_range)[1]
    inaccessible_volume = np.average(inaccessible_volumes)
    free_particles = getClustersOfSizeTimeAverage(datafilepath, 1, time_range)
        # if len(free_particles) >= 1:
            # print("and has", free_particles, "free_particles")
            # free_particles = np.average(free_particles)
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
            pass
        elif np.isfinite(x[i]) and np.isfinite(y[i]):
            pass
        else:
            del x[i]
            del y[i]
    return x,y



# get rho free from all files with given @constraints in @time_range
def getFreeParticleDensity(overviewfilepath, constraints, time_range):
    print(getFreeParticleDensity.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    # get theses values in parallel
    results = []
    for path in paths:
        results.append(pool.apply_async(__detail_getFreeParticleDensity_single_simulation_datafile,(path, time_range,)))
    rho_free_values = removeBadEntries1D([r.get() for r in results if r.get() != None])
    # print("rho_free_values", rho_free_values)
    if len(rho_free_values) > 0:
        # print("and averaged", np.average(rho_free_values))
        return np.average(rho_free_values)
    else:
        return None


#MUST be in at the EOF
pool = multiprocessing.Pool(8)