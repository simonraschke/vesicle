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
def getMatchedDirs(filepath, constraints):
    return getMatchedValues(filepath, "dir_path", constraints)



def fileOverviewHDF5(filepath):
    ph5 = pprint.PrettyPrinter(indent=4, compact=True, width=1)
    file = h5py.File(filepath, 'r')
    keys = file.keys()
    groups = list(keys)
    print("<HDF5 file", filepath+">")
    for groupname in groups:
        temp = file[groupname]
        print(temp)



def getHDF5Dataset(filepath, datasetname):
    if os.path.exists(filepath):
        dataset = h5py.File(filepath, 'r').get(datasetname).value.transpose()
        if dataset.ndim > 1:
            return dataset
        else:
            pass
    else:
        pass



def getHDF5DatasetAttributes(filepath, datasetname, attributename):
    if os.path.exists(filepath):
        return h5py.File(filepath, 'r').get(datasetname).attrs[attributename]
    else:
        pass



def getAveragedData2D(filepath, datasetname, constraints):
    print(getAveragedData2D.__name__, "  get dataset", datasetname, "  with constraints", constraints)
    dirs = getMatchedDirs(filepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    all_data = [getHDF5Dataset(path,datasetname) for path in paths]
    x_buffer = np.array([ data[0] for data in all_data ])
    y_buffer = np.array([ data[1] for data in all_data ])
    x, y = [], []
    for i in range(max([ len(array) for array in x_buffer ])):
        x_temp = []
        for array in x_buffer:
            try:
                x_temp.append(array[i])
            except IndexError:
                pass
        x.append(np.average(x_temp))
    for i in range(max([ len(array) for array in y_buffer ])):
        y_temp = []
        for array in y_buffer:
            try:
                y_temp.append(array[i])
            except IndexError:
                pass
        y.append(np.average(y_temp))
    return [x,y]



def indexOfValueRange(column,value_range):
    return [ column.index(x) for x in column if x >= min(value_range) and x <= max(value_range) ]



def extractXValueRange(x,y,value_range):
    indexes = indexOfValueRange(x,value_range)
    new_x = [x[i] for i in indexes]
    new_y = [y[i] for i in indexes]
    return [new_x, new_y]



def getVolume(filepath):
    x = float(getHDF5DatasetAttributes(filepath,"potential_energies","system.box.x"))
    y = float(getHDF5DatasetAttributes(filepath,"potential_energies","system.box.y"))
    z = float(getHDF5DatasetAttributes(filepath,"potential_energies","system.box.z"))
    return x*y*z



def getVolumes(filepath, constraints):
    print(filepath)
    dirs = getMatchedDirs(filepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    return [getVolume(path) for path in paths]



def getClustersOfSize(datafilepath, size, time_range):
    file = h5py.File(datafilepath, 'r')
    values = []
    for name, data in file['/cluster_self_assembled'].items():
        time = float(re.findall('\d+', name)[0])
        if time >= min(time_range) and time <= max(time_range):
            # values.append(data.value[0].aslist.count(size))
            unique, counts = np.unique(data.value.transpose()[0], return_counts=True)
            values.append(dict(zip(unique, counts))[size])
        else:
            pass
    return values


def getFreeParticleDensitySingleFile(datafilepath, time_range):
    # volume = np.average(getVolumes(filepath, constraints))
    volume = getVolume(datafilepath)
    data = getHDF5Dataset(datafilepath, "cluster_volumes")
    inaccessible_volumes = extractXValueRange(list(data[0]),list(data[1]),time_range)[1]
    inaccessible_volume = np.average(inaccessible_volumes)
    free_particles = getClustersOfSize(datafilepath, 1, time_range)
    print(free_particles, inaccessible_volumes)
    free_particles = np.average(free_particles)
    print(free_particles, inaccessible_volume)
    return float(free_particles)/(volume - inaccessible_volume)



def getFreeParticleDensity(overviewfilepath, constraints, time_range):
    print(getFreeParticleDensity.__name__, "  with constraints", constraints)
    dirs = getMatchedDirs(overviewfilepath,constraints)
    paths = [os.path.join(dir,"data.h5") for dir in dirs]
    return np.average([getFreeParticleDensitySingleFile(path, time_range) for path in paths])