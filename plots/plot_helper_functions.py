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



# return all dicts wich match 'matches'
# matches mus be dict
def getMatched(filepath, matches):
    chosen_dicts = []
    for paramdict in readJson(filepath):
        if all( paramdict[param] == value for param, value in matches.items() ):
            chosen_dicts.append(paramdict)
    return chosen_dicts



# return all values for 'param'
# in dicts that match all in 'matches'
# matches mus be dict
def getMatchedValues(filepath, param, matches):
    values = []
    for paramdict in getMatched(filepath, matches):
        values.append(paramdict[param])
    return sorted(list(set(values)))



# return all values for 'dir_path' 
# in dicts that match all in 'matches'
# matches mus be dict
def getMatchedDirs(filepath, matches):
    return getMatchedValues(filepath, "dir_path", matches)