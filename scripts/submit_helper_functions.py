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



WORKING_DIRECTORIES = []



def numbersListFromString(string):
    return re.findall('\d+', string)



# calculate number of threads needed when everything will be submitted
# if not agreed on, exit, else return true
def askPermission(args):
    print(askPermission.__name__)
    assert(args.t)
    assert(args.k)
    assert(args.g)
    assert(args.d)
    assert(args.it)
    jobs = len(args.t)*len(args.k)*len(args.g)*len(args.d)*args.it
    time.sleep(.1)
    cont = input("About to submit "+str(jobs) + " Jobs (" +str(jobs*args.threads) + " threads), continue? [Y/n]  ")
    if cont == 'y' or cont == 'Y': 
        print("continue submitting")
        print()
        return True
    else:
        print("not submitting")
        print()
        sys.exit()



# hirarchy of directory structure
# -temperature
# |-kappa
#  |-gamma
#   |-density
#    |-iteration
#
# function will get path and the 4 parameters
def applyFunctionToDirectoryTree(args, func):
    assert(args.dir)
    # start_dir = os.getcwd()
    root_dir = args.dir

    # temperature level
    for t in args.t:
        t_path = os.path.join(root_dir, "T"+str(t))
        
        # kappa level
        for k in args.k:
            k_path = os.path.join(t_path, "kappa"+str(k))
            
            # gamma level
            for g in args.g:
                g_path = os.path.join(k_path, "gamma"+str(g))
                
                # density level
                for d in args.d:
                    d_path = os.path.join(g_path, "dens"+str(d))

                    for it in range(args.it):
                        it_path = os.path.join(d_path,"simulation"+str(it))
                        func(args, it_path, t, k ,g, d)



# create the directory tree if not already existing
# when iteration of simulation parameters already exists
#     increase iteration counter
# safe created directories in WORKING_DIRECTORIES
def createDirectoryTree(args):
    print(createDirectoryTree.__name__)
    def createDirectory(args, path, *prms):
        # path without iteration number
        if os.path.exists(path) and args.forcenew:
            shutil.rmtree(path)
        basic_path = path.rstrip(string.digits)
        # only iteration number
        this_iteration = int(path.rsplit("simulation",maxsplit=1)[1])
        while os.path.exists(path):
            this_iteration += 1
            print(path, "already exists, ", end='')
            path = basic_path+str(this_iteration)
            print("try again with", path )
        assert(not os.path.exists(path))
        print("create", path)
        os.makedirs(path)
        WORKING_DIRECTORIES.append(path)
    applyFunctionToDirectoryTree(args, createDirectory)
    print()



# iterate a file linewise
# check line for keyword
# if keyword is found replace WHOLE line
def fileReplaceLineWithKeyword(filepath, keyword, replacement):
    replacement_counter = 0
    assert(filepath)
    for line in fileinput.input(filepath, inplace=1):
        if keyword in line:
            line = replacement
            replacement_counter += 1
            print(line, end='\n')
        else:
            print(line, end='')
    if replacement_counter == 0:
        raise Exception("unable to find keyword "+keyword+" in file "+filepath+"\n")



def stripParametersFromPath(path):
    T = re.findall(r'T.*?([0-9.-]+)',path)[-1]
    k = re.findall(r'kappa.*?([0-9.-]+)',path)[-1]
    g = re.findall(r'gamma.*?([0-9.-]+)',path)[-1]
    d = re.findall(r'dens.*?([0-9.-]+)',path)[-1]
    it = re.findall(r'simulation.*?([0-9.-]+)',path)[-1]
    return float(T), float(k), float(g), float(d), int(it)



# copy the given config file to directory tree
def copyConfigFile(args):
    print(copyConfigFile.__name__)
    assert(os.path.exists(args.config))
    for dir in WORKING_DIRECTORIES:
        assert(os.path.exists(dir))
        new_config_file = os.path.join(dir,"config.ini")
        print("copy ", args.config, " to ", new_config_file)
        shutil.copy2(args.config, new_config_file)
    print()



# copy the given config file to directory tree
def copyPrograms(args):
    print(copyPrograms.__name__)
    assert(os.path.exists(args.prog))
    assert(os.path.exists(args.analysis))
    for dir in WORKING_DIRECTORIES:
        assert(os.path.exists(dir))
        print("copy ", args.prog, " to ", dir)
        shutil.copy2(args.prog, dir)
        if len(args.analysis) > 0:
            print("copy ", args.analysis, " to ", dir)
            shutil.copy2(args.analysis, dir)
    print()



# change the relevant parameters in the config files
def updateConfigFiles(args):
    print(updateConfigFiles.__name__)
    for dir in WORKING_DIRECTORIES:
        assert(os.path.exists(dir))
        assert(os.path.exists(args.config))
        t,k,g,d,it = stripParametersFromPath(dir)
        new_config_file = os.path.join(dir,"config.ini")
        print("change  [T, kappa, gamma, dens]  to ", [t,k,g,d], " in file ", new_config_file)
        fileReplaceLineWithKeyword(new_config_file, "mobile", "mobile="+str(args.mobile))
        fileReplaceLineWithKeyword(new_config_file, "cpu_threads", "cpu_threads="+str(args.threads))
        fileReplaceLineWithKeyword(new_config_file, "time_max", "time_max="+str(args.maxtime))
        fileReplaceLineWithKeyword(new_config_file, "temperature", "temperature="+str(t))
        fileReplaceLineWithKeyword(new_config_file, "kappa", "kappa="+str(k))
        fileReplaceLineWithKeyword(new_config_file, "gamma", "gamma="+str(g))
        fileReplaceLineWithKeyword(new_config_file, "density", "density="+str(d))
    print()



# create submit script in directory
def createSubmitScripts(args):
    print(createSubmitScripts.__name__)
    for dir in WORKING_DIRECTORIES:
        filepath = os.path.join(dir, "submit.sh")
        print("create ", filepath)
        with open(filepath, "w") as slurmfile:
            print("#!/bin/bash", file=slurmfile)
            if args.group:
                print("#SBATCH -A {}".format(args.group), file=slurmfile)
            if args.mail:
                print("#SBATCH --mail-type=FAIL", file=slurmfile)
                print("#SBATCH --mail-user={}".format(args.mail), file=slurmfile)
            if args.avx:
                print("#SBATCH --constraint=\"avx\"", file=slurmfile)
            print("#SBATCH --ntasks=1", file=slurmfile)
            print("#SBATCH --nodes=1", file=slurmfile)
            print("#SBATCH --cpus-per-task={}".format(args.threads), file=slurmfile)
            print("#SBATCH --mem={}G".format(args.memory), file=slurmfile)
            hours = args.hours
            print("#SBATCH -p {}".format("long" if hours>48 else "short"), file=slurmfile)
            days, hours = divmod(hours, 24)
            print("#SBATCH --time={0:0>1}-{1:0>2}:00:00".format(days,hours), file=slurmfile)
            print("", file=slurmfile)
            print("srun $@", file=slurmfile)
    print()



# reque with dependency
def sbatchAnalysis(args,dir,jobname,jobnum,dependency="afterany"):
    old_config_file = os.path.join(dir,"config.ini")
    assert(os.path.exists(old_config_file))
    new_config_file = os.path.join(dir,"config_analysis.ini")
    shutil.copy2(old_config_file, new_config_file)
    assert(os.path.exists(new_config_file))
    program = args.analysis.rsplit("/",maxsplit=1)[1]
    command = "sbatch --dependency="+dependency+":"+jobnum+" -J \"ana_"+jobname+"\" submit.sh "+program+" --config config_analysis.ini"
    if shutil.which("sbatch") != None:
        status, jobnum = subprocess.getstatusoutput(command)
        jobnum = numbersListFromString(jobnum)[-1]
    else:
        status, jobnum = "DRYRUN", numbersListFromString("DRYRUN 1337")[-1]
    print(command)
    return jobnum



# reque with dependency
def sbatchRepeat(args,dir,jobname,jobnum,dependency="afterany"):
    old_config_file = os.path.join(dir,"config.ini")
    assert(os.path.exists(old_config_file))
    new_config_file = os.path.join(dir,"config_repeat.ini")
    shutil.copy2(old_config_file, new_config_file)
    assert(os.path.exists(new_config_file))
    try:
        fileReplaceLineWithKeyword(new_config_file, "traj=none", "traj="+str("gro"))
    except Exception as e:
        pass
    program = args.prog.rsplit("/",maxsplit=1)[1]
    for i in range(args.repeat):
        command = "sbatch --dependency="+dependency+":"+jobnum+" -J \"re"+str(i)+"_"+jobname+"\" submit.sh "+program+" --config config_repeat.ini"
        if shutil.which("sbatch") != None:
            status, jobnum = subprocess.getstatusoutput(command)
            jobnum = numbersListFromString(jobnum)[-1]
        else:
            status, jobnum = "DRYRUN", numbersListFromString("DRYRUN 137")[-1]
        print(command)
    return jobnum



# run sbatch on every simulation iteration
def sbatchAll(args):
    print(sbatchAll.__name__)
    os.chdir(args.origin)
    if shutil.which("sbatch") == None:
        print("WARNING: sbatch not available. making a dry run")
    for dir in WORKING_DIRECTORIES:
        dir = os.path.join(args.origin, dir)
        os.chdir(dir)
        program = args.prog.rsplit("/",maxsplit=1)[1]
        T,k,g,d,it = stripParametersFromPath(dir)
        name = "T"+str(T)+"kappa"+str(k)+"gamma"+str(g)+"dens"+str(d)+"it"+str(it)
        command = "sbatch -J \""+name+"\" submit.sh "+program+" --config config.ini"
        print(command)
        status, jobnum = None, None
        
        if shutil.which("sbatch") != None:
            status, jobnum = subprocess.getstatusoutput(command)
            jobnum = numbersListFromString(jobnum)[-1]
        else:
            status, jobnum = "DRYRUN", numbersListFromString("DRYRUN 17")[-1]

        if args.repeat > 0:
            jobnum = sbatchRepeat(args,dir,name,jobnum)
        if args.analysis != None:
            jobnum = sbatchAnalysis(args,dir,name,jobnum)
        os.chdir(args.origin)
    print()