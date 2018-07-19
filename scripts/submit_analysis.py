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

# SLURM submit helper script for vesicle_analysis

import os
import sys
import numpy as np
import argparse
import pprint
import subprocess
import time
import submit_helper_functions as sb


pp = pprint.PrettyPrinter(indent=4, compact=True)


parser = argparse.ArgumentParser()
parser.add_argument("--forcenew", action='store_true', help="force makeing directories new instead of increasing it")
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--dir", type=str, help="the working directory")
parser.add_argument("--filename", type=str, help="file name of trajectory file")
parser.add_argument("--depth", type=int, default=5, help="depth to walk [dir] recursively")
parser.add_argument("--config", type=str, default=None, help="program config file path from ~")
parser.add_argument("--prog", type=str, default='vesicle_analysis', help="name of program")
parser.add_argument("--analysis", type=str, default='vesicle_analysis', help="SAME AS prog! analysis program name")
parser.add_argument("--cl_ext", type=float, default=0.7, help="[CONFIG] cluster volume extension")
parser.add_argument('--threads', type=int, default=8, help="[SLURM] number of threads per job")
parser.add_argument('--memory', type=int, default=20, help="[SLURM] amount of RAM per job in GB")
parser.add_argument('--hours', type=int, default=24, help="[SLURM] number of hours for submit script")
parser.add_argument('--minhours', type=int, default=None, help="[SLURM] number of hours for minimum job runtime")
parser.add_argument('--mail', type=str, help="[SLURM] the mail to send fail news to")
parser.add_argument('--group', type=str, help="[SLURM] group membership")
parser.add_argument('--avx', action='store_true', default=True, help="[SLURM] set avx flag")
parser.add_argument('--nice',type=int, default=0, help="[SLURM] nice value of jobs")
args = parser.parse_args()
args.origin = os.getcwd() # save this directory
if args.config != None:
    args.config = os.path.join(args.origin, args.config)
args.prog = subprocess.getstatusoutput("which "+args.prog)[1]

if __name__ == "__main__":
    # print all args once
    for arg in vars(args):
        print("{0:10}".format(arg), getattr(args, arg))
    print()

    time.sleep(.1)
    prog_right = input("Is this the right program path? "+str(args.prog)+"  [Y/n]  ")
    if prog_right == 'y' or prog_right == 'Y':
        print("continue")
        print()
    else:
        print("aborting")
        print()

    sb.makeDirectoryListAnalysis(args)
    sb.askPermissionAnalysis(args)
    sb.copyConfigFileAnalysis(args)
    sb.copyProgramsAnalysis(args)
    sb.updateConfigFilesAnalysis(args)
    sb.createSubmitScripts(args)
    sb.sbatchAllAnalysis(args)