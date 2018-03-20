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

# SLURM submit helper script for vesicle

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
parser.add_argument("--repeat", type=int, default=0, help="restart args.prog x times")
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--dir", type=str, help="the working directory")
parser.add_argument("--config", type=str, default='bin/vesicle_config.ini', help="program config file path from ~")
parser.add_argument("--prog", type=str, default='vesicle', help="program name")
parser.add_argument("--analysis", type=str, default='vesicle_analysis', help="analysis program name")
parser.add_argument('--mobile', type=int, default=1000, help="mobile particles")
parser.add_argument('--maxtime', type=float, default=1e8, help="maximum overall simulation time")
parser.add_argument('-t', nargs='*', type=float, help="temperature value list")
parser.add_argument('-k', nargs='*', type=float, help="kappa value list")
parser.add_argument('-g', nargs='*', type=float, help="gamma value list")
parser.add_argument('-d', nargs='*', type=float, help="density value list")
parser.add_argument('--it', type=int, default=1, help="iterations per parameter set")
parser.add_argument('--threads', type=int, default=4, help="[SLURM] number of threads per job")
parser.add_argument('--memory', type=int, default=4, help="[SLURM] amount of RAM per job in GB")
parser.add_argument('--hours', type=int, default=48, help="[SLURM] number of hours for submit script")
parser.add_argument('--mail', type=str, help="[SLURM] the mail to send fail news to")
parser.add_argument('--group', type=str, help="[SLURM] group membership")
parser.add_argument('--avx', action='store_true', default=True, help="[SLURM] set avx flag")
args = parser.parse_args()
args.origin = os.getcwd() # save this directory
args.config = os.path.join(args.origin, args.config)
args.prog = subprocess.getstatusoutput("which "+args.prog)[1]
args.analysis = subprocess.getstatusoutput("which "+args.analysis)[1]

args.t = sorted(list(set(args.t)))
args.k = sorted(list(set(args.k)))
args.g = sorted(list(set(args.g)))
args.d = sorted(list(set(args.d)))



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

    time.sleep(.1)
    prog_right = input("Is this the right analysis program path? "+str(args.analysis)+"  [Y/n]  ")
    if prog_right == 'y' or prog_right == 'Y':
        print("continue")
        print()
    else:
        print("aborting")
        print()

    sb.askPermission(args)
    sb.createDirectoryTree(args)
    sb.copyConfigFile(args)
    sb.copyPrograms(args)
    sb.updateConfigFiles(args)
    sb.createSubmitScripts(args)
    sb.sbatchAll(args)