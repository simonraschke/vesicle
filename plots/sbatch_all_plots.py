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
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--file", type=str, help="path to config_files.json file")
parser.add_argument("--time", type=float, nargs=2, default=[0,1e10], help="path to config_files.json file")
args = parser.parse_args()

threads_ = 8
memory_ = 8
hours_ = 12
dir_of_this_script = os.path.dirname(os.path.abspath(__file__))

filepath = os.path.join(args.origin, "submit_plot.sh")
print("create ", filepath)
print()
with open(filepath, "w") as slurmfile:
    print("#!/bin/bash", file=slurmfile)
    print("#SBATCH --constraint=\"avx\"", file=slurmfile)
    print("#SBATCH --ntasks=1", file=slurmfile)
    print("#SBATCH --nodes=1", file=slurmfile)
    print("#SBATCH --cpus-per-task={}".format(threads_), file=slurmfile)
    print("#SBATCH --mem={}G".format(memory_), file=slurmfile)
    hours = hours_
    print("#SBATCH -p {}".format("long" if hours>48 else "short"), file=slurmfile)
    days, hours = divmod(hours, 24)
    print("#SBATCH --time={0:0>1}-{1:0>2}:00:00".format(days,hours), file=slurmfile)
    print("#SBATCH --signal=2@300", file=slurmfile)
    print("", file=slurmfile)
    print("srun $@", file=slurmfile)


rho_free_command = "sbatch -J \"PLOT rho_free\" submit_plot.sh " + os.path.join(dir_of_this_script,"rho_free.py") + " --file " + args.file + " --time " + str(args.time[0]) + " " + str(args.time[1])
if shutil.which("sbatch") != None:
    subprocess.getstatusoutput(rho_free_command)
print(rho_free_command)
print()

order_overall_command = "sbatch -J \"PLOT order_overall\" submit_plot.sh " + os.path.join(dir_of_this_script,"order_overall.py") + " --file " + args.file + " --time " + str(args.time[0]) + " " + str(args.time[1])
if shutil.which("sbatch") != None:
    subprocess.getstatusoutput(order_overall_command)
print(order_overall_command)
print()