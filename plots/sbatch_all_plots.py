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
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--file", type=str, help="path to config_files.json file")
parser.add_argument("--time", type=float, nargs=2, default=[0,1e10], help="path to config_files.json file")
args = parser.parse_args()

threads_ = 10
memory_ = 10
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


rho_free_args = os.path.join(dir_of_this_script,"rho_free.py") + " --file " + args.file + " --time " + str(args.time[0]) + " " + str(args.time[1])
rho_free_command = "sbatch -J \"PLOT rho_free\" submit_plot.sh " + rho_free_args

order_overall_args = os.path.join(dir_of_this_script,"order_full.py") + " --file " + args.file + " --time " + str(args.time[0]) + " " + str(args.time[1])
order_overall_command = "sbatch -J \"PLOT order_overall\" submit_plot.sh " + order_overall_args

n_avg_args = os.path.join(dir_of_this_script,"N_avg.py") + " --file " + args.file + " --time " + str(args.time[0]) + " " + str(args.time[1]) + " --min 5"
n_avg_command = "sbatch -J \"PLOT N_avg\" submit_plot.sh " + n_avg_args



rho_free_time_T026_args0 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.26 --dens 0.01 0.03 --out rho_free_time_T026_0"
rho_free_time_T026_command0 = "sbatch -J \"PLOT rho_free_time_T026\" submit_plot.sh " + rho_free_time_T026_args0

rho_free_time_T026_args1 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.26 --dens 0.04 0.06 --out rho_free_time_T026_1"
rho_free_time_T026_command1 = "sbatch -J \"PLOT rho_free_time_T026\" submit_plot.sh " + rho_free_time_T026_args1

rho_free_time_T026_args2 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.26 --dens 0.07 0.1 --out rho_free_time_T026_2"
rho_free_time_T026_command2 = "sbatch -J \"PLOT rho_free_time_T026\" submit_plot.sh " + rho_free_time_T026_args2



rho_free_time_T027_args0 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.27 --dens 0.001 0.01 --out rho_free_time_T027_0"
rho_free_time_T027_command0 = "sbatch -J \"PLOT rho_free_time_T027\" submit_plot.sh " + rho_free_time_T027_args0

rho_free_time_T027_args1 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.27 --dens 0.015 0.04 --out rho_free_time_T027_1"
rho_free_time_T027_command1 = "sbatch -J \"PLOT rho_free_time_T027\" submit_plot.sh " + rho_free_time_T027_args1

rho_free_time_T027_args2 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.27 --dens 0.05 0.1 --out rho_free_time_T027_2"
rho_free_time_T027_command2 = "sbatch -J \"PLOT rho_free_time_T027\" submit_plot.sh " + rho_free_time_T027_args2



rho_free_time_T028_args0 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.28 --dens 0.001 0.015 --out rho_free_time_T028_0"
rho_free_time_T028_command0 = "sbatch -J \"PLOT rho_free_time_T028\" submit_plot.sh " + rho_free_time_T028_args0


rho_free_time_T028_args1 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.28 --dens 0.02 0.05 --out rho_free_time_T028_1"
rho_free_time_T028_command1 = "sbatch -J \"PLOT rho_free_time_T028\" submit_plot.sh " + rho_free_time_T028_args1


rho_free_time_T028_args2 = os.path.join(dir_of_this_script,"rho_free_time.py") + " --file " + args.file + " --con temperature 0.28 --dens 0.06 0.1 --out rho_free_time_T028_2"
rho_free_time_T028_command2 = "sbatch -J \"PLOT rho_free_time_T028\" submit_plot.sh " + rho_free_time_T028_args2



N_max_time_T026_args0 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.26 --dens 0.01 0.03 --out N_max_time_T026_0"
N_max_time_T026_command0 = "sbatch -J \"PLOT N_max_time_T026\" submit_plot.sh " + N_max_time_T026_args0

N_max_time_T026_args1 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.26 --dens 0.04 0.06 --out N_max_time_T026_1"
N_max_time_T026_command1 = "sbatch -J \"PLOT N_max_time_T026\" submit_plot.sh " + N_max_time_T026_args1

N_max_time_T026_args2 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.26 --dens 0.07 0.1 --out N_max_time_T026_2"
N_max_time_T026_command2 = "sbatch -J \"PLOT N_max_time_T026\" submit_plot.sh " + N_max_time_T026_args2



N_max_time_T027_args0 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.27 --dens 0.001 0.01 --out N_max_time_T027_0"
N_max_time_T027_command0 = "sbatch -J \"PLOT N_max_time_T027\" submit_plot.sh " + N_max_time_T027_args0

N_max_time_T027_args1 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.27 --dens 0.015 0.04 --out N_max_time_T027_1"
N_max_time_T027_command1 = "sbatch -J \"PLOT N_max_time_T027\" submit_plot.sh " + N_max_time_T027_args1

N_max_time_T027_args2 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.27 --dens 0.05 0.1 --out N_max_time_T027_2"
N_max_time_T027_command2 = "sbatch -J \"PLOT N_max_time_T027\" submit_plot.sh " + N_max_time_T027_args2



N_max_time_T028_args0 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.28 --dens 0.001 0.015 --out N_max_time_T028_0"
N_max_time_T028_command0 = "sbatch -J \"PLOT N_max_time_T028\" submit_plot.sh " + N_max_time_T028_args0


N_max_time_T028_args1 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.28 --dens 0.02 0.05 --out N_max_time_T028_1"
N_max_time_T028_command1 = "sbatch -J \"PLOT N_max_time_T028\" submit_plot.sh " + N_max_time_T028_args1


N_max_time_T028_args2 = os.path.join(dir_of_this_script,"N_max_time.py") + " --file " + args.file + " --con temperature 0.28 --dens 0.06 0.1 --out N_max_time_T028_2"
N_max_time_T028_command2 = "sbatch -J \"PLOT N_max_time_T028\" submit_plot.sh " + N_max_time_T028_args2

if shutil.which("sbatch") != None:
    print(order_overall_command)
    subprocess.getstatusoutput(order_overall_command)
    print(rho_free_command)
    subprocess.getstatusoutput(rho_free_command)
    print(n_avg_command)
    subprocess.getstatusoutput(n_avg_command)
    print(rho_free_time_T026_command0)
    subprocess.getstatusoutput(rho_free_time_T026_command0)
    print(rho_free_time_T026_command1)
    subprocess.getstatusoutput(rho_free_time_T026_command1)
    print(rho_free_time_T026_command2)
    subprocess.getstatusoutput(rho_free_time_T026_command2)
    print(rho_free_time_T027_command0)
    subprocess.getstatusoutput(rho_free_time_T027_command0)
    print(rho_free_time_T027_command1)
    subprocess.getstatusoutput(rho_free_time_T027_command1)
    print(rho_free_time_T027_command2)
    subprocess.getstatusoutput(rho_free_time_T027_command2)
    print(rho_free_time_T028_command0)
    subprocess.getstatusoutput(rho_free_time_T028_command0)
    print(rho_free_time_T028_command1)
    subprocess.getstatusoutput(rho_free_time_T028_command1)
    print(rho_free_time_T028_command2)
    subprocess.getstatusoutput(rho_free_time_T028_command2)
    print(N_max_time_T026_command0)
    subprocess.getstatusoutput(N_max_time_T026_command0)
    print(N_max_time_T026_command1)
    subprocess.getstatusoutput(N_max_time_T026_command1)
    print(N_max_time_T026_command2)
    subprocess.getstatusoutput(N_max_time_T026_command2)
    print(N_max_time_T027_command0)
    subprocess.getstatusoutput(N_max_time_T027_command0)
    print(N_max_time_T027_command1)
    subprocess.getstatusoutput(N_max_time_T027_command1)
    print(N_max_time_T027_command2)
    subprocess.getstatusoutput(N_max_time_T027_command2)
    print(N_max_time_T028_command0)
    subprocess.getstatusoutput(N_max_time_T028_command0)
    print(N_max_time_T028_command1)
    subprocess.getstatusoutput(N_max_time_T028_command1)
    print(N_max_time_T028_command2)
    subprocess.getstatusoutput(N_max_time_T028_command2)
else:
    time.sleep(.1)
    prog_right = input("sbatch not found. run on this machine?   [Y/n]  ")
    if prog_right == 'y' or prog_right == 'Y':
        print("continue")
        print()
    else:
        print("aborting")
        print()
        sys.exit()
    print(rho_free_args)
    subprocess.getstatusoutput(rho_free_args)
    print(order_overall_args)
    subprocess.getstatusoutput(order_overall_args)
    print(n_avg_args)
    subprocess.getstatusoutput(n_avg_args)
    print(rho_free_time_T026_args0)
    subprocess.getstatusoutput(rho_free_time_T026_args0)
    print(rho_free_time_T026_args1)
    subprocess.getstatusoutput(rho_free_time_T026_args1)
    print(rho_free_time_T026_args2)
    subprocess.getstatusoutput(rho_free_time_T026_args2)
    print(rho_free_time_T027_args0)
    subprocess.getstatusoutput(rho_free_time_T027_args0)
    print(rho_free_time_T027_args1)
    subprocess.getstatusoutput(rho_free_time_T027_args1)
    print(rho_free_time_T027_args2)
    subprocess.getstatusoutput(rho_free_time_T027_args2)
    print(rho_free_time_T028_args0)
    subprocess.getstatusoutput(rho_free_time_T028_args0)
    print(rho_free_time_T028_args1)
    subprocess.getstatusoutput(rho_free_time_T028_args1)
    print(rho_free_time_T028_args2)
    subprocess.getstatusoutput(rho_free_time_T028_args2)
    print(N_max_time_T026_args0)
    subprocess.getstatusoutput(N_max_time_T026_args0)
    print(N_max_time_T026_args1)
    subprocess.getstatusoutput(N_max_time_T026_args1)
    print(N_max_time_T026_args2)
    subprocess.getstatusoutput(N_max_time_T026_args2)
    print(N_max_time_T027_args0)
    subprocess.getstatusoutput(N_max_time_T027_args0)
    print(N_max_time_T027_args1)
    subprocess.getstatusoutput(N_max_time_T027_args1)
    print(N_max_time_T027_args2)
    subprocess.getstatusoutput(N_max_time_T027_args2)
    print(N_max_time_T028_args0)
    subprocess.getstatusoutput(N_max_time_T028_args0)
    print(N_max_time_T028_args1)
    subprocess.getstatusoutput(N_max_time_T028_args1)
    print(N_max_time_T028_args2)
    subprocess.getstatusoutput(N_max_time_T028_args2)