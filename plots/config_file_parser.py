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
import json
import argparse
import pprint
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
import submit_helper_functions as sb
import plot_helper_functions as plthelp

pp = pprint.PrettyPrinter(indent=4, compact=True)



parser = argparse.ArgumentParser()
parser.add_argument("--origin", type=str, default=os.getcwd(), help="this directory")
parser.add_argument("--dir", type=str, help="the working directory")
parser.add_argument("--filename", type=str, default="config.ini", help="file name of config file")
parser.add_argument("--depth", type=int, default=5, help="depth to walk [dir] recursively")
args = parser.parse_args()



if __name__ == "__main__":
    # print all args once
    for arg in vars(args):
        print("{0:10}".format(arg), getattr(args, arg))
    print()

    config_files_parser_path = os.path.join(args.dir,"config_files.json")#
    list_of_config_dicts = []

    # check the directory tree for config files
    for root,dirs,files in sb.walklevel(args.dir,args.depth):
        for dir in dirs:
            config_file_path = os.path.join(root,dir,args.filename)
            if not os.path.exists(config_file_path):
                continue
            elif not os.path.getsize(config_file_path) > 0:
                continue
            else:
                print("found", args.filename, "in", os.path.join(root,dir))
                parameters = {"dir_path":os.path.abspath(os.path.join(root,dir))}
                parameters.update( {"algorithm":plthelp.fileValueFromKeyword(config_file_path, "algorithm")} )
                parameters.update( {"acceptance":plthelp.fileValueFromKeyword(config_file_path, "acceptance")} )
                parameters.update( {"thermostat":plthelp.fileValueFromKeyword(config_file_path, "thermostat")} )
                parameters.update( {"temperature":plthelp.fileValueFromKeyword(config_file_path, "temperature")} )
                parameters.update( {"mobile":plthelp.fileValueFromKeyword(config_file_path, "mobile")} )
                parameters.update( {"density":plthelp.fileValueFromKeyword(config_file_path, "density")} )
                parameters.update( {"frame_guides_grid_edge":plthelp.fileValueFromKeyword(config_file_path, "frame_guides_grid_edge")} )
                parameters.update( {"guiding_elements_each":plthelp.fileValueFromKeyword(config_file_path, "guiding_elements_each")} )
                # parameters.update( {"box.x":plthelp.fileValueFromKeyword(config_file_path, "box.x")} )
                # parameters.update( {"box.y":plthelp.fileValueFromKeyword(config_file_path, "box.y")} )
                # parameters.update( {"box.z":plthelp.fileValueFromKeyword(config_file_path, "box.z")} )
                parameters.update( {"timestep":plthelp.fileValueFromKeyword(config_file_path, "timestep")} )
                parameters.update( {"kappa":plthelp.fileValueFromKeyword(config_file_path, "kappa")} )
                parameters.update( {"gamma":plthelp.fileValueFromKeyword(config_file_path, "gamma")} )
                parameters.update( {"sw_position_min":plthelp.fileValueFromKeyword(config_file_path, "sw_position_min")} )
                parameters.update( {"sw_position_max":plthelp.fileValueFromKeyword(config_file_path, "sw_position_max")} )
                parameters.update( {"sw_position_target":plthelp.fileValueFromKeyword(config_file_path, "sw_position_target")} )
                parameters.update( {"sw_orientation_min":plthelp.fileValueFromKeyword(config_file_path, "sw_orientation_min")} )
                parameters.update( {"sw_orientation_max":plthelp.fileValueFromKeyword(config_file_path, "sw_orientation_max")} )
                parameters.update( {"sw_orientation_target":plthelp.fileValueFromKeyword(config_file_path, "sw_orientation_target")} )
                parameters.update( {"cell_min_edge":plthelp.fileValueFromKeyword(config_file_path, "cell_min_edge")} )
                parameters.update( {"max_cells_dim":plthelp.fileValueFromKeyword(config_file_path, "max_cells_dim")} )
                parameters.update( {"cluster_minimum_size":plthelp.fileValueFromKeyword(config_file_path, "cluster_minimum_size")} )
                parameters.update( {"cluster_significant_size":plthelp.fileValueFromKeyword(config_file_path, "cluster_significant_size")} )
                parameters.update( {"cluster_distance_threshold":plthelp.fileValueFromKeyword(config_file_path, "cluster_distance_threshold")} )
                parameters.update( {"cluster_volume_extension":plthelp.fileValueFromKeyword(config_file_path, "cluster_volume_extension")} )
                list_of_config_dicts += [parameters]
    plthelp.writeJson(config_files_parser_path, list_of_config_dicts)
    # pp.pprint(plthelp.readJson(config_files_parser_path))
    # pp.pprint(plthelp.getValues(config_files_parser_path, "density"))
    # pp.pprint(plthelp.getDicts(config_files_parser_path))
    
