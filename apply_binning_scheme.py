"""
apply_binning_scheme.py

This code applies the given binning scheme to a set of data, and generates root files cotaining the events in each bin.
The year of interest, size of the data, polarity and meson to be analysed must be specified using the required flags --year --size --polarity --meson. There also are the flags --input --path and --bin_path, which are not required. These are used to specify the directory where the input data is located, where the binning scheme can be found and where the output file should be written, respectively. By default it is set to be the current working directory.
It outputs the root files with the events in each individual bin, as well as a txt file with the number of events in each bin.

Author: Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk)
Last edited: 15th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #
import os
import argparse
import numpy as np
import uproot
import pandas as pd
import awkward as ak


# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size      Used to specify the amount of events the user is interested in analysing.
                The argument must be one of: [large, small, medium, 1-8]. The integers specify the number of root
                files to be read in. Large is equivalent to 8. Medium is equivalent to 4. Small takes 200000 events.
    --polarity  Used to specify the polarity of the magnet the user is interested in.
                The argument must be one of: [up, down].
    --meson     Used to specify the meson the user is interested in.
                The argument must be one of: [D0, D0bar, both].
    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --bin_path  Used to specify the directory in which the binning scheme should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    
    Returns the parsed arguments.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--year",
        type=int,
        choices=[16,17,18],
        required=True,
        help="flag to set the data taking year."
    )
    parser.add_argument(
        "--size",
        type=str,
        choices=["large", "medium", "small", "1", "2", "3", "4", "5", "6", "7", "8"],
        required=True,
        help="flag to set the data taking year."
    )
    parser.add_argument(
        "--polarity",
        type=str,
        choices=["up","down"],
        required=True,
        help="flag to set the data taking polarity."
    )
    parser.add_argument(
        "--meson",
        type=str,
        choices=["D0","D0bar","both"],
        required=True,
        help="flag to set the D0 meson flavour."
    )    
    parser.add_argument(
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--input",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found"
    )
    parser.add_argument(
        "--bin_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the binning scheme should be found"
    )
    return parser.parse_args()

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
        
# - - - - - - - MAIN BODY - - - - - - - #

args = parse_arguments()

# import data
tree_name = "D02Kpi_Tuple/DecayTree"
if args.meson=="both":
    data = uproot.concatenate(f"{args.input}/{args.polarity}_data_{args.year}_{args.size}_clean.root:{tree_name}")
else:
    data = uproot.concatenate(f"{args.input}/{args.meson}_{args.polarity}_data_{args.year}_{args.size}_clean.root:{tree_name}")

bins = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_bins.txt", delimiter=',')                          
# select particles with pT below 10 GeV/c
nevents=np.empty(0)
length = len(data["D0_PT"])

# iterate through all bins
for i in np.arange(0, 10):
    pT_mask = np.ones(length)
    pT_mask = np.logical_and(pT_mask, data["D0_PT"]>bins[0,i])
    pT_mask = np.logical_and(pT_mask, data["D0_PT"]<=bins[0,i+1])
    for j in np.arange(0,10):
        eta_mask = pT_mask
        eta_mask = np.logical_and(eta_mask, data["D0_ETA"]>bins[i+1,j])
        eta_mask = np.logical_and(eta_mask, data["D0_ETA"]<=bins[i+1,j+1])
        selected_data = data[eta_mask]
        nevents = np.append(nevents, len(selected_data["D0_PT"]))
        # Write out bin
        if args.meson=="both":
            out_file_name = f"{args.path}/{args.polarity}_{args.year}_{args.size}_bin{j}{i}.root"
        else:
            out_file_name = f"{args.path}/{args.meson}_{args.polarity}_{args.year}_{args.size}_bin{j}{i}.root"
        out_tree = "D02Kpi_Tuple/DecayTree"
        print(f"Writing to {out_file_name}...")
        out_file = uproot.recreate(out_file_name)
        branches = {column: ak.type(selected_data[column]) for column in selected_data.fields}
        out_file.mktree(out_tree, branches)
        out_file[out_tree].extend({branch: selected_data[branch] for branch in branches.keys()})
        out_file.close()

# write out number of events in each bin
np.savetxt(f"{args.bin_path}/number_of_events_{args.meson}_{args.polarity}_{args.year}_{args.size}.txt", nevents, delimiter=',')