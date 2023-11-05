"""
multiple_candidates.py

This code is used to chech which events have multiple candidates, and for those that do one of the candidates is selected at random, while the others are removed.
The year of interest, size of the data and polarity to be analysed must be specified using the required flags --year --size --polarity. There is a fourth flag --path, which is not required. This one is used to specify the directory where the input data is located, and where the output file should be written. By default it is set to be the current working directory.
It outputs the data for each into 3 seperate root files, one containing only D0 events, another containing only D0bar, and the third containing all events.
This code is heavily inspired on the work of Camille Jarvis-Stiggants and Michael England. The minor modifications to the original code have simply added flexibility to it.

Author: Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk)
Last edited: 15th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import os
import uproot as ur
import pandas as pd
import awkward as ak
import numpy as np
import argparse

# - - - - - - - FUNCTIONS - - - - - - - #

def get_multiple_candidate_selection(data_frame):
    """
    Add a boolean column of whether to reject a multiple candidate based on run & event number
    Note we randomly shuffle the columns, then the "duplicates" function returns false for the
    first occurence, then true for subsequent. So we return this as the "should we reject the event"
    decision
    """
    internal_clone = data_frame
    # get the old index so we can resort after sampling
    
    columns = ["eventNumber", "runNumber"]
    test = pd.DataFrame(data=ak.to_numpy(internal_clone[columns]),columns=columns) # for checking that the order of events is the same before shuffling and after resorting
    print("Original dataframe:\n", test.head(), "\n")
    
    # make a shuffled dataset
    # df_shuffle = internal_clone.sample(frac=1).reset_index() # frac=x means sample 100*x % of the dataset randomly, hence frac=1 means shuffle all of it
    df_shuffle = test.sample(frac=1).reset_index() # frac=x means sample 100*x % of the dataset randomly, hence frac=1 means shuffle all of it
    
    # the .duplicated() function returns a list of booleans
    # in this case it compares the eventNumber and runNumber of all entries in the tree
    df_shuffle["is_mult_cand"] = df_shuffle.duplicated(columns)
    # whatever is inside the [""] is the name of the new branch in the output file
    # the first instance of any combination of eventNumber and runNumber will have is_mult_cand=False
    # this is because it hasn"t come across this combination yet
    # all subsequent instances of the same combination of numbers will have is_mult_cand=True
    # shuffling first ensures that the first instance (that which gets False) is random
    # which is what we want (randomly pick a candidate to keep)
    # then when reading the output file for analysis, only keep events if(is_mult_cand==false)

    # resort
    df_shuffle = df_shuffle.sort_values(by=["index"]).reset_index()
    
    # shorten to only contain the following three columns (will add this new file as a friend)
    df_shuffle = df_shuffle[["eventNumber","runNumber","is_mult_cand"]]
    
    # should now put this in the original dataframe
    return df_shuffle

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
        
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
    --path      Used to specify the directory in which the output files should be written. It is not required,
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
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    
    return parser.parse_args()
        
def save_file(filename, cut_data):
    '''
    Saves a .root file containing the data in cut_data. It is written to the path specified
    by the user, and has the name given by filename.
    '''
    tree = "D02Kpi_Tuple/DecayTree"

    filename = f"{args.path}/{filename}.root" 
    print(f"Writing to {filename}...")
    outfile = ur.recreate(f"{filename}") # create the write file

    branches = {column: ak.type(cut_data[column]) for column in cut_data.fields} # here we're defining the branches we want to save to file as a dictionary of name : data type
    outfile.mktree(tree, branches)

    # now write and close
    outfile[tree].extend({branch: cut_data[branch] for branch in branches.keys()})
    outfile.close()

    return print(f'Saved file {filename}.')

def save_all(year, size):
    '''
    Iterates through the save_file function in order to write out all the data in different
    files, based on the origin of the decay (D0, D0bar or any).
    '''
    names = ["D0_", "D0bar_", ""]
    for index, dataset in enumerate([D0_DATA, D0BAR_DATA, DATA]):
        save_file(f'{names[index]}{args.polarity}_data_{year}_{size}_clean', dataset)
    
    return print(f'Saved {size} data for year 20{year}')

def split_meson(data):
    '''
    It takes data from events and splits them in 3 datasets, depending on the origin meson
    of the decay (D0, D0bar, any).
    '''
    # select positive z-momentum and remove muons
    length = len(data)
    print(f"The number of events to be analysed is {length}")
    mask = np.ones(length)
    
    # split the reconstructed particles into mesons and antimesons
    
    mask_D0 = np.logical_and(mask, data["D0_ID"]==421)
    mask_D0bar = np.logical_and(mask, data["D0_ID"]==-421)
    
    D0_data = data[mask_D0]
    D0bar_data = data[mask_D0bar]
    
    print('checkpoint: data has been masked')
        
    return D0_data, D0bar_data, data

# - - - - - - - MAIN BODY - - - - - - - #

args=parse_arguments()

np.random.seed(482022) # implemented 04/08/2022

# Import data
tree_name = "D02Kpi_Tuple/DecayTree"
data = ur.concatenate(f"{args.path}/{args.polarity}_data_{args.year}_{args.size}.root:{tree_name}")

print(f"reading file for year 20{args.year}...")

print(len(data))
if len(data) > 0:
    df = get_multiple_candidate_selection(data) # select at random candidates to remove

    print("Got dataframe:\n", df.head(), "\nCheck it matches with above!!")

    print(f"Number entries before multiple candidate cut = {len(data)}")
    is_a_multiple_candidate = df["is_mult_cand"].to_numpy()
    data = data[~is_a_multiple_candidate] # remove candidates
    print(f"Number entries after multiple candidate cut = {len(data)}")
    
    D0_DATA, D0BAR_DATA, DATA = split_meson(data) # split data based on meson
    
    save_all(args.year, args.size) # output data