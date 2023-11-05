"""
selection_of_events.py

This code is used to read in data from D0/D0bar decays to two hadrons in LHCb in the years 2016, 2017, 2018. It then proceeds to select the events that meet a set of given requirements. Finally it outputs the selected events in 2 root files. One contain data using the up polarity of the magnet, and the other the down polarity.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size. There is a third flag --path, which is not required. This one is used to specify the directory where the output files should be written. By default it is set to save the files in the current working directory.
This code is heavily inspired on the work of Camille Jarvis-Stiggants and Michael England. Here the original code has been reorganized, some comments have been added, as well as minor features to add flexibility to the code.

Author: Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk)
Last edited: 15th September 2023
"""

# - - - - - - - IMPORT STATEMENTS - - - - - - - #

import argparse
import uproot
import numpy as np
import awkward as ak
import os

# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year  Used to specify the year at which the data was taken the user is interested in.
            The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size  Used to specify the amount of events the user is interested in analysing.
            The argument must be one of: [large, small, medium, 1-8]. The integers specify the number of root
            files to be read in. Large is equivalent to 8. Medium is equivalent to 4. Small takes 200000 events.
    --path  Used to specify the directory in which the output files should be written. It is not required,
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
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    
    return parser.parse_args()

def get_data():
    '''
    Reads in the root files containing the LHCb data of D0 decays into two hadrons. It takes into
    account the year and size requested by the user, and handles the different scenarios appropiately.
    
    It only reads a set of specified variables given in the main body of the code.
    
    Returns an array of size 2 with the data using the up polarity as the first argument, and the
    down polarity as the second argument.
    '''
    tree_name = "D02Kpi_Tuple/DecayTree"

    # Data from different years have different codes in their path 
    if options.year==16:
        code_up = "00172206"
        code_down = "00172204"
    elif options.year==17:
        code_up = "00172214"
        code_down = "00172212"
    elif options.year==18:
        code_up = "00172210"
        code_down = "00172208"

    # get path to directories of the requested year
    directory_up = f"/eos/lhcb/grid/prod/lhcb/LHCb/Collision{options.year}/CHARM_D02HH_DVNTUPLE.ROOT/{code_up}/0000/"
    directory_down = f"/eos/lhcb/grid/prod/lhcb/LHCb/Collision{options.year}/CHARM_D02HH_DVNTUPLE.ROOT/{code_down}/0000/"
    
    # set the number of files to concatenate depending on the size requested
    max_events = None
    if options.size=="small":
        data_to_concatenate = np.arange(1, 2, 1)
        max_events = 200000
    elif options.size=="medium":
        data_to_concatenate = np.arange(1, 5, 1)
    elif options.size=="large":
        data_to_concatenate = np.arange(1, 9, 1)
    elif 0<int(options.size)<9:
        data_to_concatenate = np.arange(1, int(options.size)+1, 1)
    
    # reads the data from the files requested by the user and concatanates it
    data_up = uproot.concatenate((f"{directory_up}/{code_up}_0000000{i}_1.charm_d02hh_dvntuple.root:{tree_name}" for i in data_to_concatenate), expressions=read_only_these_variables, max_num_elements=max_events)
    data_down = uproot.concatenate((f"{directory_down}/{code_down}_0000000{i}_1.charm_d02hh_dvntuple.root:{tree_name}" for i in (data_to_concatenate+1)), expressions=read_only_these_variables, max_num_elements=max_events)
    
    print('checkpoint: data has been read')
    
    return data_up, data_down
    
def cut_data(data):
    '''
    It takes data from events and makes a selection based on a set of requirements.
    
    The requirements are:
    - Both hadrons must have a positive momentum in the z direction
    - Pseudorapidity must be between 0 and 6
    - Transverse momentum must be between 0 GeV and 10 GeV
    - Particle must not be muons
    - P1_PIDK > 5
    - P2_PIDK < 0
    - D0_IPCHI2_OWNPV < 1
    
    Returns the data of the events that met the requirements. It does so in an array of size 3.
    The first argument contains events with the decay of the D0 meson, the second one contains events
    with the decay of the D0bar meson, and the third one contains all events that met the selection
    criteria.
    '''
    # select positive z-momentum and remove muons
    length = len(data)
    print(f"The number of events to be analysed is {length}")
    mask = np.ones(length)
    mask = np.logical_and(mask, data["P1_PZ"]>0)
    mask = np.logical_and(mask, data["P2_PZ"]>0)
    mask = np.logical_and(mask, data["P1_ETA"]>0)
    mask = np.logical_and(mask, data["P2_ETA"]>0)
    mask = np.logical_and(mask, data["P1_ETA"]<6)
    mask = np.logical_and(mask, data["P2_ETA"]<6)
    mask = np.logical_and(mask, data["P1_PT"]>0)
    mask = np.logical_and(mask, data["P2_PT"]>0)
    mask = np.logical_and(mask, data["P1_PT"]<10000)
    mask = np.logical_and(mask, data["P2_PT"]<10000)
    mask = np.logical_and(mask, data["P1_isMuon"]==0)
    mask = np.logical_and(mask, data["P2_isMuon"]==0)
    
    # apply selection criteria for PID and remove secondary decays
    
    mask = np.logical_and(mask, data["P1_PIDK"] > 5)
    mask = np.logical_and(mask, data["P2_PIDK"] < 0)
    mask = np.logical_and(mask, np.log(data["D0_IPCHI2_OWNPV"]) < 1)
    
    # split the reconstructed particles into mesons and antimesons
    
    data = data[mask]
    print('checkpoint: data has been masked')
    
    return data

def save_file(filename, cut_data):
    '''
    Saves a .root file containing the data in cut_data. It is written to the path specified
    by the user, and has the name given by filename.
    '''
    tree = "D02Kpi_Tuple/DecayTree"

    filename = f"{options.path}/{filename}.root" 
    print(f"Writing to {filename}...")
    outfile = uproot.recreate(f"{filename}") # create the write file

    branches = {column: ak.type(cut_data[column]) for column in cut_data.fields} # here we're defining the branches we want to save to file as a dictionary of name : data type
    outfile.mktree(tree, branches)

    # now write and close
    outfile[tree].extend({branch: cut_data[branch] for branch in branches.keys()})
    outfile.close()

    return print(f'Saved file {filename}.')

def save_all(year, size):
    '''
    Iterates through the save_file function in order to write out all the data in two different
    files, based on the magnet polarity..
    '''
    names = ["up", "down"]
    for index, dataset in enumerate([DATA_UP, DATA_DOWN]):
        save_file(f'{names[index]}_data_{year}_{size}', dataset)
    
    return print(f'Saved {size} data for year 20{year}')

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
        
# - - - - - - - MAIN CODE - - - - - - - #

# Create the necessary flags
options = parse_arguments()

#Variables to be read from the original root file
read_only_these_variables = [
    'P1_PX',
    'P2_PX',
    'P1_PY',
    'P2_PY',
    'P1_PZ',
    'P2_PZ',
    'P1_ETA',
    'P2_ETA',
    'P1_PT',
    'P2_PT',
    'P1_PHI',
    'P2_PHI',
    'D0_ETA',
    'P1_ProbNNk',
    'P2_ProbNNpi',
    'P1_isMuon',
    'P2_isMuon',
    'P1_PIDK',
    'P2_PIDK',
    'P1_ID',
    'P2_ID', 
    'D0_MM',
    'D0_M',
    'D0_ID',
    'D0_PT',
    'D0_ETA',
    'D0_IPCHI2_OWNPV',
    'eventNumber',
    'runNumber'
]

# read data in
raw_data = get_data()

# selection of events
DATA_UP = cut_data(raw_data[0])
DATA_DOWN = cut_data(raw_data[1])

# save cut data 
save_all(str(options.year), options.size)
