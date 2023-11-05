"""
calculate_asymmetry.py

This code is used to process the signal normalization yields and obtain the raw asymmetries. It also uses the input from using another signal model to obtain the systematic uncertainty due to the fit model. It finally outputs the results obtained both to the secreen and to a .txt file.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size.There also are the flags --input and --path which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.
This code is  inspired on the work of Camille Jarvis-Stiggants and Michael England. The code has been completely rewritten and reorganised, and some features have been added to add flexibility to the code, but some of the original functions have been used here as well.

Author: Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk)
Last edited: 16th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import random
import os
import argparse
import numpy as np


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
    --input     Used to specify the directory in which the input data should be found. It is not required,
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
    parser.add_argument(
        "--input",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found"
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
        
def read_from_file(meson, polarity, bin_num):
    '''
    Opens a .txt files and reads the values of the signal normalization constant and its uncertainty.
    
    Returns these two values.
    '''
    with open(f'{options.input}/{bin_num}/yields_{meson}_{polarity}_{options.year}_{options.size}_bin{bin_num}.txt') as f:
        for line in f:
            currentline = line.split(",")
            Nsig = float(currentline[0])
            Nsig_err = float(currentline[1])
        f.close()
            
    return Nsig, Nsig_err

def get_yield(bin_num):
    '''
    Gets all the normalization yields, and their uncertainties, necessary to calculate the raw asymmetries.
    This takes into account both D0 and D0bar and both magnet polarities.
    
    Returns all the signal normalization constants, together with their uncertainties.
    '''
    
    yield_D0_up = read_from_file("D0", "up", bin_num)
    yield_D0bar_up = read_from_file("D0bar", "up", bin_num)
    yield_D0_down = read_from_file("D0", "down", bin_num)
    yield_D0bar_down = read_from_file("D0bar", "down", bin_num)
    
    return yield_D0_up[0], yield_D0_up[1], yield_D0bar_up[0], yield_D0bar_up[1], yield_D0_down[0], yield_D0_down[1], yield_D0bar_down[0], yield_D0bar_down[1]

def calculate_raw_asymmetry(norm_D0, norm_D0bar, bin_width, N_D0_err, N_D0bar_err):
    '''
    It takes the normalization yields for D0 and D0bar as arguments and then calculates the raw
    asymmetries from these. It also propagates the uncertainties.
    
    Returns both the asymmetry and its uncertainty as a percentage.
    '''
    
    N_D0 = abs(norm_D0)/abs(bin_width)
    N_D0bar = abs(norm_D0bar)/abs(bin_width)
    
    A = (N_D0 - N_D0bar)/(N_D0 + N_D0bar)
    A_err = 2*(((N_D0bar**2)*(N_D0_err**2) + (N_D0**2)*(N_D0bar_err**2))**0.5)*((N_D0 + N_D0bar)**(-2))
          
    return 100*A, 100*A_err

def output_results(A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err, bin_num):
    '''
    This function takes as arguments all the necessary values and outputs them to the screen in a nicely formatted way.
    It also outputs them to a .txt file, written in the directory established by the user.
    '''
    
    print('The MagUp raw asymmetry is: (', round(A_raw_up, 3), '+/-', round(A_raw_up_err, 3), ') %')
    print('The MagDown raw asymmetry is: (', round(A_raw_down, 3), '+/-', round(A_raw_down_err, 3), ') %')
    
    asymmetry = str(round(A_raw, 3)) + ' +/- ' + str(round(A_raw_err, 3)) + ' (stat) +/- '
    print(f'The 20{options.year} raw asymmetry of bin {bin_num} is:', asymmetry)
    
    array = np.array([A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err])
    np.savetxt(f"{options.path}/asymmetries_{options.year}_{options.size}_bin{bin_num}.txt", array, delimiter=',')
    
    
# - - - - - - - MAIN CODE - - - - - - - #

options = parse_arguments()

for j in range(0,10):
    for i in range (0,10):
        bin_num = str(j)+str(i)
        # get normalization yield from desired model 
        N_D0_up, N_D0_up_err, N_D0bar_up, N_D0bar_up_err, N_D0_down, N_D0_down_err, N_D0bar_down, N_D0bar_down_err = get_yield(bin_num)

        # get raw asymmetries for main model
        A_raw_up, A_raw_up_err = calculate_raw_asymmetry(N_D0_up, N_D0bar_up, 1, N_D0_up_err, N_D0bar_up_err)
        A_raw_down, A_raw_down_err = calculate_raw_asymmetry(N_D0_down, N_D0bar_down, 1, N_D0_down_err, N_D0bar_down_err)
        A_raw = (A_raw_up + A_raw_down) / 2
        A_raw_err = ((A_raw_up_err**2 + A_raw_down_err**2)**0.5) /2

        # output results
        output_results(A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err, bin_num)
        