"""
production_asymmetry.py

This code is used to calculate the systematic uncertainties associated with the choice of signal model, the choice of number of bins to use in the binned fit and the difference between the bin integrated production asymmetry and the average production asymmetry over the p_T and eta binned phase space.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size --scheme --binned_fit. There also are the flags --input --path which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.

Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 16th September 2023
"""


# - - - - - - IMPORT STATEMENTS - - - - - - #

import os
import argparse
import numpy as np

# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year              Used to specify the year at which the data was taken the user is interested in.
                        The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size              Used to specify the amount of events the user is interested in analysing.
                        The argument must be one of: [large, small, medium, 1-8]. The integers specify the number of root
                        files to be read in. Large is equivalent to 8. Medium is equivalent to 4. Small takes 200000 events.
    --path              Used to specify the directory in which the output files should be written. It is not required,
                        in the case it is not specified, the default path is the current working directory.
    --input             Used to specify the directory in which the input data should be found. It is not required,
                        in the case it is not specified, the default path is the current working directory. 
    --optimal_asymmetries_path  
                        Used to specify the directory in which the optimal asymmetry values can be found. It is not required,
                        in the case it is not specified, the default path is the current working directory.
    --model_err_asymmetries_path  
                        Used to specify the directory in which the asymmetry values for the weaker model can be found. It is not required,
                        in the case it is not specified, the default path is the current working directory.
    --bin_err_asymmetries_path  
                        Used to specify the directory in which the asymmetry values for the differing number of bins in the binned fit can be found. It is not required,
                        in the case it is not specified, the default path is the current working directory.
    --binned_fit        Used to specify if the data should be binned before performing the fit or an unbinned fit should be performed.
                        Type either y or Y for a binned fit. Type n or N for an unbinned fit.
    --scheme    Flag indicating which type of binning scheme is used. It is required.

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
    parser.add_argument(
        "--optimal_asymmetries_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the optimal asymmetry values can be found"
    )
    parser.add_argument(
        "--model_err_asymmetries_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the asymmetry values for the second model can be found"
    )
    parser.add_argument(
        "--bins_err_asymmetries_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the asymmetry values for the differing bin number can be found"
    )
    parser.add_argument(
        "--binned_fit",
        type=str,
        choices=["y", "Y", "n", "N"],
        required=True,
        help="flag to set whether a binned or an unbinned should be performed (y/n)"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        choices=["total","pT_eta","pT","eta"],
        required=True,
        help="flag to set whether a binned or an unbinned should be performed (y/n)"
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
    

# - - - - - - - MAIN CODE - - - - - - - #

options = parse_arguments()
scheme = options.scheme


### Need to read in: 
###     A_p_up and A_p_up_2 and A_p_up_3 (2 is the different model, 3 is the different number of bins)
###     Same as above but for A_p_down
###     Bin intgerated asymmetry and bin average asymmetr (A_p and A_prod_unbinned)

# How to add the individual up and down components, need to do for fit_err_model and fit_err_bins
# fit_err_up = np.abs(A_p_up - A_p_up_2)
# fit_err_down = np.abs(A_p_down - A_raw_p_2)
# fit_err = (fit_err_up**2 + fit_err_down**2)**0.5

# Then need to add in quadrature fit_err_model, fit_err_bins, fit_err_bin_integrated

# If an unbinned fit then dont need to do fit_err_bins