import os
import argparse
import numpy as np
from lhcbstyle import LHCbStyle
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import ListedColormap
from matplotlib import colormaps
import awkward as ak
import sys

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
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--bin_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the binning scheme should be found"
    )
    parser.add_argument(
        "--asymm_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the production asymmetry for each bin should be found"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        choices=["pT","eta"],
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

def read_asymmetry_values():
    with open(f'{args.path}/final_asymmetries_{args.scheme}_{args.year}_{args.size}.txt') as f:
        lines = f.readlines()
        binned_asymm = float(lines[0])
        binned_asymm_error = float(lines[1])
        unbinned_asymm = float(lines[2])
        unbinned_asymm_error = float(lines[3])
        f.close()
    return binned_asymm, binned_asymm_error , unbinned_asymm, unbinned_asymm_error

# - - - - - - - MAIN BODY - - - - - - - #
args = parse_arguments()

binned_asymm, binned_asymm_error , unbinned_asymm, unbinned_asymm_error = read_asymmetry_values()

asymmetry = []
asymmetry_error = []
for j in range(0,10):
    bin_num = str(j)
    with open(f'{args.asymm_path}/asymmetries_{args.year}_{args.size}_bin{bin_num}.txt') as f:
        lines = f.readlines()
        A_prod = float(lines[0])
        A_prod_err = float(lines[1])
    f.close()
    asymmetry.append(A_prod)
    asymmetry_error.append(A_prod_err)

x_value = []



file_path = f"{args.bin_path}/{args.year}_{args.size}_{args.scheme}_bins.txt"  # Replace with the path to your file

# Open the file in read mode
with open(file_path, 'r') as file:
    # Read all lines from the file and store them in a list
    x_value = [float(line.strip()) for line in file.readlines()]
if args.scheme == 'pT':
    x_value = [x / 1000 for x in x_value]
x_value = x_value[1:]

# Plotting
fig, ax = plt.subplots()

Data = ax.errorbar(x_value, asymmetry, yerr=asymmetry_error, fmt='o', capsize=5, color = 'black', label = 'Data')
if args.scheme == 'pT':
    ax.set_xlabel(r'$p_{T}$ [GeV$c^{-1}$]', fontsize = 16)
elif args.scheme == 'eta':
    ax.set_xlabel(r'$\eta$', fontsize = 16)

ax.set_ylabel(r'$A_{\mathrm{prod}}$ [%]', fontsize = 16)
ax.tick_params(axis='both', which='both', labelsize=12)

line1 = ax.axhline(binned_asymm, color='blue', linestyle='dashed', linewidth=1.5)
fill1 = plt.axhspan(binned_asymm-binned_asymm_error, binned_asymm+binned_asymm_error, color='blue', alpha=0.35, lw=0)

line2 = ax.axhline(unbinned_asymm, color='red', linestyle='solid', linewidth=1.5)
fill2 = plt.axhspan(unbinned_asymm-unbinned_asymm_error, unbinned_asymm+unbinned_asymm_error, color='red', alpha=0.35, lw=0)

if args.scheme == 'pT':
    ax.legend([(line1,fill1),(line2,fill2),Data],[r'Average result over $p_{T}$ bins', r'Bin integrated result','Data'])
elif args.scheme == 'eta':
    ax.legend([(line1,fill1),(line2,fill2)],[r'Average result over $\eta$ bins', r'Bin integrated result'])

if args.scheme == 'pT':
    plt.savefig(f'{args.path}/pT_Asymm_{args.year}_{args.size}.pdf', bbox_inches = "tight")
elif args.scheme == 'eta':
    plt.savefig(f'{args.path}/eta_Asymm_{args.year}_{args.size}.pdf', bbox_inches = "tight")


plt.show()