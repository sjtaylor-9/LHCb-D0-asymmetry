"""
plot_phase_space.py

This code plots relevant histograms to show the event distribution across the phase space and it also shows the edges of the bins.
The year of interest, size of the data, polarity and meson to be analysed must be specified using the required flags --year --size --polarity --meson. There also are the flags --input --path and --bin_path, which are not required. These are used to specify the directory where the input data is located, where the binning scheme can be found and where the output file should be written, respectively. By default it is set to be the current working directory.
It outputs several pdf files containing the relevant histograms.

Author: Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk)
Last edited: 15th September 2023
"""
# - - - - - - IMPORT STATEMENTS - - - - - - #
import ROOT
import os
import argparse
import numpy as np
from lhcbstyle import LHCbStyle
import uproot
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import ListedColormap

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

# Import data

tree_name = "D02Kpi_Tuple/DecayTree"
data = uproot.concatenate(f"{args.input}/{args.polarity}_data_{args.year}_{args.size}_clean.root:{tree_name}")

bins = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_bins.txt", delimiter=',')
bins[0] = bins[0]/1000
viridis = matplotlib.cm.get_cmap('YlOrRd')
newcolors = viridis(np.linspace(0, 1, 25))
newcmp = ListedColormap(newcolors)

pT = data["D0_PT"]/1000
eta = data["D0_ETA"]

# First histogram
fig = plt.figure()
ax = fig.add_subplot(111)
h2d = ax.hist2d(pT, eta, bins=200, cmap=newcmp)
ax.set_xlabel(r'$p_{T}$ [GeV/c]')
ax.set_ylabel(r'$\eta$')
ax.set_title(r'$p_{T}$ vs $\eta$')
fig.colorbar(h2d[3], ax=ax, label='Events')
plt.savefig(f'{args.path}/2D_histogram_no_bins_extended_{args.meson}_{args.polarity}_{args.year}_{args.size}.pdf')

# Second histogram
mask = np.ones(len(data["D0_PT"]))
mask = np.logical_and(mask, data["D0_PT"]<=10000)
data = data[mask]
pT = data["D0_PT"]/1000
eta = data["D0_ETA"]

viridis = matplotlib.cm.get_cmap('YlOrRd')
newcolors = viridis(np.linspace(0, 1, 10))
newcmp = ListedColormap(newcolors)

fig = plt.figure()
ax = fig.add_subplot(111)
h2d = ax.hist2d(pT, eta, bins=200, cmap=newcmp)
ax.set_xlabel(r'$p_{T}$ [GeV/c]')
ax.set_ylabel(r'$\eta$')
ax.set_title(r'$p_{T}$ vs $\eta$')
fig.colorbar(h2d[3], ax=ax, label='Events')
plt.savefig(f'{args.path}/2D_histogram_no_bins_{args.meson}_{args.polarity}_{args.year}_{args.size}.pdf')

# Third histogram
fig = plt.figure()
ax = fig.add_subplot(111)
h2d = ax.hist2d(pT, eta, bins=200, cmap=newcmp)
ax.set_xlabel(r'$p_{T}$ [GeV/c]')
ax.set_ylabel(r'$\eta$')
ax.set_title(r'$p_{T}$ vs $\eta$')
ax.set_xlim(2, 10)
ax.set_ylim(2, 5)
fig.colorbar(h2d[3], ax=ax, label='Events')

for index in np.arange(0,10):
    if index!=0:
        ax.axvline(bins[0,index], ymin=0, ymax=1, color='blue')
    for j in np.arange(0,10):
        if j!=0:
            ax.axhline(bins[index+1, j], xmin=(bins[0,index]-2)/8, xmax=(bins[0,index+1]-2)/8, color='blue')
    
plt.savefig(f'{args.path}/2D_histogram_bins_{args.meson}_{args.polarity}_{args.year}_{args.size}.pdf')