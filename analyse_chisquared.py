"""
analyse_chisquared.py

This code generetes 1D and 2d histograms representing the reduced chi squared distribution in the bins across the phase space.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size. There also are the flags --input and --path which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.
It outputs the .pdf files containing the plots generated.

Author: Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk)
Last edited: 16th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import numpy as np
import seaborn as sns

# - - - - - - - FUNCTIONS - - - - - - - #
def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size      Used to specify the amount of events the user is interested in analysing.
                The argument must be one of: [large, small, medium, 1-8]. The integers specify the number of root
                files to be read in. Large is equivalent to 8. Medium is equivalent to 4. Small takes 200000 events.
    --path      Used to specify the directory in which the output files should be written. It is not required,
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

def plot_histogram(val, string):
    """
    Plots and stores a 1D histogram, containing the distribution of chi-squared values
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.grid(lw=0.5, ls='--', alpha=0.8)
    counts, edges, bars = ax.hist(val, bins=[0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 6, 7], color='red', alpha=0.5, rwidth=0.8, lw=1.5, edgecolor='k')
    ax.set_xlabel("Reduced Chi Squared")
    ax.set_ylabel("Bins")
    ax.set_title(r"Reduced Chi Squared distribution {:}".format(string))
    plt.bar_label(bars)
    plt.savefig(f"{args.path}/1D_chisquare_distribution_{string}_{args.year}_{args.size}.pdf", dpi=600)
    plt.clf()


def plot_2Dhistogram(val, string):
    """
    Plots and stores a 2D histogram, showing the distribution of chi-squared values across the phase space
    """
    length = int(np.sqrt(len(val)))
    data = np.asarray(val).reshape(length, length)
    sns.set()
    ax2 = sns.heatmap(data, vmin=0.5, vmax=5, annot=True, annot_kws={'size': 8})
    ax2.invert_yaxis()
    fig2 = ax2.get_figure()
    plt.title(r"$\chi_R^2$ distribution {:}".format(string))
    plt.xlabel("x index")
    plt.ylabel("y index")
    fig2.savefig(f"{args.path}/2D_chisquare_distribution_{string}_{args.year}_{args.size}.pdf", dpi=600)
    plt.clf()
    
    
# - - - - - - - MAIN BODY - - - - - - - #

args = parse_arguments()
chi2 = [[],[],[],[]]

# Load data
for j in range(0, 10):
    for i in range(0, 10):
        index = f"{j}" + f"{i}"
        dataset = 0
        for meson in ["D0", "D0bar"]:
            for polarity in ["up", "down"]:
                yields = np.loadtxt(f"{args.input}/{index}/yields_{meson}_{polarity}_{args.year}_{args.size}_bin{index}.txt", delimiter=',')
                chi2[dataset] = np.append(chi2[dataset], yields[4])
                dataset += 1

# Make plots
for index, element in np.ndenumerate(["D0_up", "D0_down", "D0bar_up", "D0bar_down"]):
    plot_histogram(chi2[index[0]], element)
    plot_2Dhistogram(chi2[index[0]], element)