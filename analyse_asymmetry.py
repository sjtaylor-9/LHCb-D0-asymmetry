"""
analyse_asymmetry.py

This code generetes 1D and 2D histograms representing the asymmetry distribution in the bins across the phase space.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size. The flag --raw is used to set whether the asymmetry being analysed is raw or not, so that the output can be formatted accordingly. There also are the flags --input and --path which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.
It outputs the .pdf files containing the plots generated.

Author: Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk)
Last edited: 16th September 2023
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import numpy as np
import seaborn as sns

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
    --raw       Used to specify wheteher the asymmetry values cores pond to the raw asymmetry or not.
                Must be either y/Y or n/N
                
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
        "--raw",
        type=str,
        required=True,
        choices=["y", "Y", "n", "N"],
        help="flag to set whether asymmetry is raw"
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

def plot_histogram(val, result, uncertainty, text):
    
    initial = text[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.grid(lw=0.5, ls='--', alpha=0.2)
    ax.hist(asym_val, 15, histtype='stepfilled', ec='k', color='wheat', alpha=0.6)
    ax.set_xlabel(f"{text}Asymmetry")
    ax.set_ylabel("Bins")
    ax.set_title(f"{text}Asymmetry Distribution")
    textstr = '\n \n'.join((
        fr'Weighted $ A_{initial} =$',
        r'${:.3f} \pm {:.3f}$'.format(integrated[0], integrated[1])))

    plt.text(-4.5, 30, textstr, horizontalalignment='center', verticalalignment='center', bbox=dict(boxstyle='square,pad=1', fc='w', ec='k', alpha=1))
    fig.savefig(f"{args.path}/1D_asymmetry_distribution_{args.year}_{args.size}.pdf", dpi=600)
    plt.close(fig)

def integrated_asym(val, err):
    
    weight = err**-2
    numerator = np.sum(val*weight)
    denominator = np.sum(weight)
    weighted_mean = numerator/denominator
    uncertainty = np.sum(weight)**-0.5
    return weighted_mean, uncertainty

def plot_2Dhistogram(val, text):
    length = int(np.sqrt(len(val)))
    data = np.asarray(val).reshape(length, length)
    sns.set()
    ax2 = sns.heatmap(data, vmax=0, vmin=-2.5, annot=True, annot_kws={'size': 8}, cmap ='YlOrBr_r')
    
    ax2.invert_yaxis()
    ax2.set_title(f"{text}Asymmetry distribution")
    ax2.set_xlabel("x index")
    ax2.set_ylabel("y index")
    fig2 = ax2.get_figure()
    fig2.savefig(f"{args.path}/2D_asymmetry_distribution_{args.year}_{args.size}.pdf", dpi=600)
    plt.clf()
    plt.close(fig2)
    
# - - - - - - - MAIN BODY - - - - - - - #

args = parse_arguments()
if args.raw=="y" or args.raw=="Y":
    text="Raw "
else:
    text="Production "
    
asym_val = np.empty(0)
asym_err = np.empty(0)
for j in range(0, 10):
    for i in range(0, 10):
        index = f"{j}" + f"{i}"
        array = np.loadtxt(f"{args.input}/asymmetries_{args.year}_{args.size}_bin{index}.txt")
        asym_val = np.append(asym_val, array[0])
        asym_err = np.append(asym_err, array[1])
        dataset = 0

integrated = integrated_asym(asym_val , asym_err)
print(f"The integrated raw asymmetry is: {integrated[0]} +/- {integrated[1]}") 


plot_histogram(asym_val, integrated[0], integrated[1], text)
plot_2Dhistogram(asym_val, text)