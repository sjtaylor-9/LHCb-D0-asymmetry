
"""
production_asymmetry.py

This code is used to process the signal normalization yields and obtain the production asymmetries of the average bins and global integrated. It finally outputs the results obtained both to the secreen and to a .txt file.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size --scheme --blind. There also are the flags --input --seedval --path which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.
This code is  inspired on the work of Camille Jarvis-Stiggants and Michael England and Marc Oriol Pérez. The code has been completely rewritten and reorganised, and some features have been added to add flexibility to the code, but some of the original functions have been used here as well.

Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 16th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import random
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


# - - - - - - - FUNCTIONS - - - - - - - #

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
        
def read_from_file(meson, polarity, bin_num, scheme):
    '''
    Opens a .txt files and reads the values of the signal normalization constant and its uncertainty.
    
    Returns these two values.
    '''
    with open(f'{options.input}/{scheme}/{bin_num}/yields_{meson}_{polarity}_{options.year}_{options.size}_bin{bin_num}.txt') as f:
        for line in f:
            currentline = line.split(",")
            Nsig = float(currentline[0])
            Nsig_err = float(currentline[1])
    return Nsig, Nsig_err

def get_yield(bin_num, scheme):
    '''
    Gets all the normalization yields, and their uncertainties, necessary to calculate the raw asymmetries.
    This takes into account both D0 and D0bar and both magnet polarities.
    
    Returns all the signal normalization constants, together with their uncertainties.
    '''
    
    yield_D0_up = read_from_file("D0", "up", bin_num, scheme)
    yield_D0bar_up = read_from_file("D0bar", "up", bin_num, scheme)
    yield_D0_down = read_from_file("D0", "down", bin_num, scheme)
    yield_D0bar_down = read_from_file("D0bar", "down", bin_num, scheme)
   
    return yield_D0_up[0], yield_D0_up[1], yield_D0bar_up[0], yield_D0bar_up[1], yield_D0_down[0], yield_D0_down[1], yield_D0bar_down[0], yield_D0bar_down[1]


def A_Det():
      # Detector asymmetry for D0(Kpi) = A_raw(D+->K-pi+pi+) - A_raw(D+->Ks0pi+) - A(Ks0)
        # Running the program AsymmetryTools from GitLab over PuTTY outputs the following:
        if options.year == 16:
            A_kspi_up = -0.87534228056
            A_kspi_err_up = 0.265077797764

            A_kpipi_up = -1.35398359189
            A_kpipi_err_up = 0.13115828851

            A_kspi_down = -0.355750007642
            A_kspi_err_down = 0.247579432594

            A_kpipi_down = -0.637694362926
            A_kpipi_err_down = 0.123796822258

        elif options.year == 17:
            A_kspi_up = -0.654235263918
            A_kspi_err_up = 0.273509295656

            A_kpipi_up = -1.38335612668
            A_kpipi_err_up = 0.121595927374

            A_kspi_down = -0.126746550532
            A_kspi_err_down = 0.269402552773

            A_kpipi_down = -1.02345488078
            A_kpipi_err_down = 0.127466759335

        elif options.year == 18:
            A_kspi_up = -0.942058542057
            A_kspi_err_up = 0.270644276803

            A_kpipi_up = -1.5471233568
            A_kpipi_err_up = 0.131391285254

            A_kspi_down = -0.277769785338
            A_kspi_err_down = 0.288008079473

            A_kpipi_down = -1.27619403857
            A_kpipi_err_down = 0.129960950228

        # Production asymmetry for K0 found from a paper
        A_k0 = 0.054
        A_k0_err = 0.014

        # Calculating the detector asymmetry and the associated uncertainty    
        A_det_up = A_kpipi_up - A_kspi_up - A_k0
        A_det_up_error = np.sqrt(((A_kpipi_err_up)**2+(A_kspi_err_up)**2+(A_k0_err)**2))
        # Detector asymmetry is the same for all models
        

        A_det_down = A_kpipi_down - A_kspi_down - A_k0
        A_det_down_error =  np.sqrt(((A_kpipi_err_down)**2+(A_kspi_err_down)**2+(A_k0_err)**2))
        

        return A_det_up, A_det_down, A_det_down_error, A_det_up_error



##################################################
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

def output_results(A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err, bin_num, A_prod, A_prod_err):
    '''
    This function takes all the necessary values and outputs them to the screen in a nicely formatted way.
    It also outputs them to a .txt file, written in the directory established by the user.
    '''

    asymmetry = str(round(A_raw, 3)) + ' +/- ' + str(round(A_raw_err, 3)) + ' (stat) +/- '
    print(f'The 20{options.year} raw asymmetry of bin {bin_num} is:', asymmetry)
    print("------------------------------")
    
    array = np.array([A_prod, A_prod_err, A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err])
    np.savetxt(f"{options.path}/asymmetries_{options.year}_{options.size}_bin{bin_num}.txt", array)

def blind(A_raw, A_raw_error, string):
    ''' 
    This function blinds the values of the raw asymmetry in order to prevent bias experimenter's bias, the unintended biasing of a result in a particular direction.
    A raw is randomly multiplied by 1 or -1 then a pseudo-random constant is added. The pseudo-random constant is randomly chosen from within 3 standard deviations of a standard Gaussian distribution.
    The blinded result is reproudicble because the randomness is seeded.
    '''
    factorchoice = [-1, 1]
    factor = np.random.choice(factorchoice)
    np.random.seed(options.seedval)
    # Draw random samples from a normal (Gaussian) distribution, with a mean of 0, std of 1 and then the sample is within 3std of the mean.
    constant = np.random.normal(0, 3*A_raw_error, 1)[0]
    Ablind = factor * A_raw + constant

    if string == 'up':
        print('The blinded Magup raw asymmetry is ', round(Ablind*100, 2),'%')
    if string == 'down':
        print('The blinded MagDown raw asymmetry is ', round(Ablind*100, 2),'%')
    
    return 100*Ablind


def production_asymm(A_raw_up, A_raw_down, A_raw_up_err, A_raw_down_err, A_det_up, A_det_down, A_det_down_err, A_det_up_err):
    
    A_prod_up = A_raw_up - A_det_up
    A_prod_down = A_raw_down - A_det_down

    A_prod = (A_prod_up + A_prod_down) / 2

    
    A_prod_up_err = (A_raw_up_err**2 + A_det_up_err**2)**(0.5)
    A_prod_down_err = (A_raw_down_err**2 + A_det_down_err**2)**(0.5)
    A_prod_err = ((A_prod_up_err**2+A_prod_down_err**2)**(0.5))/2


    return A_prod, A_prod_up, A_prod_down, A_prod_up_err, A_prod_down_err, A_prod_err


def integrated_asym(val, err):

    weighted_mean, sum_weights = np.average(np.array(val), weights=np.power(np.array(err), -2), returned=True)
    uncertainty = np.power(sum_weights, -0.5)
    return weighted_mean, uncertainty



def A_prod_unbinned():

    yield_D0_up = read_from_file_global("D0", "up")
    yield_D0bar_up = read_from_file_global("D0bar", "up")
    yield_D0_down = read_from_file_global("D0", "down")
    yield_D0bar_down = read_from_file_global("D0bar", "down")
    A_det_up_local, A_det_down_local, A_det_down_error_local, A_det_up_error_local = A_Det()

    A_raw_up_local, A_raw_up_err_local = calculate_raw_asymmetry(yield_D0_up[0], yield_D0bar_up[0], 1, yield_D0_up[1], yield_D0bar_up[1])
    A_raw_down_local, A_raw_down_err_local = calculate_raw_asymmetry(yield_D0_down[0], yield_D0bar_down[0], 1, yield_D0_down[1], yield_D0bar_down[1])
    A_prod_unbinned_local= production_asymm(A_raw_up_local, A_raw_down_local, A_raw_up_err_local, A_raw_down_err_local, A_det_up_local, A_det_down_local, A_det_down_error_local, A_det_up_error_local)

    return A_prod_unbinned_local

def read_from_file_global(meson, polarity):
    with open(f'{options.input}/global/yields_{meson}_{polarity}_{options.year}_{options.size}.txt') as f:
        for line in f:
            currentline = line.split(",")
            Nsig = float(currentline[0])
            Nsig_err = float(currentline[1])
        f.close()
    return Nsig, Nsig_err




# - - - - - - - MAIN CODE - - - - - - - #

options = parse_arguments()
scheme = options.scheme


A_prod_list = []
A_prod_err_list = []
A_prod_up_list = []
A_prod_up_err_list = []
A_prod_down_list = []
A_prod_down_err_list = []

A_prod_unbinned = A_prod_unbinned()

if scheme == 'pT_eta':
    for j in range(0,10):
        for i in range (0,10):
            scheme = 'local'
            bin_num = str(j)+str(i)
            # get normalization yield from desired model 
            N_D0_up, N_D0_up_err, N_D0bar_up, N_D0bar_up_err, N_D0_down, N_D0_down_err, N_D0bar_down, N_D0bar_down_err = get_yield(bin_num,scheme)

            # get raw asymmetries for main model
            A_raw_up, A_raw_up_err = calculate_raw_asymmetry(N_D0_up, N_D0bar_up, 1, N_D0_up_err, N_D0bar_up_err)
            A_raw_down, A_raw_down_err = calculate_raw_asymmetry(N_D0_down, N_D0bar_down, 1, N_D0_down_err, N_D0bar_down_err)
            #A_raw = (A_raw_up + A_raw_down) / 2
            A_raw_err = ((A_raw_up_err**2 + A_raw_down_err**2)**0.5) /2

            A_det_up, A_det_down, A_det_down_error, A_det_up_error = A_Det()
            A_prod_bin = production_asymm(A_raw_up, A_raw_down, A_raw_up_err, A_raw_down_err, A_det_up, A_det_down, A_det_down_error, A_det_up_error)

            if options.blind == 'y' or options.blind == 'Y':
                A_unblind_up = A_raw_up
                # Asymmetry is blinded
                A_raw_up = blind(A_raw_up, A_raw_up_err, 'up')
                A_unblind_down = A_raw_down
                # Asymmetry is blinded
                A_raw_down = blind(A_raw_down, A_raw_down_err, 'down')
                A_raw = (A_raw_down + A_raw_up) /2

                # Calculating Unblind prod Asymetry of each bin
                # output results
                output_results(A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err, bin_num, A_prod_bin[0], A_prod_bin[5])
                
                A_prod_list.append(A_prod_bin[0])
                A_prod_err_list.append(A_prod_bin[5])
                A_prod_up_list.append(A_prod_bin[1])
                A_prod_up_err_list.append(A_prod_bin[3])
                A_prod_down_list.append(A_prod_bin[2])
                A_prod_down_err_list.append(A_prod_bin[4])
                


            else:
                # output results
                output_results(A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err, bin_num)
                A_prod_list.append(A_prod_bin[0])
                A_prod_err_list.append(A_prod_bin[5])
                A_prod_up_list.append(A_prod_bin[1])
                A_prod_up_err_list.append(A_prod_bin[3])
                A_prod_down_list.append(A_prod_bin[2])
                A_prod_down_err_list.append(A_prod_bin[4])
                

            # Calculates A_det
elif scheme == 'pT' or scheme == 'eta':
    for j in range(0,10):
        bin_num = str(j)
        # get normalization yield from desired model 
        N_D0_up, N_D0_up_err, N_D0bar_up, N_D0bar_up_err, N_D0_down, N_D0_down_err, N_D0bar_down, N_D0bar_down_err = get_yield(bin_num,scheme)

        # get raw asymmetries for main model
        A_raw_up, A_raw_up_err = calculate_raw_asymmetry(N_D0_up, N_D0bar_up, 1, N_D0_up_err, N_D0bar_up_err)
        A_raw_down, A_raw_down_err = calculate_raw_asymmetry(N_D0_down, N_D0bar_down, 1, N_D0_down_err, N_D0bar_down_err)
        #A_raw = (A_raw_up + A_raw_down) / 2
        A_raw_err = ((A_raw_up_err**2 + A_raw_down_err**2)**0.5) /2

        A_det_up, A_det_down, A_det_down_error, A_det_up_error = A_Det()
        A_prod_bin = production_asymm(A_raw_up, A_raw_down, A_raw_up_err, A_raw_down_err, A_det_up, A_det_down, A_det_down_error, A_det_up_error)






        if options.blind == 'y' or options.blind == 'Y':
            A_unblind_up = A_raw_up
            # Asymmetry is blinded
            A_raw_up = blind(A_raw_up, A_raw_up_err, 'up')
            A_unblind_down = A_raw_down
            # Asymmetry is blinded
            A_raw_down = blind(A_raw_down, A_raw_down_err, 'down')
            A_raw = (A_raw_down + A_raw_up) /2

            # Calculating Unblind prod Asymetry of each bin
            # output results
            output_results(A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err, bin_num, A_prod_bin[0], A_prod_bin[5])

            A_prod_list.append(A_prod_bin[0])
            A_prod_err_list.append(A_prod_bin[5])
            A_prod_up_list.append(A_prod_bin[1])
            A_prod_up_err_list.append(A_prod_bin[3])
            A_prod_down_list.append(A_prod_bin[2])
            A_prod_down_err_list.append(A_prod_bin[4])


        else:
            # output results
            output_results(A_raw, A_raw_err, A_raw_up, A_raw_up_err, A_raw_down, A_raw_down_err, bin_num)
            # A_raw_up_list_unblinded.append(A_unblind_up)
            # A_raw_down_list_unblinded.append(A_unblind_down)
            # A_raw_up_err_list.append(A_raw_up_err)
            # A_raw_down_err_list.append(A_raw_down_err)
            A_prod_list.append(A_prod_bin[0])
            A_prod_err_list.append(A_prod_bin[5])
            A_prod_up_list.append(A_prod_bin[1])
            A_prod_up_err_list.append(A_prod_bin[3])
            A_prod_down_list.append(A_prod_bin[2])
            A_prod_down_err_list.append(A_prod_bin[4])

# Calculates the integrated production asymmetry for both and unblinded asymmetries
print("------------------------------")
print("------------------------------")
print("------------------------------")

#Prod
AProd = integrated_asym(A_prod_list, A_prod_err_list)
AProd_Up = integrated_asym(A_prod_up_list, A_prod_up_err_list)
AProd_Down = integrated_asym(A_prod_down_list, A_prod_down_err_list)

print("------------------------------")
print("------------------------------")
print("------------------------------")

print('The MagUp detector asymmetry is: ', round(A_det_up, 2), '% +/-', round(A_det_up_error, 2), '%')
print('The MagDown detector asymmetry is: ', round(A_det_down, 2), '% +/-', round(A_det_down_error, 2), '%')

print("------------------------------")
print("------------------------------")
print("------------------------------")


print("------------------------------")
print("------------------------------")
print("------------------------------")
print(f"The 20{options.year} unblinded integrated MagUp prod asymmetry is: ", round(AProd_Up[0],3), "% +/-", round(AProd_Up[1], 3), '%')
print(f"The 20{options.year} unblinded integrated MagDown prod asymmetry is: ", round(AProd_Down[0],3), "% +/-", round(AProd_Down[1], 3), '%')
print(f"The 20{options.year} unblinded integrated total prod asymmetry is: ", round(AProd[0],3), "% +/-", round(AProd[1], 3), '%')



#Saving Aprod calculated from weighted average and Aprod calculated from unbinned average of D0_up, D0_down, D0bar up, D0bar down

# Aprod, Aprod error -- both from weighted mean, Aprod from unbinned average, Aprod err from unbinned average
array = np.array([AProd[0],AProd[1]
                    , A_prod_unbinned[0], A_prod_unbinned[5], AProd_Up[0],AProd_Up[1],AProd_Down[0],AProd_Down[1]])
np.savetxt(f"{options.results_path}/final_asymmetries_{options.scheme}_{options.year}_{options.size}.txt", array)

file_path = f"{options.results_path}/final_Text_asymmetries_{options.scheme}_{options.year}_{options.size}.txt"

with open(file_path, "w") as file:
    text = (
        f'A_prod_bin: {round(AProd[0], 3)}\n'
        f'A_prod_err_bin: {round(AProd[1], 3)}\n'
        f'A_prod_intergrated: {round(A_prod_unbinned[0], 3)}\n'
        f'A_prod_err_intergrated: {round(A_prod_unbinned[5], 3)}\n'
        f'A_prod_UP_bin: {round(AProd_Up[0], 3)}\n'
        f'A_prod_UP_err_bin: {round(AProd_Up[1], 3)}\n'
        f'A_prod_DOWN_bin: {round(AProd_Down[0], 3)}\n'
        f'A_prod_DOWN_err_bin: {round(AProd_Down[1], 3)}\n'
        f'A_prod_UP_intergrated: {round(A_prod_unbinned[1], 3)}\n'
        f'A_prod_UP_err_intergrated: {round(A_prod_unbinned[3], 3)}\n'
        f'A_prod_DOWN_intergrated: {round(A_prod_unbinned[2], 3)}\n'
        f'A_prod_DOWN_err_intergrated: {round(A_prod_unbinned[4], 3)}\n'
    )
    file.write(text)