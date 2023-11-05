"""
fit_global.py

This code is used to perform a global fit on the selected data. In order to do so a simulatenous fit is done on the four datasets (with different mesons and polarities). This simulatenous fit keeps all variables constant across the four fits except for the normalisation constants which are allowed to vary independently. The model used consists of a Crystal Ball function and a Gaussian distribution to model the signal and an Exponential decay to model the background.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size. It is necessary to specify if the fit should be performed on the binned data or the unbinned data using the flag --binned_fit. There is a flag --path, which is not required. This one is used to specify the directory where the input data is located, and where the output file should be written. By default it is set to be the current working directory.
It outputs the value of the constants shared in the simultaneous fit to a text file. This code is heavily inspired by Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk), however it has been redesigned so that the binned fit is succesfully performed.

Author: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 5th November 2023
"""

import ROOT
import numpy as np
import uproot
import argparse
import os
from ROOT import TChain, RooRealVar, RooDataSet, RooGaussian, RooCrystalBall, RooAddPdf, RooArgList, RooFit, RooArgSet, RooDataHist, RooExponential
import time 
start_time = time.time()
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
                in the case it is not specified, the default path is the current working directory.
    --binned_fit
                Used to specify if the data should be binned before performing the fit or an unbinned fit should be performed.
                Type either y or Y for a binned fit. Type n or N for an unbinned fit.
    
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
        "--binned_fit",
        type=str,
        choices=["y", "Y", "n", "N"],
        required=True,
        help="flag to set whether a binned or an unbinned should be performed (y/n)"
    )
    
    return parser.parse_args()

# - - - - - - - MAIN BODY - - - - - - - #
args = parse_arguments()
# Bin Parameters
numbins = 10000
lower_boundary = 1820
upper_boundary = 1910

if args.binned_fit=="y" or args.binned_fit=="Y":
    binned = True
else:
    binned = False
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR) # mute RooFit warnings
ROOT.gROOT.SetBatch(True)

# Selects invariant mass (D0_MM) of DO for MagUp
ttree_D0_up = TChain("D02Kpi_Tuple/DecayTree")
ttree_D0_up.Add(f"{args.path}/D0_up_data_{args.year}_{args.size}_clean.root")
ttree_D0_up.SetBranchStatus("*", 0)
ttree_D0_up.SetBranchStatus("D0_MM", 1)

# Selects invariant mass (D0_MM) of DO for MagDown
ttree_D0_down = TChain("D02Kpi_Tuple/DecayTree")
ttree_D0_down.Add(f"{args.path}/D0_down_data_{args.year}_{args.size}_clean.root")
ttree_D0_down.SetBranchStatus("*", 0)
ttree_D0_down.SetBranchStatus("D0_MM", 1)

# Selects invariant mass (D0_MM) of DObar for MagDown
ttree_D0bar_up = TChain("D02Kpi_Tuple/DecayTree")
ttree_D0bar_up.Add(f"{args.path}/D0bar_up_data_{args.year}_{args.size}_clean.root")
ttree_D0bar_up.SetBranchStatus("*", 0)
ttree_D0bar_up.SetBranchStatus("D0_MM", 1)

# Selects invariant mass (D0_MM) of DObar for MagDown
ttree_D0bar_down = TChain("D02Kpi_Tuple/DecayTree")
ttree_D0bar_down.Add(f"{args.path}/D0bar_down_data_{args.year}_{args.size}_clean.root")
ttree_D0bar_down.SetBranchStatus("*", 0)
ttree_D0bar_down.SetBranchStatus("D0_MM", 1)

D0_M = ROOT.RooRealVar("D0_MM", "D0 mass / [MeV/c*c]", 1810, 1910)

# Model Gaussian
mean = RooRealVar("mean", "mean", 1865, 1860, 1870)
sigma = RooRealVar("sigma", "sigma", 8.58, 5, 15)
gaussian = RooGaussian("gauss", "gauss", D0_M, mean, sigma)

# Model CrystalBall
Csig = RooRealVar("Csig", "Csig", 6.23, 5, 15)
aL = RooRealVar("aL", "aL", 2.39, -10, 10)
nL = RooRealVar("nL", "nL", 9.99, -50, 50)
aR = RooRealVar("aR", "aR", 8.50, -10, 10)
nR = RooRealVar("nR", "nR", 27.8, 0, 50)
crystal = RooCrystalBall("Crystal", "Crystal Ball", D0_M, mean, Csig, aL, nL, aR, nR)

# Model Exponential Background
a0 = RooRealVar("a0", "a0", -0.0073, -1, 0)
background = RooExponential("exponential", "exponential", D0_M, a0)

# Model Signal
frac_D0_up = RooRealVar("frac_D0_up", "frac_D0_up", 0.59, 0.4, 0.6)
frac_D0_down = RooRealVar("frac_D0_down", "frac_D0_down", 0.59, 0.4, 0.6)
frac_D0bar_up = RooRealVar("frac_D0bar_up", "frac_D0bar_up", 0.59, 0.4, 0.6)
frac_D0bar_down = RooRealVar("frac_D0bar_down", "frac_D0bar_down", 0.59, 0.4, 0.6)

if binned:
    # Creating the histograms for both polarities for D0 and D0bar by converting the TTree D0_MM data inside the TChain to a TH1(base class of ROOT histograms)
    # TTree.Draw plots a histogram with name D0_Up_Hist and given bin parameters and saves it to memory using: >>
    ttree_D0_up.Draw(f"D0_MM>>D0_Up_Hist({numbins},{lower_boundary},{upper_boundary})")
    # D0_Up_Hist recalled from memory and saved to the local variable
    D0_Up_Hist = ROOT.gPad.GetPrimitive("D0_Up_Hist")

    ttree_D0_down.Draw(f"D0_MM>>D0_Down_Hist({numbins},{lower_boundary},{upper_boundary})")
    D0_Down_Hist = ROOT.gPad.GetPrimitive("D0_Down_Hist")

    ttree_D0bar_up.Draw(f"D0_MM>>D0bar_Up_Hist({numbins},{lower_boundary},{upper_boundary})")
    D0bar_Up_Hist = ROOT.gPad.GetPrimitive("D0bar_Up_Hist")

    ttree_D0bar_down.Draw(f"D0_MM>>D0bar_Down_Hist({numbins},{lower_boundary},{upper_boundary})")
    D0bar_Down_Hist = ROOT.gPad.GetPrimitive("D0bar_Down_Hist")


    # Creating Binned container sets using RooDataHist
    Binned_D0_up = RooDataHist("Binned_D0_up", "Binned D0 Up Data", RooArgList(D0_M), D0_Up_Hist)
    Binned_D0_down = RooDataHist("Binned_D0_down", "Binned D0 Down Data", RooArgList(D0_M), D0_Down_Hist)
    Binned_D0bar_up = RooDataHist("Binned_D0bar_up", "Binned D0bar Up Data", RooArgList(D0_M), D0bar_Up_Hist)
    Binned_D0bar_down = RooDataHist("Binned_D0bar_down", "Binned D0bar Down Data", RooArgList(D0_M), D0bar_Down_Hist)

    # Creating the binned sample and simultaneous PDF
    binned_sample = ROOT.RooCategory("binned_sample", "binned_sample")
    simultaneous_pdf = ROOT.RooSimultaneous("simultaneous", "simultaneous", binned_sample)

    # Model Signal for D0 MagUp
    binned_sample.defineType("Binned_D0_up_sample")
    signal_D0_up = RooAddPdf("signal_D0_up", "signal D0 up", RooArgList(gaussian, crystal), RooArgList(frac_D0_up))
    # Generate normalization variables for D0 MagUp
    Nsig_D0_up = RooRealVar("Nsig_D0_up", "Nsig D0 up", 0.95*Binned_D0_up.numEntries(), 0.9*Binned_D0bar_up.numEntries(), Binned_D0_up.numEntries())
    Nbkg_D0_up = RooRealVar("Nbkg_D0_up", "Nbkg D0 up", 0.05*Binned_D0_up.numEntries(), 0, 0.1*Binned_D0_up.numEntries())
    # Generate model for D0 MagUp
    model_D0_up = RooAddPdf("model_D0_up", "model D0 up", [signal_D0_up, background], [Nsig_D0_up, Nbkg_D0_up])
    simultaneous_pdf.addPdf(model_D0_up, "Binned_D0_up_sample")

    # Model Signal for D0 MagDown
    binned_sample.defineType("Binned_D0_down_sample")
    signal_D0_down = RooAddPdf("signal_D0_down", "signal D0 down", RooArgList(gaussian, crystal), RooArgList(frac_D0_down))
    # Generate normalization variables for D0 MagDown
    Nsig_D0_down = RooRealVar("Nsig_D0_down", "Nsig D0 down", 0.95*Binned_D0_down.numEntries(), 0.9*Binned_D0bar_up.numEntries(), Binned_D0_down.numEntries())
    Nbkg_D0_down = RooRealVar("Nbkg_D0_down", "Nbkg D0 down", 0.05*Binned_D0_down.numEntries(), 0, 0.1*Binned_D0_down.numEntries())
    # Generate model for D0 MagDown
    model_D0_down = RooAddPdf("model_D0_down", "model D0 down", [signal_D0_down, background], [Nsig_D0_down, Nbkg_D0_down])
    simultaneous_pdf.addPdf(model_D0_down, "Binned_D0_down_sample")

    # Model Signal for D0bar MagUp
    binned_sample.defineType("Binned_D0bar_up_sample")
    signal_D0bar_up = RooAddPdf("signal_D0bar_up", "signal D0bar up", RooArgList(gaussian, crystal), RooArgList(frac_D0bar_up))
    # Generate normalization variables for D0bar MagUp
    Nsig_D0bar_up = RooRealVar("Nsig_D0bar_up", "Nsig D0bar up", 0.95*Binned_D0bar_up.numEntries(), 0.9*Binned_D0bar_up.numEntries(), Binned_D0bar_up.numEntries())
    Nbkg_D0bar_up = RooRealVar("Nbkg_D0bar_up", "Nbkg D0bar up", 0.05*Binned_D0bar_up.numEntries(), 0, 0.1*Binned_D0bar_up.numEntries())
    # Generate model for D0bar MagUp
    model_D0bar_up = RooAddPdf("model_D0bar_up", "model D0bar up", [signal_D0bar_up, background], [Nsig_D0bar_up, Nbkg_D0bar_up])
    simultaneous_pdf.addPdf(model_D0bar_up, "Binned_D0bar_up_sample")

    # Model Signal for D0bar MagDown
    binned_sample.defineType("Binned_D0bar_down_sample")
    signal_D0bar_down = RooAddPdf("signal_D0bar_down", "signal D0bar down", RooArgList(gaussian, crystal), RooArgList(frac_D0bar_down))
    # Generate normalization variables for D0bar MagDown
    Nsig_D0bar_down = RooRealVar("Nsig_D0bar_down", "Nsig D0bar down", 0.95*Binned_D0bar_down.numEntries(), 0.9*Binned_D0bar_up.numEntries(), Binned_D0bar_down.numEntries())
    Nbkg_D0bar_down = RooRealVar("Nbkg_D0bar_down", "Nbkg D0bar down", 0.05*Binned_D0bar_down.numEntries(), 0, 0.1*Binned_D0bar_down.numEntries())
    # Generate model for D0bar MagDown
    model_D0bar_down = RooAddPdf("model_D0bar_down", "model D0bar down", [signal_D0bar_down, background], [Nsig_D0bar_down, Nbkg_D0bar_down])
    simultaneous_pdf.addPdf(model_D0bar_down, "Binned_D0bar_down_sample")

    # Recombine the data into a simultaneous dataset
    imports = [ROOT.RooFit.Import("Binned_D0_up_sample", Binned_D0_up), ROOT.RooFit.Import("Binned_D0bar_up_sample", Binned_D0bar_up), ROOT.RooFit.Import("Binned_D0_down_sample", Binned_D0_down), ROOT.RooFit.Import("Binned_D0bar_down_sample", Binned_D0bar_down)]
    simultaneous_data = RooDataHist("simultaneous_data", "simultaneous data", RooArgList(D0_M), ROOT.RooFit.Index(binned_sample), *imports)

    # Performs the simultaneous fit
    fitResult = simultaneous_pdf.fitTo(simultaneous_data, PrintLevel=-1, Save=True, Extended=True)
else:
    # Creates unbinned data containers for all the meson/polarity combinations
    data_D0_up = RooDataSet("data_D0_up", "Data_D0_up", ttree_D0_up, RooArgSet(D0_M))
    data_D0bar_up = RooDataSet("data_D0bar_up", "Data_D0bar_up", ttree_D0bar_up, RooArgSet(D0_M))
    data_D0_down = RooDataSet("data_D0_down", "Data_D0_down", ttree_D0_down, RooArgSet(D0_M))
    data_D0bar_down = RooDataSet("data_D0bar_down", "Data_D0bar_down", ttree_D0bar_down, RooArgSet(D0_M))

    signal_D0_up = RooAddPdf("signal_D0_up", "signal_D0_up", RooArgList(gaussian, crystal), RooArgList(frac_D0_up))
    signal_D0_down = RooAddPdf("signal_D0_down", "signal_D0_down", RooArgList(gaussian, crystal), RooArgList(frac_D0_down))
    signal_D0bar_up = RooAddPdf("signal_D0bar_up", "signal_D0bar_up", RooArgList(gaussian, crystal), RooArgList(frac_D0bar_up))
    signal_D0bar_down = RooAddPdf("signal_D0bar_down", "signal_D0bar_down", RooArgList(gaussian, crystal), RooArgList(frac_D0bar_down))

    # Generate normalization variables
    Nsig_D0_up = ROOT.RooRealVar("Nsig_D0_up", "Nsig_D0_up", 0.95*ttree_D0_up.GetEntries(), 0, ttree_D0_up.GetEntries())
    Nsig_D0bar_up = ROOT.RooRealVar("Nsig_D0bar_up", "Nsig_D0bar_up", 0.95*ttree_D0bar_up.GetEntries(), 0, ttree_D0bar_up.GetEntries())
    Nbkg_D0_up = ROOT.RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", 0.05*ttree_D0_up.GetEntries(), 0, ttree_D0_up.GetEntries())
    Nbkg_D0bar_up = ROOT.RooRealVar("Nbkg_D0bar_up", "Nbkg_D0bar_up", 0.05*ttree_D0bar_up.GetEntries(), 0, ttree_D0bar_up.GetEntries())
    Nsig_D0_down = ROOT.RooRealVar("Nsig_D0_down", "Nsig_D0_down", 0.95*ttree_D0_down.GetEntries(), 0, ttree_D0_down.GetEntries())
    Nsig_D0bar_down = ROOT.RooRealVar("Nsig_D0bar_down", "Nsig_D0bar_down", 0.95*ttree_D0bar_down.GetEntries(), 0, ttree_D0bar_down.GetEntries())
    Nbkg_D0_down = ROOT.RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", 0.05*ttree_D0_down.GetEntries(), 0, ttree_D0_down.GetEntries())
    Nbkg_D0bar_down = ROOT.RooRealVar("Nbkg_D0bar_down", "Nbkg_D0bar_down", 0.05*ttree_D0bar_down.GetEntries(), 0, ttree_D0bar_down.GetEntries())

    # Generate models
    model_D0_up = ROOT.RooAddPdf("model_D0_up", "model_D0_up", [signal_D0_up, background], [Nsig_D0_up, Nbkg_D0_up])
    model_D0bar_up = ROOT.RooAddPdf("model_D0bar_up", "model_D0bar_up", [signal_D0bar_up, background], [Nsig_D0bar_up, Nbkg_D0bar_up])
    model_D0_down = ROOT.RooAddPdf("model_D0_down", "model_D0_down", [signal_D0_down, background], [Nsig_D0_down, Nbkg_D0_down])
    model_D0bar_down = ROOT.RooAddPdf("model_D0bar_down", "model_D0bar_down", [signal_D0bar_down, background], [Nsig_D0bar_down, Nbkg_D0bar_down])

    sample = ROOT.RooCategory("sample", "sample")
    sample.defineType("D0_up")
    sample.defineType("D0_down")
    sample.defineType("D0bar_up")
    sample.defineType("D0bar_down")

    # Combine all models in order to perform a simultaneous fit for all polarities aand mesons
    combData = ROOT.RooDataSet(
        "combData",
        "combined data",
        {D0_M},
        Index=sample,
        Import={"D0_up": data_D0_up, "D0_down": data_D0_down, "D0bar_up": data_D0bar_up, "D0bar_down": data_D0bar_down},
    )
    
    # Performs the simultaneous fit
    simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", {"D0_up": model_D0_up, "D0_down": model_D0_down, "D0bar_up": model_D0bar_up, "D0bar_down": model_D0bar_down}, sample)
    fitResult = simPdf.fitTo(combData, PrintLevel=-1, Save=True, Extended=True)

# Prints the simultaneous fit parameters
fitResult.Print()

# Get results
parameters = np.array([mean.getValV(), sigma.getValV(), Csig.getValV(), aL.getValV(), nL.getValV(), aR.getValV(), nR.getValV(), a0.getValV(), frac_D0_down.getValV(), frac_D0_up.getValV(), frac_D0bar_down.getValV(), frac_D0bar_up.getValV(), Nsig_D0_down.getValV(), Nbkg_D0_down.getValV(), Nsig_D0_up.getValV(), Nbkg_D0_up.getValV(), Nsig_D0bar_down.getValV(), Nbkg_D0bar_down.getValV(), Nsig_D0bar_up.getValV(), Nbkg_D0bar_up.getValV()])
np.savetxt(f"{args.path}/fit_parameters.txt", parameters, delimiter=',')
print("My program took", time.time() - start_time, "to run")