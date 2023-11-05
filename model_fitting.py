"""
model_fitting.py

This code us used to plot the signal Gaussian and Crystal Ball model and Exponential background model using the best fit parameters generated from fit_global.py. The code allows for a binned or unbinned fit depending on the --binned_fit parser, if the unbinned fit is requested then the plot function is called from utils.py.
It then returns the relevant plots of the best fit to the data and a .txt file containing the values and errors on the normalization constant of both signal and background, the mean and standard deviation of the pull distribution and the reduced chi squared value.
The year of interest, size of the data, meson of interest and polarity to be analysed must be specified using the required flags --year --size --meson --polarity. It is also required to specify if the fit should be done on the binned data or the unbinned data using the flag --binned_fit. There also are the flags --input --parameteers_path and --path, which are not required. These are used to specify the directory where the input data is located, where the global best-fit parameters can be found and where the output should be written, respectively. By default it is set to be the current working directory.
This code is heavily inspired by Marc Oriol Pérez, however it has been adapted to correctly plot a binned fit.

Author: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 5th November 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import ROOT
import argparse
import os
from utils import plot
import numpy as np
from ROOT import TChain, RooRealVar, RooDataSet, RooGaussian, RooCrystalBall, RooChebychev, RooAddPdf, RooArgList, RooFit, RooArgSet, RooDataHist, RooExponential
from lhcbstyle import LHCbStyle

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
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --parameters_path
                Used to specify the directory in which the global best-fit parameters should be found. It is not required,
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
        "--parameters_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the global best fit parameters are found"
    )
    parser.add_argument(
        "--input",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input files should be taken from"
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

options = parse_arguments()
# Bin parameters
numbins = 10000
lower_boundary = 1820
upper_boundary = 1910

if options.binned_fit=="y" or options.binned_fit=="Y":
    binned = True
else:
    binned = False

# Reads in the fit parameters generated by fit_global.py, these will be either for a binned/unbinned fit depending on if fit_global.py was ran as a binned fit or not
parameters = np.loadtxt(f"{options.parameters_path}/fit_parameters.txt", delimiter=',')

# Read data
ttree = TChain("D02Kpi_Tuple/DecayTree")
ttree.Add(f"{options.input}/{options.meson}_{options.polarity}_data_{options.year}_{options.size}_clean.root")

ttree.SetBranchStatus("*", 0)
ttree.SetBranchStatus("D0_MM", 1)
D0_M = RooRealVar("D0_MM", r"D0 mass / [MeVc^{-2}]", lower_boundary, upper_boundary) # D0_MM - invariant mass

# Define variables for signal model, using the best fit parameters generated from fit_global.py
mu = RooRealVar("mu", "mu", parameters[0])
Gsig = RooRealVar("sigma", "sigma", parameters[1])
Gauss = RooGaussian("Gauss", "Gaussian", D0_M, mu, Gsig)

Csig = RooRealVar("Csig", "Csig", parameters[2])
aL = RooRealVar("aL", "aL", parameters[3])
nL = RooRealVar("nL", "nL", parameters[4])
aR = RooRealVar("aR", "aR", parameters[5])
nR = RooRealVar("nR", "nR", parameters[6])
Crystal = RooCrystalBall("Crystal", "Crystal Ball", D0_M, mu, Csig, aL, nL, aR, nR)

# Model Exponential Background
a = RooRealVar("a0", "a0", parameters[7])
background = RooExponential("Exponential", "Exponential", D0_M, a)

if options.meson == "D0":
    # D0 MagDown
    if options.polarity == "down":
        frac = RooRealVar("frac_D0_down", "frac_D0_down", parameters[8])
        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[12])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[13])
    # D0 MagUp
    elif options.polarity == "up":
        frac = RooRealVar("frac_D0_up", "frac_D0_up", parameters[9])
        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[14])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[15])
elif options.meson == "D0bar":
    # D0bar MagDown
    if options.polarity == "down":
        frac = RooRealVar("frac_D0bar_down", "frac_D0bar_down", parameters[10])
        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[16])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[17])
    # D0bar MagUp
    elif options.polarity == "up":
        frac = RooRealVar("frac_D0bar_up", "frac_D0bar_up", parameters[11])
        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[18])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[19])

# Define Normalisation constants for signal and background
Nsig = RooRealVar("Nsig", "Nsig", 0.95*ttree.GetEntries(), 0, ttree.GetEntries())
Nbkg = RooRealVar("Nbkg", "Nbkg", 0.05*ttree.GetEntries(), 0, ttree.GetEntries())

# Create model
signal = RooAddPdf("signal", "signal", RooArgList(Gauss, Crystal), RooArgList(frac))
model = {
    "total": RooAddPdf("total", "Total", RooArgList(signal, background), RooArgList(Nsig, Nbkg)), # extended likelihood
    "signals": {
        Gauss.GetName(): Gauss.GetTitle(),
        Crystal.GetName(): Crystal.GetTitle(),
    },
    "backgrounds": {
        background.GetName(): background.GetTitle()
    }
}

# Fit data
if binned:
    with LHCbStyle():
        plot_type: str = "LHCb Simulation",
        # Creates the histogram for the meson by converting the TTree D0_MM data inside the TChain to a TH1(base class of ROOT histograms)
        # TTree.Draw plots a histogram with name D0_Hist and given bin parameters and saves it to memory using: >>
        ttree.Draw(f"D0_MM>>D0_Hist({numbins},{lower_boundary},{upper_boundary})")
        # D0_Hist recalled from memory and saved to the local variable
        D0_Hist = ROOT.gPad.GetPrimitive("D0_Hist")
        # Creating Binned container sets using RooDataHist
        Binned_data = RooDataHist("Binned_data", "Binned Data Set", RooArgList(D0_M), D0_Hist)

        result = model["total"].fitTo(Binned_data, RooFit.Save(True), RooFit.Extended(True))

        mD0 = 1864.84
        mD0_range = (mD0-44.84, mD0+45.16)
        mD0_bins = np.linspace(*mD0_range, numbins+1)

        frame = D0_M.frame(RooFit.Name(""))
        legend_entries = dict()

        Binned_data.plotOn(frame, ROOT.RooFit.Name("remove_me_A"))
        model["total"].plotOn(
            frame,
            RooFit.Name(model["total"].GetName()),
            RooFit.LineWidth(5),
            RooFit.LineColor(ROOT.kAzure),
        )
        pull_hist = frame.pullHist()

        legend_entries[model["total"].GetName()] = {"title": model["total"].GetTitle(), "style": "l"}

        # plot signal components
        signal_colours = [ROOT.kRed, ROOT.kSpring, ROOT.kAzure + 7, ROOT.kOrange + 7]
        signal_line_styles = [2, 7, 9, 10]
        i = 0
        for name, title in model["signals"].items():
            legend_name = f"S{i}"
            model["total"].plotOn(
                frame,
                ROOT.RooFit.Components(name),
                ROOT.RooFit.Name(legend_name),
                ROOT.RooFit.LineWidth(4),
                ROOT.RooFit.LineColor(signal_colours[i % len(signal_colours)]),
                ROOT.RooFit.LineStyle(signal_line_styles[i % len(signal_line_styles)]),
            )
            legend_entries[legend_name] = {"title": title, "style": "l"}
            i += 1

        # plot background components
        background_colours = [ROOT.kMagenta + 2, ROOT.kPink + 7, ROOT.kMagenta + 4]
        background_line_styles = [5, 8, 6]
        i = 0
        for name, title in model["backgrounds"].items():
            legend_name = f"B{i}"
            model["total"].plotOn(
                frame,
                ROOT.RooFit.Components(name),
                ROOT.RooFit.Name(legend_name),
                ROOT.RooFit.LineWidth(4),
                ROOT.RooFit.LineColor(background_colours[i % len(background_colours)]),
                ROOT.RooFit.LineStyle(background_line_styles[i % len(background_line_styles)]),
            )
            legend_entries[legend_name] = {"title": title, "style": "l"}
            i += 1
        
        # plot data points on top again
        Binned_data.plotOn(frame, ROOT.RooFit.Name("remove_me_B"))
        frame.remove("remove_me_A")
        frame.remove("remove_me_B")
        frame.addTH1(D0_Hist, "PE")
        legend_entries[D0_Hist.GetName()] = {"title": D0_Hist.GetTitle(), "style": "PE"}


        frame.SetYTitle(f"Entries MeV/c^{{2}})")

        c = ROOT.TCanvas("fit", "fit", 900, 800)
        fit_pad = ROOT.TPad("fit_pad", "fit pad", 0, 0.2, 1.0, 1.0)
        fit_pad.Draw()
        fit_pad.cd()
        frame.Draw()
        
        frame.GetXaxis().SetLabelSize(0)
        frame.GetXaxis().SetTitleSize(0)
        frame.Draw()
        title_size = frame.GetYaxis().GetTitleSize() * 2.5
        label_size = frame.GetYaxis().GetLabelSize() * 2.5

        # plot_type + total + signals + backgrounds + data
        nlines = 1 + 1 + len(model["signals"]) + len(model["backgrounds"]) + 1
        xwidth = 0.4
        ywidth = 0.03 * nlines
        legend = ROOT.TLegend(
            0.18, 0.89 - ywidth, 0.18 + xwidth, 0.89, "#bf{#it{+plot_type+}}"
        )
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(label_size*0.33)
        for key, val in legend_entries.items():
            legend.AddEntry(key, val["title"], val["style"])
        legend.Draw("same")

        # Plots the pull distribution, where bad pulls (>5 sigma away from the fit) are made to be red
        pull_frame = D0_M.frame(ROOT.RooFit.Title(" "))
        pull_TH1 = ROOT.TH1D("pull_TH1", "pull_TH1", numbins, mD0_bins)
        bad_pull_TH1 = ROOT.TH1D("bad_pull_TH1", "bad_pull_TH1", numbins, mD0_bins)
        for i in range(pull_hist.GetN()):
            if pull_hist.GetPointY(i) > 5:
                pull_TH1.SetBinContent(i + 1, 5)
                bad_pull_TH1.SetBinContent(i + 1, 5)
            elif pull_hist.GetPointY(i) < -5:
                pull_TH1.SetBinContent(i + 1, -5)
                bad_pull_TH1.SetBinContent(i + 1, -5)
            elif pull_hist.GetPointY(i) == 0:
                pull_TH1.SetBinContent(i + 1, 0)
                bad_pull_TH1.SetBinContent(i + 1, 0)
            else:
                pull_TH1.SetBinContent(i + 1, pull_hist.GetPointY(i))
                if abs(pull_hist.GetPointY(i)) >= 3:
                    bad_pull_TH1.SetBinContent(i + 1, pull_hist.GetPointY(i))

        bad_pull_TH1.SetFillColor(ROOT.kRed)
        pull_frame.addTH1(pull_TH1, "bar min0")
        pull_frame.addTH1(bad_pull_TH1, "bar min0")

        c.cd(0)
        pull_pad = ROOT.TPad("pull_pad", "pull pad", 0.0, 0.0, 1.0, 0.31)
        pull_pad.SetBottomMargin(0.4)
        pull_pad.Draw()
        pull_pad.cd()

        # Defines the axes specifications of the pull distribution plot
        pull_frame.GetXaxis().SetLabelSize(label_size)
        pull_frame.GetXaxis().SetTitleSize(title_size)
        pull_frame.GetXaxis().SetTitleOffset(1)
        pull_frame.GetYaxis().SetRangeUser(-5, 5)
        pull_frame.GetYaxis().SetNdivisions(5)
        pull_frame.GetYaxis().SetTitle("Pull [#sigma]")
        pull_frame.GetYaxis().SetLabelSize(label_size)
        pull_frame.GetYaxis().SetTitleSize(title_size)
        pull_frame.GetYaxis().SetTitleOffset(0.39)

        line = ROOT.TLine(D0_M.getMin(), 0, D0_M.getMax(), 0)
        pull_frame.Draw()
        line.Draw("same")

        three = ROOT.TLine(D0_M.getMin(), 3, D0_M.getMax(), 3)
        nthree = ROOT.TLine(D0_M.getMin(), -3, D0_M.getMax(), -3)
        three.SetLineColor(ROOT.kRed)
        three.SetLineStyle(9)
        nthree.SetLineColor(ROOT.kRed)
        nthree.SetLineStyle(9)
        three.Draw("same")
        nthree.Draw("same")

    # Saves the model
    c.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}/D0_fit_ANA.root")
    c.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}/D0_fit_ANA.C")
    c.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}/D0_fit_ANA.pdf")
    c.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}/D0_fit_ANA.jpg")

else:
    unbinned_data = RooDataSet("data", "Data", ttree, RooArgSet(D0_M))
    model["total"].fitTo(unbinned_data, RooFit.Save(), RooFit.Extended(1), RooFit.Minos(0))
    # Generate plots from the plot function in utils.py
    chi2, pull_mean, pull_std = plot(D0_M, unbinned_data, model, nbins=numbins, setlogy=False, save_to=f'{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}', plot_type=f"20{options.year} Mag{(options.polarity).title()}", meson=options.meson)
    # Write out results
    file = open(f"{options.path}/yields_{options.meson}_{options.polarity}_{options.year}_{options.size}.txt", "w")
    text = str(Nsig.getValV()) + ', ' + str(Nsig.getError()) + ', ' + str(Nbkg.getValV()) + ', ' + str(Nbkg.getError()) + ', ' + str(chi2) + ', ' + str(pull_mean) + ', ' + str(pull_std)
    file.write(text)
    file.close


print(ttree.GetEntries())