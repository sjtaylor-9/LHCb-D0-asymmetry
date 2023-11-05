import ROOT as R
from lhcbstyle import LHCbStyle

def plot(
    observable: R.RooRealVar,
    data: R.RooDataSet,
    model: dict,
    nbins: int = 50,
    nparams: int = 0,
    setlogy: bool = False,
    save_to: str = None,
    plot_type: str = "LHCb Simulation",
    meson: str = None,
) -> float:
    # everything below is now in LHCb sytle
    with LHCbStyle():
        # Draw fit with pulls beneath
        fit_frame = observable.frame(R.RooFit.Title(" "))
        legend_entries = dict()

        # plot data points
        data.plotOn(fit_frame, R.RooFit.Name(data.GetName()), R.RooFit.Binning(nbins), R.RooFit.MarkerStyle(8))
    
        # plot model
        model["total"].plotOn(
            fit_frame,
            R.RooFit.Name(model["total"].GetName()),
            R.RooFit.LineWidth(5),
            R.RooFit.LineColor(R.kAzure),
        )
        legend_entries[model["total"].GetName()] = model["total"].GetTitle()

        # plot signal components
        signal_colours = [R.kRed, R.kSpring, R.kAzure + 7, R.kOrange + 7]
        signal_line_styles = [2, 7, 9, 10]
        i = 0
        for name, title in model["signals"].items():
            legend_name = f"S{i}"
            model["total"].plotOn(
                fit_frame,
                R.RooFit.Components(name),
                R.RooFit.Name(legend_name),
                R.RooFit.LineWidth(4),
                R.RooFit.LineColor(signal_colours[i % len(signal_colours)]),
                R.RooFit.LineStyle(signal_line_styles[i % len(signal_line_styles)]),
            )
            legend_entries[legend_name] = title
            i += 1

        # plot background components
        background_colours = [R.kMagenta + 2, R.kPink + 7, R.kMagenta + 4]
        background_line_styles = [5, 8, 6]
        i = 0
        for name, title in model["backgrounds"].items():
            legend_name = f"B{i}"
            model["total"].plotOn(
                fit_frame,
                R.RooFit.Components(name),
                R.RooFit.Name(legend_name),
                R.RooFit.LineWidth(4),
                R.RooFit.LineColor(background_colours[i % len(background_colours)]),
                R.RooFit.LineStyle(background_line_styles[i % len(background_line_styles)]),
            )
            legend_entries[legend_name] = title
            i += 1

        # have to repeat it down here or else the curves disappear...
        data.plotOn(fit_frame, R.RooFit.Name(data.GetName()), R.RooFit.Binning(nbins))

        fit_canvas = R.TCanvas("fit", "fit canvas", 800, 800)
        fit_canvas.Draw()
        fit_pad = R.TPad("fit_pad", "fit pad", 0, 0.2, 1.0, 1.0)
        fit_pad.Draw()
        if setlogy:
            fit_pad.SetLogy()
        fit_pad.cd()

        # make plot for paper first
        fit_frame.GetYaxis().SetMaxDigits(3)  # TODO : a better fix than this?
        fit_frame.Draw()
        # the pull_pad object is not 1x1 but 0.31x1 so the text gets shrunk
        # 2.5 seems to multiply the lhcbstyle font size back to approx the right size
        title_size = fit_frame.GetYaxis().GetTitleSize() * 2.5
        label_size = fit_frame.GetYaxis().GetLabelSize() * 2.5
        fit_frame.GetYaxis().SetNdivisions(-406)

        # plot_type + total + signals + backgrounds
        nlines = 1 + 1 + 1 + 1 + len(model["signals"]) + len(model["backgrounds"])
        xwidth = 0.2
        ywidth = 0.04 * nlines
        legend = R.TLegend(
            0.18, 0.89 - ywidth, 0.18 + xwidth, 0.89, "#bf{#it{"+plot_type+"}}"
        )
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(label_size*0.33)
        legend.AddEntry(data.GetName(), data.GetTitle(), "pe")
        for key, val in legend_entries.items():
            legend.AddEntry(key, val, "l")
        legend.Draw("same")

        if save_to is not None:
            print(f"Saving plots to {save_to}_PAPER.root")
            fit_canvas.SaveAs(save_to + "_PAPER.root")
            fit_canvas.SaveAs(save_to + "_PAPER.pdf")

        # make plot for ana note next
        fit_frame.GetXaxis().SetLabelSize(0)
        fit_frame.GetXaxis().SetTitleSize(0)
        fit_frame.Draw()
        legend.Draw("same")

        # now pulls beneath
        pull_hist = fit_frame.pullHist(data.GetName(), model["total"].GetName())
        pull_frame = observable.frame(R.RooFit.Title(" "))
        # pull_frame.addPlotable(pull_hist, "P") # this makes points with error bars - we want bar chart
        # go about it in an awkward way:
        bin_width = pull_hist.GetPointX(1) - pull_hist.GetPointX(0)
        minimum = pull_hist.GetPointX(0) - bin_width / 2
        maximum = pull_hist.GetPointX(nbins - 1) + bin_width / 2
        pull_TH1 = R.TH1D("pull_TH1", "pull_TH1", nbins, minimum, maximum)
        for i in range(pull_hist.GetN()):
            pull_TH1.SetBinContent(i + 1, pull_hist.GetPointY(i))
        pull_frame.addTH1(pull_TH1, "bar min0")

        fit_canvas.cd(0)
        pull_pad = R.TPad("pull_pad", "pull pad", 0.0, 0.0, 1.0, 0.31)
        pull_pad.SetBottomMargin(0.4)
        pull_pad.Draw()
        pull_pad.cd()

        pull_frame.GetXaxis().SetLabelSize(label_size*0.9)
        pull_frame.GetXaxis().SetTitleSize(title_size*0.9)
        pull_frame.GetXaxis().SetTitleOffset(1)
        pull_frame.GetYaxis().SetRangeUser(-5, 5)
        pull_frame.GetYaxis().SetNdivisions(5)
        pull_frame.GetYaxis().SetTitle("Pull [#sigma]")
        pull_frame.GetYaxis().SetLabelSize(label_size*0.9)
        pull_frame.GetYaxis().SetTitleSize(title_size*0.9)
        pull_frame.GetYaxis().SetTitleOffset(0.39)

        line = R.TLine(observable.getMin(), 0, observable.getMax(), 0)
        line2 = R.TLine(observable.getMin(), 3, observable.getMax(), 3)
        line2.SetLineStyle(9)
        line2.SetLineColor(R.kRed)
        line2.SetLineWidth(2)
        line3 = R.TLine(observable.getMin(), -3, observable.getMax(), -3)
        line3.SetLineStyle(9)
        line3.SetLineColor(R.kRed)
        line3.SetLineWidth(2)
        pull_frame.Draw()
        line.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
        
        fit_canvas.cd(0)
        text = R.TPaveText(0.7, 0.8, 0.9, 0.9)
        if meson == "D0":
            text.AddText(0, 1, "#it{D^{0} #rightarrow K^{-}#pi^{+}}")
        elif meson == "D0bar":
            text.AddText(0, 1, "#it{#bar{D^{0}} #rightarrow K^{+}#pi^{-}}")
        text.SetTextSize(label_size*0.3)
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.Draw("same")

        pull = R.TH1D("Pulls", "Pulls", 50, -5, 5)
        for i in range(pull_hist.GetN()):
            pull.Fill(pull_hist.GetPointY(i))

        pull_canvas = R.TCanvas("fit_pull", "pull canvas", 600, 400)
        pull.GetXaxis().SetTitle("Distance from fit [#sigma]")
        pull.GetYaxis().SetTitle("Entries")
        pull.Draw()
        pull_mean = pull.GetMean()
        pull_std = pull.GetStdDev()

        if save_to is not None:
            print(f"Saving plots to {save_to}_ANA.root")
            fit_canvas.SaveAs(save_to + "_ANA.pdf")
            pull_canvas.SaveAs(save_to + "_pulls.pdf")

        nparams = model["total"].getParameters(data).selectByAttrib("Constant", False).getSize()
        red_chi2 = fit_frame.chiSquare(model["total"].GetName(), data.GetName(), nparams - 1)
        print(f"Reduced chi-squared for {nbins} bins and {nparams} fit parameters: {red_chi2}")

        return red_chi2, pull_mean, pull_std
    
def plot(
    observable: R.RooRealVar,
    data: R.RooDataHist,
    model: dict,
    nbins: int = 50,
    nparams: int = 0,
    setlogy: bool = False,
    save_to: str = None,
    plot_type: str = "LHCb Simulation",
    meson: str = None,
    hist : bool = True,
) -> float:
    # everything below is now in LHCb sytle
    with LHCbStyle():

        fit_frame = observable.frame(R.RooFit.Title(" "))
        legend_entries = dict()

        # plot data points
        data.plotOn(fit_frame, R.RooFit.Name(data.GetName()), R.RooFit.MarkerStyle(8))
    
        # plot model
        model["total"].plotOn(
            fit_frame,
            R.RooFit.Name(model["total"].GetName()),
            R.RooFit.LineWidth(5),
            R.RooFit.LineColor(R.kAzure),
        )
        legend_entries[model["total"].GetName()] = model["total"].GetTitle()

        # plot signal components
        signal_colours = [R.kRed, R.kSpring, R.kAzure + 7, R.kOrange + 7]
        signal_line_styles = [2, 7, 9, 10]
        i = 0
        for name, title in model["signals"].items():
            legend_name = f"S{i}"
            model["total"].plotOn(
                fit_frame,
                R.RooFit.Components(name),
                R.RooFit.Name(legend_name),
                R.RooFit.LineWidth(4),
                R.RooFit.LineColor(signal_colours[i % len(signal_colours)]),
                R.RooFit.LineStyle(signal_line_styles[i % len(signal_line_styles)]),
            )
            legend_entries[legend_name] = title
            i += 1

        # plot background components
        background_colours = [R.kMagenta + 2, R.kPink + 7, R.kMagenta + 4]
        background_line_styles = [5, 8, 6]
        i = 0
        for name, title in model["backgrounds"].items():
            legend_name = f"B{i}"
            model["total"].plotOn(
                fit_frame,
                R.RooFit.Components(name),
                R.RooFit.Name(legend_name),
                R.RooFit.LineWidth(4),
                R.RooFit.LineColor(background_colours[i % len(background_colours)]),
                R.RooFit.LineStyle(background_line_styles[i % len(background_line_styles)]),
            )
            legend_entries[legend_name] = title
            i += 1

        # have to repeat it down here or else the curves disappear...
        data.plotOn(fit_frame, R.RooFit.Name(data.GetName()))

        fit_canvas = R.TCanvas("fit", "fit canvas", 800, 800)
        fit_canvas.Draw()
        fit_pad = R.TPad("fit_pad", "fit pad", 0, 0.2, 1.0, 1.0)
        fit_pad.Draw()
        if setlogy:
            fit_pad.SetLogy()
        fit_pad.cd()

        # make plot for paper first
        fit_frame.GetYaxis().SetMaxDigits(3)  # TODO : a better fix than this?
        fit_frame.Draw()
        # the pull_pad object is not 1x1 but 0.31x1 so the text gets shrunk
        # 2.5 seems to multiply the lhcbstyle font size back to approx the right size
        title_size = fit_frame.GetYaxis().GetTitleSize() * 2.5
        label_size = fit_frame.GetYaxis().GetLabelSize() * 2.5
        fit_frame.GetYaxis().SetNdivisions(-406)

        # plot_type + total + signals + backgrounds
        nlines = 1 + 1 + 1 + 1 + len(model["signals"]) + len(model["backgrounds"])
        xwidth = 0.2
        ywidth = 0.04 * nlines
        legend = R.TLegend(
            0.18, 0.89 - ywidth, 0.18 + xwidth, 0.89, "#bf{#it{"+plot_type+"}}"
        )
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(label_size*0.33)
        legend.AddEntry(data.GetName(), data.GetTitle(), "pe")
        for key, val in legend_entries.items():
            legend.AddEntry(key, val, "l")
        legend.Draw("same")

        if save_to is not None:
            print(f"Saving plots to {save_to}_PAPER.root")
            fit_canvas.SaveAs(save_to + "_PAPER.root")
            fit_canvas.SaveAs(save_to + "_PAPER.pdf")

        # make plot for ana note next
        fit_frame.GetXaxis().SetLabelSize(0)
        fit_frame.GetXaxis().SetTitleSize(0)
        fit_frame.Draw()
        legend.Draw("same")

        # now pulls beneath
        pull_hist = fit_frame.pullHist(data.GetName(), model["total"].GetName())
        pull_frame = observable.frame(R.RooFit.Title(" "))
        # pull_frame.addPlotable(pull_hist, "P") # this makes points with error bars - we want bar chart
        # go about it in an awkward way:
        bin_width = pull_hist.GetPointX(1) - pull_hist.GetPointX(0)
        minimum = pull_hist.GetPointX(0) - bin_width / 2
        maximum = pull_hist.GetPointX(nbins - 1) + bin_width / 2
        pull_TH1 = R.TH1D("pull_TH1", "pull_TH1", nbins, minimum, maximum)
        for i in range(pull_hist.GetN()):
            pull_TH1.SetBinContent(i + 1, pull_hist.GetPointY(i))
        pull_frame.addTH1(pull_TH1, "bar min0")

        fit_canvas.cd(0)
        pull_pad = R.TPad("pull_pad", "pull pad", 0.0, 0.0, 1.0, 0.31)
        pull_pad.SetBottomMargin(0.4)
        pull_pad.Draw()
        pull_pad.cd()

        pull_frame.GetXaxis().SetLabelSize(label_size*0.9)
        pull_frame.GetXaxis().SetTitleSize(title_size*0.9)
        pull_frame.GetXaxis().SetTitleOffset(1)
        pull_frame.GetYaxis().SetRangeUser(-5, 5)
        pull_frame.GetYaxis().SetNdivisions(5)
        pull_frame.GetYaxis().SetTitle("Pull [#sigma]")
        pull_frame.GetYaxis().SetLabelSize(label_size*0.9)
        pull_frame.GetYaxis().SetTitleSize(title_size*0.9)
        pull_frame.GetYaxis().SetTitleOffset(0.39)

        line = R.TLine(observable.getMin(), 0, observable.getMax(), 0)
        line2 = R.TLine(observable.getMin(), 3, observable.getMax(), 3)
        line2.SetLineStyle(9)
        line2.SetLineColor(R.kRed)
        line2.SetLineWidth(2)
        line3 = R.TLine(observable.getMin(), -3, observable.getMax(), -3)
        line3.SetLineStyle(9)
        line3.SetLineColor(R.kRed)
        line3.SetLineWidth(2)
        pull_frame.Draw()
        line.Draw("same")
        line2.Draw("same")
        line3.Draw("same")
        
        fit_canvas.cd(0)
        text = R.TPaveText(0.7, 0.8, 0.9, 0.9)
        if meson == "D0":
            text.AddText(0, 1, "#it{D^{0} #rightarrow K^{-}#pi^{+}}")
        elif meson == "D0bar":
            text.AddText(0, 1, "#it{#bar{D^{0}} #rightarrow K^{+}#pi^{-}}")
        text.SetTextSize(label_size*0.3)
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.Draw("same")

        pull = R.TH1D("Pulls", "Pulls", 50, -5, 5)
        for i in range(pull_hist.GetN()):
            pull.Fill(pull_hist.GetPointY(i))

        pull_canvas = R.TCanvas("fit_pull", "pull canvas", 600, 400)
        pull.GetXaxis().SetTitle("Distance from fit [#sigma]")
        pull.GetYaxis().SetTitle("Entries")
        pull.Draw()
        pull_mean = pull.GetMean()
        pull_std = pull.GetStdDev()

        if save_to is not None:
            print(f"Saving plots to {save_to}_ANA.root")
            fit_canvas.SaveAs(save_to + "_ANA.pdf")
            pull_canvas.SaveAs(save_to + "_pulls.pdf")

        nparams = model["total"].getParameters(data).selectByAttrib("Constant", False).getSize()
        red_chi2 = fit_frame.chiSquare(model["total"].GetName(), data.GetName(), nparams - 1)
        print(f"Reduced chi-squared for {nbins} bins and {nparams} fit parameters: {red_chi2}")
        
        return red_chi2, pull_mean, pull_std