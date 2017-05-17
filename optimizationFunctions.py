#!/usr/bin/env python
#from sympy import *
from math import *
import ROOT
import os

class equations:
    """____________________________________________________________________
    A class that contains the equations used for the optimization
    """
    # def errorStat(self,X,DX,Y,Z):
    #   "Calculate statistical error of error"
    #   x, dx, y, dy, z, dz = symbols('x dx y dy z dz')
    #   #error = sqrt(diff(x/(y*z), x)*dx**2+diff(x/(y*z), y)*dy**2+diff(x/(y*z)
    #   #,z)*dz**2)
    #   error = sqrt(diff(x/(y*z), x)*dx**2)
    #   #errors_subs = error.subs([(x,X),(dx,DX),(y,Y),(dy,DY),(z,Z),(dz,DZ)])
    #   errors_subs = error.subs([(x,X),(dx,DX),(y,Y),(z,Z)])
    #   print "Statistical error = ", errors_subs.evalf()
    #   return errors_subs.evalf()

    # def errorSys(self,X,DX,Y,DY,A,B,C):
    #   "Calculate systematic error of error"
    #   x, dx, y, dy = symbols('x dx y dy')
    #   a, b, c = symbols('a b c')
    #   expression = sqrt( (x*a)**2 + (y*b)**2 + (y*c)**2 )/y
    #   error = sqrt((diff(expression, x)*dx)**2 + (diff(expression, y)*dy)**2)
    #   errors_subs = error.subs([(x,X),(dx,DX),(y,Y),(dy,DY),(a,A),(b,B),(c,C)])
    #   print "Systematic error = ", errors_subs.evalf()
    #   return errors_subs.evalf()

    def _sigma(self,signal, background,efficiency,lumi):
        "Calculate cross section"
        return (signal)/(efficiency*lumi)

    def _sigmaStat(self,signal,background, efficiency,lumi):
        "Calculate statistical error"
        return sqrt(signal+ background)/(efficiency*lumi)

    def _sigmaSys(self,
        sigma_,signal,  
        bkg_sys_hfake,
        bkg_sys_efake,
        bkg_sys_other,
        efficiency, 
        expSys, 
        lumiSys,
        lumi):
        "Calculate systematic error"
        sigsys = sqrt(
          ((bkg_sys_hfake["ttbar"])/(efficiency*lumi))**2 +
          ((bkg_sys_hfake["Wjets"])/(efficiency*lumi))**2 +
          ((bkg_sys_hfake["Wgamma"])/(efficiency*lumi))**2 +
          ((bkg_sys_hfake["Zjets"])/(efficiency*lumi))**2 +
          ((bkg_sys_hfake["Zgamma"])/(efficiency*lumi))**2 +
          ((bkg_sys_hfake["ST"])/(efficiency*lumi))**2 +
          ((bkg_sys_hfake["VV"])/(efficiency*lumi))**2 +

          ((bkg_sys_efake["ttbar"])/(efficiency*lumi))**2 +
          ((bkg_sys_efake["Wjets"])/(efficiency*lumi))**2 +
          ((bkg_sys_efake["Wgamma"])/(efficiency*lumi))**2 +
          ((bkg_sys_efake["Zjets"])/(efficiency*lumi))**2 +
          ((bkg_sys_efake["Zgamma"])/(efficiency*lumi))**2 +
          ((bkg_sys_efake["ST"])/(efficiency*lumi))**2 +
          ((bkg_sys_efake["VV"])/(efficiency*lumi))**2 +

          ((bkg_sys_other["ttbar"])/(efficiency*lumi))**2 +
          ((bkg_sys_other["Wjets"])/(efficiency*lumi))**2 +
          ((bkg_sys_other["Wgamma"])/(efficiency*lumi))**2 +
          ((bkg_sys_other["Zjets"])/(efficiency*lumi))**2 +
          ((bkg_sys_other["Zgamma"])/(efficiency*lumi))**2 +
          ((bkg_sys_other["ST"])/(efficiency*lumi))**2 +
          ((bkg_sys_other["VV"])/(efficiency*lumi))**2
          )
        # (sigma_*expSys)**2 +
        #     (sigma_*lumiSys)**2)
        return sigsys


class plots:
    """____________________________________________________________________
    A class to make plots
    """
    def createPlot(self,filename, channel, xlabel, lumipb, version):
        f = ROOT.TFile("optimization_outputs_histos_"+version+"/"+filename)
        tot = f.Get("cutOpt_tot_"+channel)
        sys = f.Get("cutOpt_sys_"+channel)
        stat = f.Get("cutOpt_stat_"+channel)
        sys.SetLineColor(ROOT.kRed)
        sys.SetMarkerColor(ROOT.kRed)
        sys.SetLineStyle(3)
        sys.SetLineWidth(3);
        stat.SetLineColor(ROOT.kBlue)
        stat.SetMarkerColor(ROOT.kBlue)
        stat.SetLineStyle(2)
        stat.SetLineWidth(3);
        tot.SetLineWidth(3);


        c1 = ROOT.TCanvas("c1","",800,600)
        c1.SetLeftMargin(0.15);
        c1.SetBottomMargin(0.15);
        tot.Draw("h")
        tot.SetMinimum(0)
        tot.SetMaximum(0.08)
        tot.SetTitle("")
        sys.Draw("same h")
        stat.Draw("same h")

        leg = ROOT.TLegend(0.6,0.66,0.80,0.85);
        leg.SetTextSize(0.04);
        tot.GetYaxis().SetTitle("#Delta#sigma/#sigma");
        x_label = xlabel
        if "ph_pt" in xlabel:
            x_label = "p_{T}(#gamma) [GeV]"
        if "event_njets" in xlabel:
            x_label = "n_{jets}"
        if "event_nbjets" in xlabel:
            x_label = "n_{bjets}"
        if "met_met" in xlabel:
            x_label = "MET [GeV]"
        if "mw" in xlabel:
            x_label = "m_{W_{T}} [GeV]"
        if "m_ph_el" in xlabel:
            x_label = "x [GeV]"
        if "dR_gj_all" in xlabel:
            x_label = "dR(#gamma,jet) all"
        if "dR_glep" in xlabel:
            x_label = "dR(#gamma,lep)"
        if "mll" in xlabel:
            x_label = "m(l,l)"
        if "ph_HFT_MVA" in xlabel:
            x_label = "ph_HFT_MVA"

        tot.GetXaxis().SetTitle(x_label);
        if (channel == "mujets"):
            leg.SetHeader("#mu+jets");
        elif (channel == "ee"):
            leg.SetHeader("ee");
        elif (channel == "mumu"):
            leg.SetHeader("#mu#mu");
        else:
            leg.SetHeader("e+jets");
        leg.AddEntry(tot,"#Delta#sigma_{tot}/#sigma","l");
        leg.AddEntry(sys,"#Delta#sigma_{sys}/#sigma","l");
        leg.AddEntry(stat,"#Delta#sigma_{stat}/#sigma","l");
        leg.SetBorderSize(0)
        leg.Draw()
        lumi = ROOT.TLatex();
        lumi.SetNDC();
        lumi.SetTextAlign(12);
        lumi.SetTextFont(63);
        lumi.SetTextSizePixels(15);
        lumifb = lumipb/1e3
        lumi.DrawLatex(0.35,0.85, "#it{#scale[1.2]{ATLAS}} #bf{Internal}");
        lumi.DrawLatex(0.35,0.80,"#sqrt{s}=13 TeV, " + str(round(lumifb,2)) +" fb^{-1}");
        if not os.path.exists("optimization_outputs_pngs_"+version):
            os.makedirs("optimization_outputs_pngs_"+version)

        c1.SaveAs('optimization_outputs_pngs_'+version+'/'+filename+'_'+channel+'.eps')

class systematics:
    """____________________________________________________________________
    A class to calculate up down variations of systematics
    """
    def upDownCalc(self, up, down, nominal):
        self.up = up
        self.down = down
        self.nominal = nominal
        upper = (up - nominal)/nominal
        downer = (down - nominal)/nominal
        maxValue = max(abs(upper),abs(downer))
        return maxValue



class stackedHistogram:
  """____________________________________________________________________
  A class to plot a stacked histogram as a cross check
  """
  def createStack(self,filename, region,outputName,cutName, cutValue,lumipb,version):
    f = ROOT.TFile(filename)
    h_stack = ROOT.THStack("h_stack","Post Optimisation");
    h_ttgamma = f.Get("h_ttgamma")
    h_ttgamma.SetFillColor(632)
    h_ttbar = f.Get("h_ttbar")
    h_ttbar.SetFillColor(876)
    h_Wjets = f.Get("h_Wjets")
    h_Wjets.SetFillColor(393)
    h_Wgamma = f.Get("h_Wgamma")
    h_Wgamma.SetFillColor(401)
    h_Zjets = f.Get("h_Zjets")
    h_Zjets.SetFillColor(428)
    h_Zgamma = f.Get("h_Zgamma")
    h_Zgamma.SetFillColor(62)
    h_VV = f.Get("h_VV")
    h_VV.SetFillColor(95)
    h_ST = f.Get("h_ST")
    h_ST.SetFillColor(8)
    h_stack.Add(h_ST)
    h_stack.Add(h_VV)
    h_stack.Add(h_Zgamma)
    h_stack.Add(h_Zjets)
    h_stack.Add(h_Wgamma)
    h_stack.Add(h_Wjets)
    h_stack.Add(h_ttbar)
    h_stack.Add(h_ttgamma)
    h_data_2015 = f.Get("h_data_2015")
    h_data_2015.SetFillColor(1)
    h_data_2015.SetMarkerSize(1)
    h_data_2015.SetMarkerStyle(20)

    h_all_bkg = ROOT.TH1D(h_ttgamma)
    h_all_bkg.Add(h_ttbar)
    h_all_bkg.Add(h_Wjets)
    h_all_bkg.Add(h_Wgamma)
    h_all_bkg.Add(h_Zjets)
    h_all_bkg.Add(h_Zgamma)
    h_all_bkg.Add(h_ST)
    h_all_bkg.Add(h_VV)

    h_ttgamma_cutflow = ROOT.TH1D("ttgamma_cutflow","ttgamma_cutflow", 10, -0.5, 9.5)


    canvas = ROOT.TCanvas("canvas","stacked hists",0,0,800,800);
    canvas.SetFillColor(0);

    # Upper histogram plot is pad1
   # canvas.cd(1)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    pad1.SetGridx()
    pad1.Draw()
    # # Lower ratio plot is pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()

    pad1.cd()
    h_stack.SetTitle("")

    h_stack.Draw("HIST");
    h_data_2015.Draw("same PE");

    y_top = h_stack.GetYaxis()
    y_top.SetTitle("#events/GeV ")

    leg = ROOT.TLegend(0.42,0.58,0.62,0.88);
    leg.SetTextSize(0.022);
    if (region == "mujets"):
      leg.SetHeader("#mu+jets, "+cutName+": "+str(cutValue));
    elif (region == "ee"):
       leg.SetHeader("ee, "+cutName+": "+str(cutValue));
    elif (region == "mumu"):
       leg.SetHeader("#mu#mu, "+cutName+": "+str(cutValue));
    else:
      leg.SetHeader("e+jets, "+cutName+": "+str(cutValue));
    leg.AddEntry(h_data_2015,"data 2015","p");
    leg.AddEntry(h_ttgamma,"tt#gamma","f");
    leg.AddEntry(h_ttbar,"t#bar{t}","f");
    leg.AddEntry(h_Wgamma,"W#gamma","f");
    leg.AddEntry(h_Wjets,"W+jets","f");
    leg.AddEntry(h_Zgamma,"Z#gamma","f");
    leg.AddEntry(h_Zjets,"Z+jets","f");
    leg.AddEntry(h_VV,"VV","f");
    leg.AddEntry(h_ST,"ST","f");

    leg.SetBorderSize(0)
    leg.Draw()

    lumi = ROOT.TLatex();
    lumi.SetNDC();
    lumi.SetTextAlign(12);
    lumi.SetTextFont(63);
    lumi.SetTextSizePixels(15);
    lumifb = lumipb/1e3
    lumi.DrawLatex(0.15,0.85, "#it{#scale[1.2]{ATLAS}} #bf{Internal}");
    lumi.DrawLatex(0.15,0.80,"#sqrt{s}=13 TeV, " + str(round(lumifb,2)) +" fb^{-1}");

    ratio = h_data_2015.Clone("ratio")
    ratio.SetMarkerStyle(20)
    ratio.SetMinimum(0)
    ratio.SetMaximum(2)
    #ratio.Sumw2()
    ratio.SetStats(0)
    ratio.Divide(h_all_bkg)
    ratio.SetTitle("")
    y = ratio.GetYaxis()
    y.SetTitle("data/MC ")
    y.SetNdivisions(505)
    y.SetTitleSize(20)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(15)
    x = ratio.GetXaxis()
    if "ph_pt" in cutName:
      x_label = "p_{T}(#gamma) [GeV]"
    if "event_njets" in cutName:
      x_label = "n_{jets}"
    if "event_nbjets" in cutName:
      x_label = "n_{bjets}"
    if "met_met" in cutName:
      x_label = "MET [GeV]"
    if "mw" in cutName:
      x_label = "m_{W_{T}} [GeV]"
    if "m_ph_el" in cutName:
      x_label = "m(e,#gamma) [GeV]"
    if "dR_gj_all" in cutName:
      x_label = "dR(#gamma,jet) all"
    if "dR_glep" in cutName:
      x_label = "dR(#gamma,lep)"
    if "full_cuts" in cutName:
      x_label = "p_{T}(#gamma) [GeV]"
    if "mll" in cutName:
      x_label = "m(l,l)"
    if "ph_HFT_MVA" in cutName:
      x_label = "ph_HFT_MVA"
    x.SetTitle(x_label)
    x.SetTitleSize(20)
    x.SetTitleFont(43)
    x.SetTitleOffset(3.2)
    x.SetLabelFont(43)
    x.SetLabelSize(15)
    pad2.cd()
    ratio.Draw("ep")
    line = ROOT.TF1("fa1","1",-1000,1000);
    line.Draw("same")
    line.SetLineColor(ROOT.kGray);

    
    #canvas.cd(2);
    #h_stack.Draw("nostack,e1p");
    if not os.path.exists("stacked_pngs_"+version):
        os.makedirs("stacked_pngs_"+version)
    canvas.SaveAs("stacked_pngs_"+version+"/"+region+"_"+outputName+".eps")
    # h_stack = ROOT.THStack("h_stack","Stacked 1D histograms");
    # h_ttgamma_cutflow = f.Get("h_ttgamma_cutflow")
    # h_ttbar_cutflow = f.Get("h_ttbar_cutflow")
    # h_Wjets_cutflow = f.Get("h_Wjets_cutflow")
    # h_Wgamma_cutflow = f.Get("h_Wgamma_cutflow")
    # h_Zjets_cutflow = f.Get("h_Zjets_cutflow")
    # h_Zgamma_cutflow = f.Get("h_Zgamma_cutflow")
    # h_VV_cutflow = f.Get("h_VV_cutflow")
    # h_ST_cutflow = f.Get("h_ST_cutflow")

