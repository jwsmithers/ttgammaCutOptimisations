#!/usr/bin/env python
"""____________________________________________________________________
joshua.wyatt.smith@cern.ch                  

A script to do Cut optimizations.

Setup:
none

Input:
Ntuples with any variables.

Usage:
(Uncomment needed optimization at bottom of this script)
./cutOptimizations.py

Output:
Root histograms in optimization_outputs_histos
PNG's in optimization_outputs_pngs

"""
import argparse
print "----- loading libraries -----"
def _get_args():
  parser = argparse.ArgumentParser(description='ttgamma Cut Optimizations'
    ,epilog=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  return parser.parse_args()
_get_args()

import optimizationFunctions as OF
#import array
import ROOT
import glob
import sys
import os
from math import *
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0);
########################################################################

# Set lumi scale when doing the optimizations
# This is used as a scale factor when calculating stat and syst uncert.
# Some numbers for convenience...
lumi_2015 = 3212.96;
lumi_2016 = 24799.9;
lumi_comb = 36470;
# The one that actually get's used
lumi = lumi_2015
# total_events_channel = 199760 #For AF
# total_events_channel = 199000 #For AF
total_events_channel = 590000 #For fullSim
ttgamma_sample = "ttgamma_"
bkg_samples = [ "ttbar_",
    "enugamma",
    "munugamma",
    "taunugamma",
    "Wenu",
    "Wmunu",
    "Wtaunu",
    "eegamma",
    "mumugamma",
    "tautaugamma",
    "Zee",
    "Zmumu",
    "Ztautau",
    "VV",
    "ST_other",
    "ST_Wt_inclusive"]
# Specify some groups of samples
wjets = ["Wenu","Wmunu","Wtaunu"]
wgamma = ["enugamma","munugamma","taunugamma"]
zjets = ["Zee","Zmumu","Ztautau"]
zgamma = ["eegamma","mumugamma","tautaugamma"]
st = ["ST_other", "ST_Wt_inclusive"]

# What version of our ntuples are we using?
version = "v007"
if version == "v003":
  dataWildCard = "data_2015" #V003
elif version == "v004": 
  dataWildCard = "data15" #V004
elif version == "v007":
  dataWildCard = "data15" #V006


########################################################################
# Define some cutting functions 
def ph_drLepGamma(ph_drlept,ph_drlept_cutValue):
  if ph_drlept > ph_drlept_cutValue: 
    return True

def ph_mgammleptZWindow(ph_mgammalept,ph_mgammalept_cutValue):
  if abs(ph_mgammalept - 91188) >= ph_mgammalept_cutValue: 
    return True

def ph_HFT_MVA_cut(ph_HFT_MVA, ph_HFT_MVA_cutValue):
  if ph_HFT_MVA >= ph_HFT_MVA_cutValue: 
    return True

def get_mv2c10_bin(mv2c10):
  jets = []
  for i in mv2c10:
    jets.append(abs(i))
  leading_jet_mv2c10 = max(jets)
  # sub_leading_jet_mv2c10 = sorted(jets, reverse=True)[1]
  if leading_jet_mv2c10 < 0.1758475:
    return 1
  if leading_jet_mv2c10 >= 0.1758475 and leading_jet_mv2c10 < 0.645925:
    return 2
  if leading_jet_mv2c10 >= 0.645925 and leading_jet_mv2c10 < 0.8244273:
    return 3
  if leading_jet_mv2c10 >= 0.8244273 and leading_jet_mv2c10 < 0.934906:
    return 4
  if leading_jet_mv2c10 >= 0.934906:
    return 5

########################################################################

def doOptimization(chain, cutValue, region, cutName, allCuts = False, sigOrBkgOrData = "bkg"):
  # These histograms store the offline cutflows for each sample
    h = ROOT.TH1D("temp","temp dist", 100, 0, 100)
    # Histogramms for cutflows
    cutflowBins = 11
    cutflow_start = -0.5
    cutflow_end = 10.5
    h_ttgamma_cutflow = ROOT.TH1D("ttgamma_cutflow","ttgamma_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_ttbar_cutflow = ROOT.TH1D("ttbar_cutflow","ttbar_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_Wjets_cutflow = ROOT.TH1D("Wjets_cutflow","Wjets_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_Wgamma_cutflow = ROOT.TH1D("Wgamma_cutflow","Wgamma_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_Zjets_cutflow = ROOT.TH1D("Zjets_cutflow","Zjets_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_Zgamma_cutflow = ROOT.TH1D("Zgamma_cutflow","Zgamma_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_VV_cutflow = ROOT.TH1D("VV_cutflow","VV_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_ST_cutflow = ROOT.TH1D("ST_cutflow","ST_cutflow", cutflowBins, cutflow_start, cutflow_end)
    h_data_2015_cutflow = ROOT.TH1D("data_2015_cutflow","data_2015_cutflow", cutflowBins, cutflow_start, cutflow_end)

    # h_photon_SF_UP = ROOT.TH1D("photon_UP","photon_UP", 100, 0, 100)
    # h_photon_SF_DOWN = ROOT.TH1D("photon_DOWN","photon_DOWN", 100, 0, 100)

    # Histograms for stacked plots
    h_hfake_ttbar = ROOT.TH1D("hfake_ttbar","hfake_ttbar", 100, 0, 100)
    h_hfake_ST = ROOT.TH1D("hfake_ST","hfake_ST", 100, 0, 100)
    h_hfake_VV = ROOT.TH1D("hfake_VV","hfake_VV", 100, 0, 100)
    h_hfake_Zjets = ROOT.TH1D("hfake_Zjets","hfake_Zjets", 100, 0, 100)
    h_hfake_Zgamma = ROOT.TH1D("hfake_Zgamma","hfake_Zgamma", 100, 0, 100)
    h_hfake_Wjets = ROOT.TH1D("hfake_Wjets","hfake_Wjets", 100, 0, 100)
    h_hfake_Wgamma = ROOT.TH1D("hfake_Wgamma","hfake_Wgamma", 100, 0, 100)

    h_efake_ttbar = ROOT.TH1D("efake_ttbar","efake_ttbar", 100, 0, 100)
    h_efake_ST = ROOT.TH1D("efake_ST","efake_ST", 100, 0, 100)
    h_efake_VV = ROOT.TH1D("efake_VV","efake_VV", 100, 0, 100)
    h_efake_Zjets = ROOT.TH1D("efake_Zjets","efake_Zjets", 100, 0, 100)
    h_efake_Zgamma = ROOT.TH1D("efake_Zgamma","efake_Zgamma", 100, 0, 100)
    h_efake_Wjets = ROOT.TH1D("efake_Wjets","efake_Wjets", 100, 0, 100)
    h_efake_Wgamma = ROOT.TH1D("efake_Wgamma","efake_Wgamma", 100, 0, 100)

    h_other_ttbar = ROOT.TH1D("other_ttbar","other_ttbar", 100, 0, 100)
    h_other_ST = ROOT.TH1D("other_ST","other_ST", 100, 0, 100)
    h_other_VV = ROOT.TH1D("other_VV","other_VV", 100, 0, 100)
    h_other_Zjets = ROOT.TH1D("other_Zjets","other_Zjets", 100, 0, 100)
    h_other_Zgamma = ROOT.TH1D("other_Zgamma","other_Zgamma", 100, 0, 100)
    h_other_Wjets = ROOT.TH1D("other_Wjets","other_Wjets", 100, 0, 100)
    h_other_Wgamma = ROOT.TH1D("other_Wgamma","other_Wgamma", 100, 0, 100)

    # These bins and histograms store the event count that will be used to create stacked histos
    # Unfortunately, the bins need to be adjusted depending on the veriable to optimize.
    # stacked_nbins = 50
    # nbins_start = 0
    # nbins_end = 200
    # stacked_nbins = 40
    # nbins_start = 0
    # nbins_end = 7
    stacked_nbins = 6
    nbins_start = -0.5
    nbins_end = 5.5
    # stacked_nbins = 20
    # nbins_start = 0
    # nbins_end = 1
    if not os.path.exists("stacked_"+cutName+"_"+version):
        os.makedirs("stacked_"+cutName+"_"+version)
    h_ttgamma = ROOT.TH1D("h_ttgamma","h_ttgamma", stacked_nbins, nbins_start, nbins_end)
    h_ttbar = ROOT.TH1D("h_ttbar","h_ttbar", stacked_nbins, nbins_start, nbins_end)
    h_Wjets = ROOT.TH1D("h_Wjets","h_Wjets", stacked_nbins, nbins_start, nbins_end)
    h_Zjets = ROOT.TH1D("h_Zjets","h_Zjets", stacked_nbins, nbins_start, nbins_end)
    h_Wgamma = ROOT.TH1D("h_Wgamma","h_Wgamma", stacked_nbins, nbins_start, nbins_end)
    h_Zgamma = ROOT.TH1D("h_Zgamma","h_Zgamma", stacked_nbins, nbins_start, nbins_end)
    h_VV = ROOT.TH1D("h_VV","h_VV", stacked_nbins, nbins_start, nbins_end)
    h_ST = ROOT.TH1D("h_ST","h_ST", stacked_nbins, nbins_start, nbins_end)
    h_data_2015 = ROOT.TH1D("h_data_2015","h_data_2015", stacked_nbins, nbins_start, nbins_end)


    # Define the samples that contain prompt photons
    prompt = [
    ttgamma_sample,
    "enugamma",
    "munugamma",
    "taunugamma",
    "eegamma",
    "mumugamma",
    "tautaugamma"]
    # And the samples for which the prompt will/may overlap
    overlap = [ 
    "ttbar_",
    "Wenu",
    "Wmunu",
    "Wtaunu",
    "Zee",
    "Zmumu",
    "Ztautau"]

    # Initiate counters
    total = 0
    ngoodPhDefs = 0
    overlap_removal = 0
    nbjets = 0
    HFT_MVA=0
    met = 0
    mwcut = 0
    mphel = 0
    dR_gj = 0
    dR_gl = 0

    mv2c10_85=0
    mv2c10_77=0
    mv2c10_70=0
    mv2c10_60=0
    mv2c10_le60=0

    entries = chain.GetEntries()
    print str(entries), " entries for ", str(chain.GetCurrentFile().GetName())
    # Start loop over the chain that was passed to this function.
    for i in chain:
      if i.selph_index1 < 0: continue
      filename_string=chain.GetCurrentFile().GetName()

      # Define weights, if data, all of them are 1.
      if dataWildCard in filename_string:
        totalWeight = 1
        weight_no_btag = 1
        weightToUse = 1
      else:
        # The lumi variable is used here. One could also use i.event_lumi
        totalWeight = i.weight_mc * i.weight_pileup \
          * i.weight_bTagSF_77*i.ph_SF_eff[i.selph_index1]* i.weight_leptonSF * i.weight_jvt * lumi \
          * i.event_norm;
        weight_no_btag = totalWeight/i.weight_bTagSF_77
        weightToUse=totalWeight

      total = total + 1 * weight_no_btag

      # 1 good photon
      if i.event_ngoodphotons != 1:continue
      ngoodPhDefs = ngoodPhDefs + 1 * weight_no_btag
  
      #overlap removal
      if any(x1 in filename_string for x1 in prompt):
        if (i.event_photonorigin >= 10): continue
      if any(x2 in filename_string for x2 in overlap):
        if ( i.event_photonorigin < 10 ): continue
      overlap_removal = overlap_removal + 1 * weight_no_btag

      # Photon_UP= i.weight_photonSF_ID_UP
      # Photon_DOWN= i.weight_photonSF_ID_DOWN

      #####################################################################
      # Channel specific optimization, to be set after
      # each iteration for maximized cut
      # Some easier ways to keep track:
      if region == "ejets":
        ph_mgammalept_cutValue = 5000
        ph_drlept_cutValue = 0.7
      if region == "mujets":
        met_cutValue = 0
        ph_drlept_cutValue = 0.7

      #ph_HFT_MVA_cutValue= cutValue  
      ph_HFT_MVA_cutValue = 0

      # Cuts  
      if i.event_nbjets77 < 1: continue #1 or more bjets
      nbjets = nbjets +1 * totalWeight
      ph_HFT_MVA_bool = ph_HFT_MVA_cut(i.ph_HFT_MVA[i.selph_index1],ph_HFT_MVA_cutValue)
      if ph_HFT_MVA_bool != True: continue 
      HFT_MVA = HFT_MVA + 1 * totalWeight
      if (region == "ejets"):
        ph_mgammleptZWindow_bool = ph_mgammleptZWindow(i.ph_mgammalept[i.selph_index1],ph_mgammalept_cutValue)
        if ph_mgammleptZWindow_bool != True:continue # A ph,lep mass - Zmass window cut
        mphel = mphel + 1 * totalWeight
      ph_drlept_bool = ph_drLepGamma(i.ph_drlept[i.selph_index1],ph_drlept_cutValue)
      if ph_drlept_bool != True: continue # dr_gammalept cut
      dR_gl = dR_gl +1 * totalWeight
      #####################################################################
      # We only fill at the end of the loop. This could be optimized
      h.Fill(1,weightToUse)

       # h_photon_SF_UP.Fill(1,weightToUse*Photon_UP)
       # h_photon_SF_DOWN.Fill(1,weightToUse*Photon_DOWN)

      cut_list = [total, ngoodPhDefs,overlap_removal, nbjets,HFT_MVA, mphel, dR_gl]
      cut_list_names = ["total","ngoodPhDefs", "overlap_removal", "nbjets","HFT_MVA", "mphel", "dR_gl"]

       # Purely for creating stacked histos when optimizing.
      if cutName == "ph_pt":
        plotvar = i.ph_pt[i.selph_index1]/1e3
      if cutName == "event_njets":
        plotvar = i.event_njets
      if cutName == "event_nbjets77":
        plotvar = i.event_nbjets
      if cutName == "ph_HFT_MVA":
        plotvar = i.ph_HFT_MVA[i.selph_index1]
      if cutName == "met_met":
        plotvar = i.met_met/1e3
      if cutName == "mw":
        plotvar = i.event_mwt/1e3
      if cutName == "ph_mgammalept":
        plotvar = i.ph_mgammalept[i.selph_index1]/1e3
      if cutName == "ph_dralljet":
        plotvar = i.ph_dralljet[i.selph_index1]
      if cutName == "ph_drlept":
        plotvar = i.ph_drlept[i.selph_index1]
      if cutName == "full_cuts":
        plotvar = i.ph_pt[i.selph_index1]/1e3


      highest_jet_mv2c10 = get_mv2c10_bin(i.jet_mv2c10)
      if highest_jet_mv2c10==0:
        mv2c10_85=1*totalWeight + mv2c10_85
      if highest_jet_mv2c10==1:
        mv2c10_77=1*totalWeight+ mv2c10_77
      if highest_jet_mv2c10==2:
        mv2c10_70=1*totalWeight + mv2c10_70
      if highest_jet_mv2c10==3:
        mv2c10_60=1*totalWeight + mv2c10_60
      if highest_jet_mv2c10==4:
        mv2c10_le60=1*totalWeight + mv2c10_le60

      working_points_labels=["85<Eff","85<Eff<77","77<Eff70","70<Eff<60","Eff<60"]
      working_points=[mv2c10_85,mv2c10_77,mv2c10_70,mv2c10_60,mv2c10_le60]

      # Fill the cutflows and the  histos used for stacking
      if ttgamma_sample in filename_string:
        # h_ttgamma.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_ttgamma.SetBinContent(wp,working_points[wp-1])
          h_ttgamma.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        h_ttgamma.SetBinContent(highest_jet_mv2c10,1)
        for m in range(1, len(cut_list)+1):
          h_ttgamma_cutflow.SetBinContent(m,cut_list[m-1])
          h_ttgamma_cutflow.GetXaxis().SetBinLabel(m,cut_list_names[m-1])

      if "ttbar_" in filename_string:
        # h_ttbar.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_ttbar.SetBinContent(wp,working_points[wp-1])
          h_ttbar.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for n in range(1, len(cut_list)+1):
          h_ttbar_cutflow.SetBinContent(n,cut_list[n-1])
          h_ttbar_cutflow.GetXaxis().SetBinLabel(n,cut_list_names[n-1])
        if i.event_photonorigin==10:
          h_hfake_ttbar.Fill(1,weightToUse)
        elif i.event_photonorigin==20:
          h_efake_ttbar.Fill(1,weightToUse)
        else:
          h_other_ttbar.Fill(1,weightToUse)

      if any(wj in filename_string for wj in wjets):
        # h_Wjets.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_Wjets.SetBinContent(wp,working_points[wp-1])
          h_Wjets.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for p in range(1, len(cut_list)+1):
          h_Wjets_cutflow.SetBinContent(p,cut_list[p-1])
          h_Wjets_cutflow.GetXaxis().SetBinLabel(p,cut_list_names[p-1])
        if i.event_photonorigin==10:
          h_hfake_Wjets.Fill(1,weightToUse)
        elif i.event_photonorigin==20:
          h_efake_Wjets.Fill(1,weightToUse)
        else:
          h_other_Wjets.Fill(1,weightToUse)

      if any(wg in filename_string for wg in wgamma):
        # h_Wgamma.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_Wgamma.SetBinContent(wp,working_points[wp-1])
          h_Wgamma.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for q in range(1, len(cut_list)+1):
          h_Wgamma_cutflow.SetBinContent(q,cut_list[q-1])
          h_Wgamma_cutflow.GetXaxis().SetBinLabel(q,cut_list_names[q-1])
        if i.event_photonorigin==10:
          h_hfake_Wgamma.Fill(1,weightToUse)
        elif i.event_photonorigin==20:
          h_efake_Wgamma.Fill(1,weightToUse)
        else:
          h_other_Wgamma.Fill(1,weightToUse)

      if any(zj in filename_string for zj in zjets):
        # h_Zjets.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_Zjets.SetBinContent(wp,working_points[wp-1])
          h_Zjets.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for r in range(1, len(cut_list)+1):
          h_Zjets_cutflow.SetBinContent(r,cut_list[r-1])
          h_Zjets_cutflow.GetXaxis().SetBinLabel(r,cut_list_names[r-1])
        if i.event_photonorigin==10:
          h_hfake_Zjets.Fill(1,weightToUse)
        elif i.event_photonorigin==20:
          h_efake_Zjets.Fill(1,weightToUse)
        else:
          h_other_Zjets.Fill(1,weightToUse)

      if any(zg in filename_string for zg in zgamma):
        # h_Zgamma.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_Zgamma.SetBinContent(wp,working_points[wp-1])
          h_Zgamma.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for s in range(1, len(cut_list)+1):
          h_Zgamma_cutflow.SetBinContent(s,cut_list[s-1])
          h_Zgamma_cutflow.GetXaxis().SetBinLabel(s,cut_list_names[s-1])
        if i.event_photonorigin==10:
          h_hfake_Zgamma.Fill(1,weightToUse)
        elif i.event_photonorigin==20:
          h_efake_Zgamma.Fill(1,weightToUse)
        else:
          h_other_Zgamma.Fill(1,weightToUse)

      if "VV" in filename_string:
        # h_VV.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_VV.SetBinContent(wp,working_points[wp-1])
          h_VV.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for t in range(1, len(cut_list)+1):
          h_VV_cutflow.SetBinContent(t,cut_list[t-1])
          h_VV_cutflow.GetXaxis().SetBinLabel(t,cut_list_names[t-1])
        if i.event_photonorigin==10:
          h_hfake_VV.Fill(1,weightToUse)
        elif i.event_photonorigin==20:
          h_efake_VV.Fill(1,weightToUse)
        else:
          h_other_VV.Fill(1,weightToUse)

      if "ST" in filename_string:
        # h_ST.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_ST.SetBinContent(wp,working_points[wp-1])
          h_ST.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for u in range(1, len(cut_list)+1):
          h_ST_cutflow.SetBinContent(u,cut_list[u-1])
          h_ST_cutflow.GetXaxis().SetBinLabel(u,cut_list_names[u-1])
        if i.event_photonorigin==10:
          h_hfake_ST.Fill(1,weightToUse)
        elif i.event_photonorigin==20:
          h_efake_ST.Fill(1,weightToUse)
        else:
          h_other_ST.Fill(1,weightToUse)

      if dataWildCard in filename_string:
        # h_data_2015.Fill(plotvar, weightToUse)
        for wp in range(1, len(working_points)+1):
          h_data_2015.SetBinContent(wp,working_points[wp-1])
          h_data_2015.GetXaxis().SetBinLabel(wp,working_points_labels[wp-1])

        for v in range(1, len(cut_list)+1):
          h_data_2015_cutflow.SetBinContent(v,cut_list[v-1])
          h_data_2015_cutflow.GetXaxis().SetBinLabel(v,cut_list_names[v-1])

    ########### Potential to optimize #########
    nominal = h.Integral()

    # We can also get each individual contribution
    hfakes_ttbar = h_hfake_ttbar.Integral()
    hfakes_Wjets = h_hfake_Wjets.Integral()
    hfakes_Wgamma = h_hfake_Wgamma.Integral()
    hfakes_Zjets = h_hfake_Zjets.Integral()
    hfakes_Zgamma = h_hfake_Zgamma.Integral()
    hfakes_ST = h_hfake_ST.Integral()
    hfakes_VV = h_hfake_VV.Integral()

    efakes_ttbar = h_efake_ttbar.Integral()
    efakes_Wjets = h_efake_Wjets.Integral()
    efakes_Wgamma = h_efake_Wgamma.Integral()
    efakes_Zjets = h_efake_Zjets.Integral()
    efakes_Zgamma = h_efake_Zgamma.Integral()
    efakes_ST = h_efake_ST.Integral()
    efakes_VV = h_efake_VV.Integral()

    other_ttbar = h_other_ttbar.Integral()
    other_Wjets = h_other_Wjets.Integral()
    other_Wgamma = h_other_Wgamma.Integral()
    other_Zjets = h_other_Zjets.Integral()
    other_Zgamma = h_other_Zgamma.Integral()
    other_ST = h_other_ST.Integral()
    other_VV = h_other_VV.Integral()

    # Create some dictionaries for later use
    hfakes = {"ttbar":hfakes_ttbar,"Wjets":hfakes_Wjets,"Wgamma":hfakes_Wgamma,
             "Zjets":hfakes_Zjets,"Zgamma":hfakes_Zgamma,"ST":hfakes_ST,"VV":hfakes_VV}

    efakes = {"ttbar":efakes_ttbar,"Wjets":efakes_Wjets,"Wgamma":efakes_Wgamma,
             "Zjets":efakes_Zjets,"Zgamma":efakes_Zgamma,"ST":efakes_ST,"VV":efakes_VV}

    other = {"ttbar":other_ttbar,"Wjets":other_Wjets,"Wgamma":other_Wgamma,
             "Zjets":other_Zjets,"Zgamma":other_Zgamma,"ST":other_ST,"VV":other_VV}

    # otherbkgs = h_otherbkg.Integral()

    # photon_up = h_photon_SF_UP.Integral()
    # photon_down = h_photon_SF_DOWN.Integral()
    # ph_ID_systematic = OF.systematics().upDownCalc(up=photon_up,
    #   down=photon_down,nominal=nominal) 
    while True:
      try:
        f_root = ROOT.TFile.Open("stacked_"+cutName+"_"+version+'/'+region+"_"+cutName+'_'+str(cutValue)+'_yields_and_cutflows.root','UPDATE')
        break
      except OSError, e:
        if e.errno != 17:
            raise   
        # time.sleep might help here
        pass
    f_root.cd()
    if sigOrBkgOrData == "sig":
      h_ttgamma.Write()
      h_ttgamma_cutflow.Write()
    if sigOrBkgOrData == "bkg":
      if "ttbar_" in filename_string:
        h_ttbar_cutflow.Write()
        h_ttbar.Write()

      if any(wj in filename_string for wj in wjets):
        h_Wjets.Write()
        h_Wjets_cutflow.Write()

      if any(wg in filename_string for wg in wgamma):
        h_Wgamma.Write()
        h_Wgamma_cutflow.Write()

      if any(zj in filename_string for zj in zjets):
        h_Zjets.Write()
        h_Zjets_cutflow.Write()

      if any(zg in filename_string for zg in zgamma):
        h_Zgamma.Write()
        h_Zgamma_cutflow.Write()

      if "VV" in filename_string:
        h_VV.Write()
        h_VV_cutflow.Write()

      if "ST" in filename_string:
        h_ST.Write()
        h_ST_cutflow.Write()

    if sigOrBkgOrData == "data":
      h_data_2015.Write()
      h_data_2015_cutflow.Write()
    f_root.Close()

    # Return the nominal signal count, and dictiionaries of the backgrounds
    return nominal, hfakes, efakes, other
########################################################################

def main(cutName,region,cut_range, allCuts=False):
  if version == "v003":
    path =  "/eos/atlas/user/j/jwsmith/reprocessedNtuples/v003/ntuples/SR_"+region+"/"
  elif version == "v004":
    path = "/eos/atlas/user/j/jwsmith/reprocessedNtuples/v004/ntuples/SR_"+region+"_filled/"
  elif version == "v007":
    path = "/eos/atlas/user/c/caudron/TtGamma_ntuples/v007/SR1/"+region+"/"

  if not os.path.exists(path):
      print "Uhh have you got your paths correct?"
      return

  cutOpt_stat = ROOT.TH1F("cutOpt_stat_"+region, 'cutOpt_stat_'+region, 
    len(cut_range)-1, cut_range[0], cut_range[-1])
  cutOpt_sys = ROOT.TH1F("cutOpt_sys_"+region, 'cutOpt_sys_'+region, 
    len(cut_range)-1, cut_range[0], cut_range[-1])
  cutOpt_tot = ROOT.TH1F("cutOpt_tot_"+region, 'cutOpt_tot_'+region, 
    len(cut_range)-1, cut_range[0], cut_range[-1])

  #Define bkg chains and dict
  bkg_dict={}
  bkg_chains = [ "ttbar_",
    "Wjets",
    "Wgamma",
    "Zjets",
    "Zgamma",
    "VV",
    "ST"]

  for c in range(0,len(cut_range)):
    signalChain = ROOT.TChain("nominal")
    dataChain = ROOT.TChain("nominal")
    
    for bchain in range(0,len(bkg_chains)):
      bkg_dict[bkg_chains[bchain]] = ROOT.TChain("nominal") 

    ntuples = glob.glob(path+"*")
    for i in ntuples:
      if ttgamma_sample in i:
        signalChain.Add(i)
      elif dataWildCard in i:
        dataChain.Add(i)
      else:
        for s_bkg in bkg_samples:
          if s_bkg in i:
            if s_bkg in wgamma:
              bkg_dict["Wgamma"].Add(i)
            elif s_bkg in wjets:
              bkg_dict["Wjets"].Add(i)
            elif s_bkg in zgamma:
              bkg_dict["Zgamma"].Add(i)
            elif s_bkg in zjets:
              bkg_dict["Zjets"].Add(i)
            elif s_bkg in st:
              bkg_dict["ST"].Add(i)
            else:
              bkg_dict[s_bkg].Add(i)

    print "\n+++++++++++++++++++++++++++++++++++++++"
    # Need to "optimise" to get cutflow for data. Don't actually use it for anything
    # If switching to 36 fb-1 and CR1, REMEMBER TO COMMENT THIS OUT!
    optimization_data = doOptimization(dataChain, cut_range[c], region,cutName, allCuts, sigOrBkgOrData = "data")

    optimization_signal = doOptimization(signalChain, cut_range[c], region,cutName, allCuts,sigOrBkgOrData = "sig")
    N_signal = float(optimization_signal[0])

    # Some counters
    N_background = 0
    hfake_bkg_ttbar = 0
    hfake_bkg_Wjets = 0
    hfake_bkg_Wgamma = 0
    hfake_bkg_Zjets = 0
    hfake_bkg_Zgamma = 0
    hfake_bkg_VV = 0
    hfake_bkg_ST = 0

    efake_bkg_ttbar = 0
    efake_bkg_Wjets = 0
    efake_bkg_Wgamma = 0
    efake_bkg_Zjets = 0
    efake_bkg_Zgamma = 0
    efake_bkg_VV = 0
    efake_bkg_ST = 0

    other_bkg_ttbar = 0
    other_bkg_Wjets = 0
    other_bkg_Wgamma = 0
    other_bkg_Zjets = 0
    other_bkg_Zgamma = 0
    other_bkg_VV = 0
    other_bkg_ST = 0

    for bkg in bkg_dict:
      optimization_bkg = doOptimization(bkg_dict[bkg], cut_range[c], region,cutName, allCuts,sigOrBkgOrData = "bkg")
      N_background = N_background + float(optimization_bkg[0])
      hfake_bkg_ttbar = hfake_bkg_ttbar + float(optimization_bkg[1]["ttbar"])
      hfake_bkg_Wjets = hfake_bkg_Wjets + float(optimization_bkg[1]["Wjets"])
      hfake_bkg_Wgamma = hfake_bkg_Wgamma + float(optimization_bkg[1]["Wgamma"])
      hfake_bkg_Zjets = hfake_bkg_Zjets + float(optimization_bkg[1]["Zjets"])
      hfake_bkg_Zgamma = hfake_bkg_Zgamma + float(optimization_bkg[1]["Zgamma"])
      hfake_bkg_VV = hfake_bkg_VV + float(optimization_bkg[1]["VV"])
      hfake_bkg_ST = hfake_bkg_ST + float(optimization_bkg[1]["ST"])

      efake_bkg_ttbar = efake_bkg_ttbar + float(optimization_bkg[2]["ttbar"])
      efake_bkg_Wjets = efake_bkg_Wjets + float(optimization_bkg[2]["Wjets"])
      efake_bkg_Wgamma = efake_bkg_Wgamma + float(optimization_bkg[2]["Wgamma"])
      efake_bkg_Zjets = efake_bkg_Zjets + float(optimization_bkg[2]["Zjets"])
      efake_bkg_Zgamma = efake_bkg_Zgamma + float(optimization_bkg[2]["Zgamma"])
      efake_bkg_VV = efake_bkg_VV + float(optimization_bkg[2]["VV"])
      efake_bkg_ST = efake_bkg_ST + float(optimization_bkg[2]["ST"])

      other_bkg_ttbar = other_bkg_ttbar + float(optimization_bkg[3]["ttbar"])
      other_bkg_Wjets = other_bkg_Wjets + float(optimization_bkg[3]["Wjets"])
      other_bkg_Wgamma = other_bkg_Wgamma + float(optimization_bkg[3]["Wgamma"])
      other_bkg_Zjets = other_bkg_Zjets + float(optimization_bkg[3]["Zjets"])
      other_bkg_Zgamma = other_bkg_Zgamma + float(optimization_bkg[3]["Zgamma"])
      other_bkg_VV = other_bkg_VV + float(optimization_bkg[3]["VV"])
      other_bkg_ST = other_bkg_ST + float(optimization_bkg[3]["ST"])


    # cutflowLine = region+", "+cutName+": "+str(cut_range[c]) +"\t"+str(N_signal)+ \
    #   "\t"+str(hfake_bkg)+"\t"+str(efake_bkg)+"\t"+str(other_bkg)+"\n"
    # with open("Cutflows.txt", "a") as myfile:
    #   myfile.write(cutflowLine)

    print "Slice : ",cut_range[c]
    print "Total signal events = ", N_signal 
    print "Total background events = ", N_background
    print "s/sqrt(b) = ", N_signal/sqrt(N_background)
    print "s/b = ", N_signal/N_background
    efficiency = N_signal/total_events_channel
    print "Efficiency = ", efficiency
    sigma = OF.equations()._sigma(signal=N_signal, 
      background=N_background,efficiency=efficiency,lumi=lumi)

    sigma_stat = OF.equations()._sigmaStat(signal=N_signal,
      background = N_background,efficiency=efficiency,lumi=lumi)

    print "x-section = ", sigma
    print "sigma_stat = ", sigma_stat

    # delta Sigma stat over sigma
    dsigma_stat_over_sigma = sigma_stat/sigma 
    print "sigma_STAT/sigma = ", dsigma_stat_over_sigma
    cutOpt_stat.SetBinContent(c,dsigma_stat_over_sigma)
########################################################################
    # Only used when optimizing individual variables
    # # delta Sigma sys over sigma

    # Set the bkg modelling uncertainty for each process
    bkg_sys_hfake = {"ttbar":hfake_bkg_ttbar*0.05,
    "Wjets":hfake_bkg_Wjets*0.05,
    "Wgamma":hfake_bkg_Wgamma*0.05,
    "Zjets":hfake_bkg_Zjets*0.05,
    "Zgamma":hfake_bkg_Zgamma*0.05,
    "ST":hfake_bkg_ST*0.05,
    "VV":hfake_bkg_VV*0.05}

    bkg_sys_efake = {"ttbar":efake_bkg_ttbar*0.10,
    "Wjets":efake_bkg_Wjets*0.10,
    "Wgamma":efake_bkg_Wgamma*0.10,
    "Zjets":efake_bkg_Zjets*0.10,
    "Zgamma":efake_bkg_Zgamma*0.10,
    "ST":efake_bkg_ST*0.10,
    "VV":efake_bkg_VV*0.10}

    bkg_sys_other = {"ttbar":other_bkg_ttbar*0.50,
    "Wjets":other_bkg_Wjets*0.50,
    "Wgamma":other_bkg_Wgamma*0.50,
    "Zjets":other_bkg_Zjets*0.50,
    "Zgamma":other_bkg_Zgamma*0.50,
    "ST":other_bkg_ST*0.50,
    "VV":other_bkg_VV*0.50}

    # We can also add other types of uncertainties
    if region == "ejets":
      # exp_sys = 0.051
      exp_sys = 0
    elif region == "mujets":
      # exp_sys = 0.047
      exp_sys = 0

    # Add a lumi uncertainty
    lumi_sys = 0.055

    # Pass all this to _sigmaSys function
    sigma_sys = OF.equations()._sigmaSys(sigma_=sigma, 
      signal=N_signal,
      bkg_sys_hfake=bkg_sys_hfake,
      bkg_sys_efake=bkg_sys_efake,
      bkg_sys_other=bkg_sys_other,
      efficiency=efficiency,
      expSys=exp_sys,
      lumiSys=lumi_sys, 
      lumi=lumi)


    print "sigma_sys = ", sigma_sys
    dsigma_sys_over_sigma = sigma_sys/sigma
    print "sigma_SYS/sigma = ", dsigma_sys_over_sigma 
    cutOpt_sys.SetBinContent(c,dsigma_sys_over_sigma)

    sigma_tot = sqrt(sigma_stat**2+sigma_sys**2)
    dsigma_tot_over_sigma = sigma_tot/sigma
    cutOpt_tot.SetBinContent(c,dsigma_tot_over_sigma)
    OF.stackedHistogram().createStack("stacked_"+cutName+"_"+version+'/'+region+"_"+cutName+'_'+str(cut_range[c])+'_yields_and_cutflows.root',
      region, 
      cutName+'_'+str(cut_range[c])+'_yields_and_cutflows', cutName, cut_range[c], lumi, version)

  # If all cuts is false, create the optimization plots for each cut
  if allCuts == False:
    bins = str(len(cut_range)-1)
    if not os.path.exists("optimization_outputs_histos_"+version):
        os.makedirs("optimization_outputs_histos_"+version)
    filename = cutName+"_"+bins+"_bins.root"
    cutOptFile = ROOT.TFile("optimization_outputs_histos_"+version+"/"+filename,"update")
    cutOptFile.cd()
    cutOpt_stat.Write()
    cutOpt_sys.Write()
    cutOpt_tot.Write()
    cutOptFile.Close()
    binmin = cutOpt_tot.GetMinimumBin()
    xmin = cutOpt_tot.GetXaxis().GetBinCenter(binmin)
    print '\n'
    print '=============== Minimum =============='
    print 'Region    = ',region
    print 'Cut       = ',cutName
    print "min x     = ",xmin
    print '===== Take with a pinch of salt ======'
    OF.plots().createPlot(filename,region,cutName, lumi, version)
########################################################################
regs = ["ejets","mujets"]
from joblib import Parallel, delayed
# Uncomment each region one at a time, run, adjust move on to the next.
# Remember to set cutValue in the doOptimization function!

#cutsPhHFTMVA=[0,0.05,0.10,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
# We can either run ejets and mujets in parallel:
#_ = Parallel(n_jobs=-1, verbose=5, backend="multiprocessing") \
#( delayed(main)(cutName="ph_HFT_MVA",region=i,cut_range=cutsPhHFTMVA) for i in regs) 

# Or we can run them sequentially:
# main(cutName="ph_HFT_MVA",region="ejets",cut_range=cutsPhHFTMVA)
# main(cutName="ph_HFT_MVA",region="mujets",cut_range=cutsPhHFTMVA)

# main(cutName="met_met",region="ejets",cut_range=[20,25,30,35,40,45,50,55,60,65,70])
# main(cutName="met_met",region="mujets",cut_range=[20,25,30,35,40,45,50,55,60,65,70])

# main(cutName="mw",region="ejets",cut_range=[20,25,30,35,40,45,50,55,60,65,70])
# main(cutName="mw",region="mujets",cut_range=[20,25,30,35,40,45,50,55,60,65,70])

# main(cutName="ph_mgammalept",region="ejets",cut_range=[3,4,5,6,7,8,9,10,11,12,13,14,15])


########################################################################
# These last two will run the whole chain once without optimizing. 
# Make sure all your cuts are specifed, notice allCuts = True.
# The cut_range is completely irrelavent here.
# Make sure youv'e set all your cuts correctly!

_ = Parallel(n_jobs=-1, verbose=5, backend="multiprocessing") \
( delayed(main)(cutName="full_cuts",region=i,cut_range=[1],allCuts = True) for i in regs) 
# main(cutName="full_cuts",region="ejets",cut_range=[1], allCuts = True)
# main(cutName="full_cuts",region="mujets",cut_range=[1], allCuts = True)
