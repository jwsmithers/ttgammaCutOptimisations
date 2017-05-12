#!/usr/bin/env python
"""____________________________________________________________________
joshua.wyatt.smith@cern.ch 				     			
A script to obtain offline cutflow tables.

Input:
Ntuples with any variables.

Usage:

./offline_cutFlows.py

Output:
Cutflow tables

"""
import argparse
def _get_args():
	parser = argparse.ArgumentParser(description='ttgamma offline cutFlows'
		,epilog=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
	return parser.parse_args()
_get_args()
#import array
import ROOT
import glob
import sys
import os
from math import *
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0);
########################

def main(inputs,cutLabels, channel,region,year,version):
	print "---++",version," Offline cuts"


	path = "./" 

	print "---+++ ",region,channel,year," \n"
	print "| Cuts | ttgamma | ttbar  | Wgamma | W+ jets | Zgamma | Z + jets | Single Top + gamma | Diboson+gamma | DATA (",year,") | Total (Sig+Bkg) | data/totalMC |"

	full_string = []
	f = ROOT.TFile.Open(path+"stacked_full_cuts_"+version+"/"+channel+"_full_cuts_1_yields_and_cutflows.root");
	for bin in range(0, len(cutLabels)):
		totalMC = 0
		totalDATA = 0
		for process in inputs:
	  		# open the file
			if (f == 0):
				print ("Error: cannot open file!\n");
				return;

			cutflow_process = f.Get(process);
			binContent = cutflow_process.GetBinContent(bin+1)
			if "ttgamma" in process:
				print " | " + cutLabels[bin] + " | " + str(binContent) + " | ", # TODO: Double check indices are correct
				totalMC = totalMC + binContent 

			elif "data" in process:
				totalDATA = totalDATA + binContent 
				print str(binContent) + " | ",
			else:	
				print str(binContent) + " | ",
				totalMC = totalMC + binContent 
		try:
			print totalMC," | ",totalDATA/(totalMC), " | "
		except ZeroDivisionError:
			print "- | - |"

inputs = ["ttgamma_cutflow", "ttbar_cutflow", "Wgamma_cutflow", "Wjets_cutflow", "Zgamma_cutflow", "Zjets_cutflow", "ST_cutflow", "VV_cutflow", "data_2015_cutflow"]

cut_list_names = ["total","ngoodphotons == 1","overlap_removal", "nbjets >=1","HFT_MVA", "abs(m(ph,e) - m(W))>5 GeV", "dR_gl > 0.7"]


main(inputs = inputs, cutLabels = cut_list_names, channel = "ejets", region = "SR",year ="2015", version="v007");
main(inputs = inputs, cutLabels = cut_list_names, channel = "mujets", region = "SR",year ="2015", version="v007");

