#!/usr/bin/en python

import ROOT, math
from ROOT import TFile, TH1F, TCanvas, gSystem, TAttFill, TLegend, TStyle, THStack, kViolet, kBlack, kAzure, kCyan, kGreen, kWhite, kOrange, kBlue, kRed
import CMSGraphics, CMS_lumi

# create output directory
OutputPath = 'ZZbb_plots'
gSystem.Exec("mkdir -p " + OutputPath)
print "Output directory created!"


# inputh file phath
inputFilePath = 'histos_4lbb/' 


# hist names
namelist = ['h1_M4L_CR4Lonly_4L','h1_MZ1_CR4Lonly_4L','h1_MZ2_CR4Lonly_4L','h1_eta4L_CR4Lonly_4L','h1_pt4L_CR4Lonly_4L','h1_M4L_4L','h1_MZ1_4L','h1_MZ2_4L','h1_eta4L_4L','h1_pt4L_4L','h1_methodBPT_M_4L','h1_methodBPT_PT_4L']


# read files
inFile_AllData = TFile.Open(inputFilePath + 'histos_AllData.root')
histos_AllData = []
histos_AllData.append(inFile_AllData.Get(namelist[0]))
histos_AllData.append(inFile_AllData.Get(namelist[1]))
histos_AllData.append(inFile_AllData.Get(namelist[2]))
histos_AllData.append(inFile_AllData.Get(namelist[3]))
histos_AllData.append(inFile_AllData.Get(namelist[4]))
histos_AllData.append(inFile_AllData.Get(namelist[5]))
histos_AllData.append(inFile_AllData.Get(namelist[6]))
histos_AllData.append(inFile_AllData.Get(namelist[7]))
histos_AllData.append(inFile_AllData.Get(namelist[8]))
histos_AllData.append(inFile_AllData.Get(namelist[9]))
histos_AllData.append(inFile_AllData.Get(namelist[10]))
histos_AllData.append(inFile_AllData.Get(namelist[11]))

inFile_HH4lbb = TFile.Open(inputFilePath + 'histos_HH4lbb.root')
histos_HH4lbb = []
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[0]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[1]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[2]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[3]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[4]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[5]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[6]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[7]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[8]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[9]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[10]))
histos_HH4lbb.append(inFile_HH4lbb.Get(namelist[11]))

inFile_ggH125 = TFile.Open(inputFilePath + 'histos_ggH125.root')
histos_ggH125 = []
histos_ggH125.append(inFile_ggH125.Get(namelist[0]))
histos_ggH125.append(inFile_ggH125.Get(namelist[1]))
histos_ggH125.append(inFile_ggH125.Get(namelist[2]))
histos_ggH125.append(inFile_ggH125.Get(namelist[3]))
histos_ggH125.append(inFile_ggH125.Get(namelist[4]))
histos_ggH125.append(inFile_ggH125.Get(namelist[5]))
histos_ggH125.append(inFile_ggH125.Get(namelist[6]))
histos_ggH125.append(inFile_ggH125.Get(namelist[7]))
histos_ggH125.append(inFile_ggH125.Get(namelist[8]))
histos_ggH125.append(inFile_ggH125.Get(namelist[9]))
histos_ggH125.append(inFile_ggH125.Get(namelist[10]))
histos_ggH125.append(inFile_ggH125.Get(namelist[11]))

inFile_VBFH125 = TFile.Open(inputFilePath + 'histos_VBFH125.root')
histos_VBFH125 = []
histos_VBFH125.append(inFile_VBFH125.Get(namelist[0]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[1]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[2]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[3]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[4]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[5]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[6]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[7]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[8]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[9]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[10]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[11]))

inFile_WplusH125 = TFile.Open(inputFilePath + 'histos_WplusH125.root')
histos_WplusH125 = []
histos_WplusH125.append(inFile_WplusH125.Get(namelist[0]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[1]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[2]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[3]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[4]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[5]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[6]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[7]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[8]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[9]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[10]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[11]))

inFile_WminusH125 = TFile.Open(inputFilePath + 'histos_WminusH125.root')
histos_WminusH125 = []
histos_WminusH125.append(inFile_WminusH125.Get(namelist[0]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[1]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[2]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[3]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[4]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[5]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[6]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[7]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[8]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[9]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[10]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[11]))

inFile_ZH125 = TFile.Open(inputFilePath + 'histos_ZH125.root')
histos_ZH125 = []
histos_ZH125.append(inFile_ZH125.Get(namelist[0]))
histos_ZH125.append(inFile_ZH125.Get(namelist[1]))
histos_ZH125.append(inFile_ZH125.Get(namelist[2]))
histos_ZH125.append(inFile_ZH125.Get(namelist[3]))
histos_ZH125.append(inFile_ZH125.Get(namelist[4]))
histos_ZH125.append(inFile_ZH125.Get(namelist[5]))
histos_ZH125.append(inFile_ZH125.Get(namelist[6]))
histos_ZH125.append(inFile_ZH125.Get(namelist[7]))
histos_ZH125.append(inFile_ZH125.Get(namelist[8]))
histos_ZH125.append(inFile_ZH125.Get(namelist[9]))
histos_ZH125.append(inFile_ZH125.Get(namelist[10]))
histos_ZH125.append(inFile_ZH125.Get(namelist[11]))

inFile_bbH125 = TFile.Open(inputFilePath + 'histos_bbH125.root')
histos_bbH125 = []
histos_bbH125.append(inFile_bbH125.Get(namelist[0]))
histos_bbH125.append(inFile_bbH125.Get(namelist[1]))
histos_bbH125.append(inFile_bbH125.Get(namelist[2]))
histos_bbH125.append(inFile_bbH125.Get(namelist[3]))
histos_bbH125.append(inFile_bbH125.Get(namelist[4]))
histos_bbH125.append(inFile_bbH125.Get(namelist[5]))
histos_bbH125.append(inFile_bbH125.Get(namelist[6]))
histos_bbH125.append(inFile_bbH125.Get(namelist[7]))
histos_bbH125.append(inFile_bbH125.Get(namelist[8]))
histos_bbH125.append(inFile_bbH125.Get(namelist[9]))
histos_bbH125.append(inFile_bbH125.Get(namelist[10]))
histos_bbH125.append(inFile_bbH125.Get(namelist[11]))

inFile_ttH125 = TFile.Open(inputFilePath + 'histos_ttH125.root')
histos_ttH125 = []
histos_ttH125.append(inFile_ttH125.Get(namelist[0]))
histos_ttH125.append(inFile_ttH125.Get(namelist[1]))
histos_ttH125.append(inFile_ttH125.Get(namelist[2]))
histos_ttH125.append(inFile_ttH125.Get(namelist[3]))
histos_ttH125.append(inFile_ttH125.Get(namelist[4]))
histos_ttH125.append(inFile_ttH125.Get(namelist[5]))
histos_ttH125.append(inFile_ttH125.Get(namelist[6]))
histos_ttH125.append(inFile_ttH125.Get(namelist[7]))
histos_ttH125.append(inFile_ttH125.Get(namelist[8]))
histos_ttH125.append(inFile_ttH125.Get(namelist[9]))
histos_ttH125.append(inFile_ttH125.Get(namelist[10]))
histos_ttH125.append(inFile_ttH125.Get(namelist[11]))

inFile_ggTo4e_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4e_Contin_MCFM701.root')
histos_ggTo4e_Contin_MCFM701 = []
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[0]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[1]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[2]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[3]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[4]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[5]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[6]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[7]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[8]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[9]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[10]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[11]))

inFile_ggTo4mu_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4mu_Contin_MCFM701.root')
histos_ggTo4mu_Contin_MCFM701 = []
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[0]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[1]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[2]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[3]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[4]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[5]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[6]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[7]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[8]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[9]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[10]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[11]))

inFile_ggTo4tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4tau_Contin_MCFM701.root')
histos_ggTo4tau_Contin_MCFM701 = []
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[0]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[1]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[2]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[3]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[4]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[5]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[6]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[7]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[8]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[9]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[10]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[11]))

inFile_ggTo2e2mu_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2e2mu_Contin_MCFM701.root')
histos_ggTo2e2mu_Contin_MCFM701 = []
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[0]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[1]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[2]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[3]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[4]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[5]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[6]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[7]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[8]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[9]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[10]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[11]))

inFile_ggTo2e2tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2e2tau_Contin_MCFM701.root')
histos_ggTo2e2tau_Contin_MCFM701 = []
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[0]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[1]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[2]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[3]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[4]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[5]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[6]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[7]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[8]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[9]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[10]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[11]))

inFile_ggTo2mu2tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2mu2tau_Contin_MCFM701.root')
histos_ggTo2mu2tau_Contin_MCFM701 = []
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[0]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[1]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[2]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[3]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[4]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[5]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[6]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[7]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[8]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[9]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[10]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[11]))

inFile_ZZTo4lext1 = TFile.Open(inputFilePath + 'histos_ZZTo4lext1.root')
histos_ZZTo4lext1 = []
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[0]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[1]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[2]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[3]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[4]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[5]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[6]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[7]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[8]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[9]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[10]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[11]))

inFile_TTZJets_M10_MLMext1 = TFile.Open(inputFilePath + 'histos_TTZJets_M10_MLMext1.root')
histos_TTZJets_M10_MLMext1 = []
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[0]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[1]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[2]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[3]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[4]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[5]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[6]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[7]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[8]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[9]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[10]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[11]))

inFile_TTZToLL_M1to10_MLM = TFile.Open(inputFilePath + 'histos_TTZToLL_M1to1O_MLM.root')
histos_TTZToLL_M1to10_MLM = []
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[0]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[1]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[2]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[3]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[4]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[5]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[6]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[7]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[8]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[9]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[10]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[11]))

inFile_TTWJetsToLNu = TFile.Open(inputFilePath + 'histos_TTWJetsToLNu.root')
histos_TTWJetsToLNu = []
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[0]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[1]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[2]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[3]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[4]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[5]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[6]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[7]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[8]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[9]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[10]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[11]))


print 'files read'


# --- do plots ---
i = 0 #counter for histos list
for name in namelist :

    canvas = TCanvas('canvas','canvas',800,600)

    hs = THStack('hs','')

#    integral_fondi = 0.

    # TTW
    histos_TTWJetsToLNu[i].SetFillColor(kBlue+3)
    histos_TTWJetsToLNu[i].SetLineColor(kBlack)
#    histos_TTWJetsToLNu[i].Rebin(6)
#    integral_fondi += histos_TTWJetsToLNu[i].Integral()
    hs.Add(histos_TTWJetsToLNu[i])    
    
    # gg->ZZ
    histo_ggZZ = histos_ggTo4e_Contin_MCFM701[i]
    histo_ggZZ.Add(histos_ggTo4mu_Contin_MCFM701[i])
    histo_ggZZ.Add(histos_ggTo4tau_Contin_MCFM701[i])
    histo_ggZZ.Add(histos_ggTo2e2mu_Contin_MCFM701[i])
    histo_ggZZ.Add(histos_ggTo2e2tau_Contin_MCFM701[i])
    histo_ggZZ.Add(histos_ggTo2mu2tau_Contin_MCFM701[i])
    histo_ggZZ.SetFillColor(kAzure-3)
    histo_ggZZ.SetLineColor(kBlue+2)
#    histo_ggZZ.Rebin(6)
#    integral_fondi += histo_ggZZ.Integral()
    hs.Add(histo_ggZZ)

    # TTZ
    histo_TTZ = histos_TTZJets_M10_MLMext1[i]
    histo_TTZ.Add(histos_TTZToLL_M1to10_MLM[i])
    histo_TTZ.SetFillColor(kGreen-2)
    histo_TTZ.SetLineColor(kGreen+4)
#   histo_TTZ.Rebin(6)
#   integral_fondi += histo_TTZ.Integral()
    hs.Add(histo_TTZ)


    # qq->ZZ
    histos_ZZTo4lext1[i].SetFillColor(kAzure+6)
    histos_ZZTo4lext1[i].SetLineColor(kAzure-6)
#   histos_ZZTo4lext1[i].Rebin(6)
#   integral_fondi += histos_ZZTo4lext1[i].Integral()
    hs.Add(histos_ZZTo4lext1[i])


    # SM Higgs
    histo_SMHiggs = histos_ggH125[i]
    histo_SMHiggs.Add(histos_VBFH125[i])
    histo_SMHiggs.Add(histos_WplusH125[i])
    histo_SMHiggs.Add(histos_WminusH125[i])
    histo_SMHiggs.Add(histos_ZH125[i])
    histo_SMHiggs.Add(histos_bbH125[i])
    histo_SMHiggs.Add(histos_ttH125[i])
    histo_SMHiggs.SetFillColor(kViolet+6)
    histo_SMHiggs.SetLineColor(kViolet+7)
#   histo_SMHiggs.Rebin(6)
#    integral_fondi += histo_SMHiggs.Integral()
    hs.Add(histo_SMHiggs)


    # HH->4lbb signal
    histos_HH4lbb[i].SetLineColor(kRed)
    histos_HH4lbb[i].SetLineWidth(2)
 #   histos_HH4lbb[i].Rebin(6)
    integral = histos_HH4lbb[i].Integral();
    print 'integrale ' 
    print integral
    histos_HH4lbb[i].Scale(integral * 1000.) 

    
    # ALL DATA
    histos_AllData[i].SetMarkerColor(kBlack)
    histos_AllData[i].SetMarkerStyle(20)
 #   histos_AllData[i].Rebin(6)
   

    # Draw all
    hs.SetMaximum(1.5 * max(hs.GetMaximum(),histos_AllData[i].GetMaximum()));
    hs.Draw('histo')
    histos_HH4lbb[i].Draw('histosame')
    histos_AllData[i].Draw('samep')  

    # print "-----------------"
    # print "fondi: " , integral_fondi
    # print "dati: ", histos_AllData[i].Integral()


    hs.GetXaxis().SetLabelFont(43)
    hs.GetXaxis().SetLabelSize(15)
    hs.GetXaxis().SetTitle(histos_HH4lbb[i].GetXaxis().GetTitle())
    hs.GetYaxis().SetTitleSize(20)
    hs.GetYaxis().SetTitleFont(43)
    hs.GetYaxis().SetTitleOffset(1.4)
    hs.GetYaxis().SetLabelFont(43)
    hs.GetYaxis().SetLabelSize(15)
    hs.GetYaxis().SetTitle(histos_HH4lbb[i].GetYaxis().GetTitle())
    hs.SetMinimum(0.)

    # legend
    legend = TLegend(0.74,0.64,0.94,0.87)
    legend.AddEntry(histos_HH4lbb[i],                "HH->4lbb x1000", "f")
    legend.AddEntry(histo_SMHiggs,                   "SM Higgs", "f")
    legend.AddEntry(histos_ZZTo4lext1[i],            "qq->ZZ",   "f")
    legend.AddEntry(histo_TTZ,                       "TTZ",      "f")
    legend.AddEntry(histo_ggZZ,                      "gg->ZZ",   "f")
    legend.AddEntry(histos_TTWJetsToLNu[i],          "TTW",      "f")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.SetTextFont(43)
    #legend.SetTextSize(20)
    legend.Draw()

    canvas.Update()

    #draw CMS and lumi text
    lumiText = '59.7 fb^{-1}'
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText      = "Preliminary"
    CMS_lumi.lumi_sqrtS     = lumiText + " (13 TeV)"
    CMS_lumi.cmsTextSize    = 0.6
    CMS_lumi.lumiTextSize   = 0.46
    CMS_lumi.extraOverCmsTextSize = 0.75
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(canvas, 0, 0)

    if 'M4L' in name:
        canvas.SetLogx()
 #       canvas.SetTickx();
        hs.GetXaxis().SetMoreLogLabels();
        hs.GetXaxis().SetNoExponent();

    canvas.Update()

    canvas.SaveAs(OutputPath + "/" + namelist[i] + ".png")
    
    i = i+1

    del canvas
