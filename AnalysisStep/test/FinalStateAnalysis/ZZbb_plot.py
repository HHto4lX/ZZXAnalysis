#!/usr/bin/en python

import ROOT, math
from ROOT import TFile, TH1F, TCanvas, gSystem, TAttFill, TLegend, TStyle, THStack, kViolet, kBlack, kAzure, kCyan, kGreen, kWhite, kOrange, kBlue

# create output directory
OutputPath = 'ZZbb_plots'
gSystem.Exec("mkdir -p " + OutputPath)
print "Output directory created!"


# inputh file phath
inputFilePath = '/eos/user/a/acappati/analysis_4lX/190626prod_histos/histos_4lbb_no4LmassCUT/' 


# hist names
namelist       = ['h1_M4L_4L','h1_MZ1_4L','h1_MZ2_4L']
namelist_label = ['m4l','mZ1','mZ2']

# read files
inFile_ggH125 = TFile.Open(inputFilePath + 'histos_ggH125.root')
histos_ggH125 = []
histos_ggH125.append(inFile_ggH125.Get(namelist[0]))
histos_ggH125.append(inFile_ggH125.Get(namelist[1]))
histos_ggH125.append(inFile_ggH125.Get(namelist[2]))

inFile_VBFH125 = TFile.Open(inputFilePath + 'histos_VBFH125.root')
histos_VBFH125 = []
histos_VBFH125.append(inFile_VBFH125.Get(namelist[0]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[1]))
histos_VBFH125.append(inFile_VBFH125.Get(namelist[2]))

inFile_WplusH125 = TFile.Open(inputFilePath + 'histos_WplusH125.root')
histos_WplusH125 = []
histos_WplusH125.append(inFile_WplusH125.Get(namelist[0]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[1]))
histos_WplusH125.append(inFile_WplusH125.Get(namelist[2]))

inFile_WminusH125 = TFile.Open(inputFilePath + 'histos_WminusH125.root')
histos_WminusH125 = []
histos_WminusH125.append(inFile_WminusH125.Get(namelist[0]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[1]))
histos_WminusH125.append(inFile_WminusH125.Get(namelist[2]))

inFile_ZH125 = TFile.Open(inputFilePath + 'histos_ZH125.root')
histos_ZH125 = []
histos_ZH125.append(inFile_ZH125.Get(namelist[0]))
histos_ZH125.append(inFile_ZH125.Get(namelist[1]))
histos_ZH125.append(inFile_ZH125.Get(namelist[2]))

inFile_bbH125 = TFile.Open(inputFilePath + 'histos_bbH125.root')
histos_bbH125 = []
histos_bbH125.append(inFile_bbH125.Get(namelist[0]))
histos_bbH125.append(inFile_bbH125.Get(namelist[1]))
histos_bbH125.append(inFile_bbH125.Get(namelist[2]))

inFile_ttH125 = TFile.Open(inputFilePath + 'histos_ttH125.root')
histos_ttH125 = []
histos_ttH125.append(inFile_ttH125.Get(namelist[0]))
histos_ttH125.append(inFile_ttH125.Get(namelist[1]))
histos_ttH125.append(inFile_ttH125.Get(namelist[2]))

inFile_ggTo4e_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4e_Contin_MCFM701.root')
histos_ggTo4e_Contin_MCFM701 = []
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[0]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[1]))
histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(namelist[2]))

inFile_ggTo4mu_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4mu_Contin_MCFM701.root')
histos_ggTo4mu_Contin_MCFM701 = []
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[0]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[1]))
histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(namelist[2]))

inFile_ggTo4tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4tau_Contin_MCFM701.root')
histos_ggTo4tau_Contin_MCFM701 = []
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[0]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[1]))
histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(namelist[2]))

inFile_ggTo2e2mu_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2e2mu_Contin_MCFM701.root')
histos_ggTo2e2mu_Contin_MCFM701 = []
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[0]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[1]))
histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(namelist[2]))

inFile_ggTo2e2tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2e2tau_Contin_MCFM701.root')
histos_ggTo2e2tau_Contin_MCFM701 = []
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[0]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[1]))
histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(namelist[2]))

inFile_ggTo2mu2tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2mu2tau_Contin_MCFM701.root')
histos_ggTo2mu2tau_Contin_MCFM701 = []
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[0]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[1]))
histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(namelist[2]))

inFile_ZZTo4lext1 = TFile.Open(inputFilePath + 'histos_ZZTo4lext1.root')
histos_ZZTo4lext1 = []
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[0]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[1]))
histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(namelist[2]))

inFile_TTZJets_M10_MLMext1 = TFile.Open(inputFilePath + 'histos_TTZJets_M10_MLMext1.root')
histos_TTZJets_M10_MLMext1 = []
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[0]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[1]))
histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(namelist[2]))

inFile_TTZToLL_M1to10_MLM = TFile.Open(inputFilePath + 'histos_TTZToLL_M1to1O_MLM.root')
histos_TTZToLL_M1to10_MLM = []
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[0]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[1]))
histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(namelist[2]))

inFile_TTWJetsToLNu = TFile.Open(inputFilePath + 'histos_TTWJetsToLNu.root')
histos_TTWJetsToLNu = []
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[0]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[1]))
histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(namelist[2]))


print 'files read'


# --- do plots ---
i = 0 #counter for histos list
for name in namelist :

    canvas = TCanvas('canvas','canvas',800,600)

    hs = THStack('hs','')

    # TTW
    histos_TTWJetsToLNu[i].SetFillColor(kBlue+3)
    histos_TTWJetsToLNu[i].SetLineColor(kBlack)
    hs.Add(histos_TTWJetsToLNu[i])    
    
    # gg->ZZ
    histos_ggTo4e_Contin_MCFM701[i].SetFillColor(kAzure-3)
    histos_ggTo4e_Contin_MCFM701[i].SetLineColor(kAzure-3)
    hs.Add(histos_ggTo4e_Contin_MCFM701[i])

    histos_ggTo4mu_Contin_MCFM701[i].SetFillColor(kAzure-3)
    histos_ggTo4mu_Contin_MCFM701[i].SetLineColor(kAzure-3)
    hs.Add(histos_ggTo4mu_Contin_MCFM701[i])

    histos_ggTo4tau_Contin_MCFM701[i].SetFillColor(kAzure-3)
    histos_ggTo4tau_Contin_MCFM701[i].SetLineColor(kAzure-3)
    hs.Add(histos_ggTo4tau_Contin_MCFM701[i])

    histos_ggTo2e2mu_Contin_MCFM701[i].SetFillColor(kAzure-3)
    histos_ggTo2e2mu_Contin_MCFM701[i].SetLineColor(kAzure-3)
    hs.Add(histos_ggTo2e2mu_Contin_MCFM701[i])

    histos_ggTo2e2tau_Contin_MCFM701[i].SetFillColor(kAzure-3)
    histos_ggTo2e2tau_Contin_MCFM701[i].SetLineColor(kAzure-3)
    hs.Add(histos_ggTo2e2tau_Contin_MCFM701[i])

    histos_ggTo2mu2tau_Contin_MCFM701[i].SetFillColor(kAzure-3)
    histos_ggTo2mu2tau_Contin_MCFM701[i].SetLineColor(kAzure-3)
    hs.Add(histos_ggTo2mu2tau_Contin_MCFM701[i])

    # qq->ZZ
    histos_ZZTo4lext1[i].SetFillColor(kAzure+8)
    histos_ZZTo4lext1[i].SetLineColor(kAzure+8)
    hs.Add(histos_ZZTo4lext1[i])

    # TTZ
    histos_TTZJets_M10_MLMext1[i].SetFillColor(kGreen+2)
    histos_TTZJets_M10_MLMext1[i].SetLineColor(kGreen+2)
    hs.Add(histos_TTZJets_M10_MLMext1[i])

    histos_TTZToLL_M1to10_MLM[i].SetFillColor(kGreen+2)
    histos_TTZToLL_M1to10_MLM[i].SetLineColor(kGreen+2)
    hs.Add(histos_TTZToLL_M1to10_MLM[i])

    # SM Higgs
    histos_ggH125[i].SetFillColor(kViolet+1)
    histos_ggH125[i].SetLineColor(kViolet+1)
    hs.Add(histos_ggH125[i])

    histos_VBFH125[i].SetFillColor(kViolet+1)
    histos_VBFH125[i].SetLineColor(kViolet+1)
    hs.Add(histos_VBFH125[i])
    
    histos_WplusH125[i].SetFillColor(kViolet+1)
    histos_WplusH125[i].SetLineColor(kViolet+1)
    hs.Add(histos_WplusH125[i])

    histos_WminusH125[i].SetFillColor(kViolet+1)
    histos_WminusH125[i].SetLineColor(kViolet+1)
    hs.Add(histos_WminusH125[i])

    histos_ZH125[i].SetFillColor(kViolet+1)
    histos_ZH125[i].SetLineColor(kViolet+1)
    hs.Add(histos_ZH125[i])

    histos_bbH125[i].SetFillColor(kViolet+1)
    histos_bbH125[i].SetLineColor(kViolet+1)
    hs.Add(histos_bbH125[i])

    histos_ttH125[i].SetFillColor(kViolet+1)
    histos_ttH125[i].SetLineColor(kViolet+1)
    hs.Add(histos_ttH125[i])




    hs.Draw('histo')
    
    hs.GetXaxis().SetTitle(namelist_label[i])
    hs.GetXaxis().SetLabelFont(43)
    hs.GetXaxis().SetLabelSize(15)
    hs.GetYaxis().SetTitleSize(20)
    hs.GetYaxis().SetTitleFont(43)
    hs.GetYaxis().SetTitleOffset(1.4)
    hs.GetYaxis().SetLabelFont(43)
    hs.GetYaxis().SetLabelSize(15)
    hs.GetYaxis().SetTitle("Events")

    # legend
    legend = TLegend(0.74,0.68,0.94,0.87)
    legend.AddEntry(histos_ggH125[i],                "SM Higgs", "f")
    legend.AddEntry(histos_ggTo4e_Contin_MCFM701[i], "gg->ZZ",   "f")
    legend.AddEntry(histos_ZZTo4lext1[i],            "qq->ZZ",   "f")
    legend.AddEntry(histos_TTZJets_M10_MLMext1[i],   "TTZ",      "f")
    legend.AddEntry(histos_TTWJetsToLNu[i],          "TTW",      "f")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    legend.Draw()

    canvas.Update()

    canvas.SaveAs(OutputPath + "/" + namelist[i] + ".png")
    
    i = i+1
