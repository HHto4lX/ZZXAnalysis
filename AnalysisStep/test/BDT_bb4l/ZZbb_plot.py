#!/usr/bin/en python
   
import ROOT, math
from ROOT import TFile, TH1F, TCanvas, gSystem, TAttFill, TLegend, TStyle, THStack, kViolet, kBlack, kAzure, kCyan, kGreen, kWhite, kOrange, kBlue, kRed, TPad, kCyan, TGraphErrors, kGray
import CMSGraphics, CMS_lumi

# create output directory
OutputPath = 'inputBDT_SR_plots'
gSystem.Exec("mkdir -p " + OutputPath)
print "Output directory created!"


# inputh file phath
inputFilePath = 'histos_plotsinputBDT_SR/'


# hist names
#namelist = ['h_4leptonsPt','h_MET','h_M4l']
namelist = ['h_bTagger_jet1','h_bTagger_jet2','h_pT_jet1','h_pT_jet2','h_MET','h_DeltaR','h_mjj']


# read files
inFile_AllData = TFile.Open(inputFilePath + 'histos_AllData.root')
histos_AllData = []
for iname_AllData in namelist :
    histos_AllData.append(inFile_AllData.Get(iname_AllData))


inFile_HH4lbb = TFile.Open(inputFilePath + 'histos_HH4lbb_Angela.root')
histos_HH4lbb = []
for iname_HH4lbb in namelist :
    histos_HH4lbb.append(inFile_HH4lbb.Get(iname_HH4lbb))


inFile_ggH125 = TFile.Open(inputFilePath + 'histos_ggH125.root')
histos_ggH125 = []
for iname_ggH125 in namelist :
    histos_ggH125.append(inFile_ggH125.Get(iname_ggH125))


inFile_VBFH125 = TFile.Open(inputFilePath + 'histos_VBFH125.root')
histos_VBFH125 = []
for iname_VBFH125 in namelist :
    histos_VBFH125.append(inFile_VBFH125.Get(iname_VBFH125))


inFile_WplusH125 = TFile.Open(inputFilePath + 'histos_WplusH125.root')
histos_WplusH125 = []
for iname_WplusH125 in namelist :
    histos_WplusH125.append(inFile_WplusH125.Get(iname_WplusH125))


inFile_WminusH125 = TFile.Open(inputFilePath + 'histos_WminusH125.root')
histos_WminusH125 = []
for iname_WminusH125 in namelist :
    histos_WminusH125.append(inFile_WminusH125.Get(iname_WminusH125))


inFile_ZH125 = TFile.Open(inputFilePath + 'histos_ZH125.root')
histos_ZH125 = []
for iname_ZH125 in namelist :
    histos_ZH125.append(inFile_ZH125.Get(iname_ZH125))


inFile_bbH125 = TFile.Open(inputFilePath + 'histos_bbH125.root')
histos_bbH125 = []
for iname_bbH125 in namelist :
    histos_bbH125.append(inFile_bbH125.Get(iname_bbH125))


inFile_ttH125 = TFile.Open(inputFilePath + 'histos_ttH125.root')
histos_ttH125 = []
for iname_ttH125 in namelist :
    histos_ttH125.append(inFile_ttH125.Get(iname_ttH125))


inFile_ggTo4e_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4e_Contin_MCFM701.root')
histos_ggTo4e_Contin_MCFM701 = []
for iname_ggTo4e_Contin_MCFM701 in namelist :
    histos_ggTo4e_Contin_MCFM701.append(inFile_ggTo4e_Contin_MCFM701.Get(iname_ggTo4e_Contin_MCFM701))


inFile_ggTo4mu_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4mu_Contin_MCFM701.root')
histos_ggTo4mu_Contin_MCFM701 = []
for iname_ggTo4mu_Contin_MCFM701 in namelist :
    histos_ggTo4mu_Contin_MCFM701.append(inFile_ggTo4mu_Contin_MCFM701.Get(iname_ggTo4mu_Contin_MCFM701))


inFile_ggTo4tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo4tau_Contin_MCFM701.root')
histos_ggTo4tau_Contin_MCFM701 = []
for iname_ggTo4tau_Contin_MCFM701 in namelist :
    histos_ggTo4tau_Contin_MCFM701.append(inFile_ggTo4tau_Contin_MCFM701.Get(iname_ggTo4tau_Contin_MCFM701))


inFile_ggTo2e2mu_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2e2mu_Contin_MCFM701.root')
histos_ggTo2e2mu_Contin_MCFM701 = []
for iname_ggTo2e2mu_Contin_MCFM701 in namelist :
    histos_ggTo2e2mu_Contin_MCFM701.append(inFile_ggTo2e2mu_Contin_MCFM701.Get(iname_ggTo2e2mu_Contin_MCFM701))


inFile_ggTo2e2tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2e2tau_Contin_MCFM701.root')
histos_ggTo2e2tau_Contin_MCFM701 = []
for iname_ggTo2e2tau_Contin_MCFM701 in namelist :
    histos_ggTo2e2tau_Contin_MCFM701.append(inFile_ggTo2e2tau_Contin_MCFM701.Get(iname_ggTo2e2tau_Contin_MCFM701))


inFile_ggTo2mu2tau_Contin_MCFM701 = TFile.Open(inputFilePath + 'histos_ggTo2mu2tau_Contin_MCFM701.root')
histos_ggTo2mu2tau_Contin_MCFM701 = []
for iname_ggTo2mu2tau_Contin_MCFM701 in namelist :
    histos_ggTo2mu2tau_Contin_MCFM701.append(inFile_ggTo2mu2tau_Contin_MCFM701.Get(iname_ggTo2mu2tau_Contin_MCFM701))


inFile_ZZTo4lext1 = TFile.Open(inputFilePath + 'histos_ZZTo4lamcatnlo.root')
histos_ZZTo4lext1 = []
for iname_ZZTo4lext1 in namelist :
    histos_ZZTo4lext1.append(inFile_ZZTo4lext1.Get(iname_ZZTo4lext1))


inFile_TTZJets_M10_MLMext1 = TFile.Open(inputFilePath + 'histos_TTZJets_M10_MLMext1.root')
histos_TTZJets_M10_MLMext1 = []
for iname_TTZJets_M10_MLMext1 in namelist :
    histos_TTZJets_M10_MLMext1.append(inFile_TTZJets_M10_MLMext1.Get(iname_TTZJets_M10_MLMext1))


inFile_TTZToLL_M1to10_MLM = TFile.Open(inputFilePath + 'histos_TTZToLL_M1to1O_MLM.root')
histos_TTZToLL_M1to10_MLM = []
for iname_TTZToLL_M1to10_MLM in namelist :
    histos_TTZToLL_M1to10_MLM.append(inFile_TTZToLL_M1to10_MLM.Get(iname_TTZToLL_M1to10_MLM))


inFile_TTWJetsToLNu = TFile.Open(inputFilePath + 'histos_TTWJetsToLNu.root')
histos_TTWJetsToLNu = []
for iname_TTWJetsToLNu in namelist :
    histos_TTWJetsToLNu.append(inFile_TTWJetsToLNu.Get(iname_TTWJetsToLNu))


inFile_WWZ = TFile.Open(inputFilePath + 'histos_WWZ.root')
histos_WWZ = []
for iname_WWZ in namelist :
    histos_WWZ.append(inFile_WWZ.Get(iname_WWZ))


inFile_WZZ = TFile.Open(inputFilePath + 'histos_WZZ.root')
histos_WZZ = []
for iname_WZZ in namelist :
    histos_WZZ.append(inFile_WZZ.Get(iname_WZZ))


inFile_ZZZ = TFile.Open(inputFilePath + 'histos_ZZZ.root')
histos_ZZZ = []
for iname_ZZZ in namelist :
    histos_ZZZ.append(inFile_ZZZ.Get(iname_ZZZ))


# inFile_ZX = TFile.Open(inputFilePath + 'histos_Z+X.root')
# histos_ZX = []
# for iname_ZX in namelist :
#     histos_ZX.append(inFile_ZX.Get(iname_ZX))


print 'files read'


# --- do plots ---
i = 0 #counter for histos list
for name in namelist :

    canvas = TCanvas('canvas','canvas',800,800)

    hs = THStack('hs','')


    #  tribosons
    histos_tribosons = histos_WWZ[i].Clone('histos_tribosons')
    histos_tribosons.Add(histos_WZZ[i])
    histos_tribosons.Add(histos_ZZZ[i])
    histos_tribosons.SetFillColor(kGreen-3)
    histos_tribosons.SetLineColor(kGreen-1)
    hs.Add(histos_tribosons)

    # # ZX background
    # histos_ZX[i].SetFillColor(kGreen+3)
    # histos_ZX[i].SetLineColor(kGreen+4)
    # hs.Add(histos_ZX[i])

    # TTZ
    histos_TTV = histos_TTZJets_M10_MLMext1[i].Clone('histos_TTV')
    histos_TTV.Add(histos_TTZToLL_M1to10_MLM[i])
    histos_TTV.Add(histos_TTWJetsToLNu[i])
    histos_TTV.SetFillColor(kBlue+3)
    histos_TTV.SetLineColor(kBlack)
    hs.Add(histos_TTV)

    
    # gg->ZZ
    histos_ggZZ = histos_ggTo4e_Contin_MCFM701[i].Clone('histos_ggZZ')
    histos_ggZZ.Add(histos_ggTo4mu_Contin_MCFM701[i])
    histos_ggZZ.Add(histos_ggTo4tau_Contin_MCFM701[i])
    histos_ggZZ.Add(histos_ggTo2e2mu_Contin_MCFM701[i])
    histos_ggZZ.Add(histos_ggTo2e2tau_Contin_MCFM701[i])
    histos_ggZZ.Add(histos_ggTo2mu2tau_Contin_MCFM701[i])
    histos_ggZZ.SetFillColor(kAzure-3)
    histos_ggZZ.SetLineColor(kBlue+2)
#    histo_ggZZ.Rebin(6)
#    integral_fondi += histo_ggZZ.Integral()
    hs.Add(histos_ggZZ)


    # qq->ZZ
    histos_ZZTo4lext1[i].SetFillColor(kAzure+6)
    histos_ZZTo4lext1[i].SetLineColor(kAzure-6)
#   histos_ZZTo4lext1[i].Rebin(6)
#   integral_fondi += histos_ZZTo4lext1[i].Integral()
    hs.Add(histos_ZZTo4lext1[i])


    # SM Higgs
    histos_SMHiggs = histos_ggH125[i].Clone('histos_SMHiggs')
    histos_SMHiggs.Add(histos_VBFH125[i])
    histos_SMHiggs.Add(histos_WplusH125[i])
    histos_SMHiggs.Add(histos_WminusH125[i])
    histos_SMHiggs.Add(histos_ZH125[i])
    histos_SMHiggs.Add(histos_bbH125[i])
    histos_SMHiggs.Add(histos_ttH125[i])
    histos_SMHiggs.SetFillColor(kViolet+6)
    histos_SMHiggs.SetLineColor(kViolet+7)
#   histo_SMHiggs.Rebin(6)
#    integral_fondi += histo_SMHiggs.Integral()
    hs.Add(histos_SMHiggs)


    # HH->4lbb signal
    histos_HH4lbb[i].SetLineColor(kRed)
    histos_HH4lbb[i].SetLineWidth(2)
 #   histos_HH4lbb[i].Rebin(6)
    integral = histos_HH4lbb[i].Integral();
    print 'integrale ' 
    print integral
    histos_HH4lbb[i].Scale(100.) 

    
    # ALL DATA
    histos_AllData[i].SetMarkerColor(kBlack)
    histos_AllData[i].SetLineColor(kBlack)
    histos_AllData[i].SetMarkerStyle(20)
 #   histos_AllData[i].Rebin(6)
   


    # --- upper plot pad
    pad1 = TPad("pad1","pad1", 0, 0.3, 1, 1.0)
    pad1.Draw()
    pad1.cd()


    # Draw all
    if 'M4l' in name:
#        hs.SetMaximum(1.5 * max(hs.GetMaximum(),histos_AllData[i].GetMaximum()));
        hs.SetMaximum(10e5)
        hs.SetMinimum(10e-1)

 
    #if logy (leppt e MET)
    # if '4leptonsPt' in name or 'MET' in name : 
    #     hs.SetMaximum(10e5)
    #     hs.SetMinimum(10e-1)

    else :
        hs.SetMaximum(10e4)
        hs.SetMinimum(10e-4)

    hs.Draw('histo')
    histos_HH4lbb[i].Draw('histosame')
    histos_AllData[i].Draw('samepe')  

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
#    hs.SetMinimum(0.)


    # legend
    legend = TLegend(0.78,0.61,0.94,0.87)
    legend.AddEntry(histos_HH4lbb[i],                "HH->4lbb x100", "f")
    legend.AddEntry(histos_SMHiggs,                  "SM Higgs",      "f")
    legend.AddEntry(histos_ZZTo4lext1[i],            "qq->ZZ",        "f")
    legend.AddEntry(histos_ggZZ,                     "gg->ZZ",        "f")
    legend.AddEntry(histos_TTV,                      "TTV; V=Z,W",    "f")
#    legend.AddEntry(histos_ZX[i],                    "Z+X",           "f")
    legend.AddEntry(histos_tribosons,                "VVV; V=Z,W",    "f")
    legend.AddEntry(histos_AllData[i],               "Data",          "lp")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.SetTextFont(43)
    #legend.SetTextSize(20)
    legend.Draw()

    canvas.Update()


    if '4leptonsPt' in name or 'MET' in name : 
        pad1.SetLogy()


    if 'M4l' in name:
        pad1.SetLogy()
    #     pad1.SetLogx()
    #     hs.GetXaxis().SetMoreLogLabels()
    #     hs.GetXaxis().SetNoExponent()

    else:
        pad1.SetLogy()

    canvas.Update()



    # -----------------------------------------
    # tot hist for all mc
    mc_tot = histos_SMHiggs.Clone('mc_tot') 
    mc_tot.Add(histos_ZZTo4lext1[i])
    mc_tot.Add(histos_TTV)
    mc_tot.Add(histos_ggZZ)
#    mc_tot.Add(histos_ZX[i])
    mc_tot.Add(histos_tribosons)


    # --- lower plot pad: RATIO PLOT
    canvas.cd()
    pad2 = TPad("pad2","pad2", 0, 0.05, 1, 0.3)
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()    #pad2 becomes the current pad
      


    #define ratio plot
    rp = TH1F(histos_AllData[i].Clone("rp"))
    rp.SetLineColor(kBlack)
    rp.SetMinimum(0.)
    rp.SetMaximum(2.)
    rp.SetStats(0)
    rp.Divide(mc_tot) #divide histo rp/MC
    rp.SetMarkerStyle(20)
    rp.SetMarkerColor(kBlack)
    rp.SetTitle("") 

    rp.SetYTitle("Data/MC")
    rp.GetYaxis().SetNdivisions(505)
    rp.GetYaxis().SetTitleSize(20)
    rp.GetYaxis().SetTitleFont(43)
    rp.GetYaxis().SetTitleOffset(1.4)
    rp.GetYaxis().SetLabelFont(43)
    rp.GetYaxis().SetLabelSize(15)

    rp.GetXaxis().SetTitleSize(20)
    rp.GetXaxis().SetTitleFont(43)
    rp.GetXaxis().SetTitleOffset(4.)
    rp.GetXaxis().SetLabelFont(43)
    rp.GetXaxis().SetLabelSize(15)

    

    #define mc shadow unc plot
    mc_unc_h = TH1F(mc_tot.Clone('mc_unc_h'))
    gpoint = 0
    for xbin in range(1 , mc_unc_h.GetXaxis().GetNbins() + 1) :   #estremo e' sempre escluso, quindi +1
        err = 0.
        if(mc_unc_h.GetBinContent(xbin) == 0 ) : 
            continue
        err = mc_unc_h.GetBinError(xbin) / mc_unc_h.GetBinContent(xbin)
        mc_unc_h.SetBinContent(xbin, 1)
        mc_unc_h.SetBinError(xbin, err)

    mc_unc_h.SetLineColor(1)
    mc_unc_h.SetFillStyle(3005)
    mc_unc_h.SetFillColor(kGray+3)
    mc_unc_h.SetMarkerColor(1)
    mc_unc_h.SetMarkerStyle(1)
    mc_unc_h.SetTitle('')
    mc_unc_h.SetStats(0)


    #draw
    rp.Draw("ep")
    mc_unc_h.Draw('e2 same')

    canvas.Update()



    if 'M4L' in name:
        pad2.SetLogx()
        rp.GetXaxis().SetMoreLogLabels();
        rp.GetXaxis().SetNoExponent();

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
    CMS_lumi.CMS_lumi(pad1, 0, 0)

    canvas.Update()

    canvas.SaveAs(OutputPath + "/" + namelist[i] + ".png")
    canvas.SaveAs(OutputPath + "/" + namelist[i] + ".pdf")
    canvas.SaveAs(OutputPath + "/" + namelist[i] + ".root")
    
    i = i+1

    del canvas
