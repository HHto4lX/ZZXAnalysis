// ***********************
// run with:
//
// root -l -b -q plot_mbbAndm4l_3years.C++
// ***********************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TColor.h"
#include "TFile.h"
#include "TFrame.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <vector>
#include "TLorentzVector.h"

#include "CMS_lumi.C"

using namespace std;


//*************************************************************************************
// PROCESSES
enum Process {Data=0, HH=1, ggH=2, VBF=3, VH=4, ttH=5, bbH=6, qqZZ=7, ggZZ=8, TTZ=9, TTW=10, VVV=11, ZXbkg=12}; 
const int nProcesses = 13;
TString sProcess[nProcesses] = {"Data", "HH", "ggH", "VBF", "VH", "ttH", "bbH", "qqZZ", "ggZZ", "TTZ", "TTW", "VVV", "ZXbkg"};
//*************************************************************************************



void doPlots(){

  // lumi 3 years
  TString lumiText = "137 fb^{-1}";

  // -- input names
  Int_t nPlots = 2;
  TString sPlots[] = {
    "m4l_4ljjsel",
    "mbb_4ljjsel",
  };


  // --- output dir
  TString outPath_plots = "plots_3years";
  cout<<"creating output dir "<<outPath_plots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_plots)).c_str()); // create output dir



  // --- retrieve histos from file 2016
  //  TString inFileName_2016 = "~/www/HH4lX/20200420_plots/plots_4ljjselAndYields_2016/f_histos_h1_4ljjsel_2016.root";     
  TString inFileName_2016 = "f_histos_h1_4ljjsel_2016.root";
  cout<<"Retrieving Data and MC 2016 histos from file "<<inFileName_2016<<" ..."<<endl;
  TFile* fInhistos_2016 = TFile::Open(inFileName_2016);

  TH1F* h1_2016[nPlots][nProcesses];
  for(int pl=0; pl<nPlots; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      h1_2016[pl][pr] = (TH1F*)fInhistos_2016->Get("h1_"+sPlots[pl]+"_"+sProcess[pr]+"_4l_2016");
      cout<<h1_2016[pl][pr]->GetName()<<endl;
    }
  }


  // --- retrieve histos from file 2017
  //  TString inFileName_2017 = "~/www/HH4lX/20200420_plots/plots_4ljjselAndYields_2017/f_histos_h1_4ljjsel_2017.root";     
  TString inFileName_2017 = "f_histos_h1_4ljjsel_2017.root";
  cout<<"Retrieving Data and MC 2017 histos from file "<<inFileName_2017<<" ..."<<endl;
  TFile* fInhistos_2017 = TFile::Open(inFileName_2017);

  TH1F* h1_2017[nPlots][nProcesses];
  for(int pl=0; pl<nPlots; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      h1_2017[pl][pr] = (TH1F*)fInhistos_2017->Get("h1_"+sPlots[pl]+"_"+sProcess[pr]+"_4l_2017");
      cout<<h1_2017[pl][pr]->GetName()<<endl;
    }
  }


  // --- retrieve histos from file 2018
  //  TString inFileName_2018 = "~/www/HH4lX/20200420_plots/plots_4ljjselAndYields_2018/f_histos_h1_4ljjsel_2018.root";     
  TString inFileName_2018 = "f_histos_h1_4ljjsel_2018.root";
  cout<<"Retrieving Data and MC 2018 histos from file "<<inFileName_2018<<" ..."<<endl;
  TFile* fInhistos_2018 = TFile::Open(inFileName_2018);

  TH1F* h1_2018[nPlots][nProcesses];
  for(int pl=0; pl<nPlots; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      h1_2018[pl][pr] = (TH1F*)fInhistos_2018->Get("h1_"+sPlots[pl]+"_"+sProcess[pr]+"_4l_2018");
      cout<<h1_2018[pl][pr]->GetName()<<endl;
    }
  }



  // --- sum 3 years histos
  TH1F* h1_4ljjsel[nPlots][nProcesses];
  for(int pl=0; pl<nPlots; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      h1_4ljjsel[pl][pr] = (TH1F*)h1_2016[pl][pr]->Clone("h1_4ljjsel_"+sPlots[pl]+"_"+sProcess[pr]+"_4l_3years");
      h1_4ljjsel[pl][pr]->Add(h1_2017[pl][pr]);
      h1_4ljjsel[pl][pr]->Add(h1_2018[pl][pr]);
    }
  }



  // --- define canvas, hstack and pads for plots
  TCanvas* c_4ljjsel      [nPlots];
  THStack* hs_4ljjsel     [nPlots];
  TPad*    pad1_4ljjsel   [nPlots];
  TLegend* leg_4ljjsel    [nPlots];
  TH1F*    hMCtot_4ljjsel [nPlots];
  TPad*    pad2_4ljjsel   [nPlots];
  TH1F*    rp_4ljjsel     [nPlots];
  TH1F*    hUncMC_4ljjsel [nPlots];

  //4ljjsel plots
  for(int pl=0; pl<nPlots; pl++){
    // canvas
    c_4ljjsel[pl] = new TCanvas("c_"+sPlots[pl]+"_4l_3years","c_"+sPlots[pl]+"_4l_3years",800,800);
    // hstack
    hs_4ljjsel[pl] = new THStack("hs_"+sPlots[pl],"");
    // VVV process
    h1_4ljjsel[pl][VVV]->SetFillColor(kGreen-3);
    h1_4ljjsel[pl][VVV]->SetLineColor(kGreen-1);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][VVV]); //add to hs
    // Z+X process
    h1_4ljjsel[pl][ZXbkg]->SetFillColor(kGreen+3);
    h1_4ljjsel[pl][ZXbkg]->SetLineColor(kGreen+4);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][ZXbkg]); //add to hs
    // TTV process: TTW + TTV
    h1_4ljjsel[pl][TTW]->SetFillColor(kBlue+3);
    h1_4ljjsel[pl][TTW]->SetLineColor(kBlue+3);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][TTW]); //add to hs
    h1_4ljjsel[pl][TTZ]->SetFillColor(kBlue+3);
    h1_4ljjsel[pl][TTZ]->SetLineColor(kBlue+3);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][TTZ]); //add to hs
    // ggZZ process
    h1_4ljjsel[pl][ggZZ]->SetFillColor(kAzure-3);
    h1_4ljjsel[pl][ggZZ]->SetLineColor(kBlue+2);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][ggZZ]); //add to hs
    // qqZZ process
    h1_4ljjsel[pl][qqZZ]->SetFillColor(kAzure+6);
    h1_4ljjsel[pl][qqZZ]->SetLineColor(kAzure-6);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][qqZZ]); //add to hs
    // SM Higgs processes: ggH + VBF + VH + ttH + bbH
    h1_4ljjsel[pl][ggH]->SetFillColor(kViolet+6);
    h1_4ljjsel[pl][ggH]->SetLineColor(kViolet+6);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][ggH]); //add to hs
    h1_4ljjsel[pl][VBF]->SetFillColor(kViolet+6);
    h1_4ljjsel[pl][VBF]->SetLineColor(kViolet+6);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][VBF]); //add to hs
    h1_4ljjsel[pl][VH]->SetFillColor(kViolet+6);
    h1_4ljjsel[pl][VH]->SetLineColor(kViolet+6);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][VH]); //add to hs
    h1_4ljjsel[pl][ttH]->SetFillColor(kViolet+6);
    h1_4ljjsel[pl][ttH]->SetLineColor(kViolet+6);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][ttH]); //add to hs
    h1_4ljjsel[pl][bbH]->SetFillColor(kViolet+6);
    h1_4ljjsel[pl][bbH]->SetLineColor(kViolet+6);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][bbH]); //add to hs
    // HH signal
    h1_4ljjsel[pl][HH]->SetLineStyle(7);
    h1_4ljjsel[pl][HH]->SetLineColor(kRed);
    h1_4ljjsel[pl][HH]->SetLineWidth(2);
    h1_4ljjsel[pl][HH]->Scale(100.);
    // data
    h1_4ljjsel[pl][Data]->SetMarkerColor(kBlack);
    h1_4ljjsel[pl][Data]->SetLineColor(kBlack);
    h1_4ljjsel[pl][Data]->SetMarkerStyle(20);

    // --- upper plot pad
    pad1_4ljjsel[pl] = new TPad("pad1_"+sPlots[pl],"pad1_"+sPlots[pl], 0, 0.3, 1, 1.0);
    pad1_4ljjsel[pl]->Draw();
    pad1_4ljjsel[pl]->cd();

    // hs_4ljjsel[pl]->SetMaximum(10e04);
    // hs_4ljjsel[pl]->SetMinimum(10e-04);

    if(sPlots[pl] == "mbb_4ljjsel"){
      hs_4ljjsel[pl]->SetMaximum(10.);
    }
    else if(sPlots[pl] == "m4l_4ljjsel"){
      hs_4ljjsel[pl]->SetMaximum(60.);
    }

      
    hs_4ljjsel[pl]->Draw("histo");
    h1_4ljjsel[pl][HH]->Draw("histosame");
    h1_4ljjsel[pl][Data]->Draw("samepe");

    //    hs_4ljjsel[pl]->GetXaxis()->SetLabelFont(43);
    //    hs_4ljjsel[pl]->GetXaxis()->SetLabelSize(15);
    hs_4ljjsel[pl]->GetXaxis()->SetTitle(h1_4ljjsel[pl][HH]->GetXaxis()->GetTitle());
    //    hs_4ljjsel[pl]->GetYaxis()->SetTitleSize(20);
    //    hs_4ljjsel[pl]->GetYaxis()->SetTitleFont(43);
    //    hs_4ljjsel[pl]->GetYaxis()->SetTitleOffset(1.4);
    //    hs_4ljjsel[pl]->GetYaxis()->SetLabelFont(43);
    //    hs_4ljjsel[pl]->GetYaxis()->SetLabelSize(15);
    hs_4ljjsel[pl]->GetYaxis()->SetTitle(h1_4ljjsel[pl][HH]->GetYaxis()->GetTitle());

    // --- legend
    leg_4ljjsel[pl] = new TLegend(0.47,0.57,0.67,0.89);
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][ggZZ], "gg #rightarrow ZZ #rightarrow 4l",       "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][qqZZ], "q#bar{q} #rightarrow ZZ #rightarrow 4l", "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][ggH],  "SM Higgs",                               "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][TTZ],  "ttV where V=Z,W",                        "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][VVV],  "VVV where V=Z,W",                        "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][ZXbkg],"Z+X",                                    "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][Data], "Data",                                   "lp");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][HH],   "HH #rightarrow b#bar{b}4l x100",         "l");
    leg_4ljjsel[pl]->SetTextSize(0.030);
    leg_4ljjsel[pl]->SetLineColor(0);
    leg_4ljjsel[pl]->SetLineWidth(1);
    leg_4ljjsel[pl]->SetFillColor(kWhite);
    leg_4ljjsel[pl]->SetBorderSize(0);
    leg_4ljjsel[pl]->Draw();

    c_4ljjsel[pl]->Update();

    //    pad1_4ljjsel[pl]->SetLogy();

    c_4ljjsel[pl]->Update();

    // --- tot hist for all MC
    hMCtot_4ljjsel[pl] = (TH1F*)h1_4ljjsel[pl][ggH]->Clone("hMCtot_4ljjsel_"+sPlots[pl]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][VBF]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][VH]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][ttH]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][bbH]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][qqZZ]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][ggZZ]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][TTZ]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][TTW]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][ZXbkg]);
    hMCtot_4ljjsel[pl]->Add(h1_4ljjsel[pl][VVV]);

    // --- lower pad plot
    c_4ljjsel[pl]->cd();
    pad2_4ljjsel[pl] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
    pad2_4ljjsel[pl]->SetGridy();
    pad2_4ljjsel[pl]->Draw();
    pad2_4ljjsel[pl]->cd();

    // --- define ratio plot
    rp_4ljjsel[pl] = (TH1F*)h1_4ljjsel[pl][Data]->Clone("rp_4ljjsel_"+sPlots[pl]);
    rp_4ljjsel[pl]->SetLineColor(kBlack);
    rp_4ljjsel[pl]->SetMinimum(0.);
    rp_4ljjsel[pl]->SetMaximum(2.);
    rp_4ljjsel[pl]->SetStats(0);
    rp_4ljjsel[pl]->Divide(hMCtot_4ljjsel[pl]); //divide histo rp/MC
    rp_4ljjsel[pl]->SetMarkerStyle(20);
    rp_4ljjsel[pl]->SetMarkerColor(kBlack);
    rp_4ljjsel[pl]->SetTitle("");

    rp_4ljjsel[pl]->SetYTitle("Data/#Sigma Bkg.");
    rp_4ljjsel[pl]->GetYaxis()->SetNdivisions(505);
    rp_4ljjsel[pl]->GetYaxis()->SetTitleSize(20);
    rp_4ljjsel[pl]->GetYaxis()->SetTitleFont(43);
    rp_4ljjsel[pl]->GetYaxis()->SetTitleOffset(1.4);
    rp_4ljjsel[pl]->GetYaxis()->SetLabelFont(43);
    rp_4ljjsel[pl]->GetYaxis()->SetLabelSize(15);

    rp_4ljjsel[pl]->GetXaxis()->SetTitleSize(20);
    rp_4ljjsel[pl]->GetXaxis()->SetTitleFont(43);
    rp_4ljjsel[pl]->GetXaxis()->SetTitleOffset(4.);
    rp_4ljjsel[pl]->GetXaxis()->SetLabelFont(43);
    rp_4ljjsel[pl]->GetXaxis()->SetLabelSize(15);

    // --- define mc shadow unc plot
    hUncMC_4ljjsel[pl] = (TH1F*)hMCtot_4ljjsel[pl]->Clone("hUncMC_4ljjsel_"+sPlots[pl]);
    for(int xbin=1; xbin < hUncMC_4ljjsel[pl]->GetXaxis()->GetNbins() + 1; xbin++){
      float err = 0.;
      if(hUncMC_4ljjsel[pl]->GetBinContent(xbin) == 0) continue;
      err = hUncMC_4ljjsel[pl]->GetBinError(xbin) / hUncMC_4ljjsel[pl]->GetBinContent(xbin);
      hUncMC_4ljjsel[pl]->SetBinContent(xbin, 1.);
      hUncMC_4ljjsel[pl]->SetBinError(xbin, err);
    }
    hUncMC_4ljjsel[pl]->SetLineColor(1);
    hUncMC_4ljjsel[pl]->SetFillStyle(3005);
    hUncMC_4ljjsel[pl]->SetFillColor(kGray+3);
    hUncMC_4ljjsel[pl]->SetMarkerColor(1);
    hUncMC_4ljjsel[pl]->SetMarkerStyle(1);
    hUncMC_4ljjsel[pl]->SetTitle("");
    hUncMC_4ljjsel[pl]->SetStats(0);

    // ---draw
    rp_4ljjsel[pl]->Draw("ep");
    hUncMC_4ljjsel[pl]->Draw("e2 same");

    c_4ljjsel[pl]->Update();

    // --- draw CMS and lumi text
    writeExtraText = true;
    extraText      = "Preliminary";
    lumi_sqrtS     = lumiText + " (13 TeV)";
    cmsTextSize    = 0.6;
    lumiTextSize   = 0.46;
    extraOverCmsTextSize = 0.75;
    relPosX = 0.12;
    relPosY = 0.065;
    relExtraDY = 1.1;
    //    outOfFrame = false;
    CMS_lumi(pad1_4ljjsel[pl], 0, 0);


    c_4ljjsel[pl]->SaveAs(outPath_plots + "/" + c_4ljjsel[pl]->GetName() + ".png");
    c_4ljjsel[pl]->SaveAs(outPath_plots + "/" + c_4ljjsel[pl]->GetName() + ".pdf");
    

  }


}






void plot_mbbAndm4l_3years(){
  
  doPlots();
}
