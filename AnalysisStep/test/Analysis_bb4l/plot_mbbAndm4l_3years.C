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
// FINAL STATES
enum FinalState {fs_4mu=0, fs_4e=1, fs_2e2mu=2};  // 4mu, 4e, 2e2mu
const int nFinalStates = 3;
TString sFinalState[nFinalStates+1] = {"4mu", "4e","2e2mu","4l"};
//*************************************************************************************
// PROCESSES
enum Process {Data=0, HH=1, ggH=2, VBF=3, VH=4, ttH=5, bbH=6, qqZZ=7, ggZZ=8, TTZ=9, TTW=10, VVV=11, ZXbkg=12}; 
const int nProcesses = 13;
TString sProcess[nProcesses] = {"Data", "HH", "ggH", "VBF", "VH", "ttH", "bbH", "qqZZ", "ggZZ", "TTZ", "TTW", "VVV", "ZXbkg"};
//*************************************************************************************

// lumi 3 years
TString lumiText = "138 fb^{-1}";

// update lumi for new Run2 recommendations (comb)
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
float lumi2016_old = 35.8; //fb^-1
float lumi2016_new = 36.3;
float lumi2017_old = 41.5;
float lumi2017_new = 41.5;
float lumi2018_old = 59.7;
float lumi2018_new = 59.8;

float lumiScaleComb = 1.01135794221; // ratio between new and old 2016 lumi - comb way..


// yields ZX angela
float yield_ZX_2016_4mu   = 0.790;
float yield_ZX_2016_4e    = 1.400;
float yield_ZX_2016_2e2mu = 2.640;
float yield_ZX_2017_4mu   = 1.480;
float yield_ZX_2017_4e    = 0.520;
float yield_ZX_2017_2e2mu = 2.000;
float yield_ZX_2018_4mu   = 1.600;
float yield_ZX_2018_4e    = 0.720;
float yield_ZX_2018_2e2mu = 2.220;

float yield_ZX_2016[] = {0.790, 1.400, 2.640, 0.790+1.400+2.640};
float yield_ZX_2017[] = {1.480, 0.520, 2.000, 1.480+0.520+2.000};
float yield_ZX_2018[] = {1.600, 0.720, 2.220, 1.600+0.720+2.220};



void doPlots(){

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
      if(pr != 0){ // scale histos to new lumi only for MC
        cout<<h1_2016[pl][pr]->GetName()<<endl;
	cout<<h1_2016[pl][pr]->Integral()<<" "<<h1_2016[pl][pr]->Integral() * (lumi2016_new / lumi2016_old)<<endl;
        cout<<h1_2016[pl][pr]->Integral()<<" "<<h1_2016[pl][pr]->Integral() * (lumiScaleComb)<<endl;
        //h1_2016[pl][pr]->Scale(lumi2016_new / lumi2016_old); // scale to new lumi
        h1_2016[pl][pr]->Scale(lumiScaleComb); // scale to new lumi Comb way..
        cout<<h1_2016[pl][pr]->Integral()<<endl;
      }
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
      // if(pr != 0){ // scale histos to new lumi only for MC
      //   cout<<h1_2018[pl][pr]->GetName()<<endl;
      //   cout<<h1_2018[pl][pr]->Integral()<<" "<<h1_2018[pl][pr]->Integral() * (lumi2018_new / lumi2018_old)<<endl;
      //   h1_2018[pl][pr]->Scale(lumi2018_new / lumi2018_old); // scale to new lumi
      //   cout<<h1_2018[pl][pr]->Integral()<<endl;
      // }
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
  TLegend* leg2_4ljjsel   [nPlots];

  //4ljjsel plots
  for(int pl=0; pl<nPlots; pl++){
    // canvas
    c_4ljjsel[pl] = new TCanvas("c_"+sPlots[pl]+"_4l_3years","c_"+sPlots[pl]+"_4l_3years",800,800);
    // hstack
    hs_4ljjsel[pl] = new THStack("hs_"+sPlots[pl],"");
    // VVV process
    h1_4ljjsel[pl][VVV]->SetFillColor(kOrange+1);
    h1_4ljjsel[pl][VVV]->SetLineColor(kOrange+1);
    hs_4ljjsel[pl]->Add(h1_4ljjsel[pl][VVV]); //add to hs
    // Z+X process
    h1_4ljjsel[pl][ZXbkg]->SetFillColor(kGreen-2);
    h1_4ljjsel[pl][ZXbkg]->SetLineColor(kGreen-2);
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
    //h1_4ljjsel[pl][Data]->SetBinErrorOption(TH1::kPoisson); // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
    h1_4ljjsel[pl][Data]->SetMarkerColor(kBlack);
    h1_4ljjsel[pl][Data]->SetLineColor(kBlack);
    h1_4ljjsel[pl][Data]->SetMarkerStyle(8);
    //    h1_4ljjsel[pl][Data]->SetMarkerSize(1.2);

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
    h1_4ljjsel[pl][Data]->Draw("samep e x0"); // option x0 removes orizzontal error bars

    hs_4ljjsel[pl]->GetXaxis()->SetTitleSize(27);//25
    hs_4ljjsel[pl]->GetXaxis()->SetTitleFont(43);
    hs_4ljjsel[pl]->GetXaxis()->SetTitleOffset(1.3);
    hs_4ljjsel[pl]->GetXaxis()->SetLabelFont(43);
    hs_4ljjsel[pl]->GetXaxis()->SetLabelSize(25);//23
    hs_4ljjsel[pl]->GetXaxis()->SetTitle(h1_4ljjsel[pl][HH]->GetXaxis()->GetTitle());
    if(sPlots[pl] == "m4l_4ljjsel"){
      hs_4ljjsel[pl]->GetXaxis()->SetTitle("m_{4#font[12]{l}} (GeV)");
    }
    if(sPlots[pl] == "mbb_4ljjsel"){
      hs_4ljjsel[pl]->GetXaxis()->SetTitle("m_{b#bar{b}} (GeV)");
    }
    hs_4ljjsel[pl]->GetYaxis()->SetTitleSize(27);//25
    hs_4ljjsel[pl]->GetYaxis()->SetTitleFont(43);
    hs_4ljjsel[pl]->GetYaxis()->SetTitleOffset(1.1);
    hs_4ljjsel[pl]->GetYaxis()->SetLabelFont(43);
    hs_4ljjsel[pl]->GetYaxis()->SetLabelSize(25);//23
    hs_4ljjsel[pl]->GetYaxis()->SetTitle(h1_4ljjsel[pl][HH]->GetYaxis()->GetTitle());


    // --- legend
    //leg_4ljjsel[pl] = new TLegend(0.52,0.55,0.75,0.89);
    leg_4ljjsel[pl] = new TLegend(0.60,0.48,0.82,0.89);
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][ggZZ], "gg #rightarrow ZZ #rightarrow 4#font[12]{l}",       "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][qqZZ], "q#bar{q} #rightarrow ZZ #rightarrow 4#font[12]{l}", "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][ggH],  "SM H",                                              "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][TTZ],  "t#bar{t}V, where V = Z, W",                         "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][VVV],  "VVV, where V = Z, W",                               "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][ZXbkg],"Z+X",                                               "f");
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][Data], "Data",                                              "ep"); //vertical error bar in legend
    leg_4ljjsel[pl]->AddEntry(h1_4ljjsel[pl][HH],   "HH #rightarrow b#bar{b}4#font[12]{l} x100",         "l");
    leg_4ljjsel[pl]->SetTextSize(0.037); //35
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
    pad2_4ljjsel[pl] = new TPad("pad2","pad2", 0, 0.02, 1, 0.3);//0.05
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
    rp_4ljjsel[pl]->SetMarkerStyle(8);
    rp_4ljjsel[pl]->SetMarkerColor(kBlack);
    //    rp_4ljjsel[pl]->SetMarkerSize(1.2);
    rp_4ljjsel[pl]->SetTitle("");

    rp_4ljjsel[pl]->SetYTitle("Data/#Sigma Bkg.");
    rp_4ljjsel[pl]->GetYaxis()->SetNdivisions(505);
    rp_4ljjsel[pl]->GetYaxis()->SetTitleSize(27);//25
    rp_4ljjsel[pl]->GetYaxis()->SetTitleFont(43);
    rp_4ljjsel[pl]->GetYaxis()->SetTitleOffset(1.1);
    rp_4ljjsel[pl]->GetYaxis()->SetLabelFont(43);
    rp_4ljjsel[pl]->GetYaxis()->SetLabelSize(25);//23

    rp_4ljjsel[pl]->GetXaxis()->SetTitleSize(27);//25
    rp_4ljjsel[pl]->GetXaxis()->SetTitleFont(43);
    rp_4ljjsel[pl]->GetXaxis()->SetTitleOffset(4.);
    rp_4ljjsel[pl]->GetXaxis()->SetLabelFont(43);
    rp_4ljjsel[pl]->GetXaxis()->SetLabelSize(25);//23

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
    //hUncMC_4ljjsel[pl]->SetFillStyle(3005);
    hUncMC_4ljjsel[pl]->SetFillStyle(3244);
    hUncMC_4ljjsel[pl]->SetFillColor(kGray+3);
    hUncMC_4ljjsel[pl]->SetMarkerColor(1);
    hUncMC_4ljjsel[pl]->SetMarkerStyle(1);
    hUncMC_4ljjsel[pl]->SetTitle("");
    hUncMC_4ljjsel[pl]->SetStats(0);

    // ---draw
    rp_4ljjsel[pl]->Draw("pe x0");
    hUncMC_4ljjsel[pl]->Draw("e2 same");

    // --- legend ratio plot
    leg2_4ljjsel[pl] = new TLegend(0.67,0.75,0.87,0.89);
    leg2_4ljjsel[pl]->AddEntry(hUncMC_4ljjsel[pl], "MC stat. unc.", "f");
    leg2_4ljjsel[pl]->SetTextSize(0.09);//08
    leg2_4ljjsel[pl]->SetLineColor(0);
    leg2_4ljjsel[pl]->SetLineWidth(1);
    leg2_4ljjsel[pl]->SetFillColor(kWhite);
    leg2_4ljjsel[pl]->SetBorderSize(0);

    leg2_4ljjsel[pl]->Draw("same");


    c_4ljjsel[pl]->Update();

    // --- draw CMS and lumi text
    writeExtraText = false;
    extraText      = "Preliminary";
    lumi_sqrtS     = lumiText + " (13 TeV)";
    cmsTextSize    = 0.6;
    lumiTextSize   = 0.46;
    extraOverCmsTextSize = 0.75;
    relPosX = 0.05;//0.12;
    relPosY = 0.12;//0.065;
    relExtraDY = 1.1;
    //    outOfFrame = false;
    CMS_lumi(pad1_4ljjsel[pl], 0, 0);


    c_4ljjsel[pl]->SaveAs(outPath_plots + "/" + c_4ljjsel[pl]->GetName() + ".png");
    c_4ljjsel[pl]->SaveAs(outPath_plots + "/" + c_4ljjsel[pl]->GetName() + ".pdf");
    

  }


}


void scaleYields(){

  // --- output dir
  TString outPath_yields = "yields_3years";
  cout<<"creating output dir "<<outPath_yields<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_yields)).c_str()); // create output dir


  // --- retrieve yields from file 2016
  TString inFileName_2016 = "f_yields_2016.root";
  cout<<"Retrieving 2016 yields from file "<<inFileName_2016<<" ..."<<endl;
  TFile* fInYields_2016 = TFile::Open(inFileName_2016);

  TH1F* hTemp2016;
  Float_t yield_4ljjsel_2016            [nProcesses][nFinalStates+1];
  Float_t yield_4ljjsel_2016_lumiScaled [nProcesses][nFinalStates+1];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      hTemp2016 = (TH1F*)fInYields_2016->Get("hYields_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_2016");
      yield_4ljjsel_2016[pr][fs]            = hTemp2016->GetBinContent(1);
      yield_4ljjsel_2016_lumiScaled[pr][fs] = hTemp2016->GetBinContent(1) * lumiScaleComb;
    }
  }

  // --- retrieve yields from file 2017
  TString inFileName_2017 = "f_yields_2017.root";
  cout<<"Retrieving 2017 yields from file "<<inFileName_2017<<" ..."<<endl;
  TFile* fInYields_2017 = TFile::Open(inFileName_2017);

  TH1F* hTemp2017;
  Float_t yield_4ljjsel_2017[nProcesses][nFinalStates+1];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      hTemp2017 = (TH1F*)fInYields_2017->Get("hYields_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_2017");
      yield_4ljjsel_2017[pr][fs] = hTemp2017->GetBinContent(1);
    }
  }

  // --- retrieve yields from file 2018
  TString inFileName_2018 = "f_yields_2018.root";
  cout<<"Retrieving 2018 yields from file "<<inFileName_2018<<" ..."<<endl;
  TFile* fInYields_2018 = TFile::Open(inFileName_2018);

  TH1F* hTemp2018;
  Float_t yield_4ljjsel_2018[nProcesses][nFinalStates+1];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      hTemp2018 = (TH1F*)fInYields_2018->Get("hYields_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_2018");
      yield_4ljjsel_2018[pr][fs] = hTemp2018->GetBinContent(1);
    }
  }




  // --- print yields after 4ljj sel for full Run2 
  cout<<"print yields after 4ljj sel for full Run2 ... "<<endl;
  ofstream f_yields4ljjsel;
  TString f_yields4ljjsel_name = "yields4ljjsel_FullRun2.txt";
  f_yields4ljjsel.open(string(outPath_yields) + "/" + string(f_yields4ljjsel_name));
  f_yields4ljjsel<<"--- 2016 ---"<<endl;
  f_yields4ljjsel<<"|Final state |signal HH |ttZ |ttH |bbH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |ZH+WH |others(=VVV+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel_2016[HH][fs]
                                         <<" |"<<yield_4ljjsel_2016[TTZ][fs]  
                                         <<" |"<<yield_4ljjsel_2016[ttH][fs]
                                         <<" |"<<yield_4ljjsel_2016[bbH][fs]
                                         <<" |"<<yield_4ljjsel_2016[qqZZ][fs]+yield_4ljjsel_2016[ggZZ][fs]
                                         <<" |"<<yield_4ljjsel_2016[ggH][fs]+yield_4ljjsel_2016[VBF][fs]
                                         <<" |"<<yield_4ljjsel_2016[VH][fs]
                                         <<" |"<<yield_4ljjsel_2016[VVV][fs]+yield_4ljjsel_2016[TTW][fs]
                                         <<" |"<<yield_4ljjsel_2016[ZXbkg][fs]
                                         <<" |"<<endl;
  }
  f_yields4ljjsel<<endl;
  f_yields4ljjsel<<"--- 2016 scaled to new lumi ---"<<endl;
  f_yields4ljjsel<<"|Final state |signal HH |ttZ |ttH |bbH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |ZH+WH |others(=VVV+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel_2016_lumiScaled[HH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[TTZ][fs]  
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ttH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[bbH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[qqZZ][fs]+yield_4ljjsel_2016_lumiScaled[ggZZ][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ggH][fs]+yield_4ljjsel_2016_lumiScaled[VBF][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[VH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[VVV][fs]+yield_4ljjsel_2016_lumiScaled[TTW][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ZXbkg][fs]
                                         <<" |"<<endl;
  }
  f_yields4ljjsel<<endl;
  f_yields4ljjsel<<"--- 2017 ---"<<endl;
  f_yields4ljjsel<<"|Final state |signal HH |ttZ |ttH |bbH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |ZH+WH |others(=VVV+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel_2017[HH][fs]
                                         <<" |"<<yield_4ljjsel_2017[TTZ][fs]  
                                         <<" |"<<yield_4ljjsel_2017[ttH][fs]
                                         <<" |"<<yield_4ljjsel_2017[bbH][fs]
                                         <<" |"<<yield_4ljjsel_2017[qqZZ][fs]+yield_4ljjsel_2017[ggZZ][fs]
                                         <<" |"<<yield_4ljjsel_2017[ggH][fs]+yield_4ljjsel_2017[VBF][fs]
                                         <<" |"<<yield_4ljjsel_2017[VH][fs]
                                         <<" |"<<yield_4ljjsel_2017[VVV][fs]+yield_4ljjsel_2017[TTW][fs]
                                         <<" |"<<yield_4ljjsel_2017[ZXbkg][fs]
                                         <<" |"<<endl;
  }
  f_yields4ljjsel<<endl;
  f_yields4ljjsel<<"--- 2018 ---"<<endl;
  f_yields4ljjsel<<"|Final state |signal HH |ttZ |ttH |bbH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |ZH+WH |others(=VVV+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel_2018[HH][fs]
                                         <<" |"<<yield_4ljjsel_2018[TTZ][fs]  
                                         <<" |"<<yield_4ljjsel_2018[ttH][fs]
                                         <<" |"<<yield_4ljjsel_2018[bbH][fs]
                                         <<" |"<<yield_4ljjsel_2018[qqZZ][fs]+yield_4ljjsel_2018[ggZZ][fs]
                                         <<" |"<<yield_4ljjsel_2018[ggH][fs]+yield_4ljjsel_2018[VBF][fs]
                                         <<" |"<<yield_4ljjsel_2018[VH][fs]
                                         <<" |"<<yield_4ljjsel_2018[VVV][fs]+yield_4ljjsel_2018[TTW][fs]
                                         <<" |"<<yield_4ljjsel_2018[ZXbkg][fs]
                                         <<" |"<<endl;
  }
  f_yields4ljjsel<<endl;
  f_yields4ljjsel<<"--- full Run2 ---"<<endl;
  f_yields4ljjsel<<"|Final state |signal HH |ttZ |ttH |bbH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |ZH+WH |others(=VVV+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel_2016[HH][fs]
                                                +yield_4ljjsel_2017[HH][fs]
                                                +yield_4ljjsel_2018[HH][fs]
                                         <<" |"<<yield_4ljjsel_2016[TTZ][fs]
                                                +yield_4ljjsel_2017[TTZ][fs]
                                                +yield_4ljjsel_2018[TTZ][fs]
                                         <<" |"<<yield_4ljjsel_2016[ttH][fs]
                                                +yield_4ljjsel_2017[ttH][fs]
                                                +yield_4ljjsel_2018[ttH][fs]
                                         <<" |"<<yield_4ljjsel_2016[bbH][fs]
                                                +yield_4ljjsel_2017[bbH][fs]
                                                +yield_4ljjsel_2018[bbH][fs]
                                         <<" |"<<yield_4ljjsel_2016[qqZZ][fs]+yield_4ljjsel_2016[ggZZ][fs]
                                                +yield_4ljjsel_2017[qqZZ][fs]+yield_4ljjsel_2017[ggZZ][fs]
                                                +yield_4ljjsel_2018[qqZZ][fs]+yield_4ljjsel_2018[ggZZ][fs]
                                         <<" |"<<yield_4ljjsel_2016[ggH][fs]+yield_4ljjsel_2016[VBF][fs]
                                                +yield_4ljjsel_2017[ggH][fs]+yield_4ljjsel_2017[VBF][fs]
                                                +yield_4ljjsel_2018[ggH][fs]+yield_4ljjsel_2018[VBF][fs]
                                         <<" |"<<yield_4ljjsel_2016[VH][fs]
                                                +yield_4ljjsel_2017[VH][fs]
                                                +yield_4ljjsel_2018[VH][fs]
                                         <<" |"<<yield_4ljjsel_2016[VVV][fs]+yield_4ljjsel_2016[TTW][fs]
                                                +yield_4ljjsel_2017[VVV][fs]+yield_4ljjsel_2017[TTW][fs]
                                                +yield_4ljjsel_2018[VVV][fs]+yield_4ljjsel_2018[TTW][fs]
                                         <<" |"<<yield_4ljjsel_2016[ZXbkg][fs]
                                                +yield_4ljjsel_2017[ZXbkg][fs]
                                                +yield_4ljjsel_2018[ZXbkg][fs]
                                         <<" |"<<endl;
  }
  f_yields4ljjsel<<endl;
  f_yields4ljjsel<<"--- full Run2 with 2016 scaled to new lumi ---"<<endl;
  f_yields4ljjsel<<"|Final state |signal HH |ttZ |ttH |bbH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |ZH+WH |others(=VVV+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel_2016_lumiScaled[HH][fs]
                                                +yield_4ljjsel_2017[HH][fs]
                                                +yield_4ljjsel_2018[HH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[TTZ][fs]
                                                +yield_4ljjsel_2017[TTZ][fs]
                                                +yield_4ljjsel_2018[TTZ][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ttH][fs]
                                                +yield_4ljjsel_2017[ttH][fs]
                                                +yield_4ljjsel_2018[ttH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[bbH][fs]
                                                +yield_4ljjsel_2017[bbH][fs]
                                                +yield_4ljjsel_2018[bbH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[qqZZ][fs]+yield_4ljjsel_2016_lumiScaled[ggZZ][fs]
                                                +yield_4ljjsel_2017[qqZZ][fs]+yield_4ljjsel_2017[ggZZ][fs]
                                                +yield_4ljjsel_2018[qqZZ][fs]+yield_4ljjsel_2018[ggZZ][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ggH][fs]+yield_4ljjsel_2016_lumiScaled[VBF][fs]
                                                +yield_4ljjsel_2017[ggH][fs]+yield_4ljjsel_2017[VBF][fs]
                                                +yield_4ljjsel_2018[ggH][fs]+yield_4ljjsel_2018[VBF][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[VH][fs]
                                                +yield_4ljjsel_2017[VH][fs]
                                                +yield_4ljjsel_2018[VH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[VVV][fs]+yield_4ljjsel_2016_lumiScaled[TTW][fs]
                                                +yield_4ljjsel_2017[VVV][fs]+yield_4ljjsel_2017[TTW][fs]
                                                +yield_4ljjsel_2018[VVV][fs]+yield_4ljjsel_2018[TTW][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ZXbkg][fs]
                                                +yield_4ljjsel_2017[ZXbkg][fs]
                                                +yield_4ljjsel_2018[ZXbkg][fs]
                                         <<" |"<<endl;
  }
  f_yields4ljjsel<<endl;
  f_yields4ljjsel<<"--- full Run2, scheme as plots, with 2016 scaled to new lumi ---"<<endl;
  f_yields4ljjsel<<"|Final state |signal HH |SM H |qqZZ |ggZZ |others | Z+X | total but Z+X"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel_2016_lumiScaled[HH][fs]
                                                +yield_4ljjsel_2017[HH][fs]
                                                +yield_4ljjsel_2018[HH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ggH][fs]
                                                +yield_4ljjsel_2017[ggH][fs]
                                                +yield_4ljjsel_2018[ggH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[VBF][fs]
                                                +yield_4ljjsel_2017[VBF][fs]
                                                +yield_4ljjsel_2018[VBF][fs]
                                                +yield_4ljjsel_2016_lumiScaled[VH][fs]
                                                +yield_4ljjsel_2017[VH][fs]
                                                +yield_4ljjsel_2018[VH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[ttH][fs]
                                                +yield_4ljjsel_2017[ttH][fs]
                                                +yield_4ljjsel_2018[ttH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[bbH][fs]
                                                +yield_4ljjsel_2017[bbH][fs]
                                                +yield_4ljjsel_2018[bbH][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[qqZZ][fs]
                                                +yield_4ljjsel_2017[qqZZ][fs]
                                                +yield_4ljjsel_2018[qqZZ][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ggZZ][fs]
                                                +yield_4ljjsel_2017[ggZZ][fs]
                                                +yield_4ljjsel_2018[ggZZ][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[TTZ][fs]
                                                +yield_4ljjsel_2017[TTZ][fs]
                                                +yield_4ljjsel_2018[TTZ][fs]
                                                +yield_4ljjsel_2016_lumiScaled[TTW][fs]
                                                +yield_4ljjsel_2017[TTW][fs]
                                                +yield_4ljjsel_2018[TTW][fs]
                                                +yield_4ljjsel_2016_lumiScaled[VVV][fs]
                                                +yield_4ljjsel_2017[VVV][fs]
                                                +yield_4ljjsel_2018[VVV][fs]
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[ZXbkg][fs]
                                                +yield_4ljjsel_2017[ZXbkg][fs]
                                                +yield_4ljjsel_2018[ZXbkg][fs]
                                        // total but ZX
                                         <<" |"<<yield_4ljjsel_2016_lumiScaled[HH][fs]
                                                +yield_4ljjsel_2017[HH][fs]
                                                +yield_4ljjsel_2018[HH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[ggH][fs]
                                                +yield_4ljjsel_2017[ggH][fs]
                                                +yield_4ljjsel_2018[ggH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[VBF][fs]
                                                +yield_4ljjsel_2017[VBF][fs]
                                                +yield_4ljjsel_2018[VBF][fs]
                                                +yield_4ljjsel_2016_lumiScaled[VH][fs]
                                                +yield_4ljjsel_2017[VH][fs]
                                                +yield_4ljjsel_2018[VH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[ttH][fs]
                                                +yield_4ljjsel_2017[ttH][fs]
                                                +yield_4ljjsel_2018[ttH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[bbH][fs]
                                                +yield_4ljjsel_2017[bbH][fs]
                                                +yield_4ljjsel_2018[bbH][fs]
                                                +yield_4ljjsel_2016_lumiScaled[qqZZ][fs]
                                                +yield_4ljjsel_2017[qqZZ][fs]
                                                +yield_4ljjsel_2018[qqZZ][fs]
                                                +yield_4ljjsel_2016_lumiScaled[ggZZ][fs]
                                                +yield_4ljjsel_2017[ggZZ][fs]
                                                +yield_4ljjsel_2018[ggZZ][fs]
                                                +yield_4ljjsel_2016_lumiScaled[TTZ][fs]
                                                +yield_4ljjsel_2017[TTZ][fs]
                                                +yield_4ljjsel_2018[TTZ][fs]
                                                +yield_4ljjsel_2016_lumiScaled[TTW][fs]
                                                +yield_4ljjsel_2017[TTW][fs]
                                                +yield_4ljjsel_2018[TTW][fs]
                                                +yield_4ljjsel_2016_lumiScaled[VVV][fs]
                                                +yield_4ljjsel_2017[VVV][fs]
                                                +yield_4ljjsel_2018[VVV][fs]
                                         <<" |"<<endl;
  }
  f_yields4ljjsel.close();
 
 
}



void plot_mbbAndm4l_3years(){
  
  doPlots();
  scaleYields();
}
