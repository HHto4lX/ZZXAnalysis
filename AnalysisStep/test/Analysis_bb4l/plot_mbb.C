#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TFrame.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
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


void plot_mbb(){

  TFile* inFile1 = TFile::Open("out_1bjet/f_histos_h1_2018.root");  
  TH1F* hmbb1 = (TH1F*)inFile1->Get("h1_mbb_4ljjsel_HH_4l_2018");

  TFile* inFile2 = TFile::Open("out_2bjet/f_histos_h1_2018.root");  
  TH1F* hmbb2 = (TH1F*)inFile2->Get("h1_mbb_4ljjsel_HH_4l_2018");

  TCanvas* c_mbb = new TCanvas();
  c_mbb->cd();
  hmbb1->SetLineColor(kBlue);
  hmbb1->Draw("hist");
  hmbb2->SetLineColor(kRed);
  hmbb2->Draw("hist same");
  
  TLegend* leg = new TLegend(0.65,0.7,0.94,0.9);
  leg->AddEntry(hmbb2,       "2 highest btagger", "l");
  leg->AddEntry(hmbb1,       "highest btagger +", "l");
  leg->AddEntry((TObject*)0, "highest pT","");
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kBlack);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->Draw();

  c_mbb->Update();

  c_mbb->SaveAs("mbb.png");
  



}
