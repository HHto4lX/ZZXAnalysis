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

using namespace std;


void doHistos(TString inputFilePath, TString outputFileName){

  TH1::SetDefaultSumw2(true);

  TFile* inFile;
  TTree* inTree;

  float f_ptjet1;
  float f_ptjet2;
  float f_etajet1;
  float f_etajet2;
  float f_bbmass;
  float f_JetPUValue_jet1;
  float f_JetPUValue_jet2;
  float f_JetPUID_jet1;
  float f_JetPUID_jet2;
  float f_weightsignal_nominal;
  

  TH1F* h_ptjet1          = new TH1F("h_ptjet1",          ";1st jet pT (GeV); Events/5 GeV", 60,  0.,  300.); h_ptjet1->Sumw2(true);
  TH1F* h_ptjet2          = new TH1F("h_ptjet2",          ";2nd jet pT (GeV); Events/5 GeV", 60,  0.,  300.); h_ptjet2->Sumw2(true);
  TH1F* h_etajet1         = new TH1F("h_etajet1",         ";1st jet eta; Events",            30, -2.6, 2.6);  h_etajet1->Sumw2(true);
  TH1F* h_etajet2         = new TH1F("h_etajet2",         ";2nd jet eta; Events",            30, -2.6, 2.6);  h_etajet2->Sumw2(true);
  TH1F* h_bbmass          = new TH1F("h_bbmass",          ";bb mass (GeV); Events/5 GeV",    60,  0.,  300.); h_bbmass->Sumw2(true);
  TH1F* h_JetPUValue_jet1 = new TH1F("h_JetPUValue_jet1", ";1st jet PUValue; Events",        22, -1.1, 1.1);  h_JetPUValue_jet1->Sumw2(true);
  TH1F* h_JetPUValue_jet2 = new TH1F("h_JetPUValue_jet2", ";2nd jet PUValue; Events",        22, -1.1, 1.1);  h_JetPUValue_jet2->Sumw2(true);
  TH1F* h_JetPUID_jet1    = new TH1F("h_JetPUID_jet1",    ";1st jet PUID; Events",           20, -0.2,  1.2 );  h_JetPUID_jet1->Sumw2(true);
  TH1F* h_JetPUID_jet2    = new TH1F("h_JetPUID_jet2",    ";2nd jet PUID; Events",           20, -0.2,  1.2 );  h_JetPUID_jet2->Sumw2(true);


  inFile = TFile::Open(inputFilePath);
  inTree = (TTree*)inFile->Get("reducedTree");
  inTree->SetBranchAddress("f_weightsignal_nominal", &f_weightsignal_nominal);
  inTree->SetBranchAddress("f_ptjet1",  &f_ptjet1);
  inTree->SetBranchAddress("f_ptjet2",  &f_ptjet2);
  inTree->SetBranchAddress("f_etajet1", &f_etajet1);
  inTree->SetBranchAddress("f_etajet2", &f_etajet2);
  inTree->SetBranchAddress("f_bbmass",  &f_bbmass);
  inTree->SetBranchAddress("f_JetPUValue_jet1", &f_JetPUValue_jet1);
  inTree->SetBranchAddress("f_JetPUValue_jet2", &f_JetPUValue_jet2);
  inTree->SetBranchAddress("f_JetPUID_jet1",    &f_JetPUID_jet1);
  inTree->SetBranchAddress("f_JetPUID_jet2",    &f_JetPUID_jet2);

  // loop over input tree
  Long64_t entries = inTree->GetEntries();

  for (Long64_t entry = 0; entry < entries; entry++){  

    inTree->GetEntry(entry);

    h_ptjet1         ->Fill(f_ptjet1         , f_weightsignal_nominal);  
    h_ptjet2         ->Fill(f_ptjet2         , f_weightsignal_nominal); 
    h_etajet1        ->Fill(f_etajet1        , f_weightsignal_nominal); 
    h_etajet2        ->Fill(f_etajet2        , f_weightsignal_nominal);
    h_bbmass         ->Fill(f_bbmass         , f_weightsignal_nominal);
    h_JetPUValue_jet1->Fill(f_JetPUValue_jet1, f_weightsignal_nominal);
    h_JetPUValue_jet2->Fill(f_JetPUValue_jet2, f_weightsignal_nominal);
    h_JetPUID_jet1   ->Fill(f_JetPUID_jet1   , f_weightsignal_nominal);
    h_JetPUID_jet2   ->Fill(f_JetPUID_jet2   , f_weightsignal_nominal);

  } // end loop over input tree


  // save histos into a root file
  TFile* fOutHistos = new TFile(outputFileName,"recreate");
  fOutHistos->cd();

  h_ptjet1         ->Write(); delete h_ptjet1         ;
  h_ptjet2         ->Write(); delete h_ptjet2         ;
  h_etajet1        ->Write(); delete h_etajet1        ;
  h_etajet2        ->Write(); delete h_etajet2        ;
  h_bbmass         ->Write(); delete h_bbmass         ;
  h_JetPUValue_jet1->Write(); delete h_JetPUValue_jet1;
  h_JetPUValue_jet2->Write(); delete h_JetPUValue_jet2;
  h_JetPUID_jet1   ->Write(); delete h_JetPUID_jet1   ;
  h_JetPUID_jet2   ->Write(); delete h_JetPUID_jet2   ;

  fOutHistos->Close();

} // end doHistos




void doPlots_conparison(TString nameFile1, TString nameFile2, TString outPath_comparison){

  TString lumiText = "41.5 fb^{-1}";

  

  cout<<"creating output dir "<<outPath_comparison<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_comparison)).c_str()); // create output dir


  Int_t nPlots = 9;
  TString sPlotsNames[] = {
    "ptjet1",          
    "ptjet2",         
    "etajet1",        
    "etajet2",        
    "bbmass",         
    "JetPUValue_jet1",
    "JetPUValue_jet2",
    "JetPUID_jet1",   
    "JetPUID_jet2",   
  };


  TFile* inFile1 = TFile::Open(nameFile1);  
  TH1F* h_file1[nPlots];
  for(int i=0; i<nPlots; i++){
    h_file1[i] = (TH1F*)inFile1->Get("h_"+ sPlotsNames[i]);
  }

  TFile* inFile2 = TFile::Open(nameFile2);  
  TH1F* h_file2[nPlots];
  for(int i=0; i<nPlots; i++){
    h_file2[i] = (TH1F*)inFile2->Get("h_"+ sPlotsNames[i]);
  }

  TCanvas* c_plots[nPlots];
  TLegend* leg  [nPlots];
  
  gStyle->SetOptStat(0);
 
  for(int i=0; i<nPlots; i++){
    
    c_plots[i] = new TCanvas("c_" + sPlotsNames[i], "c_" + sPlotsNames[i], 800, 600);
    c_plots[i]->cd();
    h_file2[i]->SetLineColor(kRed);
    h_file2[i]->Draw("hist");
    h_file1[i]->SetLineColor(kBlue);
    h_file1[i]->Draw("hist same");


    if(sPlotsNames[i] == "JetPUValue_jet1" || 
       sPlotsNames[i] == "JetPUValue_jet2" || 
       sPlotsNames[i] == "JetPUID_jet1"    ||
       sPlotsNames[i] == "JetPUID_jet2" ) c_plots[i]->SetLogy();

    leg[i] = new TLegend(0.63,0.8,0.77,0.89);
    leg[i]->AddEntry(h_file1[i], "old PUID", "l");
    leg[i]->AddEntry(h_file2[i], "new PUID", "l"); 
    leg[i]->SetFillColor(kWhite);
    leg[i]->SetLineColor(kWhite);
    leg[i]->SetTextFont(43);
    leg[i]->SetTextSize(20);
    leg[i]->Draw();

    c_plots[i]->Update();

    c_plots[i]->SaveAs(outPath_comparison + "/" + c_plots[i]->GetName() + ".png");


    cout<<h_file1[i]->Integral()<<" "<<h_file2[i]->Integral()<<endl;

  } // end loop ove histos


}


void plot_comparison(){

  // TString inFile1  = "ntuples_oldPUID/mvaNtuples_2bjet_2017_SR_4ljjsel_fs4mu/reduced_HH4lbb_Angela.root";
  // TString outFile1 = "fout_oldPUID_4mu.root";
  // TString inFile2  = "ntuples_newPUID/mvaNtuples_2bjet_2017_SR_4ljjsel_fs4mu/reduced_HH4lbb_Angela.root";
  // TString outFile2 = "fout_newPUID_4mu.root";
  // TString outPath_comparison = "plots_comparison_2017_4mu";

  // TString inFile1  = "ntuples_oldPUID/mvaNtuples_2bjet_2017_SR_4ljjsel_fs4e/reduced_HH4lbb_Angela.root";
  // TString outFile1 = "fout_oldPUID_4e.root";
  // TString inFile2  = "ntuples_newPUID/mvaNtuples_2bjet_2017_SR_4ljjsel_fs4e/reduced_HH4lbb_Angela.root";
  // TString outFile2 = "fout_newPUID_4e.root";
  // TString outPath_comparison = "plots_comparison_2017_4e";

  TString inFile1  = "ntuples_oldPUID/mvaNtuples_2bjet_2017_SR_4ljjsel_fs2e2mu/reduced_HH4lbb_Angela.root";
  TString outFile1 = "fout_oldPUID_2e2mu.root";
  TString inFile2  = "ntuples_newPUID/mvaNtuples_2bjet_2017_SR_4ljjsel_fs2e2mu/reduced_HH4lbb_Angela.root";
  TString outFile2 = "fout_newPUID_2e2mu.root";
  TString outPath_comparison = "plots_comparison_2017_2e2mu";



  doHistos(inFile1, outFile1);
  doHistos(inFile2, outFile2);
  doPlots_conparison(outFile1, outFile2, outPath_comparison);
}
