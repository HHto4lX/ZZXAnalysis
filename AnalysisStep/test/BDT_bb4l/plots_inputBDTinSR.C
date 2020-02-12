// ******************************************
// use: root -l -b -q plots_inputBDTinSR.C++
// ******************************************

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

#define DOBLINDHISTO 0



void doHisto(TString inputFileName, TString outputFileName){

  bool isDATA = false;
  bool isZX   = false;
  if ( inputFileName.Contains("AllData") ) isDATA = true;
  if ( inputFileName.Contains("Z+X") ) isZX   = true;
  cout<<"isDATA "<<isDATA<<endl;
  cout<<"isZX "<<isZX<<endl;


  TH1::SetDefaultSumw2(true);

  
  TFile* inFile;
  TTree* inTree;
  
  Float_t f_weight;
  Float_t f_bdiscjet1;  
  Float_t f_bdiscjet2;  
  Float_t f_ptjet1;     
  Float_t f_ptjet2;     
  Float_t f_MET_norm;   
  Float_t f_deltar_norm;
  Float_t f_massjetjet;


  TH1F* h_bTagger_jet1 = new TH1F("h_bTagger_jet1", "; DeepCSV; Events/0.04",             25, 0., 1.);   h_bTagger_jet1->Sumw2(true);
  TH1F* h_bTagger_jet2 = new TH1F("h_bTagger_jet2", "; DeepCSV; Events/0.04",             25, 0., 1.);   h_bTagger_jet2->Sumw2(true);
  TH1F* h_pT_jet1      = new TH1F("h_pT_jet1",      "; pT (GeV); Events/5 GeV",           40, 0., 200.); h_pT_jet1     ->Sumw2(true);
  TH1F* h_pT_jet2      = new TH1F("h_pT_jet2",      "; pT (GeV); Events/5 GeV",           40, 0., 200.); h_pT_jet2     ->Sumw2(true);
  TH1F* h_MET          = new TH1F("h_MET",          "; MET (GeV); Events/5 GeV",          40, 0., 200.); h_MET         ->Sumw2(true);
  TH1F* h_DeltaR       = new TH1F("h_DeltaR",       "; Delta R; Events ",                 25, 0., 10.);  h_DeltaR      ->Sumw2(true);
  TH1F* h_mjj          = new TH1F("h_mjj",          "; Mass di-jet (GeV); Events/5 GeV ", 40, 0., 200.); h_mjj         ->Sumw2(true);
  
  


  inFile = TFile::Open( inputFileName );
  inTree = (TTree*)inFile->Get("reducedTree");
  inTree->SetBranchAddress("f_weight",      &f_weight);
  inTree->SetBranchAddress("f_bdiscjet1",   &f_bdiscjet1);
  inTree->SetBranchAddress("f_bdiscjet2",   &f_bdiscjet2);
  inTree->SetBranchAddress("f_ptjet1",      &f_ptjet1);
  inTree->SetBranchAddress("f_ptjet2",      &f_ptjet2);
  inTree->SetBranchAddress("f_MET_norm",    &f_MET_norm);
  inTree->SetBranchAddress("f_deltar_norm", &f_deltar_norm);
  inTree->SetBranchAddress("f_massjetjet",  &f_massjetjet);



  // loop over input tree
  Long64_t entries = inTree->GetEntries();

  for (Long64_t entry = 0; entry < entries; entry++){  

    inTree->GetEntry(entry);

    h_bTagger_jet1->Fill(f_bdiscjet1,   f_weight);    
    h_bTagger_jet2->Fill(f_bdiscjet2,   f_weight);
    h_pT_jet1     ->Fill(f_ptjet1,      f_weight);
    h_pT_jet2     ->Fill(f_ptjet2,      f_weight);
    h_MET         ->Fill(f_MET_norm,    f_weight);
    h_DeltaR      ->Fill(f_deltar_norm, f_weight);
    h_mjj         ->Fill(f_massjetjet,  f_weight);
 

  }//end loop over input tree


  // save histos into a root file
  TFile* fOutHistos = new TFile(outputFileName,"recreate");
  fOutHistos->cd();
  
  h_bTagger_jet1->Write(); delete h_bTagger_jet1;
  h_bTagger_jet2->Write(); delete h_bTagger_jet2;
  h_pT_jet1     ->Write(); delete h_pT_jet1     ;
  h_pT_jet2     ->Write(); delete h_pT_jet2     ;
  h_MET         ->Write(); delete h_MET         ;
  h_DeltaR      ->Write(); delete h_DeltaR      ;
  h_mjj         ->Write(); delete h_mjj         ;

  fOutHistos->Close();

}



// main function
void plots_inputBDTinSR(){


  TString inputFilePath = "200212_mvaNtuples_1bjet1Pt_4l/";
  TString inputFileName[] = {"HH4lbb_Angela",
			     "ggH125",
			     "VBFH125",
                             "WplusH125",
                             "WminusH125",
                             "ZH125",
                             "bbH125",
                             "ttH125",
                             "ZZTo4lamcatnlo",
                             "TTZJets_M10_MLMext1",
                             "TTZToLL_M1to1O_MLM",
                             "TTWJetsToLNu",
                             "ggTo4e_Contin_MCFM701",
                             "ggTo4mu_Contin_MCFM701",
                             "ggTo4tau_Contin_MCFM701",
                             "ggTo2e2mu_Contin_MCFM701",
                             "ggTo2e2tau_Contin_MCFM701",
                             "ggTo2mu2tau_Contin_MCFM701",
                             "WWZ",
                             "WZZ",
                             "ZZZ",
                             //"Z+X", 
                             "AllData",
  };

  
  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<< "number of input files: " << nInputFiles<<endl;

  string outputFilePath = "histos_plotsinputBDT_SR";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir
  

  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doHisto(inputFilePath + "reduced_" + inputFileName[i] + ".root", outputFilePath + "/histos_" + inputFileName[i] + ".root");
  }

}


