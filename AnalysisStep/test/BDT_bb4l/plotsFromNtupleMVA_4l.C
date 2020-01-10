// ******************************************
// use: root -l -b -q plotsFromNtupleMVA.C++
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

  
  TFile* inFile;
  TTree* inTree;
  
  Float_t f_weight;
  Float_t f_lept1_pt;
  Float_t f_lept2_pt;
  Float_t f_lept3_pt;
  Float_t f_lept4_pt;
  Float_t f_MET_norm;
  Float_t f_ZZmass;



  TH1F* h_4leptonsPt = new TH1F("h_4leptonsPt", "; 4 leptons pT (GeV); Events / 2 GeV",  50,   0., 100.);  h_4leptonsPt->Sumw2(true);
  TH1F* h_MET        = new TH1F("h_MET",        "; MET (GeV); Events / 4 GeV",           50,   0., 200.);  h_MET->Sumw2(true);
  TH1F* h_M4l        = new TH1F("h_M4l",        "; m_{4l} (GeV); Events / 2 GeV",        65,  70., 200.);  h_M4l->Sumw2(true);


  inFile = TFile::Open( inputFileName );
  inTree = (TTree*)inFile->Get("reducedTree");
  inTree->SetBranchAddress("f_weight",   &f_weight);
  inTree->SetBranchAddress("f_lept1_pt", &f_lept1_pt);
  inTree->SetBranchAddress("f_lept2_pt", &f_lept2_pt);
  inTree->SetBranchAddress("f_lept3_pt", &f_lept3_pt);
  inTree->SetBranchAddress("f_lept4_pt", &f_lept4_pt);
  inTree->SetBranchAddress("f_MET_norm", &f_MET_norm);
  inTree->SetBranchAddress("f_ZZmass",   &f_ZZmass);


  // loop over input tree
  Long64_t entries = inTree->GetEntries();

  for (Long64_t entry = 0; entry < entries; entry++){  

    inTree->GetEntry(entry);


    // --- plots in H->4l standard sel
    // m4l plot
    if(isDATA){
      if(DOBLINDHISTO){
        if(f_ZZmass < 115 || f_ZZmass > 135){
          h_M4l->Fill(f_ZZmass, f_weight);      
        }
      } 
      else{
        h_M4l->Fill(f_ZZmass, f_weight);      
      } 
    }
    else{
      h_M4l->Fill(f_ZZmass, f_weight);      
    }
    // --- end plots in H->4l standard sel



    // --- plots in the m4l sidebands
    if(f_ZZmass < 115 || f_ZZmass > 135){
  
      h_4leptonsPt->Fill(f_lept1_pt, f_weight);  
      h_4leptonsPt->Fill(f_lept2_pt, f_weight);  
      h_4leptonsPt->Fill(f_lept3_pt, f_weight);  
      h_4leptonsPt->Fill(f_lept4_pt, f_weight);  


      h_MET->Fill(f_MET_norm, f_weight);
    
    }// end plots in the m4l sidebands

 

  }//end loop over input tree


  // save histos into a root file
  TFile* fOutHistos = new TFile(outputFileName,"recreate");
  fOutHistos->cd();
  
  h_4leptonsPt->Write(); delete h_4leptonsPt;
  h_MET->Write();        delete h_MET;
  h_M4l->Write();        delete h_M4l;

  fOutHistos->Close();

}



// main function
void plotsFromNtupleMVA_4l(){


  TString inputFilePath = "200109_mvaNtuples_4lsel_full4lmass/";
  TString inputFileName[] = {"HH4lbb",
			     "ggH125",
			     "VBFH125",
                             "WplusH125",
                             "WminusH125",
                             "ZH125",
                             "bbH125",
                             "ttH125",
                             "ZZTo4lext1",
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
                             "Z+X", 
                             "AllData",
  };

  
  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<< "number of input files: " << nInputFiles<<endl;

  string outputFilePath = "histos_4lbb_4lsel";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir
  

  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doHisto(inputFilePath + "reduced_" + inputFileName[i] + ".root", outputFilePath + "/histos_" + inputFileName[i] + ".root");
  }

}


