// -----------------------
// run with: root -b -q printTreeContent.C++
// ------------------------
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



void printTreeContent(){

  TFile* inputFile;
  TTree* inputTree;

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;

  vector<Float_t> *JetPt   = 0;
  vector<Float_t> *JetEta  = 0;
  vector<Float_t> *JetPhi  = 0;
  vector<Short_t> *JetPUID = 0;

  inputFile = TFile::Open("ZZXAnalysis.root");
  inputTree = (TTree*)inputFile->Get("ZTree/candTree");

  inputTree->SetBranchAddress("RunNumber", &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber", &nLumi);
  inputTree->SetBranchAddress("JetPt", &JetPt);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetPhi", &JetPhi);
  inputTree->SetBranchAddress("JetPUID", &JetPUID);

  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; z++){

    inputTree->GetEntry(z);

    if( (nLumi == 23216 && nEvent == 21985025) ||
        (nLumi == 23216 && nEvent == 21985385) ||
        (nLumi == 23265 && nEvent == 22031870) ||
        (nLumi == 23267 && nEvent == 22033543) ||
        (nLumi == 23267 && nEvent == 22033740) ||
        (nLumi == 23268 && nEvent == 22034336) ||
        (nLumi == 23268 && nEvent == 22034491) ||
        (nLumi == 23269 && nEvent == 22035343) ||
        (nLumi == 23272 && nEvent == 22038385) ||
        (nLumi == 23993 && nEvent == 22720467)){
      cout<<nRun<<":"<<nLumi<<":"<<nEvent<<" ";
	for(int i = 0; i<JetPt->size(); i++){
          cout<<JetPt->at(i)<<" "<<JetEta->at(i)<<" "<<JetPhi->at(i)<<" "<<JetPUID->at(i)<<" ";
        }
	cout<<endl;
    }



  }


}
