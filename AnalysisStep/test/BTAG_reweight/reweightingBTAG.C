// ******************************************
// use: root -l -b -q reweightingBTAG.C++
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

#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#include "evalEventSF.C"

using namespace std;



void reweightingBTAG(){

  float lumi = 59.7; //fb-1

  TFile* inputFile;
  TTree* inputTree;
  TH1F* hCounters;
  int NGenEvt;
  Float_t gen_sumWeights;
  Float_t partialSampleWeight;

  Float_t genBR;
  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t overallEventWeight;
  //  Float_t xsec = 0.0000447896; //pb
  Float_t xsec = 0.00001017; //pb Angela

  Short_t ZZsel;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  Float_t ZZMass;

  vector<Float_t> *JetPt     = 0;
  vector<Float_t> *JetEta    = 0;
  vector<Float_t> *JetPhi    = 0;
  vector<Float_t> *JetMass   = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Short_t> *JetHadronFlavour = 0;






  TString inFile = "/eos/user/a/acappati/samples_4lX/allsamples/HH4lbb/ZZXAnalysis.root";
  inputFile =  TFile::Open( inFile );


  hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
  NGenEvt = (Float_t)hCounters->GetBinContent(1);
  gen_sumWeights = (Float_t)hCounters->GetBinContent(40);
  partialSampleWeight = lumi * 1000 / gen_sumWeights;
  inputTree = (TTree*)inputFile->Get("ZZTree/candTree");


  inputTree->SetBranchAddress("RunNumber", &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber", &nLumi);
  inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
  //  inputTree->SetBranchAddress("xsec", &xsec);

  inputTree->SetBranchAddress("ZZsel", &ZZsel);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("ZZMass", &ZZMass); 
  
  inputTree->SetBranchAddress("JetPt", &JetPt);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetMass", &JetMass);
  inputTree->SetBranchAddress("JetPhi",  &JetPhi);
  inputTree->SetBranchAddress("JetBTagger",  &JetBTagger);
  inputTree->SetBranchAddress("JetHadronFlavour",  &JetHadronFlavour);




  Long64_t entries = inputTree->GetEntries();
  std::cout<<"Processing file: "<< inFile << "\nNumber of entries: " << entries << endl;

  for (Long64_t entry = 0; entry < entries; entry++)
  {  

    inputTree->GetEntry(entry);

    // trigger check
    if( !(ZZsel >= 0) ) continue;


    // check if only 4 lepton
    if(LepEta->size() != 4){
	cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
    }


    // mass cut: signal region
    if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV


 

    // jet in the acceptance
    vector<TLorentzVector> JetVec_; 
    vector<float> JetpT_; 
    vector<float> Jeteta_; 
    vector<float> Jethadronflavour_; 
    vector<float> JetBtagger_; 

    int d1_ = -999; // position of higest btagger jet
    int d2_ = -999; // position of highest pT jet


    for (UInt_t j = 0; j < JetPt->size(); j++)
    {
      if ( (fabs ( JetEta->at(j) ) > 2.4) || (JetPt->at(j) < 20 ) ) continue; // pt cut 20GeV from ntuplizer 
	  
      TLorentzVector temp;
      temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
      JetVec_.push_back(temp);
      JetpT_.push_back(JetPt->at(j));
      Jeteta_.push_back(JetEta->at(j));
      Jethadronflavour_.push_back(JetHadronFlavour->at(j));
      JetBtagger_.push_back(JetBTagger->at(j));
    }


    // at least 2 jets in the acceptance
    if (JetVec_.size() >= 2){

      cout<<"ciao"<<endl;

    } 



    // fill eventweight
    Double_t eventWeight = partialSampleWeight * xsec * overallEventWeight ;  //kfactor e l1prefiring non ci sono, tanto uso solo HH sample


  }// end loop over tree entries

  




}
