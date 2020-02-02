// ******************************************
// use: root -l -b -q jetMatchingStudy.C++
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



void jetMatchingStudy(){

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
  Float_t xsec;

  Short_t ZZsel;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  Float_t ZZMass;

  vector<Float_t> *JetPt     = 0;
  vector<Float_t> *JetEta    = 0;
  vector<Float_t> *JetPhi    = 0;
  vector<Float_t> *JetMass   = 0;
  vector<Float_t> *JetBTagger = 0;

  vector<Float_t> *GENjetParentID   = 0;
  vector<Float_t> *prunedGenPartID   = 0;
  vector<Float_t> *prunedGenPartPt   = 0;
  vector<Float_t> *prunedGenPartEta   = 0;
  vector<Float_t> *prunedGenPartPhi   = 0;
  vector<Float_t> *prunedGenPartMass   = 0;
  vector<Float_t> *prunedGenMotherID   = 0;

  // counters
  int count2match_Method1 = 0;
  int count1match_Method1 = 0;
  int count0match_Method1 = 0;
  int countTOTmatch_Method1 = 0;

  

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
  inputTree->SetBranchAddress("xsec", &xsec);

  inputTree->SetBranchAddress("ZZsel", &ZZsel);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("ZZMass", &ZZMass); 
  
  inputTree->SetBranchAddress("JetPt", &JetPt);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetMass", &JetMass);
  inputTree->SetBranchAddress("JetPhi",  &JetPhi);
  inputTree->SetBranchAddress("JetBTagger",  &JetBTagger);

  inputTree->SetBranchAddress("GENjetParentID",  &GENjetParentID);
  inputTree->SetBranchAddress("prunedGenPartEta", &prunedGenPartEta );
  inputTree->SetBranchAddress("prunedGenPartPhi", &prunedGenPartPhi );
  inputTree->SetBranchAddress("prunedGenPartPt", &prunedGenPartPt );
  inputTree->SetBranchAddress("prunedGenPartMass", &prunedGenPartMass );
  inputTree->SetBranchAddress("prunedGenPartID", &prunedGenPartID );
  inputTree->SetBranchAddress("prunedGenMotherID", &prunedGenMotherID );




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


    // fill eventweight
    Double_t eventWeight = partialSampleWeight * xsec * overallEventWeight ;  //kfactor e l1prefiring non ci sono, tanto uso solo HH sample


    // mass cut: signal region
    if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV



    // Delta R for Gen jet matching
    double DELTAR = 0.4;    
 
    // -------------------------------------------
    // --- jet selection METHOD1: 2 high pT jets
    vector<TLorentzVector> JetVec_Method1; // TLorentz vector with all Jets per Event
    vector<TLorentzVector> JetPair_Method1; // TLorentz vector with all Jets Pairs
    vector<float> JetpT_Method1; // vector with the pT per each Jet of the Event
    vector<float> JetBinfo_Method1; // vector with the Btagger per each Jet of the Event

    int d1_Method1 = -999; // position of higest pT jet
    int d2_Method1 = -999; // position of second-highest pT jet


    for (UInt_t j = 0; j < JetPt->size(); j++)
    {
      if ( (fabs ( JetEta->at(j) ) > 2.4) || (JetPt->at(j) < 20 ) ) continue; // pt cut 20GeV from ntuplizer 
	  
      TLorentzVector temp;
      temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
      JetVec_Method1.push_back(temp);
      JetpT_Method1.push_back(JetPt->at(j));
      JetBinfo_Method1.push_back(JetBTagger->at(j));
    }

    // at least 2 jets in the acceptance
    if (JetVec_Method1.size() >= 2){

      // count number of events with at least 2 jets
      countTOTmatch_Method1++;

      // find the position of the highest pT jet
      d1_Method1 = distance( JetpT_Method1.begin(), max_element(JetpT_Method1.begin(), JetpT_Method1.end()));

      // find the position of the second-highest pT jet
      float max = -9999.;
      for(UInt_t i=0; i<JetpT_Method1.size(); i++){
        if(i == d1_Method1) continue;
        float temp = JetpT_Method1.at(i);
        if(temp > max){
          max = temp;
          d2_Method1 = i;
        }
      }

      cout<<d1_Method1<<" "<<d2_Method1<<endl;      
    
      // find match between 1st reco jets and GEN jets
      bool bool_match1_Method1 = false;
      for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
        if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
   
          float deltaPhi1 = JetVec_Method1.at(d1_Method1).Phi() - prunedGenPartPhi->at(pr);
          if(fabs(deltaPhi1) > acos(-1)){ deltaPhi1 = (2*acos(-1)) - fabs(deltaPhi1); }
          float deltaEta1 = JetVec_Method1.at(d1_Method1).Eta() - prunedGenPartEta->at(pr);

          float deltaR1 = sqrt(deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1);

          if(deltaR1 < DELTAR){ bool_match1_Method1 = true ;}
        }
      }

      // find match between 2nd reco jets and GEN jets
      bool bool_match2_Method1 = false;
      for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
        if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
   
          float deltaPhi2 = JetVec_Method1.at(d2_Method1).Phi() - prunedGenPartPhi->at(pr);
          if(fabs(deltaPhi2) > acos(-1)){ deltaPhi2 = (2*acos(-1)) - fabs(deltaPhi2); }
          float deltaEta2 = JetVec_Method1.at(d2_Method1).Eta() - prunedGenPartEta->at(pr);

          float deltaR2 = sqrt(deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2);

          if(deltaR2 < 0.4){ bool_match2_Method1 = true ;}
        }
      }

      cout<<bool_match1_Method1<<" "<<bool_match2_Method1<<endl;

      if(bool_match1_Method1 && bool_match2_Method1) {count2match_Method1++;}
      else if( (bool_match1_Method1 && !bool_match2_Method1) || (!bool_match1_Method1 && bool_match2_Method1) ) {count1match_Method1++;}
      else {count0match_Method1++;}

      cout<<count0match_Method1<<" "<<count1match_Method1<<" "<<count2match_Method1<<" "<<countTOTmatch_Method1<<endl;

    } 
    // --- end jet selection METHOD1: 2 high pT jets
    // ----------------------------------------------


  }// end loop over tree entries




}
