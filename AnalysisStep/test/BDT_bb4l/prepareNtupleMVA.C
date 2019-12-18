// ***********************
// prepare ntuples for mva
// run with:
//
// root -l -b -q prepareNtupleMVA.C++
// ***********************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TFrame.h"
#include "TF1.h"
#include "TH1F.h"
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

#define JETSELECTION 1



void doNtuplesForMVA(TString inFile, TString outFile, float lumi)
{

  bool isDATA = false;
  bool isZX   = false;
  if ( inFile.Contains("AllData") ) isDATA = true;
  if ( inFile.Contains("Z+X") ) isZX   = true;
  cout<<"isDATA "<<isDATA<<endl;
  cout<<"isZX "<<isZX<<endl;

  // input file branches
  TFile* inputFile;
  TTree* inputTree;
  TH1F* hCounters;
  int NGenEvt;
  Float_t gen_sumWeights;
  Float_t partialSampleWeight = 1.;
  Float_t weight; //ZX weight
  
  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t overallEventWeight;
  Float_t xsec;

  Float_t KFactor_QCD_ggZZ_Nominal;
  Float_t KFactor_EW_qqZZ;
  Float_t KFactor_QCD_qqZZ_dPhi;
  Float_t KFactor_QCD_qqZZ_M;
  Float_t KFactor_QCD_qqZZ_Pt;

  Short_t ZZsel;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPt = 0;
  Float_t ZZMass;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Float_t ZZPt;
  Float_t ZZEta;
  Float_t ZZPhi;
  Short_t Z1Flav;
  Short_t Z2Flav;
  vector<Float_t> *JetPt     = 0;
  vector<Float_t> *JetEta    = 0;
  vector<Float_t> *JetPhi    = 0;
  vector<Float_t> *JetMass   = 0;
  vector<Float_t> *JetIsBtagged = 0;
  vector<Float_t> *JetBTagger = 0;
  Float_t PFMET;

  Float_t yield = 0.;


  inputFile =  TFile::Open( inFile );

  if(!isZX){
    hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
    NGenEvt = (Float_t)hCounters->GetBinContent(1);
    gen_sumWeights = (Float_t)hCounters->GetBinContent(40);
    partialSampleWeight = lumi * 1000 / gen_sumWeights;
    inputTree = (TTree*)inputFile->Get("ZZTree/candTree");
  }
  else inputTree = (TTree*)inputFile->Get("candTree");

  // set branch addresses
  if(!isZX){
    inputTree->SetBranchAddress("RunNumber", &nRun);
    inputTree->SetBranchAddress("EventNumber", &nEvent);
    inputTree->SetBranchAddress("LumiNumber", &nLumi);
    if (!isDATA) inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
    if (!isDATA) inputTree->SetBranchAddress("xsec", &xsec);
  }
  if(isZX){ inputTree->SetBranchAddress("weight", &weight); }
  inputTree->SetBranchAddress("ZZsel", &ZZsel);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("ZZMass", &ZZMass);  
  inputTree->SetBranchAddress("Z1Flav", &Z1Flav);
  inputTree->SetBranchAddress("Z2Flav", &Z2Flav);
  inputTree->SetBranchAddress("Z1Mass", &Z1Mass);
  inputTree->SetBranchAddress("Z2Mass", &Z2Mass);
  inputTree->SetBranchAddress("ZZPt", &ZZPt);
  inputTree->SetBranchAddress("ZZEta", &ZZEta);
  inputTree->SetBranchAddress("ZZPhi", &ZZPhi);
  inputTree->SetBranchAddress("JetPt", &JetPt);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetMass", &JetMass);
  inputTree->SetBranchAddress("JetPhi",  &JetPhi);
  inputTree->SetBranchAddress("JetBTagger",  &JetBTagger);
  inputTree->SetBranchAddress("JetIsBtagged",  &JetIsBtagged);  
  inputTree->SetBranchAddress("PFMET",  &PFMET);  
  // ggZZ samples
  if(inFile.Contains("ggTo"))
  {
    inputTree->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal); 
  }
  // qqZZ samples  
  if(inFile.Contains("ZZTo4l"))       
  {
    inputTree->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
    inputTree->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi);
    inputTree->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
    inputTree->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt);
  }

  
  //output file 
  float f_lept1_ptsignal_4mu;
  float f_lept2_ptsignal_4mu;
  float f_lept3_ptsignal_4mu;
  float f_lept4_ptsignal_4mu;
  float f_bdiscjet1signal_4mu;
  float f_bdiscjet2signal_4mu;
  float f_deltarsignal_4mu;
  float f_METsignal_4mu;
  float f_weightsignal_4mu;

  float f_ZZmass;
  float f_bbmass;
  float f_HHmass;


  TFile *f = new TFile(outFile,"recreate");
  TTree *tnew = new TTree("reducedTree","");

  tnew->Branch("f_lept1_pt",    &f_lept1_ptsignal_4mu);
  tnew->Branch("f_lept2_pt",    &f_lept2_ptsignal_4mu);
  tnew->Branch("f_lept3_pt",    &f_lept3_ptsignal_4mu);
  tnew->Branch("f_lept4_pt",    &f_lept4_ptsignal_4mu);
  tnew->Branch("f_bdiscjet1",   &f_bdiscjet1signal_4mu);
  tnew->Branch("f_bdiscjet2",   &f_bdiscjet2signal_4mu);
  tnew->Branch("f_deltar_norm", &f_deltarsignal_4mu); 
  tnew->Branch("f_MET_norm",    &f_METsignal_4mu); 
  tnew->Branch("f_weight",      &f_weightsignal_4mu); 
  tnew->Branch("f_ZZmass",      &f_ZZmass); 
  tnew->Branch("f_bbmass",      &f_bbmass); 
  tnew->Branch("f_HHmass",      &f_HHmass); 


  // loop over input tree
  Long64_t entries = inputTree->GetEntries();
  std::cout<<"Processing file: "<< inFile << "\nNumber of entries: " << entries << endl;

  for (Long64_t entry = 0; entry < entries; entry++)
  {  
    inputTree->GetEntry(entry);
    

    if( !(ZZsel>=0) ) continue;

    // 4l selection
    if(LepEta->size() != 4)
      {
	cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
      }


    // mass cut
    if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV

 


    if(JETSELECTION){
      // jet quantities
      vector<TLorentzVector> JetVec; // TLorentz vector with all Jets per Event
      vector<TLorentzVector> JetPair; // TLorentz vector with all Jets Pairs
      vector<float> JetBinfo; // vector with b-tag info per each Jet of the Event
      
  
      for (UInt_t j = 0; j < JetPt->size(); j++)
        {
  	if ( (fabs ( JetEta->at(j) ) > 2.4) || (JetPt->at(j) < 30 ) ) continue; // pt cut 20GeV from ntuplizer reduced to 30
  	  
  	TLorentzVector temp;
  	temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
  	JetVec.push_back(temp);
  	JetBinfo.push_back(JetBTagger->at(j));
        }
  
  
  
      // at least 2 jets in the acceptance
      if (JetVec.size() < 2) continue;   
          
      
  
      // get and save vector with max btagger value
      f_bdiscjet1signal_4mu = *max_element(JetBinfo.begin(), JetBinfo.end());
  
      // get and save btagger value of the second jet (the one with max pt)
      int d1_maxbtag = distance( JetBinfo.begin(), max_element(JetBinfo.begin(), JetBinfo.end()));
  
      float maxJetPt = -9999.;
      int d2_maxJetPt = -9999;
      for(UInt_t i=0; i<JetVec.size(); i++){
  
        if(i == d1_maxbtag) continue;
        float temp = JetVec.at(i).Pt();
        if (temp > maxJetPt){
          maxJetPt = temp;
          d2_maxJetPt = i;
        }
      }
      // save btagger value of the second jet (the one with max pt)
      f_bdiscjet2signal_4mu = JetBinfo.at(d2_maxJetPt);
  
  
      // build H->bb tlorentzvector
      TLorentzVector Hbb_Vec = JetVec.at(d1_maxbtag) + JetVec.at(d2_maxJetPt);
  
      // build H-H DeltaR
      float DeltaR = sqrt( (( ZZPhi - Hbb_Vec.Phi() )*( ZZPhi - Hbb_Vec.Phi() )) + (( ZZEta - Hbb_Vec.Eta() )*( ZZEta - Hbb_Vec.Eta() )) );
  
      // save DeltaR
      f_deltarsignal_4mu = DeltaR;
    


      // build HZZ tlorentzvector
      TLorentzVector HZZ_Vec;
      HZZ_Vec.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, ZZMass);
    
      // build HH tlorentzvector
      TLorentzVector HH_Vec = HZZ_Vec + Hbb_Vec;


      // save bb and HH masses
      f_bbmass = Hbb_Vec.M();
      f_HHmass = HH_Vec.M();

    } // end if JETSELECTION





    // fill k factors and event weights
    Float_t kfactor = 1.;
    // qqZZ sample
    if(inFile.Contains("ZZTo4l")) { kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; }
    //ggZZ samples                       
    else if(inFile.Contains("ggTo")) { kfactor = KFactor_QCD_ggZZ_Nominal; }

    Double_t eventWeight = 1.;
    if(!isDATA && !isZX) eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight;
    if(isZX) eventWeight = weight; //ZX weight




    // save ZZ mass
    f_ZZmass = ZZMass;

    // save lepton pt
    f_lept1_ptsignal_4mu = LepPt->at(0);
    f_lept2_ptsignal_4mu = LepPt->at(1);
    f_lept3_ptsignal_4mu = LepPt->at(2);
    f_lept4_ptsignal_4mu = LepPt->at(3);

    // save MET
    f_METsignal_4mu = PFMET;    
 

    // save event weight
    f_weightsignal_4mu = eventWeight; 

    yield += eventWeight; 


    

    tnew->Fill();
  }

  cout<<"yield "<<yield<<endl;
  
  f->cd();
  tnew->Write();
  f->Close();

}


void prepareNtupleMVA()
{

  float lumi = 59.74; //fb-1

  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/allsamples/";
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
                             //"Zzto4lamcatnlo",
                             // "DY3JetsToLL_M50",
                             // "DY2JetsToLL_M50" 
                             };


  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<< "number of input files: " << nInputFiles<<endl;


  string outputFilePath = "191214_mvaNtuples";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir


  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doNtuplesForMVA(inputFilePath + inputFileName[i] + "/ZZXAnalysis.root", outputFilePath + "/reduced_" + inputFileName[i] + ".root", lumi);
  }


}
