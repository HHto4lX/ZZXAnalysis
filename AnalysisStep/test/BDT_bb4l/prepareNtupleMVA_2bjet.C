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
#define MERGE2E2MU 1

enum FinalState {fs_4mu=0, fs_4e=1, fs_2e2mu=2, fs_2mu2e=3};  // 4mu, 4e, 2e2mu, 2mu2e
const int nFinalState = 4;
string FinalState[nFinalState+1] = {"4mu", "4e","2e2mu","2mu2e", "4L"};



void doNtuplesForMVA(TString inFile, TString outFile, float lumi)
{

  bool isDATA = false;
  bool isZX   = false;
  if ( inFile.Contains("AllData") ) isDATA = true;
  if ( inFile.Contains("ZXbkg") )   isZX   = true;
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
  Float_t L1prefiringWeight;

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

  vector<Float_t> *JetPt_JESUp = 0;
  vector<Float_t> *JetPt_JESDown = 0;
  vector<Float_t> *JetPt_JERUp = 0;
  vector<Float_t> *JetPt_JERDown = 0;

  Float_t yield_4e = 0.;
  Float_t yield_4mu = 0.;
  Float_t yield_2e2mu = 0.;
  Float_t yield_4l = 0.;


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
    if( !isDATA) inputTree->SetBranchAddress("L1prefiringWeight", &L1prefiringWeight);
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
  inputTree->SetBranchAddress("JetPt_JESUp", &JetPt_JESUp);
  inputTree->SetBranchAddress("JetPt_JESDown", &JetPt_JESDown);
  inputTree->SetBranchAddress("JetPt_JERUp", &JetPt_JERUp);
  inputTree->SetBranchAddress("JetPt_JERDown", &JetPt_JERDown);
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
  float f_lept1_ptsignal  = -9999.;
  float f_lept2_ptsignal  = -9999.;
  float f_lept3_ptsignal  = -9999.;
  float f_lept4_ptsignal  = -9999.;
  float f_massjetjet      = -9999.;
  float f_ptjet1          = -9999.;
  float f_ptjet2          = -9999.;
  float f_bdiscjet1signal = -9999.;
  float f_bdiscjet2signal = -9999.;
  float f_deltarsignal    = -9999.;
  float f_METsignal       = -9999.;
  float f_weightsignal    = -9999.;

  float f_Z1Mass = -9999.;
  float f_Z2Mass = -9999.;
  float f_ZZmass = -9999.;
  float f_bbmass = -9999.;
  float f_HHmass = -9999.;


  TFile *f = new TFile(outFile,"recreate");
  TTree *tnew = new TTree("reducedTree","");

  tnew->Branch("f_lept1_pt",    &f_lept1_ptsignal);
  tnew->Branch("f_lept2_pt",    &f_lept2_ptsignal);
  tnew->Branch("f_lept3_pt",    &f_lept3_ptsignal);
  tnew->Branch("f_lept4_pt",    &f_lept4_ptsignal);
  tnew->Branch("f_massjetjet",  &f_massjetjet);
  tnew->Branch("f_ptjet1",      &f_ptjet1);
  tnew->Branch("f_ptjet2",      &f_ptjet2);
  tnew->Branch("f_bdiscjet1",   &f_bdiscjet1signal);
  tnew->Branch("f_bdiscjet2",   &f_bdiscjet2signal);
  tnew->Branch("f_deltar_norm", &f_deltarsignal); 
  tnew->Branch("f_MET_norm",    &f_METsignal); 
  tnew->Branch("f_weight",      &f_weightsignal); 
  tnew->Branch("f_Z1Mass",      &f_Z1Mass); 
  tnew->Branch("f_Z2Mass",      &f_Z2Mass); 
  tnew->Branch("f_ZZmass",      &f_ZZmass); 
  tnew->Branch("f_bbmass",      &f_bbmass); 
  tnew->Branch("f_HHmass",      &f_HHmass); 



  int currentFinalState;


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



    // select final state
    currentFinalState = -1;
    if (!isZX){
	if (Z1Flav == -121){
	  if (Z2Flav == -121) {currentFinalState = fs_4e;}
          else if ( Z2Flav == -169) {currentFinalState = fs_2e2mu;}
          else { cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z2Flav = " << Z2Flav << endl;}
        }
	else if (Z1Flav == -169){
	  if (Z2Flav == -121) { currentFinalState = fs_2mu2e; }
	  else if (Z2Flav == -169) { currentFinalState = fs_4mu; }
          else { cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z2Flav = " << Z2Flav << endl; }
	}
	else{ cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z1Flav = " << Z1Flav << endl; }

	if(MERGE2E2MU && ( currentFinalState == fs_2mu2e )) { currentFinalState = fs_2e2mu; }
    }

    else if (isZX){
	if (Z1Flav == -121){
	  if (Z2Flav == +121) { currentFinalState = fs_4e; }
          else if (Z2Flav == +169) { currentFinalState = fs_2e2mu; }
          else { cerr << "error, Z2Flav: " << endl; }
        }
	else if (Z1Flav == -169){
	  if (Z2Flav == +121) { currentFinalState = fs_2mu2e; }
	  else if (Z2Flav == +169) { currentFinalState = fs_4mu; }
	  else { cerr << "error, Z2Flav: " << endl; }
	}
	else{ cerr << "error, Z2Flav: " << endl; }

	if(MERGE2E2MU && ( currentFinalState == fs_2mu2e )) currentFinalState = fs_2e2mu;
      }



    // fill k factors and event weights
    Float_t kfactor = 1.;
    // qqZZ sample
    if(inFile.Contains("ZZTo4l")) { kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; }
    //ggZZ samples                       
    else if(inFile.Contains("ggTo")) { kfactor = KFactor_QCD_ggZZ_Nominal; }

    Double_t eventWeight = 1.;
    if(!isDATA && !isZX) eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight * L1prefiringWeight;
    if(isZX) eventWeight = weight; //ZX weight



    //save only events for 1 final state at the time
    if(currentFinalState != fs_4mu)   continue;  // save 4mu only
    //    if(currentFinalState != fs_4e)    continue;  // save 4e only
    //    if(currentFinalState != fs_2e2mu) continue;  // save 2e2mu only
    //    cout<<currentFinalState<<endl;



    // //***SYNC*** out on file for sync1
    // ofstream fout1;
    // fout1.open("sync1_ggH_4mu_4lsel.txt",ios::app);
    // fout1<<nRun<<":"<<nLumi<<":"<<nEvent<<":"<<Z1Mass<<":"<<Z2Mass<<":"<<ZZMass<<":"<<JetPt->size()<<":"<<eventWeight<<endl;
    // fout1.close();

    


    // // mass cut: signal region
    if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV
    // mass cut: side bands
    //    if(ZZMass >= 115 && ZZMass <= 135) continue; // ZZMass < 115 or ZZMass > 135 GeV

 



    // //***SYNC*** sync, at least 1 jet
    // if (JetPt->size() >= 1){  
    //   // out on file for sync2
    //   ofstream fout2;
    //   fout2.open("sync2_ggH_4mu_4lsel1jet.txt",ios::app);
    //   fout2<<nRun<<":"<<nLumi<<":"<<nEvent<<":"<<Z1Mass<<":"<<Z2Mass<<":"<<ZZMass<<":"<<JetPt->size()<<":"<<eventWeight<<endl;
    //   fout2.close();
    // }

    //JETSELECTION---------------------------------------------------
    // at least 2 jets in the acceptance
    if (JetPt->size() < 2) continue; 
    // cout<<JetPt->size()<<" "<<JetPt_JESUp->size()<<" "<<JetPt_JESDown->size()<<" "<<JetPt_JERUp->size()<<" "<<JetPt_JERDown->size()<<endl;
    // for(int i=0; i<JetPt->size(); i++){
    //   cout<<JetPt->at(i)<<" "<<JetPt_JESUp->at(i)<<" "<<JetPt_JESDown->at(i)<<" "<<JetPt_JERUp->at(i)<<" "<<JetPt_JERDown->at(i)<<endl;
    // }



    // //***SYNC*** out on file for sync3
    // ofstream fout3;
    // fout3.open("sync3_ggH_4mu_4lsel2jet.txt",ios::app);
    // fout3<<nRun<<":"<<nLumi<<":"<<nEvent<<":"<<Z1Mass<<":"<<Z2Mass<<":"<<ZZMass<<":"<<JetPt->size()<<":"<<eventWeight<<endl;
    // fout3.close();

        

    // get and save vector with max btagger value
    f_bdiscjet1signal = *max_element(JetBTagger->begin(), JetBTagger->end());

    // get and save btagger value of the second jet (the one with second max btag)
    int d1_maxbtag = distance( JetBTagger->begin(), max_element(JetBTagger->begin(), JetBTagger->end()));

    float maxJetBtag = -9999.;
    int d2_maxbtag = -9999;
    for(UInt_t i=0; i<JetBTagger->size(); i++){

      if(i == d1_maxbtag) continue;
      float temp = JetBTagger->at(i);
      if (temp > maxJetBtag){
        maxJetBtag = temp;
        d2_maxbtag = i;
      }
    }
    // save btagger value of the second jet (the one with second max btag)
    f_bdiscjet2signal = JetBTagger->at(d2_maxbtag);

    // save jet pT
    f_ptjet1 = JetPt->at(d1_maxbtag);
    f_ptjet2 = JetPt->at(d2_maxbtag);


    // build 2 jets tlorentzvectors
    TLorentzVector tlzvec_j1_;
    tlzvec_j1_.SetPtEtaPhiM(JetPt->at(d1_maxbtag), JetEta->at(d1_maxbtag), JetPhi->at(d1_maxbtag), JetMass->at(d1_maxbtag));
    TLorentzVector tlzvec_j2_;
    tlzvec_j2_.SetPtEtaPhiM(JetPt->at(d2_maxbtag), JetEta->at(d2_maxbtag), JetPhi->at(d2_maxbtag), JetMass->at(d2_maxbtag));


    // build H->bb tlorentzvector
    TLorentzVector Hbb_Vec = tlzvec_j1_ + tlzvec_j2_;

    // build H-H DeltaR
    float DeltaPhi = ZZPhi - Hbb_Vec.Phi();
    if( fabs(DeltaPhi) > acos(-1) ) { DeltaPhi = (2*acos(-1)) - fabs(DeltaPhi); }
    float DeltaEta = ZZEta - Hbb_Vec.Eta();
    float DeltaR = sqrt( DeltaPhi*DeltaPhi + DeltaEta*DeltaEta );

    // save DeltaR
    f_deltarsignal = DeltaR;
  


    // build HZZ tlorentzvector
    TLorentzVector HZZ_Vec;
    HZZ_Vec.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, ZZMass);
  
    // build HH tlorentzvector
    TLorentzVector HH_Vec = HZZ_Vec + Hbb_Vec;


    // save bb and HH masses
    f_bbmass = Hbb_Vec.M();
    f_HHmass = HH_Vec.M();

    // save jet jet inv mass
    f_massjetjet = Hbb_Vec.M();

    //JETSELECTION---------------------------------------------------



    // save Z1 and Z2 mass
    f_Z1Mass = Z1Mass;
    f_Z2Mass = Z2Mass;

    // save ZZ mass
    f_ZZmass = ZZMass;

    // save lepton pt
    f_lept1_ptsignal = LepPt->at(0);
    f_lept2_ptsignal = LepPt->at(1);
    f_lept3_ptsignal = LepPt->at(2);
    f_lept4_ptsignal = LepPt->at(3);

    // save MET
    f_METsignal = PFMET;    
 

    // save event weight
    f_weightsignal = eventWeight; 



    // yields
    if(currentFinalState == fs_4e)    { yield_4e += eventWeight; }
    if(currentFinalState == fs_4mu)   { yield_4mu += eventWeight; }
    if(currentFinalState == fs_2e2mu) { yield_2e2mu += eventWeight; }
    yield_4l += eventWeight; 


    

    tnew->Fill();
  }


  cout<<"yields "<<endl;
  cout<<"4e    :"<<yield_4e<<endl;
  cout<<"4mu   :"<<yield_4mu<<endl;
  cout<<"2e2mu :"<<yield_2e2mu<<endl;
  cout<<"4l    :"<<yield_4l<<endl;
  cout<<"***************"<<endl;
  
  f->cd();
  tnew->Write();
  f->Close();

}


void prepareNtupleMVA_2bjet()
{

  //  float lumi = 59.7; //fb-1 2018
  float lumi = 41.5; //fb-1 2017


  //  TString inputFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/";
  TString inputFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2017/";
  TString inputFileName[] = {
    //       "AllData", 
    "HH4lbb_Angela",
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "bbH125",
    "ttH125",
    //    "ZZTo4lext2",
    "ZZTo4l",
    "ggTo4e_Contin_MCFM701",
    "ggTo4mu_Contin_MCFM701",
    "ggTo4tau_Contin_MCFM701",
    "ggTo2e2mu_Contin_MCFM701",
    "ggTo2e2tau_Contin_MCFM701",
    "ggTo2mu2tau_Contin_MCFM701",
    "TTZToLLNuNu_M10",
    "TTWJetsToLNu",
    "WWZ",
    "WZZ",
    "ZZZ",
    // "HToWW125_ggH",
    // "HToWW125_VBFH",
    // "HToWW125_HWplusJ",
    // "HToWW125_HWminusJ",
    // "HToWW125_HZJ",
    // "HToWW125_bbH",
    // "ZXbkg",
  };


  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<< "number of input files: " << nInputFiles<<endl;


  string outputFilePath = "mvaNtuples_2bjet_2017_SR_4ljjsel_fs4mu";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir


  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doNtuplesForMVA(inputFilePath + inputFileName[i] + "/ZZXAnalysis.root", outputFilePath + "/reduced_" + inputFileName[i] + ".root", lumi);
  }


  //***SYNC*** call function for sync only
  //  doNtuplesForMVA("../ZZXAnalysis.root", outputFilePath + "/reduced_sync.root", lumi);
}
