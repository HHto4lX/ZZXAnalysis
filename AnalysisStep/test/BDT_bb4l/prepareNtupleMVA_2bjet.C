
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

// includes for btag SF
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#include "evalEventSF.C"


using namespace std;

#define MERGE2E2MU 1

enum FinalState {fs_4mu=0, fs_4e=1, fs_2e2mu=2, fs_2mu2e=3};  // 4mu, 4e, 2e2mu, 2mu2e
const int nFinalState = 4;
string FinalState[nFinalState+1] = {"4mu", "4e","2e2mu","2mu2e", "4L"};



void doNtuplesForMVA(TString inFile, TString outFile, float lumi, TString syear, TString massRegion, TString selection, TString finalstate)
{

  bool isDATA = false;
  bool isZX   = false;
  if ( inFile.Contains("AllData") ) isDATA = true;
  if ( inFile.Contains("ZXbkg") )   isZX   = true;
  cout<<"isDATA "<<isDATA<<endl;
  cout<<"isZX "<<isZX<<endl;

  
  // arrays for BTagSF norm
  float sum_events = 0.;            
  float sum_BTagSF = 0.;            

  ///////////// SET UP B-TAG CALIBRATION ///////////////
    
  // set up calibration + reader
  cout << "Loading the .csv file..." << endl;
    
  std::string inputCSVfile = "";
  if(syear=="2016")     { inputCSVfile = "../../data/BTagging/DeepCSV_2016LegacySF_V1.csv"; }
  else if(syear=="2017"){ inputCSVfile = "../../data/BTagging/DeepCSV_94XSF_V5_B_F.csv"; }
  else if(syear=="2018"){ inputCSVfile = "../../data/BTagging/DeepCSV_102XSF_V1.csv"; }
  std::string measType            = "iterativefit";
  std::string sysType             = "central";
  std::string sysTypeJESUp        = "up_jes";
  std::string sysTypeJESDown      = "down_jes";
  std::string sysTypeHFUp         = "up_hf";
  std::string sysTypeHFDown       = "down_hf";
  std::string sysTypeLFUp         = "up_lf";
  std::string sysTypeLFDown       = "down_lf";
  std::string sysTypehfstats1Up   = "up_hfstats1";
  std::string sysTypehfstats1Down = "down_hfstats1";
  std::string sysTypehfstats2Up   = "up_hfstats2";
  std::string sysTypehfstats2Down = "down_hfstats2";
  std::string sysTypelfstats1Up   = "up_lfstats1";
  std::string sysTypelfstats1Down = "down_lfstats1";
  std::string sysTypelfstats2Up   = "up_lfstats2";
  std::string sysTypelfstats2Down = "down_lfstats2";
  std::string sysTypecfErr1Up     = "up_cferr1";
  std::string sysTypecfErr1Down   = "down_cferr1";
  std::string sysTypecfErr2Up     = "up_cferr2";
  std::string sysTypecfErr2Down   = "down_cferr2";
    
  BTagCalibration calib("csvv2", inputCSVfile);
  //nominal
  BTagCalibrationReader CSVreader(BTagEntry::OP_RESHAPING, sysType);       
  CSVreader.load(calib, BTagEntry::FLAV_B, measType);
  CSVreader.load(calib, BTagEntry::FLAV_C, measType);
  CSVreader.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_jes shift
  BTagCalibrationReader CSVreaderJESUp(BTagEntry::OP_RESHAPING, sysTypeJESUp);       
  CSVreaderJESUp.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderJESUp.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderJESUp.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_jes shift
  BTagCalibrationReader CSVreaderJESDown(BTagEntry::OP_RESHAPING, sysTypeJESDown);       
  CSVreaderJESDown.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderJESDown.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderJESDown.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_hf shift
  BTagCalibrationReader CSVreaderHFUp(BTagEntry::OP_RESHAPING, sysTypeHFUp);       
  CSVreaderHFUp.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderHFUp.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderHFUp.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_hf shift
  BTagCalibrationReader CSVreaderHFDown(BTagEntry::OP_RESHAPING, sysTypeHFDown);       
  CSVreaderHFDown.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderHFDown.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderHFDown.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_lf shift
  BTagCalibrationReader CSVreaderLFUp(BTagEntry::OP_RESHAPING, sysTypeLFUp);       
  CSVreaderLFUp.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderLFUp.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderLFUp.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_lf shift
  BTagCalibrationReader CSVreaderLFDown(BTagEntry::OP_RESHAPING, sysTypeLFDown);       
  CSVreaderLFDown.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderLFDown.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderLFDown.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_hfstats1 shift
  BTagCalibrationReader CSVreaderhfstats1Up(BTagEntry::OP_RESHAPING, sysTypehfstats1Up);       
  CSVreaderhfstats1Up.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderhfstats1Up.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderhfstats1Up.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_hfstats1 shift
  BTagCalibrationReader CSVreaderhfstats1Down(BTagEntry::OP_RESHAPING, sysTypehfstats1Down);       
  CSVreaderhfstats1Down.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderhfstats1Down.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderhfstats1Down.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_lfstats2 shift
  BTagCalibrationReader CSVreaderhfstats2Up(BTagEntry::OP_RESHAPING, sysTypehfstats2Up);       
  CSVreaderhfstats2Up.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderhfstats2Up.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderhfstats2Up.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_lfstats2 shift
  BTagCalibrationReader CSVreaderhfstats2Down(BTagEntry::OP_RESHAPING, sysTypehfstats2Down);       
  CSVreaderhfstats2Down.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderhfstats2Down.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderhfstats2Down.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_lfstats1 shift
  BTagCalibrationReader CSVreaderlfstats1Up(BTagEntry::OP_RESHAPING, sysTypelfstats1Up);       
  CSVreaderlfstats1Up.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderlfstats1Up.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderlfstats1Up.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_lfstats1 shift
  BTagCalibrationReader CSVreaderlfstats1Down(BTagEntry::OP_RESHAPING, sysTypelfstats1Down);       
  CSVreaderlfstats1Down.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderlfstats1Down.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderlfstats1Down.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_lfstats2 shift
  BTagCalibrationReader CSVreaderlfstats2Up(BTagEntry::OP_RESHAPING, sysTypelfstats2Up);       
  CSVreaderlfstats2Up.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderlfstats2Up.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderlfstats2Up.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_lfstats2 shift
  BTagCalibrationReader CSVreaderlfstats2Down(BTagEntry::OP_RESHAPING, sysTypelfstats2Down);       
  CSVreaderlfstats2Down.load(calib, BTagEntry::FLAV_B, measType);
  CSVreaderlfstats2Down.load(calib, BTagEntry::FLAV_C, measType);
  CSVreaderlfstats2Down.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_cferr1 shift
  BTagCalibrationReader CSVreadercfErr1Up(BTagEntry::OP_RESHAPING, sysTypecfErr1Up);       
  CSVreadercfErr1Up.load(calib, BTagEntry::FLAV_B, measType);
  CSVreadercfErr1Up.load(calib, BTagEntry::FLAV_C, measType);
  CSVreadercfErr1Up.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_cferr1 shift
  BTagCalibrationReader CSVreadercfErr1Down(BTagEntry::OP_RESHAPING, sysTypecfErr1Down);       
  CSVreadercfErr1Down.load(calib, BTagEntry::FLAV_B, measType);
  CSVreadercfErr1Down.load(calib, BTagEntry::FLAV_C, measType);
  CSVreadercfErr1Down.load(calib, BTagEntry::FLAV_UDSG, measType);
  //up_cferr2 shift
  BTagCalibrationReader CSVreadercfErr2Up(BTagEntry::OP_RESHAPING, sysTypecfErr2Up);       
  CSVreadercfErr2Up.load(calib, BTagEntry::FLAV_B, measType);
  CSVreadercfErr2Up.load(calib, BTagEntry::FLAV_C, measType);
  CSVreadercfErr2Up.load(calib, BTagEntry::FLAV_UDSG, measType);
  //down_cferr2 shift
  BTagCalibrationReader CSVreadercfErr2Down(BTagEntry::OP_RESHAPING, sysTypecfErr2Down);       
  CSVreadercfErr2Down.load(calib, BTagEntry::FLAV_B, measType);
  CSVreadercfErr2Down.load(calib, BTagEntry::FLAV_C, measType);
  CSVreadercfErr2Down.load(calib, BTagEntry::FLAV_UDSG, measType);
    
  cout << "Input CSV weight file = " << inputCSVfile << "; measurementType = " << measType << ";" << endl;

  ///////////////////////////

  

  // input file branches
  TFile* inputFile;
  TTree* inputTree;
  TH1F* hCounters;
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
  //  vector<Float_t> *LepPhi = 0;
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
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetHadronFlavour = 0;
  Float_t PFMET;

  Float_t LHEweight_QCDscale_muR1_muF1;
  Float_t LHEweight_QCDscale_muR1_muF2;
  Float_t LHEweight_QCDscale_muR1_muF0p5;
  Float_t LHEweight_QCDscale_muR2_muF1;
  Float_t LHEweight_QCDscale_muR2_muF2;
  Float_t LHEweight_QCDscale_muR2_muF0p5;
  Float_t LHEweight_QCDscale_muR0p5_muF1;
  Float_t LHEweight_QCDscale_muR0p5_muF2;
  Float_t LHEweight_QCDscale_muR0p5_muF0p5;
  Float_t LHEweight_PDFVariation_Up;
  Float_t LHEweight_PDFVariation_Dn;

  // vector<Float_t> *JetPt_JESUp = 0;
  // vector<Float_t> *JetPt_JESDown = 0;
  // vector<Float_t> *JetPt_JERUp = 0;
  // vector<Float_t> *JetPt_JERDown = 0;

  Float_t yield_4e = 0.;
  Float_t yield_4mu = 0.;
  Float_t yield_2e2mu = 0.;
  Float_t yield_4l = 0.;


  inputFile =  TFile::Open( inFile );

  if(!isZX){
    hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
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
  //  inputTree->SetBranchAddress("LepPhi", &LepPhi);
  inputTree->SetBranchAddress("ZZMass", &ZZMass);  
  inputTree->SetBranchAddress("Z1Flav", &Z1Flav);
  inputTree->SetBranchAddress("Z2Flav", &Z2Flav);
  inputTree->SetBranchAddress("Z1Mass", &Z1Mass);
  inputTree->SetBranchAddress("Z2Mass", &Z2Mass);
  inputTree->SetBranchAddress("ZZPt", &ZZPt);
  inputTree->SetBranchAddress("ZZEta", &ZZEta);
  inputTree->SetBranchAddress("ZZPhi", &ZZPhi);
  inputTree->SetBranchAddress("JetPt", &JetPt);
  // inputTree->SetBranchAddress("JetPt_JESUp", &JetPt_JESUp);
  // inputTree->SetBranchAddress("JetPt_JESDown", &JetPt_JESDown);
  // inputTree->SetBranchAddress("JetPt_JERUp", &JetPt_JERUp);
  // inputTree->SetBranchAddress("JetPt_JERDown", &JetPt_JERDown);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetMass", &JetMass);
  inputTree->SetBranchAddress("JetPhi",  &JetPhi);
  inputTree->SetBranchAddress("JetBTagger",  &JetBTagger);
  inputTree->SetBranchAddress("JetHadronFlavour",  &JetHadronFlavour);
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
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF1",     &LHEweight_QCDscale_muR1_muF1);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF2",     &LHEweight_QCDscale_muR1_muF2);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF0p5",   &LHEweight_QCDscale_muR1_muF0p5);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF1",     &LHEweight_QCDscale_muR2_muF1);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF2",     &LHEweight_QCDscale_muR2_muF2);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF0p5",   &LHEweight_QCDscale_muR2_muF0p5);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF1",   &LHEweight_QCDscale_muR0p5_muF1);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF2",   &LHEweight_QCDscale_muR0p5_muF2);
  inputTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF0p5", &LHEweight_QCDscale_muR0p5_muF0p5);
  inputTree->SetBranchAddress("LHEweight_PDFVariation_Up",        &LHEweight_PDFVariation_Up);
  inputTree->SetBranchAddress("LHEweight_PDFVariation_Dn",        &LHEweight_PDFVariation_Dn);
  


  
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
  // branches for tests
  float f_lept1_eta  = -9999.;
  float f_lept2_eta  = -9999.;
  float f_lept3_eta  = -9999.;
  float f_lept4_eta  = -9999.;
  float f_lept1_phi  = -9999.;
  float f_lept2_phi  = -9999.;
  float f_lept3_phi  = -9999.;
  float f_lept4_phi  = -9999.;
  float f_etajet1    = -9999.;
  float f_etajet2    = -9999.;
  float f_phijet1    = -9999.;
  float f_phijet2    = -9999.;
  float f_deltaPhiHH = -9999.;
  // branches for weights
  float f_weightsignal_nominal     = -9999.;
  float f_weightsignal_JESUp       = -9999.;
  float f_weightsignal_JESDown     = -9999.;
  float f_weightsignal_HFUp        = -9999.;
  float f_weightsignal_HFDown      = -9999.;
  float f_weightsignal_LFUp        = -9999.;
  float f_weightsignal_LFDown      = -9999.;
  float f_weightsignal_hfstats1Up  = -9999.;
  float f_weightsignal_hfstats1Down= -9999.;
  float f_weightsignal_hfstats2Up  = -9999.;
  float f_weightsignal_hfstats2Down= -9999.;
  float f_weightsignal_lfstats1Up  = -9999.;
  float f_weightsignal_lfstats1Down= -9999.;
  float f_weightsignal_lfstats2Up  = -9999.;
  float f_weightsignal_lfstats2Down= -9999.;
  float f_weightsignal_cfErr1Up    = -9999.;
  float f_weightsignal_cfErr1Down  = -9999.;
  float f_weightsignal_cfErr2Up    = -9999.;
  float f_weightsignal_cfErr2Down  = -9999.;
  // branches for QCD and PDF scale
  Float_t f_LHEweight_QCDscale_muR1_muF1     = -9999.;
  Float_t f_LHEweight_QCDscale_muR1_muF2     = -9999.;
  Float_t f_LHEweight_QCDscale_muR1_muF0p5   = -9999.;
  Float_t f_LHEweight_QCDscale_muR2_muF1     = -9999.;
  Float_t f_LHEweight_QCDscale_muR2_muF2     = -9999.;
  Float_t f_LHEweight_QCDscale_muR2_muF0p5   = -9999.;
  Float_t f_LHEweight_QCDscale_muR0p5_muF1   = -9999.;
  Float_t f_LHEweight_QCDscale_muR0p5_muF2   = -9999.;
  Float_t f_LHEweight_QCDscale_muR0p5_muF0p5 = -9999.;
  Float_t f_LHEweight_PDFVariation_Up        = -9999.;
  Float_t f_LHEweight_PDFVariation_Dn        = -9999.;
  // branches in piu'
  float f_Z1Mass = -9999.;
  float f_Z2Mass = -9999.;
  float f_ZZmass = -9999.;
  float f_bbmass = -9999.;
  float f_HHmass = -9999.;
  float f_HadronFlavourjet1 = -9999.;
  float f_HadronFlavourjet2 = -9999.;


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
  tnew->Branch("f_Z1Mass",      &f_Z1Mass); 
  tnew->Branch("f_Z2Mass",      &f_Z2Mass); 
  tnew->Branch("f_ZZmass",      &f_ZZmass); 
  tnew->Branch("f_bbmass",      &f_bbmass); 
  tnew->Branch("f_HHmass",      &f_HHmass); 
  tnew->Branch("f_lept1_eta",   &f_lept1_eta);
  tnew->Branch("f_lept2_eta",   &f_lept2_eta);
  tnew->Branch("f_lept3_eta",   &f_lept3_eta);
  tnew->Branch("f_lept4_eta",   &f_lept4_eta);
  tnew->Branch("f_lept1_phi",   &f_lept1_phi);
  tnew->Branch("f_lept2_phi",   &f_lept2_phi);
  tnew->Branch("f_lept3_phi",   &f_lept3_phi);
  tnew->Branch("f_lept4_phi",   &f_lept4_phi);
  tnew->Branch("f_etajet1",     &f_etajet1);
  tnew->Branch("f_etajet2",     &f_etajet2);
  tnew->Branch("f_phijet1",     &f_phijet1);
  tnew->Branch("f_phijet2",     &f_phijet2);
  tnew->Branch("f_deltaPhiHH",  &f_deltaPhiHH);
  tnew->Branch("f_HadronFlavourjet1", &f_HadronFlavourjet1);
  tnew->Branch("f_HadronFlavourjet2", &f_HadronFlavourjet2);
  tnew->Branch("f_weightsignal_nominal",      &f_weightsignal_nominal     ); 
  tnew->Branch("f_weightsignal_JESUp",        &f_weightsignal_JESUp       ); 
  tnew->Branch("f_weightsignal_JESDown",      &f_weightsignal_JESDown     ); 
  tnew->Branch("f_weightsignal_HFUp",         &f_weightsignal_HFUp        ); 
  tnew->Branch("f_weightsignal_HFDown",       &f_weightsignal_HFDown      ); 
  tnew->Branch("f_weightsignal_LFUp",         &f_weightsignal_LFUp        ); 
  tnew->Branch("f_weightsignal_LFDown",       &f_weightsignal_LFDown      ); 
  tnew->Branch("f_weightsignal_hfstats1Up",   &f_weightsignal_hfstats1Up  ); 
  tnew->Branch("f_weightsignal_hfstats1Down", &f_weightsignal_hfstats1Down); 
  tnew->Branch("f_weightsignal_hfstats2Up",   &f_weightsignal_hfstats2Up  ); 
  tnew->Branch("f_weightsignal_hfstats2Down", &f_weightsignal_hfstats2Down); 
  tnew->Branch("f_weightsignal_lfstats1Up",   &f_weightsignal_lfstats1Up  ); 
  tnew->Branch("f_weightsignal_lfstats1Down", &f_weightsignal_lfstats1Down); 
  tnew->Branch("f_weightsignal_lfstats2Up",   &f_weightsignal_lfstats2Up  ); 
  tnew->Branch("f_weightsignal_lfstats2Down", &f_weightsignal_lfstats2Down); 
  tnew->Branch("f_weightsignal_cfErr1Up",     &f_weightsignal_cfErr1Up    ); 
  tnew->Branch("f_weightsignal_cfErr1Down",   &f_weightsignal_cfErr1Down  ); 
  tnew->Branch("f_weightsignal_cfErr2Up",     &f_weightsignal_cfErr2Up    ); 
  tnew->Branch("f_weightsignal_cfErr2Down",   &f_weightsignal_cfErr2Down  ); 
  tnew->Branch("f_LHEweight_QCDscale_muR1_muF1",     &f_LHEweight_QCDscale_muR1_muF1);
  tnew->Branch("f_LHEweight_QCDscale_muR1_muF2",     &f_LHEweight_QCDscale_muR1_muF2);    
  tnew->Branch("f_LHEweight_QCDscale_muR1_muF0p5",   &f_LHEweight_QCDscale_muR1_muF0p5);  
  tnew->Branch("f_LHEweight_QCDscale_muR2_muF1",     &f_LHEweight_QCDscale_muR2_muF1);
  tnew->Branch("f_LHEweight_QCDscale_muR2_muF2",     &f_LHEweight_QCDscale_muR2_muF2);
  tnew->Branch("f_LHEweight_QCDscale_muR2_muF0p5",   &f_LHEweight_QCDscale_muR2_muF0p5);
  tnew->Branch("f_LHEweight_QCDscale_muR0p5_muF1",   &f_LHEweight_QCDscale_muR0p5_muF1);
  tnew->Branch("f_LHEweight_QCDscale_muR0p5_muF2",   &f_LHEweight_QCDscale_muR0p5_muF2);
  tnew->Branch("f_LHEweight_QCDscale_muR0p5_muF0p5", &f_LHEweight_QCDscale_muR0p5_muF0p5);
  tnew->Branch("f_LHEweight_PDFVariation_Up",        &f_LHEweight_PDFVariation_Up);
  tnew->Branch("f_LHEweight_PDFVariation_Dn",        &f_LHEweight_PDFVariation_Dn);


  int currentFinalState;

  
  // --------------------------------------------------------
  // --- first loop over input tree to get norm for BtagSF
  if(inFile.Contains("AllData") || inFile.Contains("ZXbkg")){
      sum_events = 1.;      
      sum_BTagSF = 1.;
  }
  else{
    Long64_t entries1 = inputTree->GetEntries();    
    cout<<"First loop over input files to get norm for BTagSF ..."<<endl;
    cout<<"Processing file: "<< inFile << " (" << entries1 <<" entries) ..."<< endl;    
    for(Long64_t z=0; z<entries1; z++){

      inputTree->GetEntry(z);
  
      if( !(ZZsel>=0) ) continue;

      // 4l selection
      if(LepEta->size() != 4){
        cout << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
        continue;
      }


      // compute SF
      double * scaleFactors;
      scaleFactors = evalEventSF( int(JetPt->size()), JetHadronFlavour, JetEta, JetPt, JetBTagger, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down );

      // total counters for BTagSF norm 
      sum_events += 1.; 
      sum_BTagSF += scaleFactors[0]; 
    
    } // end first loop over entries 
  }// end else

  cout<<inFile<<" "<<sum_events<<" "<<sum_BTagSF<<endl;

  // --- control for norm
  if( sum_events == 0. || std::isnan(sum_events) ){ sum_events = 1.; }
  if( sum_BTagSF == 0. || std::isnan(sum_BTagSF) ){ sum_BTagSF = 1.; }

  cout<<inFile<<" "<<sum_events<<" "<<sum_BTagSF<<endl;
  // --------------------------------------------------------
  



  // --- second loop over input tree
  Long64_t entries2 = inputTree->GetEntries();
  cout<<"Second loop over input files to do all the rest ..."<<endl;
  cout<<"Processing file: "<< inFile << " (" << entries2 <<" entries) ..."<< endl;   
  for (Long64_t entry = 0; entry < entries2; entry++)
  {  
    inputTree->GetEntry(entry);
    

    if( !(ZZsel>=0) ) continue;

    // 4l selection
    if(LepEta->size() != 4)
      {
	cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
      }


    // --- fill k factors 
    Float_t kfactor = 1.;
    if(inFile.Contains("ZZTo4l")) { kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; } // qqZZ sample
    else if(inFile.Contains("ggTo")) { kfactor = KFactor_QCD_ggZZ_Nominal; }  //ggZZ samples 

    // compute SF
    double * scaleFactors;
    if(!isDATA && !isZX){
      scaleFactors = evalEventSF( int(JetPt->size()), JetHadronFlavour, JetEta, JetPt, JetBTagger, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down );
    }

    // --- event weights
    Double_t eventWeight = 1.;
    if(!isDATA && !isZX){ eventWeight = partialSampleWeight *xsec *kfactor *overallEventWeight *L1prefiringWeight; }
    if(isZX) eventWeight = weight; //ZX weight





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




    // --- save only events for 1 final state at the time
    if(finalstate == "fs4mu"){ 
      if(currentFinalState != fs_4mu){ continue; } // save 4mu only
    } 
    else if(finalstate == "fs4e"){
      if(currentFinalState != fs_4e){ continue; } // save 4e only
    }
    else if(finalstate == "fs2e2mu"){
      if(currentFinalState != fs_2e2mu) continue;  // save 2e2mu only
    }
    else cerr<<"wrong final state"<<endl;


    

    // --- mass region
    if(massRegion == "SR"){
      // mass cut: signal region
      if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV
    }
    else if(massRegion == "sidebands"){
      // mass cut: side bands
      if(ZZMass >= 115 && ZZMass <= 135) continue; // ZZMass < 115 or ZZMass > 135 GeV
    }
    else if(massRegion == "fullmass"){
      // full mass range
    }
    else cerr<<"wrong mass region!"<<endl;
 



    //JETSELECTION---------------------------------------------------
    // at least 2 jets in the acceptance
    if (JetPt->size() < 2) continue; 
    // cout<<JetPt->size()<<" "<<JetPt_JESUp->size()<<" "<<JetPt_JESDown->size()<<" "<<JetPt_JERUp->size()<<" "<<JetPt_JERDown->size()<<endl;
    // for(int i=0; i<JetPt->size(); i++){
    //   cout<<JetPt->at(i)<<" "<<JetPt_JESUp->at(i)<<" "<<JetPt_JESDown->at(i)<<" "<<JetPt_JERUp->at(i)<<" "<<JetPt_JERDown->at(i)<<endl;
    // }


        

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

    // save jet hadronflavour
    f_HadronFlavourjet1 = JetHadronFlavour->at(d1_maxbtag);
    f_HadronFlavourjet2 = JetHadronFlavour->at(d2_maxbtag);


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

    // JETSELECTION---------------------------------------------------



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
 
    // save branches for tests
    f_lept1_eta  = LepEta->at(0);
    f_lept2_eta  = LepEta->at(1);
    f_lept3_eta  = LepEta->at(2);
    f_lept4_eta  = LepEta->at(3);
    // f_lept1_phi  = LepPhi->at(0);
    // f_lept2_phi  = LepPhi->at(1); 
    // f_lept3_phi  = LepPhi->at(2); 
    // f_lept4_phi  = LepPhi->at(3); 
    f_etajet1    = JetEta->at(d1_maxbtag);
    f_etajet2    = JetEta->at(d2_maxbtag);
    f_phijet1    = JetPhi->at(d1_maxbtag);
    f_phijet2    = JetPhi->at(d2_maxbtag);
    f_deltaPhiHH = DeltaPhi;



    // save event weight
    if(isDATA || isZX){ 
      f_weightsignal_nominal = eventWeight;
    }
    else{
      f_weightsignal_nominal = eventWeight * scaleFactors[0] * sum_events/sum_BTagSF;
    }
 
    // save event weight for BTagSF syst
    f_weightsignal_JESUp        = eventWeight * scaleFactors[1]  * sum_events/sum_BTagSF;
    f_weightsignal_JESDown      = eventWeight * scaleFactors[2]  * sum_events/sum_BTagSF;
    f_weightsignal_HFUp         = eventWeight * scaleFactors[3]  * sum_events/sum_BTagSF;
    f_weightsignal_HFDown       = eventWeight * scaleFactors[4]  * sum_events/sum_BTagSF;
    f_weightsignal_LFUp         = eventWeight * scaleFactors[5]  * sum_events/sum_BTagSF;
    f_weightsignal_LFDown       = eventWeight * scaleFactors[6]  * sum_events/sum_BTagSF;
    f_weightsignal_hfstats1Up   = eventWeight * scaleFactors[7]  * sum_events/sum_BTagSF;
    f_weightsignal_hfstats1Down = eventWeight * scaleFactors[8]  * sum_events/sum_BTagSF;
    f_weightsignal_hfstats2Up   = eventWeight * scaleFactors[9]  * sum_events/sum_BTagSF;
    f_weightsignal_hfstats2Down = eventWeight * scaleFactors[10] * sum_events/sum_BTagSF;
    f_weightsignal_lfstats1Up   = eventWeight * scaleFactors[11] * sum_events/sum_BTagSF;
    f_weightsignal_lfstats1Down = eventWeight * scaleFactors[12] * sum_events/sum_BTagSF;
    f_weightsignal_lfstats2Up   = eventWeight * scaleFactors[13] * sum_events/sum_BTagSF;
    f_weightsignal_lfstats2Down = eventWeight * scaleFactors[14] * sum_events/sum_BTagSF;
    f_weightsignal_cfErr1Up     = eventWeight * scaleFactors[15] * sum_events/sum_BTagSF;
    f_weightsignal_cfErr1Down   = eventWeight * scaleFactors[16] * sum_events/sum_BTagSF;
    f_weightsignal_cfErr2Up     = eventWeight * scaleFactors[17] * sum_events/sum_BTagSF;
    f_weightsignal_cfErr2Down   = eventWeight * scaleFactors[18] * sum_events/sum_BTagSF;

    // save branches for QCD and PDF scale syst
    f_LHEweight_QCDscale_muR1_muF1     = LHEweight_QCDscale_muR1_muF1;    
    f_LHEweight_QCDscale_muR1_muF2     = LHEweight_QCDscale_muR1_muF2; 
    f_LHEweight_QCDscale_muR1_muF0p5   = LHEweight_QCDscale_muR1_muF0p5;
    f_LHEweight_QCDscale_muR2_muF1     = LHEweight_QCDscale_muR2_muF1;
    f_LHEweight_QCDscale_muR2_muF2     = LHEweight_QCDscale_muR2_muF2;
    f_LHEweight_QCDscale_muR2_muF0p5   = LHEweight_QCDscale_muR2_muF0p5;
    f_LHEweight_QCDscale_muR0p5_muF1   = LHEweight_QCDscale_muR0p5_muF1;
    f_LHEweight_QCDscale_muR0p5_muF2   = LHEweight_QCDscale_muR0p5_muF2;
    f_LHEweight_QCDscale_muR0p5_muF0p5 = LHEweight_QCDscale_muR0p5_muF0p5;
    f_LHEweight_PDFVariation_Up        = LHEweight_PDFVariation_Up;
    f_LHEweight_PDFVariation_Dn        = LHEweight_PDFVariation_Dn;



    // yields
    if(currentFinalState == fs_4e)    { yield_4e    += f_weightsignal_nominal; }
    if(currentFinalState == fs_4mu)   { yield_4mu   += f_weightsignal_nominal; }
    if(currentFinalState == fs_2e2mu) { yield_2e2mu += f_weightsignal_nominal; }
    yield_4l += f_weightsignal_nominal; 


    

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
  // --- mass region
  TString massRegion = "SR";
  //  TString massRegion = "sidebands";
  //  TString massRegion = "fullmass";

  // --- selection
  TString selection = "4ljjsel";
  //  TString selection = "4lsel";

  // --- finalstate
  TString finalstate = "fs4mu";
  //  TString finalstate = "fs4e";
  //  TString finalstate = "fs2e2mu";
  
  // -- year
  TString syear = "2016";
  //  TString syear = "2017";
  //  TString syear = "2018";

  // --- lumi
  float lumi = 35.8; //fb-1 2016
  //  float lumi = 41.5; //fb-1 2017
  //  float lumi = 59.7; //fb-1 2018
  cout<<lumi<<endl;


  TString inputFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2016/";
  //  TString inputFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2017/";
  //  TString inputFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/";
  TString inputFileName[] = {
    // "AllData", 
    // "HH4lbb_Angela",
    // "HH4lbb_Ilirjan",
    "ggH125",
    // "ggH125minlo",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "bbH125",
    "ttH125",
    // "ZZTo4lext2",
    "ZZTo4l",
    "ggTo4e_Contin_MCFM701",
    "ggTo4mu_Contin_MCFM701",
    "ggTo4tau_Contin_MCFM701",
    "ggTo2e2mu_Contin_MCFM701",
    "ggTo2e2tau_Contin_MCFM701",
    "ggTo2mu2tau_Contin_MCFM701",
    // "TTZToLLNuNu_M10",
    // "TTWJetsToLNu",
    // "WWZ",
    // "WZZ",
    // "ZZZ",
    // "ZXbkg_4ljjsel",
    // "DYJetsToLL_M50",
    // "TTTo2L2Nu",
  };


  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<< "number of input files: " << nInputFiles<<endl;

  

  TString outputFilePath = "mvaNtuples_2bjet_" + syear + "_" + massRegion + "_" + selection + "_" + finalstate;
  cout<<"Output dir: "<<outputFilePath<<endl;
  gSystem->Exec(("mkdir -p "+string(outputFilePath)).c_str()); // create output dir


  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doNtuplesForMVA(inputFilePath + inputFileName[i]+"/ZZXAnalysis.root",outputFilePath +"/reduced_"+ inputFileName[i]+".root",lumi,syear,massRegion,selection,finalstate);
  }


  //***SYNC*** call function for sync only
  //  doNtuplesForMVA("../ZZXAnalysis.root", outputFilePath + "/reduced_sync.root", lumi);
}
