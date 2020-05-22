// ***********************
// run with:
//
// root -l -b -q analysis_4lbb_2bjet.C++
// ***********************

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


// includes for btag SF
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#include "evalEventSF.C"


#include "CMS_lumi.C"



using namespace std;

#define REDOHISTOS 1

//******************
//int year = 2016;
//int year = 2017;
int year = 2018;
//******************



//*************************************************************************************
// FINAL STATES
enum FinalState {fs_4mu=0, fs_4e=1, fs_2e2mu=2};  // 4mu, 4e, 2e2mu
const int nFinalStates = 3;
TString sFinalState[nFinalStates+1] = {"4mu", "4e","2e2mu","4l"};
//*************************************************************************************
// PROCESSES
enum Process {Data=0, HH=1, ggH=2, VBF=3, VH=4, ttH=5, bbH=6, qqZZ=7, ggZZ=8, TTZ=9, TTW=10, VVV=11, ZXbkg=12}; 
const int nProcesses = 13;
TString sProcess[nProcesses] = {"Data", "HH", "ggH", "VBF", "VH", "ttH", "bbH", "qqZZ", "ggZZ", "TTZ", "TTW", "VVV", "ZXbkg"};
//*************************************************************************************




//*************************
//*** doHistos function ***
//*************************
void doHistos()
{

  TH1::SetDefaultSumw2(true);

  //---lumi and input path
  float lumi = 0.;
  TString inFilePath;
  TString inDataPath;
  TString sYear;
  Float_t rescale_ZX[nFinalStates];
  if(year==2016){
    lumi       = 35.8; //fb-1 2016
    sYear      = "2016";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2016/";
    inDataPath = "/eos/user/a/acappati/samples_HH4lbb/samples_2016/";
    rescale_ZX[fs_4mu]   = 0.79; //FIXME
    rescale_ZX[fs_4e]    = 1.40; //FIXME
    rescale_ZX[fs_2e2mu] = 2.64; //FIXME:da levare
  }
  else if(year==2017){
    lumi       = 41.5; //fb-1 2017
    sYear      = "2017";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2017/";
    inDataPath = "/eos/user/a/acappati/samples_HH4lbb/samples_2017/";
    rescale_ZX[fs_4mu]   = 1.48; //FIXME
    rescale_ZX[fs_4e]    = 0.52; //FIXME
    rescale_ZX[fs_2e2mu] = 2.00; //FIXME:da levare
  }
  else if(year==2018){
    lumi       = 59.7; //fb-1 2018
    sYear      = "2018";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/";
    inDataPath = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/";
    rescale_ZX[fs_4mu]   = 1.60; //FIXME
    rescale_ZX[fs_4e]    = 0.72; //FIXME
    rescale_ZX[fs_2e2mu] = 2.22; //FIXME:da levare
  }
  else{ 
    cout<<"wrong year selected!"<<endl;
  }
  cout<<"Year chosen: "<<year<<endl;


  //---datasets
  //  static int nDatasets = 22;
  TString datasets[] = {
    "AllData", 
    "HH4lbb_Angela",
    //"HH4lbb_Ilirjan",
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "bbH125",
    "ttH125",
    //"ZZTo4lamcatnlo",
    "ZZTo4lext2",
    //"ZZTo4l",
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
    //"ZXbkg_4ljjsel",
    "ZXbkg_4lsel",
  };
  static size_t nDatasets = sizeof(datasets)/sizeof(datasets[0]);
  cout<< "number of input files: " << nDatasets<<endl;

  for(int i=0; i<nDatasets; i++){
    cout<<datasets[i]<<" ";
  }
  cout<<endl;

  // arrays for BTagSF norm
  float sum_events[nDatasets];            
  float sum_BTagSF[nDatasets];
  for(int i=0; i<nDatasets; i++){
    sum_events[i] = 0.;           
    sum_BTagSF[i] = 0.;           
  }



  ///////////// SET UP B-TAG CALIBRATION ///////////////
    
  // set up calibration + reader
  cout << "Loading the .csv file..." << endl;
    
  std::string inputCSVfile = "";
  if(year==2016)     { inputCSVfile = "../../data/BTagging/DeepCSV_2016LegacySF_V1.csv"; }
  else if(year==2017){ inputCSVfile = "../../data/BTagging/DeepCSV_94XSF_V5_B_F.csv"; }
  else if(year==2018){ inputCSVfile = "../../data/BTagging/DeepCSV_102XSF_V1.csv"; }
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




  // --- input file branches
  TFile* inputFile[nDatasets];
  TTree* inputTree[nDatasets];
  TH1F* hCounters[nDatasets];
  Float_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];
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
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetHadronFlavour = 0;
  Float_t PFMET;

  // define yields histos and nEvents histos (no weight)
  TH1F* hYields_4lsel  [nProcesses][nFinalStates+1];
  TH1F* hYields_4ljjsel[nProcesses][nFinalStates+1];
  TH1F* hEvents_4lsel  [nProcesses][nFinalStates+1];
  TH1F* hEvents_4ljjsel[nProcesses][nFinalStates+1];
  TH1F* hYields_4ljjsel_sidebands[nProcesses][nFinalStates+1];
  TH1F* hYields_4ljjsel_sidebandSX[nProcesses][nFinalStates+1];
  TH1F* hYields_4ljjsel_sidebandDX[nProcesses][nFinalStates+1];
  TH1F* hYields_4ljjsel_sidebandSXDX[nProcesses][nFinalStates+1];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      hYields_4lsel  [pr][fs] = new TH1F("hYields_4lsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hYields_4lsel  [pr][fs]->Sumw2(true);
      hYields_4ljjsel[pr][fs] = new TH1F("hYields_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hYields_4ljjsel[pr][fs]->Sumw2(true);
      hEvents_4lsel  [pr][fs] = new TH1F("hEvents_4lsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hEvents_4lsel  [pr][fs]->Sumw2(true);
      hEvents_4ljjsel[pr][fs] = new TH1F("hEvents_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hEvents_4ljjsel[pr][fs]->Sumw2(true);
      hYields_4ljjsel_sidebands[pr][fs] = new TH1F("hYields_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hYields_4ljjsel_sidebands[pr][fs]->Sumw2(true);
      hYields_4ljjsel_sidebandSX[pr][fs] = new TH1F("hYields_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hYields_4ljjsel_sidebandSX[pr][fs]->Sumw2(true);
      hYields_4ljjsel_sidebandDX[pr][fs] = new TH1F("hYields_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hYields_4ljjsel_sidebandDX[pr][fs]->Sumw2(true);
      hYields_4ljjsel_sidebandSXDX[pr][fs] = new TH1F("hYields_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,"",1,0.,1.);
      hYields_4ljjsel_sidebandSXDX[pr][fs]->Sumw2(true);
    }
  }

  // define 1D histos
  // --- 4ljjsel
  // fullmass range
  TH1F* h1_m4l_4ljjsel  [nProcesses][nFinalStates+1]; 
  // sidebands
  TH1F* h1_pT4l_4ljjsel_sidebands    [nProcesses][nFinalStates+1];
  TH1F* h1_j1btag_4ljjsel_sidebands  [nProcesses][nFinalStates+1];
  TH1F* h1_j2btag_4ljjsel_sidebands  [nProcesses][nFinalStates+1];
  TH1F* h1_j1pT_4ljjsel_sidebands    [nProcesses][nFinalStates+1];
  TH1F* h1_j2pT_4ljjsel_sidebands    [nProcesses][nFinalStates+1];
  TH1F* h1_MET_4ljjsel_sidebands     [nProcesses][nFinalStates+1];
  TH1F* h1_DeltaRhh_4ljjsel_sidebands[nProcesses][nFinalStates+1];
  TH1F* h1_mbb_4ljjsel_sidebands     [nProcesses][nFinalStates+1];
  // sidebandSX: m4l in 95-115 GeV
  TH1F* h1_m4l_4ljjsel_sidebandSX     [nProcesses][nFinalStates+1]; 
  TH1F* h1_pT4l_4ljjsel_sidebandSX    [nProcesses][nFinalStates+1];
  TH1F* h1_j1btag_4ljjsel_sidebandSX  [nProcesses][nFinalStates+1];
  TH1F* h1_j2btag_4ljjsel_sidebandSX  [nProcesses][nFinalStates+1];
  TH1F* h1_j1pT_4ljjsel_sidebandSX    [nProcesses][nFinalStates+1];
  TH1F* h1_j2pT_4ljjsel_sidebandSX    [nProcesses][nFinalStates+1];
  TH1F* h1_MET_4ljjsel_sidebandSX     [nProcesses][nFinalStates+1];
  TH1F* h1_DeltaRhh_4ljjsel_sidebandSX[nProcesses][nFinalStates+1];
  TH1F* h1_mbb_4ljjsel_sidebandSX     [nProcesses][nFinalStates+1];
  // sidebandDX: m4l in 135-170 GeV
  TH1F* h1_m4l_4ljjsel_sidebandDX     [nProcesses][nFinalStates+1];
  TH1F* h1_pT4l_4ljjsel_sidebandDX    [nProcesses][nFinalStates+1];
  TH1F* h1_j1btag_4ljjsel_sidebandDX  [nProcesses][nFinalStates+1];
  TH1F* h1_j2btag_4ljjsel_sidebandDX  [nProcesses][nFinalStates+1];
  TH1F* h1_j1pT_4ljjsel_sidebandDX    [nProcesses][nFinalStates+1];
  TH1F* h1_j2pT_4ljjsel_sidebandDX    [nProcesses][nFinalStates+1];
  TH1F* h1_MET_4ljjsel_sidebandDX     [nProcesses][nFinalStates+1];
  TH1F* h1_DeltaRhh_4ljjsel_sidebandDX[nProcesses][nFinalStates+1];
  TH1F* h1_mbb_4ljjsel_sidebandDX     [nProcesses][nFinalStates+1];
  // sidebandSXDX: m4l in 95-115 & 135-170 GeV
  TH1F* h1_m4l_4ljjsel_sidebandSXDX     [nProcesses][nFinalStates+1];
  TH1F* h1_pT4l_4ljjsel_sidebandSXDX    [nProcesses][nFinalStates+1];
  TH1F* h1_j1btag_4ljjsel_sidebandSXDX  [nProcesses][nFinalStates+1];
  TH1F* h1_j2btag_4ljjsel_sidebandSXDX  [nProcesses][nFinalStates+1];
  TH1F* h1_j1pT_4ljjsel_sidebandSXDX    [nProcesses][nFinalStates+1];
  TH1F* h1_j2pT_4ljjsel_sidebandSXDX    [nProcesses][nFinalStates+1];
  TH1F* h1_MET_4ljjsel_sidebandSXDX     [nProcesses][nFinalStates+1];
  TH1F* h1_DeltaRhh_4ljjsel_sidebandSXDX[nProcesses][nFinalStates+1];
  TH1F* h1_mbb_4ljjsel_sidebandSXDX     [nProcesses][nFinalStates+1];
  // mass cut plots: (BDT input histos)
  TH1F* h1_pT4l_4ljjsel    [nProcesses][nFinalStates+1];
  TH1F* h1_j1btag_4ljjsel  [nProcesses][nFinalStates+1];
  TH1F* h1_j2btag_4ljjsel  [nProcesses][nFinalStates+1];
  TH1F* h1_j1pT_4ljjsel    [nProcesses][nFinalStates+1];
  TH1F* h1_j2pT_4ljjsel    [nProcesses][nFinalStates+1];
  TH1F* h1_MET_4ljjsel     [nProcesses][nFinalStates+1];
  TH1F* h1_DeltaRhh_4ljjsel[nProcesses][nFinalStates+1];
  TH1F* h1_mbb_4ljjsel     [nProcesses][nFinalStates+1];
  
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // --- 4ljjsel
      // fullmass range
      h1_m4l_4ljjsel  [pr][fs] = new TH1F("h1_m4l_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events/4 GeV", 33, 70., 202.);
      h1_m4l_4ljjsel  [pr][fs]->Sumw2(true);
      // sidebands
      h1_pT4l_4ljjsel_sidebands    [pr][fs] = new TH1F("h1_pT4l_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";4 leptons pT (GeV); Events/2 GeV", 50, 0., 100.);
      h1_pT4l_4ljjsel_sidebands    [pr][fs]->Sumw2(true);
      h1_j1btag_4ljjsel_sidebands  [pr][fs] = new TH1F("h1_j1btag_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j1btag_4ljjsel_sidebands  [pr][fs]->Sumw2(true);
      h1_j2btag_4ljjsel_sidebands  [pr][fs] = new TH1F("h1_j2btag_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j2btag_4ljjsel_sidebands  [pr][fs]->Sumw2(true);
      h1_j1pT_4ljjsel_sidebands    [pr][fs] = new TH1F("h1_j1pT_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j1pT_4ljjsel_sidebands    [pr][fs]->Sumw2(true);
      h1_j2pT_4ljjsel_sidebands    [pr][fs] = new TH1F("h1_j2pT_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j2pT_4ljjsel_sidebands    [pr][fs]->Sumw2(true);
      h1_MET_4ljjsel_sidebands     [pr][fs] = new TH1F("h1_MET_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";MET (GeV); Events/5 GeV", 40, 0., 200.);
      h1_MET_4ljjsel_sidebands     [pr][fs]->Sumw2(true);
      h1_DeltaRhh_4ljjsel_sidebands[pr][fs] = new TH1F("h1_DeltaRhh_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";Delta R H4l-Hbb; Events/0.4 ", 25, 0., 10.);
      h1_DeltaRhh_4ljjsel_sidebands[pr][fs]->Sumw2(true);
      h1_mbb_4ljjsel_sidebands     [pr][fs] = new TH1F("h1_mbb_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{jj} (GeV); Events/5 GeV", 40, 0., 200.);
      h1_mbb_4ljjsel_sidebands     [pr][fs]->Sumw2(true);

      // sidebandSX
      h1_m4l_4ljjsel_sidebandSX     [pr][fs] = new TH1F("h1_m4l_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events/4 GeV", 5, 95., 115.);
      h1_m4l_4ljjsel_sidebandSX     [pr][fs]->Sumw2(true);
      h1_pT4l_4ljjsel_sidebandSX    [pr][fs] = new TH1F("h1_pT4l_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";4 leptons pT (GeV); Events/2 GeV", 50, 0., 100.);
      h1_pT4l_4ljjsel_sidebandSX    [pr][fs]->Sumw2(true);
      h1_j1btag_4ljjsel_sidebandSX  [pr][fs] = new TH1F("h1_j1btag_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j1btag_4ljjsel_sidebandSX  [pr][fs]->Sumw2(true);
      h1_j2btag_4ljjsel_sidebandSX  [pr][fs] = new TH1F("h1_j2btag_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j2btag_4ljjsel_sidebandSX  [pr][fs]->Sumw2(true);
      h1_j1pT_4ljjsel_sidebandSX    [pr][fs] = new TH1F("h1_j1pT_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j1pT_4ljjsel_sidebandSX    [pr][fs]->Sumw2(true);
      h1_j2pT_4ljjsel_sidebandSX    [pr][fs] = new TH1F("h1_j2pT_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j2pT_4ljjsel_sidebandSX    [pr][fs]->Sumw2(true);
      h1_MET_4ljjsel_sidebandSX     [pr][fs] = new TH1F("h1_MET_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";MET (GeV); Events/5 GeV", 40, 0., 200.);
      h1_MET_4ljjsel_sidebandSX     [pr][fs]->Sumw2(true);
      h1_DeltaRhh_4ljjsel_sidebandSX[pr][fs] = new TH1F("h1_DeltaRhh_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";Delta R H4l-Hbb; Events/0.4 ", 25, 0., 10.);
      h1_DeltaRhh_4ljjsel_sidebandSX[pr][fs]->Sumw2(true);
      h1_mbb_4ljjsel_sidebandSX     [pr][fs] = new TH1F("h1_mbb_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{jj} (GeV); Events/5 GeV", 40, 0., 200.);
      h1_mbb_4ljjsel_sidebandSX     [pr][fs]->Sumw2(true);

      // sidebandDX
      h1_m4l_4ljjsel_sidebandDX     [pr][fs] = new TH1F("h1_m4l_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events/4 GeV", 9, 135., 171.);
      h1_m4l_4ljjsel_sidebandDX     [pr][fs]->Sumw2(true);
      h1_pT4l_4ljjsel_sidebandDX    [pr][fs] = new TH1F("h1_pT4l_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";4 leptons pT (GeV); Events/2 GeV", 50, 0., 100.);
      h1_pT4l_4ljjsel_sidebandDX    [pr][fs]->Sumw2(true);
      h1_j1btag_4ljjsel_sidebandDX  [pr][fs] = new TH1F("h1_j1btag_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j1btag_4ljjsel_sidebandDX  [pr][fs]->Sumw2(true);
      h1_j2btag_4ljjsel_sidebandDX  [pr][fs] = new TH1F("h1_j2btag_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j2btag_4ljjsel_sidebandDX  [pr][fs]->Sumw2(true);
      h1_j1pT_4ljjsel_sidebandDX    [pr][fs] = new TH1F("h1_j1pT_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j1pT_4ljjsel_sidebandDX    [pr][fs]->Sumw2(true);
      h1_j2pT_4ljjsel_sidebandDX    [pr][fs] = new TH1F("h1_j2pT_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j2pT_4ljjsel_sidebandDX    [pr][fs]->Sumw2(true);
      h1_MET_4ljjsel_sidebandDX     [pr][fs] = new TH1F("h1_MET_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";MET (GeV); Events/5 GeV", 40, 0., 200.);
      h1_MET_4ljjsel_sidebandDX     [pr][fs]->Sumw2(true);
      h1_DeltaRhh_4ljjsel_sidebandDX[pr][fs] = new TH1F("h1_DeltaRhh_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";Delta R H4l-Hbb; Events/0.4 ", 25, 0., 10.);
      h1_DeltaRhh_4ljjsel_sidebandDX[pr][fs]->Sumw2(true);
      h1_mbb_4ljjsel_sidebandDX     [pr][fs] = new TH1F("h1_mbb_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{jj} (GeV); Events/5 GeV", 40, 0., 200.);
      h1_mbb_4ljjsel_sidebandDX     [pr][fs]->Sumw2(true);

      // sidebandSXDX
      h1_m4l_4ljjsel_sidebandSXDX     [pr][fs] = new TH1F("h1_m4l_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events/4 GeV", 19, 95., 171.);
      h1_m4l_4ljjsel_sidebandSXDX     [pr][fs]->Sumw2(true);
      h1_pT4l_4ljjsel_sidebandSXDX    [pr][fs] = new TH1F("h1_pT4l_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";4 leptons pT (GeV); Events/2 GeV", 50, 0., 100.);
      h1_pT4l_4ljjsel_sidebandSXDX    [pr][fs]->Sumw2(true);
      h1_j1btag_4ljjsel_sidebandSXDX  [pr][fs] = new TH1F("h1_j1btag_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j1btag_4ljjsel_sidebandSXDX  [pr][fs]->Sumw2(true);
      h1_j2btag_4ljjsel_sidebandSXDX  [pr][fs] = new TH1F("h1_j2btag_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j2btag_4ljjsel_sidebandSXDX  [pr][fs]->Sumw2(true);
      h1_j1pT_4ljjsel_sidebandSXDX    [pr][fs] = new TH1F("h1_j1pT_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j1pT_4ljjsel_sidebandSXDX    [pr][fs]->Sumw2(true);
      h1_j2pT_4ljjsel_sidebandSXDX    [pr][fs] = new TH1F("h1_j2pT_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j2pT_4ljjsel_sidebandSXDX    [pr][fs]->Sumw2(true);
      h1_MET_4ljjsel_sidebandSXDX     [pr][fs] = new TH1F("h1_MET_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";MET (GeV); Events/5 GeV", 40, 0., 200.);
      h1_MET_4ljjsel_sidebandSXDX     [pr][fs]->Sumw2(true);
      h1_DeltaRhh_4ljjsel_sidebandSXDX[pr][fs] = new TH1F("h1_DeltaRhh_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";Delta R H4l-Hbb; Events/0.4 ", 25, 0., 10.);
      h1_DeltaRhh_4ljjsel_sidebandSXDX[pr][fs]->Sumw2(true);
      h1_mbb_4ljjsel_sidebandSXDX     [pr][fs] = new TH1F("h1_mbb_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{jj} (GeV); Events/5 GeV", 40, 0., 200.);
      h1_mbb_4ljjsel_sidebandSXDX     [pr][fs]->Sumw2(true);


      //BDT input histos
      h1_pT4l_4ljjsel    [pr][fs] = new TH1F("h1_pT4l_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";4 leptons pT (GeV); Events/2 GeV", 50, 0., 100.);
      h1_pT4l_4ljjsel    [pr][fs]->Sumw2(true);
      h1_j1btag_4ljjsel  [pr][fs] = new TH1F("h1_j1btag_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j1btag_4ljjsel  [pr][fs]->Sumw2(true);
      h1_j2btag_4ljjsel  [pr][fs] = new TH1F("h1_j2btag_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 DeepCSV; Events/0.04", 25, 0., 1.);
      h1_j2btag_4ljjsel  [pr][fs]->Sumw2(true);
      h1_j1pT_4ljjsel    [pr][fs] = new TH1F("h1_j1pT_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j1pT_4ljjsel    [pr][fs]->Sumw2(true);
      h1_j2pT_4ljjsel    [pr][fs] = new TH1F("h1_j2pT_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 pT(GeV); Events/5 GeV", 40, 0., 200.);
      h1_j2pT_4ljjsel    [pr][fs]->Sumw2(true);
      h1_MET_4ljjsel     [pr][fs] = new TH1F("h1_MET_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";MET (GeV); Events/5 GeV", 40, 0., 200.);
      h1_MET_4ljjsel     [pr][fs]->Sumw2(true);
      h1_DeltaRhh_4ljjsel[pr][fs] = new TH1F("h1_DeltaRhh_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";Delta R H4l-Hbb; Events/0.4 ", 25, 0., 10.);
      h1_DeltaRhh_4ljjsel[pr][fs]->Sumw2(true);
      h1_mbb_4ljjsel     [pr][fs] = new TH1F("h1_mbb_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{jj} (GeV); Events/5 GeV", 40, 0., 200.); 
      h1_mbb_4ljjsel     [pr][fs]->Sumw2(true);
    }
  }

  // define 2D histos
  TH2F* h2_m4lvsmbb_4ljjsel[nProcesses][nFinalStates];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates; fs++){
      h2_m4lvsmbb_4ljjsel[pr][fs] = new TH2F("h2_m4lvsmbb_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} ;m_{bb}", 20, 115., 135., 240, 60., 180.);
    }
  }
  

  int currentFinalState;
  int currentProcess;  


  //--- loop over all datasets
  for(int d=0; d<nDatasets; d++){

    currentProcess = -1;

    if(datasets[d]=="AllData") currentProcess = Data;
    if(datasets[d]=="HH4lbb_Angela") currentProcess = HH;
    //    if(datasets[d]=="HH4lbb_Ilirjan") currentProcess = HH;
    if(datasets[d]=="ggH125") currentProcess = ggH;
    if(datasets[d]=="VBFH125") currentProcess = VBF;
    if(datasets[d]=="WplusH125" ||
       datasets[d]=="WminusH125" || 
       datasets[d]=="ZH125") currentProcess = VH; 
    if(datasets[d]=="ttH125") currentProcess = ttH;
    if(datasets[d]=="bbH125") currentProcess = bbH;
    //    if(datasets[d]=="ZZTo4lamcatnlo") currentProcess = qqZZ;
    if(datasets[d]=="ZZTo4lext2") currentProcess = qqZZ;
    //if(datasets[d]=="ZZTo4l") currentProcess = qqZZ;
    if(datasets[d]=="ggTo4e_Contin_MCFM701" ||
       datasets[d]=="ggTo4mu_Contin_MCFM701" ||
       datasets[d]=="ggTo4tau_Contin_MCFM701" ||
       datasets[d]=="ggTo2e2mu_Contin_MCFM701" ||
       datasets[d]=="ggTo2e2tau_Contin_MCFM701" ||
       datasets[d]=="ggTo2mu2tau_Contin_MCFM701") currentProcess = ggZZ;
    if(datasets[d]=="TTZToLLNuNu_M10") currentProcess = TTZ; 
    if(datasets[d]=="TTWJetsToLNu") currentProcess = TTW;
    if(datasets[d]=="WWZ" ||
       datasets[d]=="WZZ" ||
       datasets[d]=="ZZZ") currentProcess = VVV;
    //    if(datasets[d]=="ZXbkg_4ljjsel") currentProcess = ZXbkg;
    if(datasets[d]=="ZXbkg_4lsel") currentProcess = ZXbkg;
    

    // select input file
    TString inputFileName;
    if(currentProcess == Data) inputFileName = inDataPath + datasets[d] + "/ZZXAnalysis.root";
    else inputFileName = inFilePath + datasets[d] + "/ZZXAnalysis.root";
    cout<<"Opening file "<<inputFileName<<" ..."<<endl;
    inputFile[d] = TFile::Open(inputFileName);

    if(currentProcess == ZXbkg){
      hCounters[d] = 0;
      gen_sumWeights[d] = 0.;
      partialSampleWeight[d] = 0;
      inputTree[d] = (TTree*)inputFile[d]->Get("candTree");
    }
    else{
      hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
      gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);
      partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d];
      inputTree[d] = (TTree*)inputFile[d]->Get("ZZTree/candTree");
    }
    cout<<"debug "<<datasets[d]<<" "<<currentProcess<<" "<<sProcess[currentProcess]<<endl;

    // set branch addresses
    if(currentProcess != ZXbkg){
      inputTree[d]->SetBranchAddress("RunNumber", &nRun);
      inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
      inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
      if(currentProcess != Data) inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
      if(currentProcess != Data) inputTree[d]->SetBranchAddress("xsec", &xsec);
      if(currentProcess != Data) inputTree[d]->SetBranchAddress("L1prefiringWeight", &L1prefiringWeight);
    } 
    if(currentProcess == ZXbkg) inputTree[d]->SetBranchAddress("weight", &weight);
    inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);  
    inputTree[d]->SetBranchAddress("Z1Flav", &Z1Flav);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("ZZEta", &ZZEta);
    inputTree[d]->SetBranchAddress("ZZPhi", &ZZPhi);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("JetPhi",  &JetPhi);
    inputTree[d]->SetBranchAddress("JetBTagger",  &JetBTagger);
    inputTree[d]->SetBranchAddress("JetHadronFlavour",  &JetHadronFlavour);
    inputTree[d]->SetBranchAddress("PFMET",  &PFMET);  
    if(currentProcess == ggZZ){
      inputTree[d]->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
    }
    if(currentProcess == qqZZ){
      inputTree[d]->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
      inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi);
      inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
      inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt);
    }


    // --------------------------------------------------------
    // --- first loop over input tree to get norm for BtagSF
    cout<<"First loop over input files to get norm for BTagSF ..."<<endl;
    if(currentProcess == Data || currentProcess == ZXbkg){
      sum_events[d] = 1.;      
      sum_BTagSF[d] = 1.;
    }
    else{
      Long64_t entries1 = inputTree[d]->GetEntries();    
      cout<<"Processing file: "<< datasets[d] << " (" << entries1 <<" entries) ..."<< endl;    
      for(Long64_t z=0; z<entries1; z++){

        inputTree[d]->GetEntry(z);
  
        if( !(ZZsel>=0) ) continue;

        // 4l selection
        if(LepEta->size() != 4){
          cout << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	  continue;
        }


        // compute SF
        double * scaleFactors;
        scaleFactors = evalEventSF( int(JetPt->size()), JetHadronFlavour, JetEta, JetPt, JetBTagger, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down );

        // total counters for BTagSF norm --- 4lsel
        sum_events[d] += 1.; 
        sum_BTagSF[d] += scaleFactors[0]; 

    
      } // end first loop over entries 
    }// end else

    cout<<datasets[d]<<" "<<sum_events[d]<<" "<<sum_BTagSF[d]<<endl;

    // --- control for norm
    if( sum_events[d] == 0. || std::isnan(sum_events[d]) ){ sum_events[d] = 1.; }
    if( sum_BTagSF[d] == 0. || std::isnan(sum_BTagSF[d]) ){ sum_BTagSF[d] = 1.; }

    cout<<datasets[d]<<" "<<sum_events[d]<<" "<<sum_BTagSF[d]<<endl;
    // --------------------------------------------------------


    // --- second loop over input tree to do all the rest
    Long64_t entries2 = inputTree[d]->GetEntries();    
    cout<<"Second loop over input files to do all the rest ..."<<endl;
    cout<<"Processing file: "<< datasets[d] << " (" << entries2 <<" entries) ..."<< endl;    
    for(Long64_t z=0; z<entries2; z++){

      inputTree[d]->GetEntry(z);

  
      if( !(ZZsel>=0) ) continue;

      // 4l selection
      if(LepEta->size() != 4){
	cout << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
      }

      
      // --- fill k factors
      Float_t kfactor = 1.;
      if(currentProcess == qqZZ) { kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; } // qqZZ sample                      
      else if(currentProcess == ggZZ) { kfactor = KFactor_QCD_ggZZ_Nominal; } //ggZZ samples 

      // compute SF
      double * scaleFactors;
      if(currentProcess != Data && currentProcess != ZXbkg){
        scaleFactors = evalEventSF( int(JetPt->size()), JetHadronFlavour, JetEta, JetPt, JetBTagger, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down );
      }

      // --- event weights
      Double_t eventWeight = 1.;
      if(currentProcess != Data && currentProcess != ZXbkg){
        eventWeight = partialSampleWeight[d] *xsec *kfactor *overallEventWeight *L1prefiringWeight *scaleFactors[0] *sum_events[d]/sum_BTagSF[d];
      }
      if(currentProcess == ZXbkg) eventWeight = weight; //ZX weight



      // --- select final state
      currentFinalState = -1;
      if (currentProcess != ZXbkg){
	if (Z1Flav == -121){
	  if (Z2Flav == -121) {currentFinalState = fs_4e;}
          else if ( Z2Flav == -169) {currentFinalState = fs_2e2mu;}
          else { cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z2Flav = " << Z2Flav << endl;}
        }
	else if (Z1Flav == -169){
	  if (Z2Flav == -121) { currentFinalState = fs_2e2mu; }
	  else if (Z2Flav == -169) { currentFinalState = fs_4mu; }
          else { cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z2Flav = " << Z2Flav << endl; }
	}
	else{ cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z1Flav = " << Z1Flav << endl; }
      }

      else if (currentProcess == ZXbkg){
	if (Z1Flav == -121){
	  if (Z2Flav == +121) { currentFinalState = fs_4e; }
          else if (Z2Flav == +169) { currentFinalState = fs_2e2mu; }
          else { cerr << "error, Z2Flav: " << endl; }
        }
	else if (Z1Flav == -169){
	  if (Z2Flav == +121) { currentFinalState = fs_2e2mu; }
	  else if (Z2Flav == +169) { currentFinalState = fs_4mu; }
	  else { cerr << "error, Z2Flav: " << endl; }
	}
	else{ cerr << "error, Z2Flav: " << endl; }
      }





      // --- fill yields after 4l sel
      hYields_4lsel[currentProcess][currentFinalState]->Fill(0.5, eventWeight);

      // --- fill number of non weighted events selected after 4l sel
      hEvents_4lsel[currentProcess][currentFinalState]->Fill(0.5, 1.);


 


      //JETSELECTION--------------------------------------------------  
      // at least 2 jets in the acceptance
      if (JetPt->size() < 2) continue;   
          
      // index of jet with max btagger value (j1)
      int dj1_ = distance( JetBTagger->begin(), max_element(JetBTagger->begin(), JetBTagger->end()));
      
      // index of jet with 2nd highest btagger value (j2)  
      float max2JetBtag = -9999.;
      int dj2_ = -9999;
      for(UInt_t i=0; i<JetBTagger->size(); i++){
  
        if(i == dj1_) continue;
        float temp = JetBTagger->at(i);
        if (temp > max2JetBtag){
          max2JetBtag = temp;
          dj2_ = i;
        }
      }


      
      // build 2 jets tlorentzvectors
      TLorentzVector tlzvec_j1_;
      tlzvec_j1_.SetPtEtaPhiM(JetPt->at(dj1_), JetEta->at(dj1_), JetPhi->at(dj1_), JetMass->at(dj1_));
      TLorentzVector tlzvec_j2_;
      tlzvec_j2_.SetPtEtaPhiM(JetPt->at(dj2_), JetEta->at(dj2_), JetPhi->at(dj2_), JetMass->at(dj2_));
      
      
      // build H->bb tlorentzvector
      TLorentzVector Hbb_tlzvec = tlzvec_j1_ + tlzvec_j2_;
      Float_t bbMass = Hbb_tlzvec.M();

      // build H-H DeltaR
      float DeltaPhi = ZZPhi - Hbb_tlzvec.Phi();
      if( fabs(DeltaPhi) > acos(-1) ) { DeltaPhi = (2*acos(-1)) - fabs(DeltaPhi); }
      float DeltaEta = ZZEta - Hbb_tlzvec.Eta();
      float DeltaR = sqrt( DeltaPhi*DeltaPhi + DeltaEta*DeltaEta );






      // --- fill histos after 4ljj sel (no mass cut)
      h1_m4l_4ljjsel  [currentProcess][currentFinalState]->Fill(ZZMass, eventWeight);



      // --- MASSCUT: mass cut: signal region
      if(ZZMass >= 115. && ZZMass <= 135.){  // 115 < ZZMass < 135 GeV

        // --- fill yields after 4ljj sel
        hYields_4ljjsel[currentProcess][currentFinalState]->Fill(0.5, eventWeight);
        // --- fill number of non weighted events selected after 4ljj sel
        hEvents_4ljjsel[currentProcess][currentFinalState]->Fill(0.5, 1.);
        
        // --- fill histos after 4ljj sel
        // (BDT input histos)
        for(int i=0; i<LepPt->size(); i++){
          h1_pT4l_4ljjsel[currentProcess][currentFinalState]->Fill(LepPt->at(i), eventWeight);
        }
        h1_j1btag_4ljjsel  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj1_), eventWeight);
        h1_j2btag_4ljjsel  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj2_), eventWeight);
        h1_j1pT_4ljjsel    [currentProcess][currentFinalState]->Fill(JetPt->at(dj1_),      eventWeight);  
        h1_j2pT_4ljjsel    [currentProcess][currentFinalState]->Fill(JetPt->at(dj2_),      eventWeight);
        h1_MET_4ljjsel     [currentProcess][currentFinalState]->Fill(PFMET,                eventWeight);
        h1_DeltaRhh_4ljjsel[currentProcess][currentFinalState]->Fill(DeltaR,               eventWeight);
        h1_mbb_4ljjsel     [currentProcess][currentFinalState]->Fill(bbMass,               eventWeight);
    
        //2D histo
        h2_m4lvsmbb_4ljjsel[currentProcess][currentFinalState]->Fill(ZZMass, bbMass, eventWeight);

      }
      else{

        // --- fill yields in the sidebands
        hYields_4ljjsel_sidebands[currentProcess][currentFinalState]->Fill(0.5, eventWeight);

        // --- fill histos after 4ljj sel: sidebands
        for(int i=0; i<LepPt->size(); i++){
          h1_pT4l_4ljjsel_sidebands[currentProcess][currentFinalState]->Fill(LepPt->at(i), eventWeight);
        }
        h1_j1btag_4ljjsel_sidebands  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj1_), eventWeight);
        h1_j2btag_4ljjsel_sidebands  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj2_), eventWeight);
        h1_j1pT_4ljjsel_sidebands    [currentProcess][currentFinalState]->Fill(JetPt->at(dj1_),      eventWeight);  
        h1_j2pT_4ljjsel_sidebands    [currentProcess][currentFinalState]->Fill(JetPt->at(dj2_),      eventWeight);
        h1_MET_4ljjsel_sidebands     [currentProcess][currentFinalState]->Fill(PFMET,                eventWeight);
        h1_DeltaRhh_4ljjsel_sidebands[currentProcess][currentFinalState]->Fill(DeltaR,               eventWeight);
        h1_mbb_4ljjsel_sidebands     [currentProcess][currentFinalState]->Fill(bbMass,               eventWeight);

      } //end else (sidebands)

      // --- sidebandSX
      if(ZZMass >= 95. && ZZMass < 115.){

	// --- fill yields in sidebandSX
        hYields_4ljjsel_sidebandSX[currentProcess][currentFinalState]->Fill(0.5, eventWeight);

        // --- fill histos after 4ljj sel: sidebandSX
        for(int i=0; i<LepPt->size(); i++){
          h1_pT4l_4ljjsel_sidebandSX[currentProcess][currentFinalState]->Fill(LepPt->at(i), eventWeight);
        }
        h1_m4l_4ljjsel_sidebandSX     [currentProcess][currentFinalState]->Fill(ZZMass,               eventWeight);
        h1_j1btag_4ljjsel_sidebandSX  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj1_), eventWeight);
        h1_j2btag_4ljjsel_sidebandSX  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj2_), eventWeight);
        h1_j1pT_4ljjsel_sidebandSX    [currentProcess][currentFinalState]->Fill(JetPt->at(dj1_),      eventWeight);  
        h1_j2pT_4ljjsel_sidebandSX    [currentProcess][currentFinalState]->Fill(JetPt->at(dj2_),      eventWeight);
        h1_MET_4ljjsel_sidebandSX     [currentProcess][currentFinalState]->Fill(PFMET,                eventWeight);
        h1_DeltaRhh_4ljjsel_sidebandSX[currentProcess][currentFinalState]->Fill(DeltaR,               eventWeight);
        h1_mbb_4ljjsel_sidebandSX     [currentProcess][currentFinalState]->Fill(bbMass,               eventWeight);

      } // end sidebandSX

      // --- sidebandDX
      if(ZZMass > 135. && ZZMass <= 170.){

	// --- fill yields in sidebandDX
        hYields_4ljjsel_sidebandDX[currentProcess][currentFinalState]->Fill(0.5, eventWeight);

        // --- fill histos after 4ljj sel: sidebandDX
        for(int i=0; i<LepPt->size(); i++){
          h1_pT4l_4ljjsel_sidebandDX[currentProcess][currentFinalState]->Fill(LepPt->at(i), eventWeight);
        }
        h1_m4l_4ljjsel_sidebandDX     [currentProcess][currentFinalState]->Fill(ZZMass,               eventWeight);
        h1_j1btag_4ljjsel_sidebandDX  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj1_), eventWeight);
        h1_j2btag_4ljjsel_sidebandDX  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj2_), eventWeight);
        h1_j1pT_4ljjsel_sidebandDX    [currentProcess][currentFinalState]->Fill(JetPt->at(dj1_),      eventWeight);  
        h1_j2pT_4ljjsel_sidebandDX    [currentProcess][currentFinalState]->Fill(JetPt->at(dj2_),      eventWeight);
        h1_MET_4ljjsel_sidebandDX     [currentProcess][currentFinalState]->Fill(PFMET,                eventWeight);
        h1_DeltaRhh_4ljjsel_sidebandDX[currentProcess][currentFinalState]->Fill(DeltaR,               eventWeight);
        h1_mbb_4ljjsel_sidebandDX     [currentProcess][currentFinalState]->Fill(bbMass,               eventWeight);

      } // end sidebandDX


      // --- sidebandSXDX
      if((ZZMass >= 95. && ZZMass < 115.) || (ZZMass > 135. && ZZMass <= 170.)){

	// --- fill yields in sidebandSXDX
        hYields_4ljjsel_sidebandSXDX[currentProcess][currentFinalState]->Fill(0.5, eventWeight);

        // --- fill histos after 4ljj sel: sidebandSXDX
        for(int i=0; i<LepPt->size(); i++){
          h1_pT4l_4ljjsel_sidebandSXDX[currentProcess][currentFinalState]->Fill(LepPt->at(i), eventWeight);
        }
        h1_m4l_4ljjsel_sidebandSXDX     [currentProcess][currentFinalState]->Fill(ZZMass,               eventWeight);
        h1_j1btag_4ljjsel_sidebandSXDX  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj1_), eventWeight);
        h1_j2btag_4ljjsel_sidebandSXDX  [currentProcess][currentFinalState]->Fill(JetBTagger->at(dj2_), eventWeight);
        h1_j1pT_4ljjsel_sidebandSXDX    [currentProcess][currentFinalState]->Fill(JetPt->at(dj1_),      eventWeight);  
        h1_j2pT_4ljjsel_sidebandSXDX    [currentProcess][currentFinalState]->Fill(JetPt->at(dj2_),      eventWeight);
        h1_MET_4ljjsel_sidebandSXDX     [currentProcess][currentFinalState]->Fill(PFMET,                eventWeight);
        h1_DeltaRhh_4ljjsel_sidebandSXDX[currentProcess][currentFinalState]->Fill(DeltaR,               eventWeight);
        h1_mbb_4ljjsel_sidebandSXDX     [currentProcess][currentFinalState]->Fill(bbMass,               eventWeight);

      } // end sidebandSXDX


    }//end loop over tree events

  }//end loop over datasets


  // --- RESCALE ZX HISTOS ---
  // --- riscalo solo histo che disegno
  for(int fs=0; fs<nFinalStates; fs++){
    // // fullmass range
    // h1_m4l_4ljjsel           [ZXbkg][fs]->Scale(rescale_ZX[fs]);
    // sidebands
    h1_pT4l_4ljjsel_sidebands    [ZXbkg][fs]->Scale(2.);
    h1_j1btag_4ljjsel_sidebands  [ZXbkg][fs]->Scale(2.);
    h1_j2btag_4ljjsel_sidebands  [ZXbkg][fs]->Scale(2.);
    h1_j1pT_4ljjsel_sidebands    [ZXbkg][fs]->Scale(2.);
    h1_j2pT_4ljjsel_sidebands    [ZXbkg][fs]->Scale(2.);
    h1_MET_4ljjsel_sidebands     [ZXbkg][fs]->Scale(2.);
    h1_DeltaRhh_4ljjsel_sidebands[ZXbkg][fs]->Scale(2.);
    h1_mbb_4ljjsel_sidebands     [ZXbkg][fs]->Scale(2.);
    // // sidebandSX
    // h1_pT4l_4ljjsel_sidebandSX    [ZXbkg][fs]->Scale(2.);
    // h1_m4l_4ljjsel_sidebandSX     [ZXbkg][fs]->Scale(2.);
    // h1_j1btag_4ljjsel_sidebandSX  [ZXbkg][fs]->Scale(2.);
    // h1_j2btag_4ljjsel_sidebandSX  [ZXbkg][fs]->Scale(2.);
    // h1_j1pT_4ljjsel_sidebandSX    [ZXbkg][fs]->Scale(2.);
    // h1_j2pT_4ljjsel_sidebandSX    [ZXbkg][fs]->Scale(2.);
    // h1_MET_4ljjsel_sidebandSX     [ZXbkg][fs]->Scale(2.);
    // h1_DeltaRhh_4ljjsel_sidebandSX[ZXbkg][fs]->Scale(2.);
    // h1_mbb_4ljjsel_sidebandSX     [ZXbkg][fs]->Scale(2.);
    // // sidebandDX
    // h1_pT4l_4ljjsel_sidebandDX    [ZXbkg][fs]->Scale(2.);
    // h1_m4l_4ljjsel_sidebandDX     [ZXbkg][fs]->Scale(2.);
    // h1_j1btag_4ljjsel_sidebandDX  [ZXbkg][fs]->Scale(2.);
    // h1_j2btag_4ljjsel_sidebandDX  [ZXbkg][fs]->Scale(2.);
    // h1_j1pT_4ljjsel_sidebandDX    [ZXbkg][fs]->Scale(2.);
    // h1_j2pT_4ljjsel_sidebandDX    [ZXbkg][fs]->Scale(2.);
    // h1_MET_4ljjsel_sidebandDX     [ZXbkg][fs]->Scale(2.);
    // h1_DeltaRhh_4ljjsel_sidebandDX[ZXbkg][fs]->Scale(2.);
    // h1_mbb_4ljjsel_sidebandDX     [ZXbkg][fs]->Scale(2.);
    // (BDT input)
    h1_pT4l_4ljjsel    [ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_pT4l_4ljjsel    [ZXbkg][fs]->Integral());
    h1_j1btag_4ljjsel  [ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_j1btag_4ljjsel  [ZXbkg][fs]->Integral());
    h1_j2btag_4ljjsel  [ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_j2btag_4ljjsel  [ZXbkg][fs]->Integral());
    h1_j1pT_4ljjsel    [ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_j1pT_4ljjsel    [ZXbkg][fs]->Integral());
    h1_j2pT_4ljjsel    [ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_j2pT_4ljjsel    [ZXbkg][fs]->Integral());
    h1_MET_4ljjsel     [ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_MET_4ljjsel     [ZXbkg][fs]->Integral());
    h1_DeltaRhh_4ljjsel[ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_DeltaRhh_4ljjsel[ZXbkg][fs]->Integral());
    h1_mbb_4ljjsel     [ZXbkg][fs]->Scale(rescale_ZX[fs] / h1_mbb_4ljjsel     [ZXbkg][fs]->Integral());
  }


  //---fill inclusive yields and histos
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates; fs++){
      hYields_4lsel   [pr][nFinalStates]->Add(hYields_4lsel   [pr][fs]);
      hYields_4ljjsel [pr][nFinalStates]->Add(hYields_4ljjsel [pr][fs]);
      hEvents_4lsel   [pr][nFinalStates]->Add(hEvents_4lsel   [pr][fs]);
      hEvents_4ljjsel [pr][nFinalStates]->Add(hEvents_4ljjsel [pr][fs]);
      hYields_4ljjsel_sidebands [pr][nFinalStates]->Add(hYields_4ljjsel_sidebands [pr][fs]);
      hYields_4ljjsel_sidebandSX[pr][nFinalStates]->Add(hYields_4ljjsel_sidebandSX[pr][fs]);
      hYields_4ljjsel_sidebandDX[pr][nFinalStates]->Add(hYields_4ljjsel_sidebandDX[pr][fs]);
      hYields_4ljjsel_sidebandSXDX[pr][nFinalStates]->Add(hYields_4ljjsel_sidebandSXDX[pr][fs]);

      // (h1 histos)
      // 4ljjsel
      // fullmass range
      h1_m4l_4ljjsel  [pr][nFinalStates]->Add(h1_m4l_4ljjsel  [pr][fs]); 
      // sidebands
      h1_pT4l_4ljjsel_sidebands    [pr][nFinalStates]->Add(h1_pT4l_4ljjsel_sidebands    [pr][fs]);
      h1_j1btag_4ljjsel_sidebands  [pr][nFinalStates]->Add(h1_j1btag_4ljjsel_sidebands  [pr][fs]);
      h1_j2btag_4ljjsel_sidebands  [pr][nFinalStates]->Add(h1_j2btag_4ljjsel_sidebands  [pr][fs]);
      h1_j1pT_4ljjsel_sidebands    [pr][nFinalStates]->Add(h1_j1pT_4ljjsel_sidebands    [pr][fs]);
      h1_j2pT_4ljjsel_sidebands    [pr][nFinalStates]->Add(h1_j2pT_4ljjsel_sidebands    [pr][fs]);
      h1_MET_4ljjsel_sidebands     [pr][nFinalStates]->Add(h1_MET_4ljjsel_sidebands     [pr][fs]);
      h1_DeltaRhh_4ljjsel_sidebands[pr][nFinalStates]->Add(h1_DeltaRhh_4ljjsel_sidebands[pr][fs]);
      h1_mbb_4ljjsel_sidebands     [pr][nFinalStates]->Add(h1_mbb_4ljjsel_sidebands     [pr][fs]);
      // sidebandSX
      h1_m4l_4ljjsel_sidebandSX     [pr][nFinalStates]->Add(h1_m4l_4ljjsel_sidebandSX     [pr][fs]);
      h1_pT4l_4ljjsel_sidebandSX    [pr][nFinalStates]->Add(h1_pT4l_4ljjsel_sidebandSX    [pr][fs]);
      h1_j1btag_4ljjsel_sidebandSX  [pr][nFinalStates]->Add(h1_j1btag_4ljjsel_sidebandSX  [pr][fs]);
      h1_j2btag_4ljjsel_sidebandSX  [pr][nFinalStates]->Add(h1_j2btag_4ljjsel_sidebandSX  [pr][fs]);
      h1_j1pT_4ljjsel_sidebandSX    [pr][nFinalStates]->Add(h1_j1pT_4ljjsel_sidebandSX    [pr][fs]);
      h1_j2pT_4ljjsel_sidebandSX    [pr][nFinalStates]->Add(h1_j2pT_4ljjsel_sidebandSX    [pr][fs]);
      h1_MET_4ljjsel_sidebandSX     [pr][nFinalStates]->Add(h1_MET_4ljjsel_sidebandSX     [pr][fs]);
      h1_DeltaRhh_4ljjsel_sidebandSX[pr][nFinalStates]->Add(h1_DeltaRhh_4ljjsel_sidebandSX[pr][fs]);
      h1_mbb_4ljjsel_sidebandSX     [pr][nFinalStates]->Add(h1_mbb_4ljjsel_sidebandSX     [pr][fs]);
      // sidebandDX
      h1_m4l_4ljjsel_sidebandDX     [pr][nFinalStates]->Add(h1_m4l_4ljjsel_sidebandDX     [pr][fs]);
      h1_pT4l_4ljjsel_sidebandDX    [pr][nFinalStates]->Add(h1_pT4l_4ljjsel_sidebandDX    [pr][fs]);
      h1_j1btag_4ljjsel_sidebandDX  [pr][nFinalStates]->Add(h1_j1btag_4ljjsel_sidebandDX  [pr][fs]);
      h1_j2btag_4ljjsel_sidebandDX  [pr][nFinalStates]->Add(h1_j2btag_4ljjsel_sidebandDX  [pr][fs]);
      h1_j1pT_4ljjsel_sidebandDX    [pr][nFinalStates]->Add(h1_j1pT_4ljjsel_sidebandDX    [pr][fs]);
      h1_j2pT_4ljjsel_sidebandDX    [pr][nFinalStates]->Add(h1_j2pT_4ljjsel_sidebandDX    [pr][fs]);
      h1_MET_4ljjsel_sidebandDX     [pr][nFinalStates]->Add(h1_MET_4ljjsel_sidebandDX     [pr][fs]);
      h1_DeltaRhh_4ljjsel_sidebandDX[pr][nFinalStates]->Add(h1_DeltaRhh_4ljjsel_sidebandDX[pr][fs]);
      h1_mbb_4ljjsel_sidebandDX     [pr][nFinalStates]->Add(h1_mbb_4ljjsel_sidebandDX     [pr][fs]);
      // sidebandSXDX
      h1_m4l_4ljjsel_sidebandSXDX     [pr][nFinalStates]->Add(h1_m4l_4ljjsel_sidebandSXDX     [pr][fs]);
      h1_pT4l_4ljjsel_sidebandSXDX    [pr][nFinalStates]->Add(h1_pT4l_4ljjsel_sidebandSXDX    [pr][fs]);
      h1_j1btag_4ljjsel_sidebandSXDX  [pr][nFinalStates]->Add(h1_j1btag_4ljjsel_sidebandSXDX  [pr][fs]);
      h1_j2btag_4ljjsel_sidebandSXDX  [pr][nFinalStates]->Add(h1_j2btag_4ljjsel_sidebandSXDX  [pr][fs]);
      h1_j1pT_4ljjsel_sidebandSXDX    [pr][nFinalStates]->Add(h1_j1pT_4ljjsel_sidebandSXDX    [pr][fs]);
      h1_j2pT_4ljjsel_sidebandSXDX    [pr][nFinalStates]->Add(h1_j2pT_4ljjsel_sidebandSXDX    [pr][fs]);
      h1_MET_4ljjsel_sidebandSXDX     [pr][nFinalStates]->Add(h1_MET_4ljjsel_sidebandSXDX     [pr][fs]);
      h1_DeltaRhh_4ljjsel_sidebandSXDX[pr][nFinalStates]->Add(h1_DeltaRhh_4ljjsel_sidebandSXDX[pr][fs]);
      h1_mbb_4ljjsel_sidebandSXDX     [pr][nFinalStates]->Add(h1_mbb_4ljjsel_sidebandSXDX     [pr][fs]);
      // (BDT input histos)
      h1_pT4l_4ljjsel    [pr][nFinalStates]->Add(h1_pT4l_4ljjsel    [pr][fs]);
      h1_j1btag_4ljjsel  [pr][nFinalStates]->Add(h1_j1btag_4ljjsel  [pr][fs]);
      h1_j2btag_4ljjsel  [pr][nFinalStates]->Add(h1_j2btag_4ljjsel  [pr][fs]);
      h1_j1pT_4ljjsel    [pr][nFinalStates]->Add(h1_j1pT_4ljjsel    [pr][fs]);
      h1_j2pT_4ljjsel    [pr][nFinalStates]->Add(h1_j2pT_4ljjsel    [pr][fs]);
      h1_MET_4ljjsel     [pr][nFinalStates]->Add(h1_MET_4ljjsel     [pr][fs]);
      h1_DeltaRhh_4ljjsel[pr][nFinalStates]->Add(h1_DeltaRhh_4ljjsel[pr][fs]);
      h1_mbb_4ljjsel     [pr][nFinalStates]->Add(h1_mbb_4ljjsel     [pr][fs]);

    }
  }

  //---save yields in a root file
  TString fout_yields_name = "f_yields_"+ sYear + ".root";
  TFile* fout_yields = new TFile(fout_yields_name, "recreate");
  fout_yields->cd();
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      hYields_4lsel  [pr][fs]->Write();
      hYields_4ljjsel[pr][fs]->Write();
      hEvents_4lsel  [pr][fs]->Write();
      hEvents_4ljjsel[pr][fs]->Write();
      hYields_4ljjsel_sidebands [pr][fs]->Write();
      hYields_4ljjsel_sidebandSX[pr][fs]->Write();
      hYields_4ljjsel_sidebandDX[pr][fs]->Write();
      hYields_4ljjsel_sidebandSXDX[pr][fs]->Write();
    }
  }
  fout_yields->Close();

  //---save 1D histos in a root file
  TString fout_1Dhistos_name = "f_histos_h1_4ljjsel_" + sYear +".root";
  TFile* fout_1Dhistos = new TFile(fout_1Dhistos_name, "recreate");
  fout_1Dhistos->cd();
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      //4ljjsel
      // fullmass range
      h1_m4l_4ljjsel  [pr][fs]->Write();
      // sidebands
      h1_pT4l_4ljjsel_sidebands    [pr][fs]->Write();
      h1_j1btag_4ljjsel_sidebands  [pr][fs]->Write();     
      h1_j2btag_4ljjsel_sidebands  [pr][fs]->Write();
      h1_j1pT_4ljjsel_sidebands    [pr][fs]->Write();
      h1_j2pT_4ljjsel_sidebands    [pr][fs]->Write();
      h1_MET_4ljjsel_sidebands     [pr][fs]->Write();
      h1_DeltaRhh_4ljjsel_sidebands[pr][fs]->Write();
      h1_mbb_4ljjsel_sidebands     [pr][fs]->Write();     
      // sidebandSX
      h1_m4l_4ljjsel_sidebandSX     [pr][fs]->Write();
      h1_pT4l_4ljjsel_sidebandSX    [pr][fs]->Write();
      h1_j1btag_4ljjsel_sidebandSX  [pr][fs]->Write();     
      h1_j2btag_4ljjsel_sidebandSX  [pr][fs]->Write();
      h1_j1pT_4ljjsel_sidebandSX    [pr][fs]->Write();
      h1_j2pT_4ljjsel_sidebandSX    [pr][fs]->Write();
      h1_MET_4ljjsel_sidebandSX     [pr][fs]->Write();
      h1_DeltaRhh_4ljjsel_sidebandSX[pr][fs]->Write();
      h1_mbb_4ljjsel_sidebandSX     [pr][fs]->Write();     
      // sidebandDX
      h1_m4l_4ljjsel_sidebandDX     [pr][fs]->Write();
      h1_pT4l_4ljjsel_sidebandDX    [pr][fs]->Write();
      h1_j1btag_4ljjsel_sidebandDX  [pr][fs]->Write();     
      h1_j2btag_4ljjsel_sidebandDX  [pr][fs]->Write();
      h1_j1pT_4ljjsel_sidebandDX    [pr][fs]->Write();
      h1_j2pT_4ljjsel_sidebandDX    [pr][fs]->Write();
      h1_MET_4ljjsel_sidebandDX     [pr][fs]->Write();
      h1_DeltaRhh_4ljjsel_sidebandDX[pr][fs]->Write();
      h1_mbb_4ljjsel_sidebandDX     [pr][fs]->Write();     
      // sidebandSXDX
      h1_m4l_4ljjsel_sidebandSXDX     [pr][fs]->Write();
      h1_pT4l_4ljjsel_sidebandSXDX    [pr][fs]->Write();
      h1_j1btag_4ljjsel_sidebandSXDX  [pr][fs]->Write();     
      h1_j2btag_4ljjsel_sidebandSXDX  [pr][fs]->Write();
      h1_j1pT_4ljjsel_sidebandSXDX    [pr][fs]->Write();
      h1_j2pT_4ljjsel_sidebandSXDX    [pr][fs]->Write();
      h1_MET_4ljjsel_sidebandSXDX     [pr][fs]->Write();
      h1_DeltaRhh_4ljjsel_sidebandSXDX[pr][fs]->Write();
      h1_mbb_4ljjsel_sidebandSXDX     [pr][fs]->Write();     
      // mass cut: BDT input
      h1_pT4l_4ljjsel    [pr][fs]->Write();
      h1_j1btag_4ljjsel  [pr][fs]->Write();     
      h1_j2btag_4ljjsel  [pr][fs]->Write();
      h1_j1pT_4ljjsel    [pr][fs]->Write();
      h1_j2pT_4ljjsel    [pr][fs]->Write();
      h1_MET_4ljjsel     [pr][fs]->Write();
      h1_DeltaRhh_4ljjsel[pr][fs]->Write();
      h1_mbb_4ljjsel     [pr][fs]->Write();     
    }
  }
  fout_1Dhistos->Close();
    
  //---save 2D histos in a root file
  TString fout_2Dhistos_name = "f_histos_h2_4ljjsel_" + sYear + ".root";
  TFile* fout_2Dhistos = new TFile(fout_2Dhistos_name, "recreate");  
  fout_2Dhistos->cd();
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates; fs++){
      h2_m4lvsmbb_4ljjsel[pr][fs]->Write();
    }
  }
  fout_2Dhistos->Close();


}//end doHistos function


//************************************
//*** printYields_forSync function ***
//************************************
void printYields_forSync(){

  //---input path
  TString sYear;
  if(year==2016)      sYear = "2016";
  else if(year==2017) sYear = "2017";
  else if(year==2018) sYear = "2018";
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;

 
  // retrieve yields histos from file
  TString inFileName = "f_yields_" + sYear + ".root";     
  cout<<"Retrieving Data and MC yields histos from file "<<inFileName<<" ..."<<endl;
  TFile* fInYields = TFile::Open(inFileName);

  TH1F* hTemp0;
  TH1F* hTemp1;
  TH1F* hTemp2;
  TH1F* hTemp3;
  Float_t yield_4lsel   [nProcesses][nFinalStates+1];
  Float_t yield_4ljjsel [nProcesses][nFinalStates+1];
  Float_t nEvent_4lsel  [nProcesses][nFinalStates+1];
  Float_t nEvent_4ljjsel[nProcesses][nFinalStates+1];
  TH1F* yield_4ljjsel_sidebands[nProcesses][nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSX[nProcesses][nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandDX[nProcesses][nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSXDX[nProcesses][nFinalStates+1];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      hTemp0 = (TH1F*)fInYields->Get("hYields_4lsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      yield_4lsel[pr][fs] = hTemp0->GetBinContent(1);
      hTemp1 = (TH1F*)fInYields->Get("hYields_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      yield_4ljjsel[pr][fs] = hTemp1->GetBinContent(1);
      hTemp2 = (TH1F*)fInYields->Get("hEvents_4lsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      nEvent_4lsel[pr][fs] = hTemp2->GetBinContent(1);  
      hTemp3 = (TH1F*)fInYields->Get("hEvents_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      nEvent_4ljjsel[pr][fs] = hTemp3->GetBinContent(1);  

      yield_4ljjsel_sidebands[pr][fs]  = (TH1F*)fInYields->Get("hYields_4ljjsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      yield_4ljjsel_sidebandSX[pr][fs] = (TH1F*)fInYields->Get("hYields_4ljjsel_sidebandSX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      yield_4ljjsel_sidebandDX[pr][fs] = (TH1F*)fInYields->Get("hYields_4ljjsel_sidebandDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      yield_4ljjsel_sidebandSXDX[pr][fs] = (TH1F*)fInYields->Get("hYields_4ljjsel_sidebandSXDX_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
    }
  }


  // // **********************************************
  // // *** PRINT NUMBER OF EVENTS SEL (no weight) ***
  // // **********************************************
  // // --- print number of event selected after 4l sel for sync (no weight)
  // cout<<"print number of event selected after 4l sel for sync (no weight) ... "<<endl; 
  // ofstream f_events4lsel_sync;
  // TString f_events4lsel_sync_name = "nEvents_4lsel_perSync_"+ sYear + ".txt";
  // f_events4lsel_sync.open(f_events4lsel_sync_name);
  // f_events4lsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  // for(int fs=0; fs<nFinalStates+1; fs++){
  //   f_events4lsel_sync<<"|"<<sFinalState[fs]<<" |"<<nEvent_4lsel[HH][fs]<<" |"<<nEvent_4lsel[TTZ][fs]<<" |"<<nEvent_4lsel[ttH][fs]<<" |"<<nEvent_4lsel[qqZZ][fs]+nEvent_4lsel[ggZZ][fs]<<" |"<<nEvent_4lsel[ggH][fs]+nEvent_4lsel[VBF][fs]<<" |"<<nEvent_4lsel[VVV][fs]+nEvent_4lsel[VH][fs]+nEvent_4lsel[TTW][fs]<<" |"<<nEvent_4lsel[ZXbkg][fs]<<" |"<<endl;
  // }
  // f_events4lsel_sync.close();

  // // --- print number of event selected after 4ljj sel for sync (no weight)
  // cout<<"print number of event selected after 4ljj sel for sync (no weight) ... "<<endl; 
  // ofstream f_events4ljjsel_sync;
  // TString f_events4ljjsel_sync_name = "nEvents_4ljjsel_perSync_"+ sYear + ".txt";
  // f_events4ljjsel_sync.open(f_events4ljjsel_sync_name);
  // f_events4ljjsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  // for(int fs=0; fs<nFinalStates+1; fs++){
  //   f_events4ljjsel_sync<<"|"<<sFinalState[fs]<<" |"<<nEvent_4ljjsel[HH][fs]<<" |"<<nEvent_4ljjsel[TTZ][fs]<<" |"<<nEvent_4ljjsel[ttH][fs]<<" |"<<nEvent_4ljjsel[qqZZ][fs]+nEvent_4ljjsel[ggZZ][fs]<<" |"<<nEvent_4ljjsel[ggH][fs]+nEvent_4ljjsel[VBF][fs]<<" |"<<nEvent_4ljjsel[VVV][fs]+nEvent_4ljjsel[VH][fs]+nEvent_4ljjsel[TTW][fs]<<" |"<<nEvent_4ljjsel[ZXbkg][fs]<<" |"<<endl;
  // }
  // f_events4ljjsel_sync.close();



  // ********************
  // *** PRINT YIELDS ***
  // ********************
  // --- print yields after 4l sel for sync 
  cout<<"print yields after 4l sel for sync ... "<<endl;
  ofstream f_yields4lsel_sync;
  TString f_yields4lsel_sync_name = "yields4lsel_perSync_"+ sYear + ".txt";
  f_yields4lsel_sync.open(f_yields4lsel_sync_name);
  f_yields4lsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4lsel_sync<<"|"<<sFinalState[fs]<<" |"<<yield_4lsel[HH][fs]<<" |"<<yield_4lsel[TTZ][fs]<<" |"<<yield_4lsel[ttH][fs]<<" |"<<yield_4lsel[qqZZ][fs]+yield_4lsel[ggZZ][fs]<<" |"<<yield_4lsel[ggH][fs]+yield_4lsel[VBF][fs]<<" |"<<yield_4lsel[VVV][fs]+yield_4lsel[VH][fs]+yield_4lsel[TTW][fs]<<" |"<<yield_4lsel[ZXbkg][fs]<<" |"<<endl;
  }
  f_yields4lsel_sync.close();

  // --- print yields after 4l sel per process
  cout<<"print yields after 4l sel per process ... "<<endl;
  ofstream f_yields4lsel_process;
  TString f_yields4lsel_process_name = "yields4lsel_perProcess_"+ sYear + ".txt";
  f_yields4lsel_process.open(f_yields4lsel_process_name);
  f_yields4lsel_process<<"|Process |"<<sFinalState[fs_4mu]<<" |"<<sFinalState[fs_4e]<<" |"<<sFinalState[fs_2e2mu]<<" |"<<endl;
  for(int pr=0; pr<nProcesses; pr++){
    f_yields4lsel_process<<"|"<<sProcess[pr]<<" |"<<yield_4lsel[pr][fs_4mu]<<" |"<<yield_4lsel[pr][fs_4e]<<" |"<<yield_4lsel[pr][fs_2e2mu]<<" |"<<endl;
  }
  f_yields4lsel_process.close();


  // --- print yields after 4ljj sel for sync 
  cout<<"print yields after 4ljj sel for sync ... "<<endl;
  ofstream f_yields4ljjsel_sync;
  TString f_yields4ljjsel_sync_name = "yields4ljjsel_perSync_"+ sYear + ".txt";
  f_yields4ljjsel_sync.open(f_yields4ljjsel_sync_name);
  f_yields4ljjsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel_sync<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel[HH][fs]<<" |"<<yield_4ljjsel[TTZ][fs]<<" |"<<yield_4ljjsel[ttH][fs]<<" |"<<yield_4ljjsel[qqZZ][fs]+yield_4ljjsel[ggZZ][fs]<<" |"<<yield_4ljjsel[ggH][fs]+yield_4ljjsel[VBF][fs]<<" |"<<yield_4ljjsel[VVV][fs]+yield_4ljjsel[VH][fs]+yield_4ljjsel[TTW][fs]<<" |"<<yield_4ljjsel[ZXbkg][fs]<<" |"<<endl;
  }
  f_yields4ljjsel_sync.close();

  // --- print yields after 4ljj sel per process
  cout<<"print yields after 4ljj sel per process ... "<<endl;
  ofstream f_yields4ljjsel_process;
  TString f_yields4ljjsel_process_name = "yields4ljjsel_perProcess_"+ sYear + ".txt";
  f_yields4ljjsel_process.open(f_yields4ljjsel_process_name);
  f_yields4ljjsel_process<<"|Process |"<<sFinalState[fs_4mu]<<" |"<<sFinalState[fs_4e]<<" |"<<sFinalState[fs_2e2mu]<<" |"<<endl;
  for(int pr=0; pr<nProcesses; pr++){
    f_yields4ljjsel_process<<"|"<<sProcess[pr]<<" |"<<yield_4ljjsel[pr][fs_4mu]<<" |"<<yield_4ljjsel[pr][fs_4e]<<" |"<<yield_4ljjsel[pr][fs_2e2mu]<<" |"<<endl;
  }
  f_yields4ljjsel_process.close();



  // ******************************
  // *** PRINT YIELDS SIDEBANDS ***
  // ******************************
  // --- print yields after 4ljj sel in sidebands 
  cout<<"print yields after 4ljj sel in sidebands ... "<<endl;
  ofstream f_yields4ljjsel_sidebands;
  TString f_yields4ljjsel_sidebands_name = "yields4ljjsel_sidebands_"+ sYear + ".txt";

  TH1F* yield_4ljjsel_sidebands_ZZ[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebands_Higgs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebands_others[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebands_tot_Ange[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebands_SMHiggs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebands_TTV[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebands_tot_Ale[nFinalStates+1];

  f_yields4ljjsel_sidebands.open(f_yields4ljjsel_sidebands_name); 

  // print yields in Angela style
  f_yields4ljjsel_sidebands<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |Sum bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // ZZ
    yield_4ljjsel_sidebands_ZZ[fs] = (TH1F*)yield_4ljjsel_sidebands[qqZZ][fs]->Clone("hYields_4ljjsel_sidebands_ZZ_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebands_ZZ[fs]->Add(yield_4ljjsel_sidebands[ggZZ][fs]); 
    // Higgs
    yield_4ljjsel_sidebands_Higgs[fs] = (TH1F*)yield_4ljjsel_sidebands[ggH][fs]->Clone("hYields_4ljjsel_sidebands_Higgs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebands_Higgs[fs]->Add(yield_4ljjsel_sidebands[VBF][fs]);
    // others
    yield_4ljjsel_sidebands_others[fs] = (TH1F*)yield_4ljjsel_sidebands[VVV][fs]->Clone("hYields_4ljjsel_sidebands_others_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebands_others[fs]->Add(yield_4ljjsel_sidebands[VH][fs]);
    yield_4ljjsel_sidebands_others[fs]->Add(yield_4ljjsel_sidebands[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebands_tot_Ange[fs] = (TH1F*)yield_4ljjsel_sidebands[TTZ][fs]->Clone("hYields_4ljjsel_sidebands_totAnge_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[ttH][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[qqZZ][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[ggZZ][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[ggH][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[VBF][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[VVV][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[VH][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[TTW][fs]);
    yield_4ljjsel_sidebands_tot_Ange[fs]->Add(yield_4ljjsel_sidebands[ZXbkg][fs]);

    // print
    f_yields4ljjsel_sidebands<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebands[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands[TTZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[TTZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands[ttH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[ttH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands_ZZ[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands_ZZ[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands_Higgs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands_Higgs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands_others[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands_others[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[ZXbkg][fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebands_tot_Ange[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands_tot_Ange[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebands[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }




  f_yields4ljjsel_sidebands<<" "<<endl;
  f_yields4ljjsel_sidebands<<" "<<endl;
  f_yields4ljjsel_sidebands<<" "<<endl;
  f_yields4ljjsel_sidebands<<" "<<endl;

  // print yields in Ale style
  f_yields4ljjsel_sidebands<<"|Final state |signal HH |SM Higgs |qqZZ |ggZZ |TTV |Z+X |VVV |sum all bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // SM Higgs
    yield_4ljjsel_sidebands_SMHiggs[fs] = (TH1F*)yield_4ljjsel_sidebands[ggH][fs]->Clone("hYields_4ljjsel_sidebands_SMHiggs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebands_SMHiggs[fs]->Add(yield_4ljjsel_sidebands[VBF][fs]);
    yield_4ljjsel_sidebands_SMHiggs[fs]->Add(yield_4ljjsel_sidebands[VH][fs]);
    yield_4ljjsel_sidebands_SMHiggs[fs]->Add(yield_4ljjsel_sidebands[ttH][fs]);
    yield_4ljjsel_sidebands_SMHiggs[fs]->Add(yield_4ljjsel_sidebands[bbH][fs]);
    // TTV
    yield_4ljjsel_sidebands_TTV[fs] = (TH1F*)yield_4ljjsel_sidebands[TTZ][fs]->Clone("hYields_4ljjsel_sidebands_TTV_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebands_TTV[fs]->Add(yield_4ljjsel_sidebands[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebands_tot_Ale[fs] = (TH1F*)yield_4ljjsel_sidebands[ggH][fs]->Clone("hYields_4ljjsel_sidebands_totAle_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[VBF][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[VH][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[ttH][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[bbH][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[qqZZ][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[ggZZ][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[TTZ][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[TTW][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[ZXbkg][fs]);
    yield_4ljjsel_sidebands_tot_Ale[fs]->Add(yield_4ljjsel_sidebands[VVV][fs]);


    f_yields4ljjsel_sidebands<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebands[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands_SMHiggs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands_SMHiggs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands[qqZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[qqZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands[ggZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[ggZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands_TTV[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands_TTV[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[ZXbkg][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands[VVV][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[VVV][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebands_tot_Ale[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands_tot_Ale[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebands[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebands[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }
  f_yields4ljjsel_sidebands.close();
  



  // *******************************
  // *** PRINT YIELDS SIDEBANDSX ***
  // *******************************
  // --- print yields after 4ljj sel in sidebandSX 
  cout<<"print yields after 4ljj sel in sidebandSX ... "<<endl;
  ofstream f_yields4ljjsel_sidebandSX;
  TString f_yields4ljjsel_sidebandSX_name = "yields4ljjsel_sidebandSX_"+ sYear + ".txt";

  TH1F* yield_4ljjsel_sidebandSX_ZZ[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSX_Higgs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSX_others[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSX_tot_Ange[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSX_SMHiggs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSX_TTV[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSX_tot_Ale[nFinalStates+1];

  f_yields4ljjsel_sidebandSX.open(f_yields4ljjsel_sidebandSX_name); 

  // print yields in Angela style
  f_yields4ljjsel_sidebandSX<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |Sum bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // ZZ
    yield_4ljjsel_sidebandSX_ZZ[fs] = (TH1F*)yield_4ljjsel_sidebandSX[qqZZ][fs]->Clone("hYields_4ljjsel_sidebandSX_ZZ_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSX_ZZ[fs]->Add(yield_4ljjsel_sidebandSX[ggZZ][fs]); 
    // Higgs
    yield_4ljjsel_sidebandSX_Higgs[fs] = (TH1F*)yield_4ljjsel_sidebandSX[ggH][fs]->Clone("hYields_4ljjsel_sidebandSX_Higgs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSX_Higgs[fs]->Add(yield_4ljjsel_sidebandSX[VBF][fs]);
    // others
    yield_4ljjsel_sidebandSX_others[fs] = (TH1F*)yield_4ljjsel_sidebandSX[VVV][fs]->Clone("hYields_4ljjsel_sidebandSX_others_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSX_others[fs]->Add(yield_4ljjsel_sidebandSX[VH][fs]);
    yield_4ljjsel_sidebandSX_others[fs]->Add(yield_4ljjsel_sidebandSX[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebandSX_tot_Ange[fs] = (TH1F*)yield_4ljjsel_sidebandSX[TTZ][fs]->Clone("hYields_4ljjsel_sidebandSX_totAnge_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[ttH][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[qqZZ][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[ggZZ][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[ggH][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[VBF][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[VVV][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[VH][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[TTW][fs]);
    yield_4ljjsel_sidebandSX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSX[ZXbkg][fs]);

    // print
    f_yields4ljjsel_sidebandSX<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebandSX[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX[TTZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[TTZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX[ttH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[ttH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX_ZZ[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX_ZZ[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX_Higgs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX_Higgs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX_others[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX_others[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[ZXbkg][fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandSX_tot_Ange[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX_tot_Ange[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandSX[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }




  f_yields4ljjsel_sidebandSX<<" "<<endl;
  f_yields4ljjsel_sidebandSX<<" "<<endl;
  f_yields4ljjsel_sidebandSX<<" "<<endl;
  f_yields4ljjsel_sidebandSX<<" "<<endl;

  // print yields in Ale style
  f_yields4ljjsel_sidebandSX<<"|Final state |signal HH |SM Higgs |qqZZ |ggZZ |TTV |Z+X |VVV |sum all bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // SM Higgs
    yield_4ljjsel_sidebandSX_SMHiggs[fs] = (TH1F*)yield_4ljjsel_sidebandSX[ggH][fs]->Clone("hYields_4ljjsel_sidebandSX_SMHiggs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSX[VBF][fs]);
    yield_4ljjsel_sidebandSX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSX[VH][fs]);
    yield_4ljjsel_sidebandSX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSX[ttH][fs]);
    yield_4ljjsel_sidebandSX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSX[bbH][fs]);
    // TTV
    yield_4ljjsel_sidebandSX_TTV[fs] = (TH1F*)yield_4ljjsel_sidebandSX[TTZ][fs]->Clone("hYields_4ljjsel_sidebandSX_TTV_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSX_TTV[fs]->Add(yield_4ljjsel_sidebandSX[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebandSX_tot_Ale[fs] = (TH1F*)yield_4ljjsel_sidebandSX[ggH][fs]->Clone("hYields_4ljjsel_sidebandSX_totAle_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[VBF][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[VH][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[ttH][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[bbH][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[qqZZ][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[ggZZ][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[TTZ][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[TTW][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[ZXbkg][fs]);
    yield_4ljjsel_sidebandSX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSX[VVV][fs]);


    f_yields4ljjsel_sidebandSX<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebandSX[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX_SMHiggs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX_SMHiggs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX[qqZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[qqZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX[ggZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[ggZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX_TTV[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX_TTV[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[ZXbkg][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX[VVV][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[VVV][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSX_tot_Ale[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX_tot_Ale[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandSX[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSX[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }
  f_yields4ljjsel_sidebandSX.close();




  // *******************************
  // *** PRINT YIELDS SIDEBANDDX ***
  // *******************************
  // --- print yields after 4ljj sel in sidebandDX 
  cout<<"print yields after 4ljj sel in sidebandDX ... "<<endl;
  ofstream f_yields4ljjsel_sidebandDX;
  TString f_yields4ljjsel_sidebandDX_name = "yields4ljjsel_sidebandDX_"+ sYear + ".txt";

  TH1F* yield_4ljjsel_sidebandDX_ZZ[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandDX_Higgs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandDX_others[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandDX_tot_Ange[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandDX_SMHiggs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandDX_TTV[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandDX_tot_Ale[nFinalStates+1];

  f_yields4ljjsel_sidebandDX.open(f_yields4ljjsel_sidebandDX_name); 

  // print yields in Angela style
  f_yields4ljjsel_sidebandDX<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |Sum bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // ZZ
    yield_4ljjsel_sidebandDX_ZZ[fs] = (TH1F*)yield_4ljjsel_sidebandDX[qqZZ][fs]->Clone("hYields_4ljjsel_sidebandDX_ZZ_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandDX_ZZ[fs]->Add(yield_4ljjsel_sidebandDX[ggZZ][fs]); 
    // Higgs
    yield_4ljjsel_sidebandDX_Higgs[fs] = (TH1F*)yield_4ljjsel_sidebandDX[ggH][fs]->Clone("hYields_4ljjsel_sidebandDX_Higgs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandDX_Higgs[fs]->Add(yield_4ljjsel_sidebandDX[VBF][fs]);
    // others
    yield_4ljjsel_sidebandDX_others[fs] = (TH1F*)yield_4ljjsel_sidebandDX[VVV][fs]->Clone("hYields_4ljjsel_sidebandDX_others_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandDX_others[fs]->Add(yield_4ljjsel_sidebandDX[VH][fs]);
    yield_4ljjsel_sidebandDX_others[fs]->Add(yield_4ljjsel_sidebandDX[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebandDX_tot_Ange[fs] = (TH1F*)yield_4ljjsel_sidebandDX[TTZ][fs]->Clone("hYields_4ljjsel_sidebandDX_totAnge_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[ttH][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[qqZZ][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[ggZZ][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[ggH][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[VBF][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[VVV][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[VH][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[TTW][fs]);
    yield_4ljjsel_sidebandDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandDX[ZXbkg][fs]);

    // print
    f_yields4ljjsel_sidebandDX<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebandDX[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX[TTZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[TTZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX[ttH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[ttH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX_ZZ[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX_ZZ[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX_Higgs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX_Higgs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX_others[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX_others[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[ZXbkg][fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandDX_tot_Ange[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX_tot_Ange[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandDX[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }




  f_yields4ljjsel_sidebandDX<<" "<<endl;
  f_yields4ljjsel_sidebandDX<<" "<<endl;
  f_yields4ljjsel_sidebandDX<<" "<<endl;
  f_yields4ljjsel_sidebandDX<<" "<<endl;

  // print yields in Ale style
  f_yields4ljjsel_sidebandDX<<"|Final state |signal HH |SM Higgs |qqZZ |ggZZ |TTV |Z+X |VVV |sum all bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // SM Higgs
    yield_4ljjsel_sidebandDX_SMHiggs[fs] = (TH1F*)yield_4ljjsel_sidebandDX[ggH][fs]->Clone("hYields_4ljjsel_sidebandDX_SMHiggs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandDX[VBF][fs]);
    yield_4ljjsel_sidebandDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandDX[VH][fs]);
    yield_4ljjsel_sidebandDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandDX[ttH][fs]);
    yield_4ljjsel_sidebandDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandDX[bbH][fs]);
    // TTV
    yield_4ljjsel_sidebandDX_TTV[fs] = (TH1F*)yield_4ljjsel_sidebandDX[TTZ][fs]->Clone("hYields_4ljjsel_sidebandDX_TTV_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandDX_TTV[fs]->Add(yield_4ljjsel_sidebandDX[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebandDX_tot_Ale[fs] = (TH1F*)yield_4ljjsel_sidebandDX[ggH][fs]->Clone("hYields_4ljjsel_sidebandDX_totAle_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[VBF][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[VH][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[ttH][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[bbH][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[qqZZ][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[ggZZ][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[TTZ][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[TTW][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[ZXbkg][fs]);
    yield_4ljjsel_sidebandDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandDX[VVV][fs]);


    f_yields4ljjsel_sidebandDX<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebandDX[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX_SMHiggs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX_SMHiggs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX[qqZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[qqZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX[ggZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[ggZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX_TTV[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX_TTV[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[ZXbkg][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX[VVV][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[VVV][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandDX_tot_Ale[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX_tot_Ale[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandDX[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandDX[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }
  f_yields4ljjsel_sidebandDX.close();



  // *******************************
  // *** PRINT YIELDS SIDEBANDDX ***
  // *******************************
  // --- print yields after 4ljj sel in sidebandSXDX 
  cout<<"print yields after 4ljj sel in sidebandSXDX ... "<<endl;
  ofstream f_yields4ljjsel_sidebandSXDX;
  TString f_yields4ljjsel_sidebandSXDX_name = "yields4ljjsel_sidebandSXDX_"+ sYear + ".txt";

  TH1F* yield_4ljjsel_sidebandSXDX_ZZ[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSXDX_Higgs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSXDX_others[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSXDX_tot_Ange[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSXDX_SMHiggs[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSXDX_TTV[nFinalStates+1];
  TH1F* yield_4ljjsel_sidebandSXDX_tot_Ale[nFinalStates+1];

  f_yields4ljjsel_sidebandSXDX.open(f_yields4ljjsel_sidebandSXDX_name); 

  // print yields in Angela style
  f_yields4ljjsel_sidebandSXDX<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF) |others(=VVV+VH+TTW) |Z+X |Sum bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // ZZ
    yield_4ljjsel_sidebandSXDX_ZZ[fs] = (TH1F*)yield_4ljjsel_sidebandSXDX[qqZZ][fs]->Clone("hYields_4ljjsel_sidebandSXDX_ZZ_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSXDX_ZZ[fs]->Add(yield_4ljjsel_sidebandSXDX[ggZZ][fs]); 
    // Higgs
    yield_4ljjsel_sidebandSXDX_Higgs[fs] = (TH1F*)yield_4ljjsel_sidebandSXDX[ggH][fs]->Clone("hYields_4ljjsel_sidebandSXDX_Higgs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSXDX_Higgs[fs]->Add(yield_4ljjsel_sidebandSXDX[VBF][fs]);
    // others
    yield_4ljjsel_sidebandSXDX_others[fs] = (TH1F*)yield_4ljjsel_sidebandSXDX[VVV][fs]->Clone("hYields_4ljjsel_sidebandSXDX_others_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSXDX_others[fs]->Add(yield_4ljjsel_sidebandSXDX[VH][fs]);
    yield_4ljjsel_sidebandSXDX_others[fs]->Add(yield_4ljjsel_sidebandSXDX[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs] = (TH1F*)yield_4ljjsel_sidebandSXDX[TTZ][fs]->Clone("hYields_4ljjsel_sidebandSXDX_totAnge_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[ttH][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[qqZZ][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[ggZZ][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[ggH][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[VBF][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[VVV][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[VH][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[TTW][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->Add(yield_4ljjsel_sidebandSXDX[ZXbkg][fs]);

    // print
    f_yields4ljjsel_sidebandSXDX<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebandSXDX[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX[TTZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[TTZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX[ttH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[ttH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX_ZZ[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX_ZZ[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX_Higgs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX_Higgs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX_others[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX_others[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[ZXbkg][fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX_tot_Ange[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandSXDX[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }


  f_yields4ljjsel_sidebandSXDX<<" "<<endl;
  f_yields4ljjsel_sidebandSXDX<<" "<<endl;
  f_yields4ljjsel_sidebandSXDX<<" "<<endl;
  f_yields4ljjsel_sidebandSXDX<<" "<<endl;

  // print yields in Ale style
  f_yields4ljjsel_sidebandSXDX<<"|Final state |signal HH |SM Higgs |qqZZ |ggZZ |TTV |Z+X |VVV |sum all bkg |Data |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    // add histos
    // SM Higgs
    yield_4ljjsel_sidebandSXDX_SMHiggs[fs] = (TH1F*)yield_4ljjsel_sidebandSXDX[ggH][fs]->Clone("hYields_4ljjsel_sidebandSXDX_SMHiggs_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSXDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSXDX[VBF][fs]);
    yield_4ljjsel_sidebandSXDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSXDX[VH][fs]);
    yield_4ljjsel_sidebandSXDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSXDX[ttH][fs]);
    yield_4ljjsel_sidebandSXDX_SMHiggs[fs]->Add(yield_4ljjsel_sidebandSXDX[bbH][fs]);
    // TTV
    yield_4ljjsel_sidebandSXDX_TTV[fs] = (TH1F*)yield_4ljjsel_sidebandSXDX[TTZ][fs]->Clone("hYields_4ljjsel_sidebandSXDX_TTV_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSXDX_TTV[fs]->Add(yield_4ljjsel_sidebandSXDX[TTW][fs]);
    // tot histo
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs] = (TH1F*)yield_4ljjsel_sidebandSXDX[ggH][fs]->Clone("hYields_4ljjsel_sidebandSXDX_totAle_"+sFinalState[fs]+"_"+sYear);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[VBF][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[VH][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[ttH][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[bbH][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[qqZZ][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[ggZZ][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[TTZ][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[TTW][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[ZXbkg][fs]);
    yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->Add(yield_4ljjsel_sidebandSXDX[VVV][fs]);


    f_yields4ljjsel_sidebandSXDX<<" |"<<sFinalState[fs]
			     <<" |"<<yield_4ljjsel_sidebandSXDX[HH][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[HH][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX_SMHiggs[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX_SMHiggs[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX[qqZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[qqZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX[ggZZ][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[ggZZ][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX_TTV[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX_TTV[fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX[ZXbkg][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[ZXbkg][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX[VVV][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[VVV][fs]->GetBinError(1)
			     <<" |"<<yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX_tot_Ale[fs]->GetBinError(1)
                             <<" |"<<yield_4ljjsel_sidebandSXDX[Data][fs]->GetBinContent(1)<<" +- "<<yield_4ljjsel_sidebandSXDX[Data][fs]->GetBinError(1)
                             <<" |"<<endl;
  }
  f_yields4ljjsel_sidebandSXDX.close();
    


}// end function printyields_forSync



//*************************************
//*** printYields_forCards function ***
//*************************************
void  printYields_forCards(){

  cout<<"print yields for cards ..."<<endl;

  //---input path
  TString sYear;
  if(year==2016)      sYear = "2016";
  else if(year==2017) sYear = "2017";
  else if(year==2018) sYear = "2018";
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;

 
  // retrieve yields histos from file
  TString inFileName = "f_yields_" + sYear + ".root";     
  cout<<"Retrieving Data and MC yields histos from file "<<inFileName<<" ..."<<endl;
  TFile* fInYields = TFile::Open(inFileName);

  TH1F* hTemp1;
  Float_t yield_4ljjsel [nProcesses][nFinalStates+1];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      hTemp1 = (TH1F*)fInYields->Get("hYields_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
      yield_4ljjsel[pr][fs] = hTemp1->GetBinContent(1);
    }
  }

  //da fare per FS (only 4mu, 4e, 2e2mu)
  cout<<"print yields after 4ljj sel per cards ... "<<endl;
  for(int fs=0; fs<nFinalStates; fs++){
    ofstream f_yields4ljjsel_cards;
    TString f_yields4ljjsel_cards_name = "yields4ljjsel_perCards_"+ sYear + "_" + sFinalState[fs] + ".yaml";
    f_yields4ljjsel_cards.open(f_yields4ljjsel_cards_name);
    f_yields4ljjsel_cards<<"year"<<sYear<<": "<<endl;
    for(int pr=1; pr<nProcesses; pr++){   //pr starts from 1=HH, not to print Data yield (0)
      f_yields4ljjsel_cards<<"    "<<sProcess[pr]<<": '"<<yield_4ljjsel[pr][fs]<<"'"<<endl;
    }
    f_yields4ljjsel_cards.close();
  }

  


}// end function printYields_forCards



//*********************************
//*** doPlots_inputBDT function ***
//*********************************
void doPlots_inputBDT(){

 cout<<"do BDT input plots ..."<<endl;

  //---input path
  TString sYear;
  TString lumiText;
  if(year==2016){
      sYear    = "2016";
      lumiText = "35.8 fb^{-1}";
  }
  else if(year==2017){
      sYear    = "2017";
      lumiText = "41.5 fb^{-1}";
  }
  else if(year==2018){
      sYear    = "2018";
      lumiText = "59.7 fb^{-1}";
  }
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;


  TString outPath_inputBDTplots = "plots_inputBDT_" + sYear;
  cout<<"creating output dir "<<outPath_inputBDTplots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_inputBDTplots)).c_str()); // create output dir



  // retrieve yields histos from file
  TString inFileName = "f_histos_h1_4ljjsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  // (BDT input histos)
  Int_t nBDTinputHistos = 8;
  TString sBDTInputNames[] = {
    "pT4l_4ljjsel",
    "j1btag_4ljjsel",
    "j2btag_4ljjsel",
    "j1pT_4ljjsel",
    "j2pT_4ljjsel",
    "MET_4ljjsel",    
    "DeltaRhh_4ljjsel",
    "mbb_4ljjsel",
  };
  TH1F* h1_BDTinput_4ljjsel[nBDTinputHistos][nProcesses][nFinalStates+1];
  for(int bdtIn=0; bdtIn<nBDTinputHistos; bdtIn++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates+1; fs++){
        h1_BDTinput_4ljjsel[bdtIn][pr][fs] = (TH1F*)fInhistos->Get("h1_"+sBDTInputNames[bdtIn]+"_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
        cout<<h1_BDTinput_4ljjsel[bdtIn][pr][fs]->GetName()<<endl;
      }
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_BDTinput_4ljjsel     [nBDTinputHistos][nFinalStates+1];
  THStack* hs_BDTinput_4ljjsel    [nBDTinputHistos][nFinalStates+1];
  TPad*    pad1_BDTinput_4ljjsel  [nBDTinputHistos][nFinalStates+1];
  TLegend* leg_BDTinput_4ljjsel   [nBDTinputHistos][nFinalStates+1];
  TH1F*    hMCtot_BDTinput_4ljjsel[nBDTinputHistos][nFinalStates+1];
  TPad*    pad2_BDTinput_4ljjsel  [nBDTinputHistos][nFinalStates+1];
  TH1F*    rp_BDTinput_4ljjsel    [nBDTinputHistos][nFinalStates+1];
  TH1F*    hUncMC_BDTinput_4ljjsel[nBDTinputHistos][nFinalStates+1];
  
  //BDT input plots
  for(int bdtIn=0; bdtIn<nBDTinputHistos; bdtIn++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // canvas
      c_BDTinput_4ljjsel[bdtIn][fs] = new TCanvas("c_InputBDT_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]+"_"+sYear,"c_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]+"_"+sYear,800,800);
      // hstack
      hs_BDTinput_4ljjsel[bdtIn][fs] = new THStack("hs_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs],"");
      // VVV process
      h1_BDTinput_4ljjsel[bdtIn][VVV][fs]->SetFillColor(kGreen-3);
      h1_BDTinput_4ljjsel[bdtIn][VVV][fs]->SetLineColor(kGreen-1);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VVV][fs]); //add to hs
      // TTV process: TTW + TTV
      h1_BDTinput_4ljjsel[bdtIn][TTW][fs]->SetFillColor(kBlue+3);
      h1_BDTinput_4ljjsel[bdtIn][TTW][fs]->SetLineColor(kBlue+3);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTW][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]->SetFillColor(kBlue+3);
      h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]->SetLineColor(kBlue+3);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]); //add to hs
      // ggZZ process
      h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]->SetFillColor(kAzure-3);
      h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]->SetLineColor(kBlue+2);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]); //add to hs
      // qqZZ process
      h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]->SetFillColor(kAzure+6);
      h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]->SetLineColor(kAzure-6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]); //add to hs
      // SM Higgs processes: ggH + VBF + VH + ttH + bbH
      h1_BDTinput_4ljjsel[bdtIn][ggH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][ggH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ggH][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][VBF][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][VBF][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VBF][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][VH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][VH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VH][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][ttH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][ttH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ttH][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][bbH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][bbH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][bbH][fs]); //add to hs
      // HH signal
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->SetLineColor(kRed);
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->SetLineWidth(2);
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->Scale(100.);
      // data
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->SetMarkerColor(kBlack);
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->SetLineColor(kBlack);
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->SetMarkerStyle(20);
    
      // --- upper plot pad
      pad1_BDTinput_4ljjsel[bdtIn][fs] = new TPad("pad1_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs],"pad1_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs], 0, 0.3, 1, 1.0);
      pad1_BDTinput_4ljjsel[bdtIn][fs]->Draw();
      pad1_BDTinput_4ljjsel[bdtIn][fs]->cd();
    
      hs_BDTinput_4ljjsel[bdtIn][fs]->SetMaximum(10e04);
      hs_BDTinput_4ljjsel[bdtIn][fs]->SetMinimum(10e-04);
    
      hs_BDTinput_4ljjsel[bdtIn][fs]->Draw("histo");
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->Draw("histosame");
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->Draw("samepe");

      hs_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelFont(43);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelSize(15);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitle(h1_BDTinput_4ljjsel[bdtIn][HH][fs]->GetXaxis()->GetTitle());
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleSize(20);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleFont(43);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleOffset(1.4);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelFont(43);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelSize(15);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitle(h1_BDTinput_4ljjsel[bdtIn][HH][fs]->GetYaxis()->GetTitle());

      // --- legend
      leg_BDTinput_4ljjsel[bdtIn][fs] = new TLegend(0.78,0.61,0.94,0.87);
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][Data][fs], "Data",          "lp");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][HH][fs],   "HH->4lbb x100", "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][ggH][fs],  "SM Higgs",      "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs], "qq->ZZ",        "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs], "gg->ZZ",        "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][TTZ][fs],  "TTV; V=Z,W",    "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][VVV][fs],  "VVV; V=Z,W",    "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->SetFillColor(kWhite);
      leg_BDTinput_4ljjsel[bdtIn][fs]->SetLineColor(kBlack);
      leg_BDTinput_4ljjsel[bdtIn][fs]->SetTextFont(43);
      leg_BDTinput_4ljjsel[bdtIn][fs]->Draw();

      c_BDTinput_4ljjsel[bdtIn][fs]->Update();

      pad1_BDTinput_4ljjsel[bdtIn][fs]->SetLogy();

      c_BDTinput_4ljjsel[bdtIn][fs]->Update();

      // --- tot hist for all MC
      hMCtot_BDTinput_4ljjsel[bdtIn][fs] = (TH1F*)h1_BDTinput_4ljjsel[bdtIn][ggH][fs]->Clone("hMCtot_BDTinput_4ljjsel_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VBF][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VH][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ttH][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][bbH][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTW][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VVV][fs]);
      
      // --- lower pad plot
      c_BDTinput_4ljjsel[bdtIn][fs]->cd();
      pad2_BDTinput_4ljjsel[bdtIn][fs] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
      pad2_BDTinput_4ljjsel[bdtIn][fs]->SetGridy();
      pad2_BDTinput_4ljjsel[bdtIn][fs]->Draw();
      pad2_BDTinput_4ljjsel[bdtIn][fs]->cd();

      // --- define ratio plot
      rp_BDTinput_4ljjsel[bdtIn][fs] = (TH1F*)h1_BDTinput_4ljjsel[bdtIn][Data][fs]->Clone("rp_BDTinput_4ljjsel_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetLineColor(kBlack);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMinimum(0.);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMaximum(2.);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetStats(0);
      rp_BDTinput_4ljjsel[bdtIn][fs]->Divide(hMCtot_BDTinput_4ljjsel[bdtIn][fs]); //divide histo rp/MC
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerStyle(20);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerColor(kBlack);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetTitle("");

      rp_BDTinput_4ljjsel[bdtIn][fs]->SetYTitle("Data/#Sigma bkg");
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetNdivisions(505);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleSize(20);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleOffset(1.4);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelSize(15);

      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitleSize(20);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitleFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitleOffset(4.);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelSize(15);

      // --- define mc shadow unc plot
      hUncMC_BDTinput_4ljjsel[bdtIn][fs] = (TH1F*)hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Clone("hUncMC_BDTinput_4ljjsel_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]);
      for(int xbin=1; xbin < hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->GetNbins() + 1; xbin++){
        float err = 0.;
        if(hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetBinContent(xbin) == 0) continue;
        err = hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetBinError(xbin) / hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetBinContent(xbin);
        hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetBinContent(xbin, 1.);
        hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetBinError(xbin, err);
      }
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetLineColor(1);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetFillStyle(3005);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetFillColor(kGray+3);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerColor(1);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerStyle(1);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetTitle("");
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetStats(0);

      // ---draw
      rp_BDTinput_4ljjsel[bdtIn][fs]->Draw("ep");
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->Draw("e2 same");

      c_BDTinput_4ljjsel[bdtIn][fs]->Update();

      // --- draw CMS and lumi text
      writeExtraText = true;
      extraText      = "Preliminary";
      lumi_sqrtS     = lumiText + " (13 TeV)";
      cmsTextSize    = 0.6;
      lumiTextSize   = 0.46;
      extraOverCmsTextSize = 0.75;
      relPosX = 0.12;
      CMS_lumi(pad1_BDTinput_4ljjsel[bdtIn][fs], 0, 0);


      c_BDTinput_4ljjsel[bdtIn][fs]->SaveAs(outPath_inputBDTplots + "/" + c_BDTinput_4ljjsel[bdtIn][fs]->GetName() + ".png");
      c_BDTinput_4ljjsel[bdtIn][fs]->SaveAs(outPath_inputBDTplots + "/" + c_BDTinput_4ljjsel[bdtIn][fs]->GetName() + ".pdf");
      

    }
  }

} // end doPlots_inputBDT function





//****************************************
//*** doPlots_inputBDT_withZX function ***
//****************************************
void doPlots_inputBDT_withZX(){

 cout<<"do BDT input plots ..."<<endl;

  //---input path
  TString sYear;
  TString lumiText;
  if(year==2016){
      sYear    = "2016";
      lumiText = "35.8 fb^{-1}";
  }
  else if(year==2017){
      sYear    = "2017";
      lumiText = "41.5 fb^{-1}";
  }
  else if(year==2018){
      sYear    = "2018";
      lumiText = "59.7 fb^{-1}";
  }
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;


  TString outPath_inputBDTwithZXplots = "plots_inputBDT_withZX_" + sYear;
  cout<<"creating output dir "<<outPath_inputBDTwithZXplots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_inputBDTwithZXplots)).c_str()); // create output dir



  // retrieve yields histos from file
  TString inFileName = "f_histos_h1_4ljjsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  // (BDT input histos)
  Int_t nBDTinputHistos = 8;
  TString sBDTInputNames[] = {
    "pT4l_4ljjsel",
    "j1btag_4ljjsel",
    "j2btag_4ljjsel",
    "j1pT_4ljjsel",
    "j2pT_4ljjsel",
    "MET_4ljjsel",    
    "DeltaRhh_4ljjsel",
    "mbb_4ljjsel",
  };
  TH1F* h1_BDTinput_4ljjsel[nBDTinputHistos][nProcesses][nFinalStates+1];
  for(int bdtIn=0; bdtIn<nBDTinputHistos; bdtIn++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates+1; fs++){
        h1_BDTinput_4ljjsel[bdtIn][pr][fs] = (TH1F*)fInhistos->Get("h1_"+sBDTInputNames[bdtIn]+"_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
        cout<<h1_BDTinput_4ljjsel[bdtIn][pr][fs]->GetName()<<endl;
      }
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_BDTinput_4ljjsel     [nBDTinputHistos][nFinalStates+1];
  THStack* hs_BDTinput_4ljjsel    [nBDTinputHistos][nFinalStates+1];
  TPad*    pad1_BDTinput_4ljjsel  [nBDTinputHistos][nFinalStates+1];
  TLegend* leg_BDTinput_4ljjsel   [nBDTinputHistos][nFinalStates+1];
  TH1F*    hMCtot_BDTinput_4ljjsel[nBDTinputHistos][nFinalStates+1];
  TPad*    pad2_BDTinput_4ljjsel  [nBDTinputHistos][nFinalStates+1];
  TH1F*    rp_BDTinput_4ljjsel    [nBDTinputHistos][nFinalStates+1];
  TH1F*    hUncMC_BDTinput_4ljjsel[nBDTinputHistos][nFinalStates+1];
  
  //BDT input plots
  for(int bdtIn=0; bdtIn<nBDTinputHistos; bdtIn++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // canvas
      c_BDTinput_4ljjsel[bdtIn][fs] = new TCanvas("c_InputBDT_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]+"_"+sYear,"c_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]+"_"+sYear,800,800);
      // hstack
      hs_BDTinput_4ljjsel[bdtIn][fs] = new THStack("hs_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs],"");
      // VVV process
      h1_BDTinput_4ljjsel[bdtIn][VVV][fs]->SetFillColor(kGreen-3);
      h1_BDTinput_4ljjsel[bdtIn][VVV][fs]->SetLineColor(kGreen-1);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VVV][fs]); //add to hs
      // Z+X process
      h1_BDTinput_4ljjsel[bdtIn][ZXbkg][fs]->SetFillColor(kGreen+3);
      h1_BDTinput_4ljjsel[bdtIn][ZXbkg][fs]->SetLineColor(kGreen+4);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ZXbkg][fs]); //add to hs
      // TTV process: TTW + TTZ
      h1_BDTinput_4ljjsel[bdtIn][TTW][fs]->SetFillColor(kBlue+3);
      h1_BDTinput_4ljjsel[bdtIn][TTW][fs]->SetLineColor(kBlue+3);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTW][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]->SetFillColor(kBlue+3);
      h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]->SetLineColor(kBlue+3);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]); //add to hs
      // ggZZ process
      h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]->SetFillColor(kAzure-3);
      h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]->SetLineColor(kBlue+2);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]); //add to hs
      // qqZZ process
      h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]->SetFillColor(kAzure+6);
      h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]->SetLineColor(kAzure-6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]); //add to hs
      // SM Higgs processes: ggH + VBF + VH + ttH + bbH
      h1_BDTinput_4ljjsel[bdtIn][ggH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][ggH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ggH][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][VBF][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][VBF][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VBF][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][VH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][VH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VH][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][ttH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][ttH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ttH][fs]); //add to hs
      h1_BDTinput_4ljjsel[bdtIn][bbH][fs]->SetFillColor(kViolet+6);
      h1_BDTinput_4ljjsel[bdtIn][bbH][fs]->SetLineColor(kViolet+6);
      hs_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][bbH][fs]); //add to hs
      // HH signal
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->SetLineColor(kRed);
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->SetLineWidth(2);
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->Scale(100.);
      // data
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->SetMarkerColor(kBlack);
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->SetLineColor(kBlack);
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->SetMarkerStyle(20);
    
      // --- upper plot pad
      pad1_BDTinput_4ljjsel[bdtIn][fs] = new TPad("pad1_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs],"pad1_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs], 0, 0.3, 1, 1.0);
      pad1_BDTinput_4ljjsel[bdtIn][fs]->Draw();
      pad1_BDTinput_4ljjsel[bdtIn][fs]->cd();
    
      hs_BDTinput_4ljjsel[bdtIn][fs]->SetMaximum(10e04);
      hs_BDTinput_4ljjsel[bdtIn][fs]->SetMinimum(10e-04);
    
      hs_BDTinput_4ljjsel[bdtIn][fs]->Draw("histo");
      h1_BDTinput_4ljjsel[bdtIn][HH][fs]->Draw("histosame");
      h1_BDTinput_4ljjsel[bdtIn][Data][fs]->Draw("samepe");

      hs_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelFont(43);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelSize(15);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitle(h1_BDTinput_4ljjsel[bdtIn][HH][fs]->GetXaxis()->GetTitle());
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleSize(20);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleFont(43);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleOffset(1.4);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelFont(43);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelSize(15);
      hs_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitle(h1_BDTinput_4ljjsel[bdtIn][HH][fs]->GetYaxis()->GetTitle());

      // --- legend
      leg_BDTinput_4ljjsel[bdtIn][fs] = new TLegend(0.78,0.61,0.94,0.87);
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][Data][fs], "Data",          "lp");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][HH][fs],   "HH->4lbb x100", "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][ggH][fs],  "SM Higgs",      "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs], "qq->ZZ",        "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs], "gg->ZZ",        "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][TTZ][fs],  "TTV; V=Z,W",    "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][ZXbkg][fs],"Z+X",           "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->AddEntry(h1_BDTinput_4ljjsel[bdtIn][VVV][fs],  "VVV; V=Z,W",    "f");
      leg_BDTinput_4ljjsel[bdtIn][fs]->SetFillColor(kWhite);
      leg_BDTinput_4ljjsel[bdtIn][fs]->SetLineColor(kBlack);
      leg_BDTinput_4ljjsel[bdtIn][fs]->SetTextFont(43);
      leg_BDTinput_4ljjsel[bdtIn][fs]->Draw();

      c_BDTinput_4ljjsel[bdtIn][fs]->Update();

      pad1_BDTinput_4ljjsel[bdtIn][fs]->SetLogy();

      c_BDTinput_4ljjsel[bdtIn][fs]->Update();

      // --- tot hist for all MC
      hMCtot_BDTinput_4ljjsel[bdtIn][fs] = (TH1F*)h1_BDTinput_4ljjsel[bdtIn][ggH][fs]->Clone("hMCtot_BDTinput_4ljjsel_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VBF][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VH][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ttH][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][bbH][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][qqZZ][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ggZZ][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTZ][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][TTW][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][ZXbkg][fs]);
      hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Add(h1_BDTinput_4ljjsel[bdtIn][VVV][fs]);
      
      // --- lower pad plot
      c_BDTinput_4ljjsel[bdtIn][fs]->cd();
      pad2_BDTinput_4ljjsel[bdtIn][fs] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
      pad2_BDTinput_4ljjsel[bdtIn][fs]->SetGridy();
      pad2_BDTinput_4ljjsel[bdtIn][fs]->Draw();
      pad2_BDTinput_4ljjsel[bdtIn][fs]->cd();

      // --- define ratio plot
      rp_BDTinput_4ljjsel[bdtIn][fs] = (TH1F*)h1_BDTinput_4ljjsel[bdtIn][Data][fs]->Clone("rp_BDTinput_4ljjsel_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetLineColor(kBlack);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMinimum(0.);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMaximum(2.);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetStats(0);
      rp_BDTinput_4ljjsel[bdtIn][fs]->Divide(hMCtot_BDTinput_4ljjsel[bdtIn][fs]); //divide histo rp/MC
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerStyle(20);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerColor(kBlack);
      rp_BDTinput_4ljjsel[bdtIn][fs]->SetTitle("");

      rp_BDTinput_4ljjsel[bdtIn][fs]->SetYTitle("Data/#Sigma bkg");
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetNdivisions(505);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleSize(20);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetTitleOffset(1.4);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetYaxis()->SetLabelSize(15);

      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitleSize(20);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitleFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetTitleOffset(4.);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelFont(43);
      rp_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->SetLabelSize(15);

      // --- define mc shadow unc plot
      hUncMC_BDTinput_4ljjsel[bdtIn][fs] = (TH1F*)hMCtot_BDTinput_4ljjsel[bdtIn][fs]->Clone("hUncMC_BDTinput_4ljjsel_"+sBDTInputNames[bdtIn]+"_"+sFinalState[fs]);
      for(int xbin=1; xbin < hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetXaxis()->GetNbins() + 1; xbin++){
        float err = 0.;
        if(hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetBinContent(xbin) == 0) continue;
        err = hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetBinError(xbin) / hUncMC_BDTinput_4ljjsel[bdtIn][fs]->GetBinContent(xbin);
        hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetBinContent(xbin, 1.);
        hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetBinError(xbin, err);
      }
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetLineColor(1);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetFillStyle(3005);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetFillColor(kGray+3);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerColor(1);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetMarkerStyle(1);
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetTitle("");
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->SetStats(0);

      // ---draw
      rp_BDTinput_4ljjsel[bdtIn][fs]->Draw("ep");
      hUncMC_BDTinput_4ljjsel[bdtIn][fs]->Draw("e2 same");

      c_BDTinput_4ljjsel[bdtIn][fs]->Update();

      // --- draw CMS and lumi text
      writeExtraText = true;
      extraText      = "Preliminary";
      lumi_sqrtS     = lumiText + " (13 TeV)";
      cmsTextSize    = 0.6;
      lumiTextSize   = 0.46;
      extraOverCmsTextSize = 0.75;
      relPosX = 0.12;
      CMS_lumi(pad1_BDTinput_4ljjsel[bdtIn][fs], 0, 0);


      c_BDTinput_4ljjsel[bdtIn][fs]->SaveAs(outPath_inputBDTwithZXplots + "/" + c_BDTinput_4ljjsel[bdtIn][fs]->GetName() + ".png");
      c_BDTinput_4ljjsel[bdtIn][fs]->SaveAs(outPath_inputBDTwithZXplots + "/" + c_BDTinput_4ljjsel[bdtIn][fs]->GetName() + ".pdf");
      

    }
  }

} // end doPlots_inputBDT_withZX function




//*********************************
//*** doPlots_4ljjsel function ***
//*********************************
void doPlots_4ljjsel(){

 cout<<"do plots after 4ljj sel..."<<endl;

  //---input path
  TString sYear;
  TString lumiText;
  if(year==2016){
      sYear    = "2016";
      lumiText = "35.8 fb^{-1}";
  }
  else if(year==2017){
      sYear    = "2017";
      lumiText = "41.5 fb^{-1}";
  }
  else if(year==2018){
      sYear    = "2018";
      lumiText = "59.7 fb^{-1}";
  }
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;


  TString outPath_4ljjselplots = "plots_4ljjsel_" + sYear;
  cout<<"creating output dir "<<outPath_4ljjselplots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_4ljjselplots)).c_str()); // create output dir



  // retrieve histos from file
  TString inFileName = "f_histos_h1_4ljjsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  Int_t nPlots = 10;
  TString sPlots[] = {
    "m4l_4ljjsel",
    "mbb_4ljjsel",
    "pT4l_4ljjsel_sidebands",
    "j1btag_4ljjsel_sidebands",
    "j2btag_4ljjsel_sidebands",
    "j1pT_4ljjsel_sidebands",
    "j2pT_4ljjsel_sidebands",
    "MET_4ljjsel_sidebands",    
    "DeltaRhh_4ljjsel_sidebands",
    "mbb_4ljjsel_sidebands",
  };
  TH1F* h1_4ljjsel[nPlots][nProcesses][nFinalStates+1];
  for(int pl=0; pl<nPlots; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates+1; fs++){
        h1_4ljjsel[pl][pr][fs] = (TH1F*)fInhistos->Get("h1_"+sPlots[pl]+"_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
        cout<<h1_4ljjsel[pl][pr][fs]->GetName()<<endl;
      }
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_4ljjsel      [nPlots][nFinalStates+1];
  THStack* hs_4ljjsel     [nPlots][nFinalStates+1];
  TPad*    pad1_4ljjsel   [nPlots][nFinalStates+1];
  TLegend* leg_4ljjsel    [nPlots][nFinalStates+1];
  TH1F*    hMCtot_4ljjsel [nPlots][nFinalStates+1];
  TPad*    pad2_4ljjsel   [nPlots][nFinalStates+1];
  TH1F*    rp_4ljjsel     [nPlots][nFinalStates+1];
  TH1F*    hUncMC_4ljjsel [nPlots][nFinalStates+1];

  
  //4ljjsel plots
  for(int pl=0; pl<nPlots; pl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // canvas
      c_4ljjsel[pl][fs] = new TCanvas("c_"+sPlots[pl]+"_"+sFinalState[fs]+"_"+sYear,"c_"+sPlots[pl]+"_"+sFinalState[fs]+"_"+sYear,800,800);
      // hstack
      hs_4ljjsel[pl][fs] = new THStack("hs_"+sPlots[pl]+"_"+sFinalState[fs],"");
      // VVV process
      h1_4ljjsel[pl][VVV][fs]->SetFillColor(kGreen-3);
      h1_4ljjsel[pl][VVV][fs]->SetLineColor(kGreen-1);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][VVV][fs]); //add to hs
      // Z+X process
      h1_4ljjsel[pl][ZXbkg][fs]->SetFillColor(kGreen+3);
      h1_4ljjsel[pl][ZXbkg][fs]->SetLineColor(kGreen+4);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][ZXbkg][fs]); //add to hs
      // TTV process: TTW + TTV
      h1_4ljjsel[pl][TTW][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel[pl][TTW][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][TTW][fs]); //add to hs
      h1_4ljjsel[pl][TTZ][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel[pl][TTZ][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][TTZ][fs]); //add to hs
      // ggZZ process
      h1_4ljjsel[pl][ggZZ][fs]->SetFillColor(kAzure-3);
      h1_4ljjsel[pl][ggZZ][fs]->SetLineColor(kBlue+2);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][ggZZ][fs]); //add to hs
      // qqZZ process
      h1_4ljjsel[pl][qqZZ][fs]->SetFillColor(kAzure+6);
      h1_4ljjsel[pl][qqZZ][fs]->SetLineColor(kAzure-6);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][qqZZ][fs]); //add to hs
      // SM Higgs processes: ggH + VBF + VH + ttH + bbH
      h1_4ljjsel[pl][ggH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel[pl][ggH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][ggH][fs]); //add to hs
      h1_4ljjsel[pl][VBF][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel[pl][VBF][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][VBF][fs]); //add to hs
      h1_4ljjsel[pl][VH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel[pl][VH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][VH][fs]); //add to hs
      h1_4ljjsel[pl][ttH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel[pl][ttH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][ttH][fs]); //add to hs
      h1_4ljjsel[pl][bbH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel[pl][bbH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][bbH][fs]); //add to hs
      // HH signal
      h1_4ljjsel[pl][HH][fs]->SetLineColor(kRed);
      h1_4ljjsel[pl][HH][fs]->SetLineWidth(2);
      h1_4ljjsel[pl][HH][fs]->Scale(100.);
      // data
      h1_4ljjsel[pl][Data][fs]->SetMarkerColor(kBlack);
      h1_4ljjsel[pl][Data][fs]->SetLineColor(kBlack);
      h1_4ljjsel[pl][Data][fs]->SetMarkerStyle(20);

      // --- upper plot pad
      pad1_4ljjsel[pl][fs] = new TPad("pad1_"+sPlots[pl]+"_"+sFinalState[fs],"pad1_"+sPlots[pl]+"_"+sFinalState[fs], 0, 0.3, 1, 1.0);
      pad1_4ljjsel[pl][fs]->Draw();
      pad1_4ljjsel[pl][fs]->cd();

      hs_4ljjsel[pl][fs]->SetMaximum(10e04);
      hs_4ljjsel[pl][fs]->SetMinimum(10e-04);
        
      hs_4ljjsel[pl][fs]->Draw("histo");
      h1_4ljjsel[pl][HH][fs]->Draw("histosame");
      h1_4ljjsel[pl][Data][fs]->Draw("samepe");

      hs_4ljjsel[pl][fs]->GetXaxis()->SetLabelFont(43);
      hs_4ljjsel[pl][fs]->GetXaxis()->SetLabelSize(15);
      hs_4ljjsel[pl][fs]->GetXaxis()->SetTitle(h1_4ljjsel[pl][HH][fs]->GetXaxis()->GetTitle());
      hs_4ljjsel[pl][fs]->GetYaxis()->SetTitleSize(20);
      hs_4ljjsel[pl][fs]->GetYaxis()->SetTitleFont(43);
      hs_4ljjsel[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      hs_4ljjsel[pl][fs]->GetYaxis()->SetLabelFont(43);
      hs_4ljjsel[pl][fs]->GetYaxis()->SetLabelSize(15);
      hs_4ljjsel[pl][fs]->GetYaxis()->SetTitle(h1_4ljjsel[pl][HH][fs]->GetYaxis()->GetTitle());

      // --- legend
      leg_4ljjsel[pl][fs] = new TLegend(0.78,0.61,0.94,0.87);
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][Data][fs], "Data",          "lp");
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][HH][fs],   "HH->4lbb x100", "f");
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][ggH][fs],  "SM Higgs",      "f");
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][qqZZ][fs], "qq->ZZ",        "f");
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][ggZZ][fs], "gg->ZZ",        "f");
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][TTZ][fs],  "TTV; V=Z,W",    "f");
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][ZXbkg][fs],"Z+X",           "f");
      leg_4ljjsel[pl][fs]->AddEntry(h1_4ljjsel[pl][VVV][fs],  "VVV; V=Z,W",    "f");
      leg_4ljjsel[pl][fs]->SetFillColor(kWhite);
      leg_4ljjsel[pl][fs]->SetLineColor(kBlack);
      leg_4ljjsel[pl][fs]->SetTextFont(43);
      leg_4ljjsel[pl][fs]->Draw();

      c_4ljjsel[pl][fs]->Update();

      pad1_4ljjsel[pl][fs]->SetLogy();

      c_4ljjsel[pl][fs]->Update();

      // --- tot hist for all MC
      hMCtot_4ljjsel[pl][fs] = (TH1F*)h1_4ljjsel[pl][ggH][fs]->Clone("hMCtot_4ljjsel_"+sPlots[pl]+"_"+sFinalState[fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][VBF][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][VH][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][ttH][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][bbH][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][qqZZ][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][ggZZ][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][TTZ][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][TTW][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][ZXbkg][fs]);
      hMCtot_4ljjsel[pl][fs]->Add(h1_4ljjsel[pl][VVV][fs]);

      // --- lower pad plot
      c_4ljjsel[pl][fs]->cd();
      pad2_4ljjsel[pl][fs] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
      pad2_4ljjsel[pl][fs]->SetGridy();
      pad2_4ljjsel[pl][fs]->Draw();
      pad2_4ljjsel[pl][fs]->cd();

      // --- define ratio plot
      rp_4ljjsel[pl][fs] = (TH1F*)h1_4ljjsel[pl][Data][fs]->Clone("rp_4ljjsel_"+sPlots[pl]+"_"+sFinalState[fs]);
      rp_4ljjsel[pl][fs]->SetLineColor(kBlack);
      rp_4ljjsel[pl][fs]->SetMinimum(0.);
      rp_4ljjsel[pl][fs]->SetMaximum(2.);
      rp_4ljjsel[pl][fs]->SetStats(0);
      rp_4ljjsel[pl][fs]->Divide(hMCtot_4ljjsel[pl][fs]); //divide histo rp/MC
      rp_4ljjsel[pl][fs]->SetMarkerStyle(20);
      rp_4ljjsel[pl][fs]->SetMarkerColor(kBlack);
      rp_4ljjsel[pl][fs]->SetTitle("");

      rp_4ljjsel[pl][fs]->SetYTitle("Data/#Sigma bkg");
      rp_4ljjsel[pl][fs]->GetYaxis()->SetNdivisions(505);
      rp_4ljjsel[pl][fs]->GetYaxis()->SetTitleSize(20);
      rp_4ljjsel[pl][fs]->GetYaxis()->SetTitleFont(43);
      rp_4ljjsel[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      rp_4ljjsel[pl][fs]->GetYaxis()->SetLabelFont(43);
      rp_4ljjsel[pl][fs]->GetYaxis()->SetLabelSize(15);

      rp_4ljjsel[pl][fs]->GetXaxis()->SetTitleSize(20);
      rp_4ljjsel[pl][fs]->GetXaxis()->SetTitleFont(43);
      rp_4ljjsel[pl][fs]->GetXaxis()->SetTitleOffset(4.);
      rp_4ljjsel[pl][fs]->GetXaxis()->SetLabelFont(43);
      rp_4ljjsel[pl][fs]->GetXaxis()->SetLabelSize(15);

      // --- define mc shadow unc plot
      hUncMC_4ljjsel[pl][fs] = (TH1F*)hMCtot_4ljjsel[pl][fs]->Clone("hUncMC_4ljjsel_"+sPlots[pl]+"_"+sFinalState[fs]);
      for(int xbin=1; xbin < hUncMC_4ljjsel[pl][fs]->GetXaxis()->GetNbins() + 1; xbin++){
        float err = 0.;
        if(hUncMC_4ljjsel[pl][fs]->GetBinContent(xbin) == 0) continue;
        err = hUncMC_4ljjsel[pl][fs]->GetBinError(xbin) / hUncMC_4ljjsel[pl][fs]->GetBinContent(xbin);
        hUncMC_4ljjsel[pl][fs]->SetBinContent(xbin, 1.);
        hUncMC_4ljjsel[pl][fs]->SetBinError(xbin, err);
      }
      hUncMC_4ljjsel[pl][fs]->SetLineColor(1);
      hUncMC_4ljjsel[pl][fs]->SetFillStyle(3005);
      hUncMC_4ljjsel[pl][fs]->SetFillColor(kGray+3);
      hUncMC_4ljjsel[pl][fs]->SetMarkerColor(1);
      hUncMC_4ljjsel[pl][fs]->SetMarkerStyle(1);
      hUncMC_4ljjsel[pl][fs]->SetTitle("");
      hUncMC_4ljjsel[pl][fs]->SetStats(0);

      // ---draw
      rp_4ljjsel[pl][fs]->Draw("ep");
      hUncMC_4ljjsel[pl][fs]->Draw("e2 same");

      c_4ljjsel[pl][fs]->Update();

      // --- draw CMS and lumi text
      writeExtraText = true;
      extraText      = "Preliminary";
      lumi_sqrtS     = lumiText + " (13 TeV)";
      cmsTextSize    = 0.6;
      lumiTextSize   = 0.46;
      extraOverCmsTextSize = 0.75;
      relPosX = 0.12;
      CMS_lumi(pad1_4ljjsel[pl][fs], 0, 0);


      c_4ljjsel[pl][fs]->SaveAs(outPath_4ljjselplots + "/" + c_4ljjsel[pl][fs]->GetName() + ".png");
      c_4ljjsel[pl][fs]->SaveAs(outPath_4ljjselplots + "/" + c_4ljjsel[pl][fs]->GetName() + ".pdf");
      

    }
  }



} // end function doPlots_4ljjsel




//*******************************************
//*** doPlots_4ljjsel_sidebandSX function ***
//*******************************************
void doPlots_4ljjsel_sidebandSX(){

 cout<<"do plots after 4ljj sel in sidebandSX..."<<endl;

  //---input path
  TString sYear;
  TString lumiText;
  if(year==2016){
      sYear    = "2016";
      lumiText = "35.8 fb^{-1}";
  }
  else if(year==2017){
      sYear    = "2017";
      lumiText = "41.5 fb^{-1}";
  }
  else if(year==2018){
      sYear    = "2018";
      lumiText = "59.7 fb^{-1}";
  }
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;


  TString outPath_4ljjsel_sidebandSX_plots = "plots_4ljjsel_sidebandSX_" + sYear;
  cout<<"creating output dir "<<outPath_4ljjsel_sidebandSX_plots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_4ljjsel_sidebandSX_plots)).c_str()); // create output dir



  // retrieve histos from file
  TString inFileName = "f_histos_h1_4ljjsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  Int_t nPlots_sidebandSX = 9;
  TString sPlots_sidebandSX[] = {
    "m4l_4ljjsel_sidebandSX",
    "pT4l_4ljjsel_sidebandSX",
    "j1btag_4ljjsel_sidebandSX",
    "j2btag_4ljjsel_sidebandSX",
    "j1pT_4ljjsel_sidebandSX",
    "j2pT_4ljjsel_sidebandSX",
    "MET_4ljjsel_sidebandSX",    
    "DeltaRhh_4ljjsel_sidebandSX",
    "mbb_4ljjsel_sidebandSX",
  };
  TH1F* h1_4ljjsel_sidebandSX[nPlots_sidebandSX][nProcesses][nFinalStates+1];
  for(int pl=0; pl<nPlots_sidebandSX; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates+1; fs++){
        h1_4ljjsel_sidebandSX[pl][pr][fs] = (TH1F*)fInhistos->Get("h1_"+sPlots_sidebandSX[pl]+"_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
        cout<<h1_4ljjsel_sidebandSX[pl][pr][fs]->GetName()<<endl;
      }
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_4ljjsel_sidebandSX      [nPlots_sidebandSX][nFinalStates+1];
  THStack* hs_4ljjsel_sidebandSX     [nPlots_sidebandSX][nFinalStates+1];
  TPad*    pad1_4ljjsel_sidebandSX   [nPlots_sidebandSX][nFinalStates+1];
  TLegend* leg_4ljjsel_sidebandSX    [nPlots_sidebandSX][nFinalStates+1];
  TH1F*    hMCtot_4ljjsel_sidebandSX [nPlots_sidebandSX][nFinalStates+1];
  TPad*    pad2_4ljjsel_sidebandSX   [nPlots_sidebandSX][nFinalStates+1];
  TH1F*    rp_4ljjsel_sidebandSX     [nPlots_sidebandSX][nFinalStates+1];
  TH1F*    hUncMC_4ljjsel_sidebandSX [nPlots_sidebandSX][nFinalStates+1];

  
  //4ljjsel plots
  for(int pl=0; pl<nPlots_sidebandSX; pl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // canvas
      c_4ljjsel_sidebandSX[pl][fs] = new TCanvas("c_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs]+"_"+sYear,"c_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs]+"_"+sYear,800,800);
      // hstack
      hs_4ljjsel_sidebandSX[pl][fs] = new THStack("hs_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs],"");
      // VVV process
      h1_4ljjsel_sidebandSX[pl][VVV][fs]->SetFillColor(kGreen-3);
      h1_4ljjsel_sidebandSX[pl][VVV][fs]->SetLineColor(kGreen-1);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][VVV][fs]); //add to hs
      // Z+X process
      h1_4ljjsel_sidebandSX[pl][ZXbkg][fs]->SetFillColor(kGreen+3);
      h1_4ljjsel_sidebandSX[pl][ZXbkg][fs]->SetLineColor(kGreen+4);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][ZXbkg][fs]); //add to hs
      // TTV process: TTW + TTV
      h1_4ljjsel_sidebandSX[pl][TTW][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel_sidebandSX[pl][TTW][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][TTW][fs]); //add to hs
      h1_4ljjsel_sidebandSX[pl][TTZ][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel_sidebandSX[pl][TTZ][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][TTZ][fs]); //add to hs
      // ggZZ process
      h1_4ljjsel_sidebandSX[pl][ggZZ][fs]->SetFillColor(kAzure-3);
      h1_4ljjsel_sidebandSX[pl][ggZZ][fs]->SetLineColor(kBlue+2);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][ggZZ][fs]); //add to hs
      // qqZZ process
      h1_4ljjsel_sidebandSX[pl][qqZZ][fs]->SetFillColor(kAzure+6);
      h1_4ljjsel_sidebandSX[pl][qqZZ][fs]->SetLineColor(kAzure-6);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][qqZZ][fs]); //add to hs
      // SM Higgs processes: ggH + VBF + VH + ttH + bbH
      h1_4ljjsel_sidebandSX[pl][ggH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSX[pl][ggH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][ggH][fs]); //add to hs
      h1_4ljjsel_sidebandSX[pl][VBF][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSX[pl][VBF][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][VBF][fs]); //add to hs
      h1_4ljjsel_sidebandSX[pl][VH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSX[pl][VH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][VH][fs]); //add to hs
      h1_4ljjsel_sidebandSX[pl][ttH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSX[pl][ttH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][ttH][fs]); //add to hs
      h1_4ljjsel_sidebandSX[pl][bbH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSX[pl][bbH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][bbH][fs]); //add to hs
      // HH signal
      h1_4ljjsel_sidebandSX[pl][HH][fs]->SetLineColor(kRed);
      h1_4ljjsel_sidebandSX[pl][HH][fs]->SetLineWidth(2);
      h1_4ljjsel_sidebandSX[pl][HH][fs]->Scale(100.);
      // data
      h1_4ljjsel_sidebandSX[pl][Data][fs]->SetMarkerColor(kBlack);
      h1_4ljjsel_sidebandSX[pl][Data][fs]->SetLineColor(kBlack);
      h1_4ljjsel_sidebandSX[pl][Data][fs]->SetMarkerStyle(20);

      // --- upper plot pad
      pad1_4ljjsel_sidebandSX[pl][fs] = new TPad("pad1_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs],"pad1_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs], 0, 0.3, 1, 1.0);
      pad1_4ljjsel_sidebandSX[pl][fs]->Draw();
      pad1_4ljjsel_sidebandSX[pl][fs]->cd();

      hs_4ljjsel_sidebandSX[pl][fs]->SetMaximum(10e04);
      hs_4ljjsel_sidebandSX[pl][fs]->SetMinimum(10e-04);
        
      hs_4ljjsel_sidebandSX[pl][fs]->Draw("histo");
      h1_4ljjsel_sidebandSX[pl][HH][fs]->Draw("histosame");
      h1_4ljjsel_sidebandSX[pl][Data][fs]->Draw("samepe");

      hs_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetLabelFont(43);
      hs_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetLabelSize(15);
      hs_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetTitle(h1_4ljjsel_sidebandSX[pl][HH][fs]->GetXaxis()->GetTitle());
      hs_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetTitleSize(20);
      hs_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetTitleFont(43);
      hs_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      hs_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetLabelFont(43);
      hs_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetLabelSize(15);
      hs_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetTitle(h1_4ljjsel_sidebandSX[pl][HH][fs]->GetYaxis()->GetTitle());

      // --- legend
      leg_4ljjsel_sidebandSX[pl][fs] = new TLegend(0.78,0.61,0.94,0.87);
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][Data][fs], "Data",          "lp");
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][HH][fs],   "HH->4lbb x100", "f");
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][ggH][fs],  "SM Higgs",      "f");
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][qqZZ][fs], "qq->ZZ",        "f");
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][ggZZ][fs], "gg->ZZ",        "f");
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][TTZ][fs],  "TTV; V=Z,W",    "f");
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][ZXbkg][fs],"Z+X",           "f");
      leg_4ljjsel_sidebandSX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSX[pl][VVV][fs],  "VVV; V=Z,W",    "f");
      leg_4ljjsel_sidebandSX[pl][fs]->SetFillColor(kWhite);
      leg_4ljjsel_sidebandSX[pl][fs]->SetLineColor(kBlack);
      leg_4ljjsel_sidebandSX[pl][fs]->SetTextFont(43);
      leg_4ljjsel_sidebandSX[pl][fs]->Draw();

      c_4ljjsel_sidebandSX[pl][fs]->Update();

      pad1_4ljjsel_sidebandSX[pl][fs]->SetLogy();

      c_4ljjsel_sidebandSX[pl][fs]->Update();

      // --- tot hist for all MC
      hMCtot_4ljjsel_sidebandSX[pl][fs] = (TH1F*)h1_4ljjsel_sidebandSX[pl][ggH][fs]->Clone("hMCtot_4ljjsel_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][VBF][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][VH][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][ttH][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][bbH][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][qqZZ][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][ggZZ][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][TTZ][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][TTW][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][ZXbkg][fs]);
      hMCtot_4ljjsel_sidebandSX[pl][fs]->Add(h1_4ljjsel_sidebandSX[pl][VVV][fs]);

      // --- lower pad plot
      c_4ljjsel_sidebandSX[pl][fs]->cd();
      pad2_4ljjsel_sidebandSX[pl][fs] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
      pad2_4ljjsel_sidebandSX[pl][fs]->SetGridy();
      pad2_4ljjsel_sidebandSX[pl][fs]->Draw();
      pad2_4ljjsel_sidebandSX[pl][fs]->cd();

      // --- define ratio plot
      rp_4ljjsel_sidebandSX[pl][fs] = (TH1F*)h1_4ljjsel_sidebandSX[pl][Data][fs]->Clone("rp_4ljjsel_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs]);
      rp_4ljjsel_sidebandSX[pl][fs]->SetLineColor(kBlack);
      rp_4ljjsel_sidebandSX[pl][fs]->SetMinimum(0.);
      rp_4ljjsel_sidebandSX[pl][fs]->SetMaximum(2.);
      rp_4ljjsel_sidebandSX[pl][fs]->SetStats(0);
      rp_4ljjsel_sidebandSX[pl][fs]->Divide(hMCtot_4ljjsel_sidebandSX[pl][fs]); //divide histo rp/MC
      rp_4ljjsel_sidebandSX[pl][fs]->SetMarkerStyle(20);
      rp_4ljjsel_sidebandSX[pl][fs]->SetMarkerColor(kBlack);
      rp_4ljjsel_sidebandSX[pl][fs]->SetTitle("");

      rp_4ljjsel_sidebandSX[pl][fs]->SetYTitle("Data/#Sigma bkg");
      rp_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetNdivisions(505);
      rp_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetTitleSize(20);
      rp_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetTitleFont(43);
      rp_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      rp_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetLabelFont(43);
      rp_4ljjsel_sidebandSX[pl][fs]->GetYaxis()->SetLabelSize(15);

      rp_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetTitleSize(20);
      rp_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetTitleFont(43);
      rp_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetTitleOffset(4.);
      rp_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetLabelFont(43);
      rp_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->SetLabelSize(15);

      // --- define mc shadow unc plot
      hUncMC_4ljjsel_sidebandSX[pl][fs] = (TH1F*)hMCtot_4ljjsel_sidebandSX[pl][fs]->Clone("hUncMC_4ljjsel_"+sPlots_sidebandSX[pl]+"_"+sFinalState[fs]);
      for(int xbin=1; xbin < hUncMC_4ljjsel_sidebandSX[pl][fs]->GetXaxis()->GetNbins() + 1; xbin++){
        float err = 0.;
        if(hUncMC_4ljjsel_sidebandSX[pl][fs]->GetBinContent(xbin) == 0) continue;
        err = hUncMC_4ljjsel_sidebandSX[pl][fs]->GetBinError(xbin) / hUncMC_4ljjsel_sidebandSX[pl][fs]->GetBinContent(xbin);
        hUncMC_4ljjsel_sidebandSX[pl][fs]->SetBinContent(xbin, 1.);
        hUncMC_4ljjsel_sidebandSX[pl][fs]->SetBinError(xbin, err);
      }
      hUncMC_4ljjsel_sidebandSX[pl][fs]->SetLineColor(1);
      hUncMC_4ljjsel_sidebandSX[pl][fs]->SetFillStyle(3005);
      hUncMC_4ljjsel_sidebandSX[pl][fs]->SetFillColor(kGray+3);
      hUncMC_4ljjsel_sidebandSX[pl][fs]->SetMarkerColor(1);
      hUncMC_4ljjsel_sidebandSX[pl][fs]->SetMarkerStyle(1);
      hUncMC_4ljjsel_sidebandSX[pl][fs]->SetTitle("");
      hUncMC_4ljjsel_sidebandSX[pl][fs]->SetStats(0);

      // ---draw
      rp_4ljjsel_sidebandSX[pl][fs]->Draw("ep");
      hUncMC_4ljjsel_sidebandSX[pl][fs]->Draw("e2 same");

      c_4ljjsel_sidebandSX[pl][fs]->Update();

      // --- draw CMS and lumi text
      writeExtraText = true;
      extraText      = "Preliminary";
      lumi_sqrtS     = lumiText + " (13 TeV)";
      cmsTextSize    = 0.6;
      lumiTextSize   = 0.46;
      extraOverCmsTextSize = 0.75;
      relPosX = 0.12;
      CMS_lumi(pad1_4ljjsel_sidebandSX[pl][fs], 0, 0);


      c_4ljjsel_sidebandSX[pl][fs]->SaveAs(outPath_4ljjsel_sidebandSX_plots + "/" + c_4ljjsel_sidebandSX[pl][fs]->GetName() + ".png");
      c_4ljjsel_sidebandSX[pl][fs]->SaveAs(outPath_4ljjsel_sidebandSX_plots + "/" + c_4ljjsel_sidebandSX[pl][fs]->GetName() + ".pdf");
      

    }
  }



} // end function doPlots_4ljjsel_sidebandSX





//*******************************************
//*** doPlots_4ljjsel_sidebandDX function ***
//*******************************************
void doPlots_4ljjsel_sidebandDX(){

 cout<<"do plots after 4ljj sel in sidebandDX..."<<endl;

  //---input path
  TString sYear;
  TString lumiText;
  if(year==2016){
      sYear    = "2016";
      lumiText = "35.8 fb^{-1}";
  }
  else if(year==2017){
      sYear    = "2017";
      lumiText = "41.5 fb^{-1}";
  }
  else if(year==2018){
      sYear    = "2018";
      lumiText = "59.7 fb^{-1}";
  }
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;


  TString outPath_4ljjsel_sidebandDX_plots = "plots_4ljjsel_sidebandDX_" + sYear;
  cout<<"creating output dir "<<outPath_4ljjsel_sidebandDX_plots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_4ljjsel_sidebandDX_plots)).c_str()); // create output dir



  // retrieve histos from file
  TString inFileName = "f_histos_h1_4ljjsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  Int_t nPlots_sidebandDX = 9;
  TString sPlots_sidebandDX[] = {
    "m4l_4ljjsel_sidebandDX",
    "pT4l_4ljjsel_sidebandDX",
    "j1btag_4ljjsel_sidebandDX",
    "j2btag_4ljjsel_sidebandDX",
    "j1pT_4ljjsel_sidebandDX",
    "j2pT_4ljjsel_sidebandDX",
    "MET_4ljjsel_sidebandDX",    
    "DeltaRhh_4ljjsel_sidebandDX",
    "mbb_4ljjsel_sidebandDX",
  };
  TH1F* h1_4ljjsel_sidebandDX[nPlots_sidebandDX][nProcesses][nFinalStates+1];
  for(int pl=0; pl<nPlots_sidebandDX; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates+1; fs++){
        h1_4ljjsel_sidebandDX[pl][pr][fs] = (TH1F*)fInhistos->Get("h1_"+sPlots_sidebandDX[pl]+"_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
        cout<<h1_4ljjsel_sidebandDX[pl][pr][fs]->GetName()<<endl;
      }
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_4ljjsel_sidebandDX      [nPlots_sidebandDX][nFinalStates+1];
  THStack* hs_4ljjsel_sidebandDX     [nPlots_sidebandDX][nFinalStates+1];
  TPad*    pad1_4ljjsel_sidebandDX   [nPlots_sidebandDX][nFinalStates+1];
  TLegend* leg_4ljjsel_sidebandDX    [nPlots_sidebandDX][nFinalStates+1];
  TH1F*    hMCtot_4ljjsel_sidebandDX [nPlots_sidebandDX][nFinalStates+1];
  TPad*    pad2_4ljjsel_sidebandDX   [nPlots_sidebandDX][nFinalStates+1];
  TH1F*    rp_4ljjsel_sidebandDX     [nPlots_sidebandDX][nFinalStates+1];
  TH1F*    hUncMC_4ljjsel_sidebandDX [nPlots_sidebandDX][nFinalStates+1];

  
  //4ljjsel plots
  for(int pl=0; pl<nPlots_sidebandDX; pl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // canvas
      c_4ljjsel_sidebandDX[pl][fs] = new TCanvas("c_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs]+"_"+sYear,"c_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs]+"_"+sYear,800,800);
      // hstack
      hs_4ljjsel_sidebandDX[pl][fs] = new THStack("hs_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs],"");
      // VVV process
      h1_4ljjsel_sidebandDX[pl][VVV][fs]->SetFillColor(kGreen-3);
      h1_4ljjsel_sidebandDX[pl][VVV][fs]->SetLineColor(kGreen-1);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][VVV][fs]); //add to hs
      // Z+X process
      h1_4ljjsel_sidebandDX[pl][ZXbkg][fs]->SetFillColor(kGreen+3);
      h1_4ljjsel_sidebandDX[pl][ZXbkg][fs]->SetLineColor(kGreen+4);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][ZXbkg][fs]); //add to hs
      // TTV process: TTW + TTV
      h1_4ljjsel_sidebandDX[pl][TTW][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel_sidebandDX[pl][TTW][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][TTW][fs]); //add to hs
      h1_4ljjsel_sidebandDX[pl][TTZ][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel_sidebandDX[pl][TTZ][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][TTZ][fs]); //add to hs
      // ggZZ process
      h1_4ljjsel_sidebandDX[pl][ggZZ][fs]->SetFillColor(kAzure-3);
      h1_4ljjsel_sidebandDX[pl][ggZZ][fs]->SetLineColor(kBlue+2);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][ggZZ][fs]); //add to hs
      // qqZZ process
      h1_4ljjsel_sidebandDX[pl][qqZZ][fs]->SetFillColor(kAzure+6);
      h1_4ljjsel_sidebandDX[pl][qqZZ][fs]->SetLineColor(kAzure-6);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][qqZZ][fs]); //add to hs
      // SM Higgs processes: ggH + VBF + VH + ttH + bbH
      h1_4ljjsel_sidebandDX[pl][ggH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandDX[pl][ggH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][ggH][fs]); //add to hs
      h1_4ljjsel_sidebandDX[pl][VBF][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandDX[pl][VBF][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][VBF][fs]); //add to hs
      h1_4ljjsel_sidebandDX[pl][VH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandDX[pl][VH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][VH][fs]); //add to hs
      h1_4ljjsel_sidebandDX[pl][ttH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandDX[pl][ttH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][ttH][fs]); //add to hs
      h1_4ljjsel_sidebandDX[pl][bbH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandDX[pl][bbH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][bbH][fs]); //add to hs
      // HH signal
      h1_4ljjsel_sidebandDX[pl][HH][fs]->SetLineColor(kRed);
      h1_4ljjsel_sidebandDX[pl][HH][fs]->SetLineWidth(2);
      h1_4ljjsel_sidebandDX[pl][HH][fs]->Scale(100.);
      // data
      h1_4ljjsel_sidebandDX[pl][Data][fs]->SetMarkerColor(kBlack);
      h1_4ljjsel_sidebandDX[pl][Data][fs]->SetLineColor(kBlack);
      h1_4ljjsel_sidebandDX[pl][Data][fs]->SetMarkerStyle(20);

      // --- upper plot pad
      pad1_4ljjsel_sidebandDX[pl][fs] = new TPad("pad1_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs],"pad1_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs], 0, 0.3, 1, 1.0);
      pad1_4ljjsel_sidebandDX[pl][fs]->Draw();
      pad1_4ljjsel_sidebandDX[pl][fs]->cd();

      hs_4ljjsel_sidebandDX[pl][fs]->SetMaximum(10e04);
      hs_4ljjsel_sidebandDX[pl][fs]->SetMinimum(10e-04);
        
      hs_4ljjsel_sidebandDX[pl][fs]->Draw("histo");
      h1_4ljjsel_sidebandDX[pl][HH][fs]->Draw("histosame");
      h1_4ljjsel_sidebandDX[pl][Data][fs]->Draw("samepe");

      hs_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetLabelFont(43);
      hs_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetLabelSize(15);
      hs_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetTitle(h1_4ljjsel_sidebandDX[pl][HH][fs]->GetXaxis()->GetTitle());
      hs_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetTitleSize(20);
      hs_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetTitleFont(43);
      hs_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      hs_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetLabelFont(43);
      hs_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetLabelSize(15);
      hs_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetTitle(h1_4ljjsel_sidebandDX[pl][HH][fs]->GetYaxis()->GetTitle());

      // --- legend
      leg_4ljjsel_sidebandDX[pl][fs] = new TLegend(0.78,0.61,0.94,0.87);
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][Data][fs], "Data",          "lp");
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][HH][fs],   "HH->4lbb x100", "f");
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][ggH][fs],  "SM Higgs",      "f");
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][qqZZ][fs], "qq->ZZ",        "f");
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][ggZZ][fs], "gg->ZZ",        "f");
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][TTZ][fs],  "TTV; V=Z,W",    "f");
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][ZXbkg][fs],"Z+X",           "f");
      leg_4ljjsel_sidebandDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandDX[pl][VVV][fs],  "VVV; V=Z,W",    "f");
      leg_4ljjsel_sidebandDX[pl][fs]->SetFillColor(kWhite);
      leg_4ljjsel_sidebandDX[pl][fs]->SetLineColor(kBlack);
      leg_4ljjsel_sidebandDX[pl][fs]->SetTextFont(43);
      leg_4ljjsel_sidebandDX[pl][fs]->Draw();

      c_4ljjsel_sidebandDX[pl][fs]->Update();

      pad1_4ljjsel_sidebandDX[pl][fs]->SetLogy();

      c_4ljjsel_sidebandDX[pl][fs]->Update();

      // --- tot hist for all MC
      hMCtot_4ljjsel_sidebandDX[pl][fs] = (TH1F*)h1_4ljjsel_sidebandDX[pl][ggH][fs]->Clone("hMCtot_4ljjsel_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][VBF][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][VH][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][ttH][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][bbH][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][qqZZ][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][ggZZ][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][TTZ][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][TTW][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][ZXbkg][fs]);
      hMCtot_4ljjsel_sidebandDX[pl][fs]->Add(h1_4ljjsel_sidebandDX[pl][VVV][fs]);

      // --- lower pad plot
      c_4ljjsel_sidebandDX[pl][fs]->cd();
      pad2_4ljjsel_sidebandDX[pl][fs] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
      pad2_4ljjsel_sidebandDX[pl][fs]->SetGridy();
      pad2_4ljjsel_sidebandDX[pl][fs]->Draw();
      pad2_4ljjsel_sidebandDX[pl][fs]->cd();

      // --- define ratio plot
      rp_4ljjsel_sidebandDX[pl][fs] = (TH1F*)h1_4ljjsel_sidebandDX[pl][Data][fs]->Clone("rp_4ljjsel_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs]);
      rp_4ljjsel_sidebandDX[pl][fs]->SetLineColor(kBlack);
      rp_4ljjsel_sidebandDX[pl][fs]->SetMinimum(0.);
      rp_4ljjsel_sidebandDX[pl][fs]->SetMaximum(2.);
      rp_4ljjsel_sidebandDX[pl][fs]->SetStats(0);
      rp_4ljjsel_sidebandDX[pl][fs]->Divide(hMCtot_4ljjsel_sidebandDX[pl][fs]); //divide histo rp/MC
      rp_4ljjsel_sidebandDX[pl][fs]->SetMarkerStyle(20);
      rp_4ljjsel_sidebandDX[pl][fs]->SetMarkerColor(kBlack);
      rp_4ljjsel_sidebandDX[pl][fs]->SetTitle("");

      rp_4ljjsel_sidebandDX[pl][fs]->SetYTitle("Data/#Sigma bkg");
      rp_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetNdivisions(505);
      rp_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetTitleSize(20);
      rp_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetTitleFont(43);
      rp_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      rp_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetLabelFont(43);
      rp_4ljjsel_sidebandDX[pl][fs]->GetYaxis()->SetLabelSize(15);

      rp_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetTitleSize(20);
      rp_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetTitleFont(43);
      rp_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetTitleOffset(4.);
      rp_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetLabelFont(43);
      rp_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->SetLabelSize(15);

      // --- define mc shadow unc plot
      hUncMC_4ljjsel_sidebandDX[pl][fs] = (TH1F*)hMCtot_4ljjsel_sidebandDX[pl][fs]->Clone("hUncMC_4ljjsel_"+sPlots_sidebandDX[pl]+"_"+sFinalState[fs]);
      for(int xbin=1; xbin < hUncMC_4ljjsel_sidebandDX[pl][fs]->GetXaxis()->GetNbins() + 1; xbin++){
        float err = 0.;
        if(hUncMC_4ljjsel_sidebandDX[pl][fs]->GetBinContent(xbin) == 0) continue;
        err = hUncMC_4ljjsel_sidebandDX[pl][fs]->GetBinError(xbin) / hUncMC_4ljjsel_sidebandDX[pl][fs]->GetBinContent(xbin);
        hUncMC_4ljjsel_sidebandDX[pl][fs]->SetBinContent(xbin, 1.);
        hUncMC_4ljjsel_sidebandDX[pl][fs]->SetBinError(xbin, err);
      }
      hUncMC_4ljjsel_sidebandDX[pl][fs]->SetLineColor(1);
      hUncMC_4ljjsel_sidebandDX[pl][fs]->SetFillStyle(3005);
      hUncMC_4ljjsel_sidebandDX[pl][fs]->SetFillColor(kGray+3);
      hUncMC_4ljjsel_sidebandDX[pl][fs]->SetMarkerColor(1);
      hUncMC_4ljjsel_sidebandDX[pl][fs]->SetMarkerStyle(1);
      hUncMC_4ljjsel_sidebandDX[pl][fs]->SetTitle("");
      hUncMC_4ljjsel_sidebandDX[pl][fs]->SetStats(0);

      // ---draw
      rp_4ljjsel_sidebandDX[pl][fs]->Draw("ep");
      hUncMC_4ljjsel_sidebandDX[pl][fs]->Draw("e2 same");

      c_4ljjsel_sidebandDX[pl][fs]->Update();

      // --- draw CMS and lumi text
      writeExtraText = true;
      extraText      = "Preliminary";
      lumi_sqrtS     = lumiText + " (13 TeV)";
      cmsTextSize    = 0.6;
      lumiTextSize   = 0.46;
      extraOverCmsTextSize = 0.75;
      relPosX = 0.12;
      CMS_lumi(pad1_4ljjsel_sidebandDX[pl][fs], 0, 0);


      c_4ljjsel_sidebandDX[pl][fs]->SaveAs(outPath_4ljjsel_sidebandDX_plots + "/" + c_4ljjsel_sidebandDX[pl][fs]->GetName() + ".png");
      c_4ljjsel_sidebandDX[pl][fs]->SaveAs(outPath_4ljjsel_sidebandDX_plots + "/" + c_4ljjsel_sidebandDX[pl][fs]->GetName() + ".pdf");
      

    }
  }



} // end function doPlots_4ljjsel_sidebandDX





//*******************************************
//*** doPlots_4ljjsel_sidebandSXDX function ***
//*******************************************
void doPlots_4ljjsel_sidebandSXDX(){

 cout<<"do plots after 4ljj sel in sidebandSXDX..."<<endl;

  //---input path
  TString sYear;
  TString lumiText;
  if(year==2016){
      sYear    = "2016";
      lumiText = "35.8 fb^{-1}";
  }
  else if(year==2017){
      sYear    = "2017";
      lumiText = "41.5 fb^{-1}";
  }
  else if(year==2018){
      sYear    = "2018";
      lumiText = "59.7 fb^{-1}";
  }
  else cout<<"wrong year selected!"<<endl;
  cout<<"Year chosen: "<<year<<endl;


  TString outPath_4ljjsel_sidebandSXDX_plots = "plots_4ljjsel_sidebandSXDX_" + sYear;
  cout<<"creating output dir "<<outPath_4ljjsel_sidebandSXDX_plots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_4ljjsel_sidebandSXDX_plots)).c_str()); // create output dir



  // retrieve histos from file
  TString inFileName = "f_histos_h1_4ljjsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  Int_t nPlots_sidebandSXDX = 9;
  TString sPlots_sidebandSXDX[] = {
    "m4l_4ljjsel_sidebandSXDX",
    "pT4l_4ljjsel_sidebandSXDX",
    "j1btag_4ljjsel_sidebandSXDX",
    "j2btag_4ljjsel_sidebandSXDX",
    "j1pT_4ljjsel_sidebandSXDX",
    "j2pT_4ljjsel_sidebandSXDX",
    "MET_4ljjsel_sidebandSXDX",    
    "DeltaRhh_4ljjsel_sidebandSXDX",
    "mbb_4ljjsel_sidebandSXDX",
  };
  TH1F* h1_4ljjsel_sidebandSXDX[nPlots_sidebandSXDX][nProcesses][nFinalStates+1];
  for(int pl=0; pl<nPlots_sidebandSXDX; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates+1; fs++){
        h1_4ljjsel_sidebandSXDX[pl][pr][fs] = (TH1F*)fInhistos->Get("h1_"+sPlots_sidebandSXDX[pl]+"_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
        cout<<h1_4ljjsel_sidebandSXDX[pl][pr][fs]->GetName()<<endl;
      }
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_4ljjsel_sidebandSXDX      [nPlots_sidebandSXDX][nFinalStates+1];
  THStack* hs_4ljjsel_sidebandSXDX     [nPlots_sidebandSXDX][nFinalStates+1];
  TPad*    pad1_4ljjsel_sidebandSXDX   [nPlots_sidebandSXDX][nFinalStates+1];
  TLegend* leg_4ljjsel_sidebandSXDX    [nPlots_sidebandSXDX][nFinalStates+1];
  TH1F*    hMCtot_4ljjsel_sidebandSXDX [nPlots_sidebandSXDX][nFinalStates+1];
  TPad*    pad2_4ljjsel_sidebandSXDX   [nPlots_sidebandSXDX][nFinalStates+1];
  TH1F*    rp_4ljjsel_sidebandSXDX     [nPlots_sidebandSXDX][nFinalStates+1];
  TH1F*    hUncMC_4ljjsel_sidebandSXDX [nPlots_sidebandSXDX][nFinalStates+1];

  
  //4ljjsel plots
  for(int pl=0; pl<nPlots_sidebandSXDX; pl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // canvas
      c_4ljjsel_sidebandSXDX[pl][fs] = new TCanvas("c_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs]+"_"+sYear,"c_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs]+"_"+sYear,800,800);
      // hstack
      hs_4ljjsel_sidebandSXDX[pl][fs] = new THStack("hs_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs],"");
      // VVV process
      h1_4ljjsel_sidebandSXDX[pl][VVV][fs]->SetFillColor(kGreen-3);
      h1_4ljjsel_sidebandSXDX[pl][VVV][fs]->SetLineColor(kGreen-1);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][VVV][fs]); //add to hs
      // Z+X process
      h1_4ljjsel_sidebandSXDX[pl][ZXbkg][fs]->SetFillColor(kGreen+3);
      h1_4ljjsel_sidebandSXDX[pl][ZXbkg][fs]->SetLineColor(kGreen+4);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][ZXbkg][fs]); //add to hs
      // TTV process: TTW + TTV
      h1_4ljjsel_sidebandSXDX[pl][TTW][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel_sidebandSXDX[pl][TTW][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][TTW][fs]); //add to hs
      h1_4ljjsel_sidebandSXDX[pl][TTZ][fs]->SetFillColor(kBlue+3);
      h1_4ljjsel_sidebandSXDX[pl][TTZ][fs]->SetLineColor(kBlue+3);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][TTZ][fs]); //add to hs
      // ggZZ process
      h1_4ljjsel_sidebandSXDX[pl][ggZZ][fs]->SetFillColor(kAzure-3);
      h1_4ljjsel_sidebandSXDX[pl][ggZZ][fs]->SetLineColor(kBlue+2);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][ggZZ][fs]); //add to hs
      // qqZZ process
      h1_4ljjsel_sidebandSXDX[pl][qqZZ][fs]->SetFillColor(kAzure+6);
      h1_4ljjsel_sidebandSXDX[pl][qqZZ][fs]->SetLineColor(kAzure-6);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][qqZZ][fs]); //add to hs
      // SM Higgs processes: ggH + VBF + VH + ttH + bbH
      h1_4ljjsel_sidebandSXDX[pl][ggH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSXDX[pl][ggH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][ggH][fs]); //add to hs
      h1_4ljjsel_sidebandSXDX[pl][VBF][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSXDX[pl][VBF][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][VBF][fs]); //add to hs
      h1_4ljjsel_sidebandSXDX[pl][VH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSXDX[pl][VH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][VH][fs]); //add to hs
      h1_4ljjsel_sidebandSXDX[pl][ttH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSXDX[pl][ttH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][ttH][fs]); //add to hs
      h1_4ljjsel_sidebandSXDX[pl][bbH][fs]->SetFillColor(kViolet+6);
      h1_4ljjsel_sidebandSXDX[pl][bbH][fs]->SetLineColor(kViolet+6);
      hs_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][bbH][fs]); //add to hs
      // HH signal
      h1_4ljjsel_sidebandSXDX[pl][HH][fs]->SetLineColor(kRed);
      h1_4ljjsel_sidebandSXDX[pl][HH][fs]->SetLineWidth(2);
      h1_4ljjsel_sidebandSXDX[pl][HH][fs]->Scale(100.);
      // data
      h1_4ljjsel_sidebandSXDX[pl][Data][fs]->SetMarkerColor(kBlack);
      h1_4ljjsel_sidebandSXDX[pl][Data][fs]->SetLineColor(kBlack);
      h1_4ljjsel_sidebandSXDX[pl][Data][fs]->SetMarkerStyle(20);

      // --- upper plot pad
      pad1_4ljjsel_sidebandSXDX[pl][fs] = new TPad("pad1_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs],"pad1_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs], 0, 0.3, 1, 1.0);
      pad1_4ljjsel_sidebandSXDX[pl][fs]->Draw();
      pad1_4ljjsel_sidebandSXDX[pl][fs]->cd();

      hs_4ljjsel_sidebandSXDX[pl][fs]->SetMaximum(10e04);
      hs_4ljjsel_sidebandSXDX[pl][fs]->SetMinimum(10e-04);
        
      hs_4ljjsel_sidebandSXDX[pl][fs]->Draw("histo");
      h1_4ljjsel_sidebandSXDX[pl][HH][fs]->Draw("histosame");
      h1_4ljjsel_sidebandSXDX[pl][Data][fs]->Draw("samepe");

      hs_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetLabelFont(43);
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetLabelSize(15);
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetTitle(h1_4ljjsel_sidebandSXDX[pl][HH][fs]->GetXaxis()->GetTitle());
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetTitleSize(20);
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetTitleFont(43);
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetLabelFont(43);
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetLabelSize(15);
      hs_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetTitle(h1_4ljjsel_sidebandSXDX[pl][HH][fs]->GetYaxis()->GetTitle());

      // --- legend
      leg_4ljjsel_sidebandSXDX[pl][fs] = new TLegend(0.78,0.61,0.94,0.87);
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][Data][fs], "Data",          "lp");
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][HH][fs],   "HH->4lbb x100", "f");
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][ggH][fs],  "SM Higgs",      "f");
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][qqZZ][fs], "qq->ZZ",        "f");
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][ggZZ][fs], "gg->ZZ",        "f");
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][TTZ][fs],  "TTV; V=Z,W",    "f");
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][ZXbkg][fs],"Z+X",           "f");
      leg_4ljjsel_sidebandSXDX[pl][fs]->AddEntry(h1_4ljjsel_sidebandSXDX[pl][VVV][fs],  "VVV; V=Z,W",    "f");
      leg_4ljjsel_sidebandSXDX[pl][fs]->SetFillColor(kWhite);
      leg_4ljjsel_sidebandSXDX[pl][fs]->SetLineColor(kBlack);
      leg_4ljjsel_sidebandSXDX[pl][fs]->SetTextFont(43);
      leg_4ljjsel_sidebandSXDX[pl][fs]->Draw();

      c_4ljjsel_sidebandSXDX[pl][fs]->Update();

      pad1_4ljjsel_sidebandSXDX[pl][fs]->SetLogy();

      c_4ljjsel_sidebandSXDX[pl][fs]->Update();

      // --- tot hist for all MC
      hMCtot_4ljjsel_sidebandSXDX[pl][fs] = (TH1F*)h1_4ljjsel_sidebandSXDX[pl][ggH][fs]->Clone("hMCtot_4ljjsel_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][VBF][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][VH][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][ttH][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][bbH][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][qqZZ][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][ggZZ][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][TTZ][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][TTW][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][ZXbkg][fs]);
      hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Add(h1_4ljjsel_sidebandSXDX[pl][VVV][fs]);

      // --- lower pad plot
      c_4ljjsel_sidebandSXDX[pl][fs]->cd();
      pad2_4ljjsel_sidebandSXDX[pl][fs] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
      pad2_4ljjsel_sidebandSXDX[pl][fs]->SetGridy();
      pad2_4ljjsel_sidebandSXDX[pl][fs]->Draw();
      pad2_4ljjsel_sidebandSXDX[pl][fs]->cd();

      // --- define ratio plot
      rp_4ljjsel_sidebandSXDX[pl][fs] = (TH1F*)h1_4ljjsel_sidebandSXDX[pl][Data][fs]->Clone("rp_4ljjsel_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs]);
      rp_4ljjsel_sidebandSXDX[pl][fs]->SetLineColor(kBlack);
      rp_4ljjsel_sidebandSXDX[pl][fs]->SetMinimum(0.);
      rp_4ljjsel_sidebandSXDX[pl][fs]->SetMaximum(2.);
      rp_4ljjsel_sidebandSXDX[pl][fs]->SetStats(0);
      rp_4ljjsel_sidebandSXDX[pl][fs]->Divide(hMCtot_4ljjsel_sidebandSXDX[pl][fs]); //divide histo rp/MC
      rp_4ljjsel_sidebandSXDX[pl][fs]->SetMarkerStyle(20);
      rp_4ljjsel_sidebandSXDX[pl][fs]->SetMarkerColor(kBlack);
      rp_4ljjsel_sidebandSXDX[pl][fs]->SetTitle("");

      rp_4ljjsel_sidebandSXDX[pl][fs]->SetYTitle("Data/#Sigma bkg");
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetNdivisions(505);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetTitleSize(20);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetTitleFont(43);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetLabelFont(43);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetYaxis()->SetLabelSize(15);

      rp_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetTitleSize(20);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetTitleFont(43);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetTitleOffset(4.);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetLabelFont(43);
      rp_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->SetLabelSize(15);

      // --- define mc shadow unc plot
      hUncMC_4ljjsel_sidebandSXDX[pl][fs] = (TH1F*)hMCtot_4ljjsel_sidebandSXDX[pl][fs]->Clone("hUncMC_4ljjsel_"+sPlots_sidebandSXDX[pl]+"_"+sFinalState[fs]);
      for(int xbin=1; xbin < hUncMC_4ljjsel_sidebandSXDX[pl][fs]->GetXaxis()->GetNbins() + 1; xbin++){
        float err = 0.;
        if(hUncMC_4ljjsel_sidebandSXDX[pl][fs]->GetBinContent(xbin) == 0) continue;
        err = hUncMC_4ljjsel_sidebandSXDX[pl][fs]->GetBinError(xbin) / hUncMC_4ljjsel_sidebandSXDX[pl][fs]->GetBinContent(xbin);
        hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetBinContent(xbin, 1.);
        hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetBinError(xbin, err);
      }
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetLineColor(1);
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetFillStyle(3005);
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetFillColor(kGray+3);
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetMarkerColor(1);
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetMarkerStyle(1);
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetTitle("");
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->SetStats(0);

      // ---draw
      rp_4ljjsel_sidebandSXDX[pl][fs]->Draw("ep");
      hUncMC_4ljjsel_sidebandSXDX[pl][fs]->Draw("e2 same");

      c_4ljjsel_sidebandSXDX[pl][fs]->Update();

      // --- draw CMS and lumi text
      writeExtraText = true;
      extraText      = "Preliminary";
      lumi_sqrtS     = lumiText + " (13 TeV)";
      cmsTextSize    = 0.6;
      lumiTextSize   = 0.46;
      extraOverCmsTextSize = 0.75;
      relPosX = 0.12;
      CMS_lumi(pad1_4ljjsel_sidebandSXDX[pl][fs], 0, 0);


      c_4ljjsel_sidebandSXDX[pl][fs]->SaveAs(outPath_4ljjsel_sidebandSXDX_plots + "/" + c_4ljjsel_sidebandSXDX[pl][fs]->GetName() + ".png");
      c_4ljjsel_sidebandSXDX[pl][fs]->SaveAs(outPath_4ljjsel_sidebandSXDX_plots + "/" + c_4ljjsel_sidebandSXDX[pl][fs]->GetName() + ".pdf");
      

    }
  }



} // end function doPlots_4ljjsel_sidebandSXDX






//*********************
//*** main function ***
//*********************
void analysis_4lbb_2bjet_yieldsAnd4ljjselPlots()
{


  if(REDOHISTOS) doHistos();

  printYields_forSync();

  //  printYields_forCards();

  doPlots_inputBDT();

  doPlots_inputBDT_withZX();

  doPlots_4ljjsel();

  doPlots_4ljjsel_sidebandSX();

  doPlots_4ljjsel_sidebandDX();

  doPlots_4ljjsel_sidebandSXDX();

}
