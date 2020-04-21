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
int year = 2016;
//int year = 2017;
//int year = 2018;
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
  if(year==2016){
    lumi       = 35.8; //fb-1 2016
    sYear      = "2016";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2016/";
    inDataPath = "/eos/user/a/acappati/samples_HH4lbb/samples_2016/";
  }
  else if(year==2017){
    lumi       = 41.5; //fb-1 2017
    sYear      = "2017";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2017/";
    inDataPath = "/eos/user/a/acappati/samples_HH4lbb/samples_2017/";
  }
  else if(year==2018){
    lumi       = 59.7; //fb-1 2018
    sYear      = "2018";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/";
    inDataPath = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/";
  }
  else{ 
    cout<<"wrong year selected!"<<endl;
  }
  cout<<"Year chosen: "<<year<<endl;


  //---datasets
  //  static int nDatasets = 22;
  TString datasets[] = {
    "AllData", 
    //"HH4lbb_Angela",
    "HH4lbb_Ilirjan",
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "bbH125",
    "ttH125",
    //"ZZTo4lamcatnlo",
    //"ZZTo4lext2",
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


  // define 1D histos
  // 4lsel
  TH1F* h1_m4l_4lsel           [nProcesses][nFinalStates+1];
  TH1F* h1_MET_4lsel_sidebands [nProcesses][nFinalStates+1];
  TH1F* h1_pT4l_4lsel_sidebands[nProcesses][nFinalStates+1];
  
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // 4lsel
      h1_m4l_4lsel    [pr][fs] = new TH1F("h1_m4l_4lsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events/4 GeV", 33, 70., 202.);
      h1_m4l_4lsel    [pr][fs]->Sumw2(true);
      h1_MET_4lsel_sidebands [pr][fs] = new TH1F("h1_MET_4lsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";MET (GeV); Events/5 GeV", 40, 0., 200.);
      h1_MET_4lsel_sidebands [pr][fs]->Sumw2(true);
      h1_pT4l_4lsel_sidebands[pr][fs] = new TH1F("h1_pT4l_4lsel_sidebands_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";4 leptons pT (GeV); Events/4 GeV", 25, 0., 100.);
      h1_pT4l_4lsel_sidebands[pr][fs]->Sumw2(true);
    }
  }

  

  int currentFinalState;
  int currentProcess;  


  //--- loop over all datasets
  for(int d=0; d<nDatasets; d++){

    currentProcess = -1;

    if(datasets[d]=="AllData") currentProcess = Data;
    //if(datasets[d]=="HH4lbb_Angela") currentProcess = HH;
    if(datasets[d]=="HH4lbb_Ilirjan") currentProcess = HH;
    if(datasets[d]=="ggH125") currentProcess = ggH;
    if(datasets[d]=="VBFH125") currentProcess = VBF;
    if(datasets[d]=="WplusH125" ||
       datasets[d]=="WminusH125" || 
       datasets[d]=="ZH125") currentProcess = VH; 
    if(datasets[d]=="ttH125") currentProcess = ttH;
    if(datasets[d]=="bbH125") currentProcess = bbH;
    //    if(datasets[d]=="ZZTo4lamcatnlo") currentProcess = qqZZ;
    //    if(datasets[d]=="ZZTo4lext2") currentProcess = qqZZ;
    if(datasets[d]=="ZZTo4l") currentProcess = qqZZ;
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
    if(currentProcess == Data || currentProcess == ZXbkg){
      sum_events[d] = 1.;
      sum_BTagSF[d] = 1.;
    }
    else{
      Long64_t entries1 = inputTree[d]->GetEntries();    
      cout<<"First loop over input files to get norm for BTagSF ..."<<endl;
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



      // --- fill histos 4l sel
      h1_m4l_4lsel[currentProcess][currentFinalState]->Fill(ZZMass, eventWeight);
      // (in 4l sidebands)
      if(ZZMass<115 || ZZMass>135){ 
        h1_MET_4lsel_sidebands[currentProcess][currentFinalState]->Fill(PFMET, eventWeight);
        for(int i=0; i<LepPt->size(); i++){
          h1_pT4l_4lsel_sidebands[currentProcess][currentFinalState]->Fill(LepPt->at(i), eventWeight);
        }
      }

    

    }//end loop over tree events

  }//end loop over datasets




  //---fill inclusive yields and histos
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates; fs++){

      // (h1 histos)
      // 4lsel
      h1_m4l_4lsel           [pr][nFinalStates]->Add(h1_m4l_4lsel    [pr][fs]);
      h1_MET_4lsel_sidebands [pr][nFinalStates]->Add(h1_MET_4lsel_sidebands[pr][fs]);
      h1_pT4l_4lsel_sidebands[pr][nFinalStates]->Add(h1_pT4l_4lsel_sidebands[pr][fs]);

    }
  }


  //---save 1D histos in a root file
  TString fout_1Dhistos_name = "f_histos_h1_4lsel_" + sYear +".root";
  TFile* fout_1Dhistos = new TFile(fout_1Dhistos_name, "recreate");
  fout_1Dhistos->cd();
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      //4lsel
      h1_m4l_4lsel           [pr][fs]->Write();
      h1_MET_4lsel_sidebands [pr][fs]->Write();
      h1_pT4l_4lsel_sidebands[pr][fs]->Write();
    }
  }
  fout_1Dhistos->Close();
    

}//end doHistos function




//*********************************
//*** doPlots_4lsel function ***
//*********************************
void doPlots_4lsel(){

 cout<<"do plots after 4l sel..."<<endl;

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


  TString outPath_4lselplots = "plots_4lsel_" + sYear;
  cout<<"creating output dir "<<outPath_4lselplots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_4lselplots)).c_str()); // create output dir



  // retrieve histos from file
  TString inFileName = "f_histos_h1_4lsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  Int_t nPlots = 3;
  TString sPlots[] = {
    "m4l_4lsel",
    "MET_4lsel_sidebands",
    "pT4l_4lsel_sidebands",
  };
  TH1F* h1_4lsel[nPlots][nProcesses][nFinalStates+1];
  for(int pl=0; pl<nPlots; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates+1; fs++){
        h1_4lsel[pl][pr][fs] = (TH1F*)fInhistos->Get("h1_"+sPlots[pl]+"_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear);
        cout<<h1_4lsel[pl][pr][fs]->GetName()<<endl;
      }
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_4lsel      [nPlots][nFinalStates+1];
  THStack* hs_4lsel     [nPlots][nFinalStates+1];
  TPad*    pad1_4lsel   [nPlots][nFinalStates+1];
  TLegend* leg_4lsel    [nPlots][nFinalStates+1];
  TH1F*    hMCtot_4lsel [nPlots][nFinalStates+1];
  TPad*    pad2_4lsel   [nPlots][nFinalStates+1];
  TH1F*    rp_4lsel     [nPlots][nFinalStates+1];
  TH1F*    hUncMC_4lsel [nPlots][nFinalStates+1];

  
  //4lsel plots
  for(int pl=0; pl<nPlots; pl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      // canvas
      c_4lsel[pl][fs] = new TCanvas("c_"+sPlots[pl]+"_"+sFinalState[fs]+"_"+sYear,"c_"+sPlots[pl]+"_"+sFinalState[fs]+"_"+sYear,800,800);
      // hstack
      hs_4lsel[pl][fs] = new THStack("hs_"+sPlots[pl]+"_"+sFinalState[fs],"");
      // VVV process
      h1_4lsel[pl][VVV][fs]->SetFillColor(kGreen-3);
      h1_4lsel[pl][VVV][fs]->SetLineColor(kGreen-1);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][VVV][fs]); //add to hs
      // Z+X process
      h1_4lsel[pl][ZXbkg][fs]->SetFillColor(kGreen+3);
      h1_4lsel[pl][ZXbkg][fs]->SetLineColor(kGreen+4);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][ZXbkg][fs]); //add to hs
      // TTV process: TTW + TTV
      h1_4lsel[pl][TTW][fs]->SetFillColor(kBlue+3);
      h1_4lsel[pl][TTW][fs]->SetLineColor(kBlue+3);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][TTW][fs]); //add to hs
      h1_4lsel[pl][TTZ][fs]->SetFillColor(kBlue+3);
      h1_4lsel[pl][TTZ][fs]->SetLineColor(kBlue+3);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][TTZ][fs]); //add to hs
      // ggZZ process
      h1_4lsel[pl][ggZZ][fs]->SetFillColor(kAzure-3);
      h1_4lsel[pl][ggZZ][fs]->SetLineColor(kBlue+2);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][ggZZ][fs]); //add to hs
      // qqZZ process
      h1_4lsel[pl][qqZZ][fs]->SetFillColor(kAzure+6);
      h1_4lsel[pl][qqZZ][fs]->SetLineColor(kAzure-6);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][qqZZ][fs]); //add to hs
      // SM Higgs processes: ggH + VBF + VH + ttH + bbH
      h1_4lsel[pl][ggH][fs]->SetFillColor(kViolet+6);
      h1_4lsel[pl][ggH][fs]->SetLineColor(kViolet+6);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][ggH][fs]); //add to hs
      h1_4lsel[pl][VBF][fs]->SetFillColor(kViolet+6);
      h1_4lsel[pl][VBF][fs]->SetLineColor(kViolet+6);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][VBF][fs]); //add to hs
      h1_4lsel[pl][VH][fs]->SetFillColor(kViolet+6);
      h1_4lsel[pl][VH][fs]->SetLineColor(kViolet+6);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][VH][fs]); //add to hs
      h1_4lsel[pl][ttH][fs]->SetFillColor(kViolet+6);
      h1_4lsel[pl][ttH][fs]->SetLineColor(kViolet+6);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][ttH][fs]); //add to hs
      h1_4lsel[pl][bbH][fs]->SetFillColor(kViolet+6);
      h1_4lsel[pl][bbH][fs]->SetLineColor(kViolet+6);
      hs_4lsel[pl][fs]->Add(h1_4lsel[pl][bbH][fs]); //add to hs
      // HH signal
      h1_4lsel[pl][HH][fs]->SetLineColor(kRed);
      h1_4lsel[pl][HH][fs]->SetLineWidth(2);
      h1_4lsel[pl][HH][fs]->Scale(100.);
      // data
      h1_4lsel[pl][Data][fs]->SetMarkerColor(kBlack);
      h1_4lsel[pl][Data][fs]->SetLineColor(kBlack);
      h1_4lsel[pl][Data][fs]->SetMarkerStyle(20);

      // --- upper plot pad
      pad1_4lsel[pl][fs] = new TPad("pad1_"+sPlots[pl]+"_"+sFinalState[fs],"pad1_"+sPlots[pl]+"_"+sFinalState[fs], 0, 0.3, 1, 1.0);
      pad1_4lsel[pl][fs]->Draw();
      pad1_4lsel[pl][fs]->cd();

      hs_4lsel[pl][fs]->SetMaximum(10e05);
      hs_4lsel[pl][fs]->SetMinimum(10e-01);
        
      hs_4lsel[pl][fs]->Draw("histo");
      h1_4lsel[pl][HH][fs]->Draw("histosame");
      h1_4lsel[pl][Data][fs]->Draw("samepe");

      hs_4lsel[pl][fs]->GetXaxis()->SetLabelFont(43);
      hs_4lsel[pl][fs]->GetXaxis()->SetLabelSize(15);
      hs_4lsel[pl][fs]->GetXaxis()->SetTitle(h1_4lsel[pl][HH][fs]->GetXaxis()->GetTitle());
      hs_4lsel[pl][fs]->GetYaxis()->SetTitleSize(20);
      hs_4lsel[pl][fs]->GetYaxis()->SetTitleFont(43);
      hs_4lsel[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      hs_4lsel[pl][fs]->GetYaxis()->SetLabelFont(43);
      hs_4lsel[pl][fs]->GetYaxis()->SetLabelSize(15);
      hs_4lsel[pl][fs]->GetYaxis()->SetTitle(h1_4lsel[pl][HH][fs]->GetYaxis()->GetTitle());

      // --- legend
      leg_4lsel[pl][fs] = new TLegend(0.78,0.61,0.94,0.87);
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][Data][fs], "Data",          "lp");
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][HH][fs],   "HH->4lbb x100", "f");
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][ggH][fs],  "SM Higgs",      "f");
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][qqZZ][fs], "qq->ZZ",        "f");
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][ggZZ][fs], "gg->ZZ",        "f");
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][TTZ][fs],  "TTV; V=Z,W",    "f");
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][ZXbkg][fs],"Z+X",           "f");
      leg_4lsel[pl][fs]->AddEntry(h1_4lsel[pl][VVV][fs],  "VVV; V=Z,W",    "f");
      leg_4lsel[pl][fs]->SetFillColor(kWhite);
      leg_4lsel[pl][fs]->SetLineColor(kBlack);
      leg_4lsel[pl][fs]->SetTextFont(43);
      leg_4lsel[pl][fs]->Draw();

      c_4lsel[pl][fs]->Update();

      pad1_4lsel[pl][fs]->SetLogy();

      c_4lsel[pl][fs]->Update();

      // --- tot hist for all MC
      hMCtot_4lsel[pl][fs] = (TH1F*)h1_4lsel[pl][ggH][fs]->Clone("hMCtot_4lsel_"+sPlots[pl]+"_"+sFinalState[fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][VBF][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][VH][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][ttH][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][bbH][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][qqZZ][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][ggZZ][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][TTZ][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][TTW][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][ZXbkg][fs]);
      hMCtot_4lsel[pl][fs]->Add(h1_4lsel[pl][VVV][fs]);

      // --- lower pad plot
      c_4lsel[pl][fs]->cd();
      pad2_4lsel[pl][fs] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
      pad2_4lsel[pl][fs]->SetGridy();
      pad2_4lsel[pl][fs]->Draw();
      pad2_4lsel[pl][fs]->cd();

      // --- define ratio plot
      rp_4lsel[pl][fs] = (TH1F*)h1_4lsel[pl][Data][fs]->Clone("rp_4lsel_"+sPlots[pl]+"_"+sFinalState[fs]);
      rp_4lsel[pl][fs]->SetLineColor(kBlack);
      rp_4lsel[pl][fs]->SetMinimum(0.);
      rp_4lsel[pl][fs]->SetMaximum(2.);
      rp_4lsel[pl][fs]->SetStats(0);
      rp_4lsel[pl][fs]->Divide(hMCtot_4lsel[pl][fs]); //divide histo rp/MC
      rp_4lsel[pl][fs]->SetMarkerStyle(20);
      rp_4lsel[pl][fs]->SetMarkerColor(kBlack);
      rp_4lsel[pl][fs]->SetTitle("");

      rp_4lsel[pl][fs]->SetYTitle("Data/MC");
      rp_4lsel[pl][fs]->GetYaxis()->SetNdivisions(505);
      rp_4lsel[pl][fs]->GetYaxis()->SetTitleSize(20);
      rp_4lsel[pl][fs]->GetYaxis()->SetTitleFont(43);
      rp_4lsel[pl][fs]->GetYaxis()->SetTitleOffset(1.4);
      rp_4lsel[pl][fs]->GetYaxis()->SetLabelFont(43);
      rp_4lsel[pl][fs]->GetYaxis()->SetLabelSize(15);

      rp_4lsel[pl][fs]->GetXaxis()->SetTitleSize(20);
      rp_4lsel[pl][fs]->GetXaxis()->SetTitleFont(43);
      rp_4lsel[pl][fs]->GetXaxis()->SetTitleOffset(4.);
      rp_4lsel[pl][fs]->GetXaxis()->SetLabelFont(43);
      rp_4lsel[pl][fs]->GetXaxis()->SetLabelSize(15);

      // --- define mc shadow unc plot
      hUncMC_4lsel[pl][fs] = (TH1F*)hMCtot_4lsel[pl][fs]->Clone("hUncMC_4lsel_"+sPlots[pl]+"_"+sFinalState[fs]);
      for(int xbin=1; xbin < hUncMC_4lsel[pl][fs]->GetXaxis()->GetNbins() + 1; xbin++){
        float err = 0.;
        if(hUncMC_4lsel[pl][fs]->GetBinContent(xbin) == 0) continue;
        err = hUncMC_4lsel[pl][fs]->GetBinError(xbin) / hUncMC_4lsel[pl][fs]->GetBinContent(xbin);
        hUncMC_4lsel[pl][fs]->SetBinContent(xbin, 1.);
        hUncMC_4lsel[pl][fs]->SetBinError(xbin, err);
      }
      hUncMC_4lsel[pl][fs]->SetLineColor(1);
      hUncMC_4lsel[pl][fs]->SetFillStyle(3005);
      hUncMC_4lsel[pl][fs]->SetFillColor(kGray+3);
      hUncMC_4lsel[pl][fs]->SetMarkerColor(1);
      hUncMC_4lsel[pl][fs]->SetMarkerStyle(1);
      hUncMC_4lsel[pl][fs]->SetTitle("");
      hUncMC_4lsel[pl][fs]->SetStats(0);

      // ---draw
      rp_4lsel[pl][fs]->Draw("ep");
      hUncMC_4lsel[pl][fs]->Draw("e2 same");

      c_4lsel[pl][fs]->Update();

      // --- draw CMS and lumi text
      writeExtraText = true;
      extraText      = "Preliminary";
      lumi_sqrtS     = lumiText + " (13 TeV)";
      cmsTextSize    = 0.6;
      lumiTextSize   = 0.46;
      extraOverCmsTextSize = 0.75;
      relPosX = 0.12;
      CMS_lumi(pad1_4lsel[pl][fs], 0, 0);


      c_4lsel[pl][fs]->SaveAs(outPath_4lselplots + "/" + c_4lsel[pl][fs]->GetName() + ".png");
      c_4lsel[pl][fs]->SaveAs(outPath_4lselplots + "/" + c_4lsel[pl][fs]->GetName() + ".pdf");
      

    }
  }



} // end function doPlots_4lsel




//*********************
//*** main function ***
//*********************
void analysis_4lbb_2bjet_4lselPlots()
{


  if(REDOHISTOS) doHistos();

  doPlots_4lsel();

}
