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
// PROCESSES
enum Process {Data=0, TTV=1, VVV=2, DY=3, TTbar=4}; 
const int nProcesses = 5;
TString sProcess[nProcesses] = {"Data", "TTV", "VVV", "DY", "TTbar"};
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
  TString datasets[] = {
    "AllData", 
    "TTZToLLNuNu_M10",
    "TTWJetsToLNu",
    "WWZ",
    "WZZ",
    "ZZZ",
    "DYJetsToLL_M50",
    "TTTo2L2Nu",
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

  Short_t Zsel;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPt = 0;
  Float_t ZMass;
  vector<Float_t> *JetPt     = 0;
  vector<Float_t> *JetEta    = 0;
  vector<Float_t> *JetPhi    = 0;
  vector<Float_t> *JetMass   = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetHadronFlavour = 0;


  // define 1D histos
  // 2l2jsel
  TH1F* h1_mll_2l2jsel         [nProcesses];
  TH1F* h1_mjj_2l2jsel         [nProcesses];
  TH1F* h1_j1btag_2l2jsel      [nProcesses];
  TH1F* h1_j2btag_2l2jsel      [nProcesses];
  TH1F* h1_j1pT_20-200_2l2jsel [nProcesses];
  TH1F* h1_j2pT_20-200_2l2jsel [nProcesses];
  TH1F* h1_j1pT_20-1200_2l2jsel[nProcesses];
  TH1F* h1_j2pT_20-1200_2l2jsel[nProcesses];


  
  for(int pr=0; pr<nProcesses; pr++){

    h1_mll_2l2jsel         [pr] = new TH1F("h1_mll_2l2jsel_"+sProcess[pr]+"_"+sYear,";m_{ll}; Events/1 GeV", 60, 60., 120.);
    h1_mll_2l2jsel         [pr]->Sumw2(true);
    h1_mjj_2l2jsel         [pr] = new TH1F("h1_mjj_2l2jsel_"+sProcess[pr]+"_"+sYear,";m_{jj} (GeV); Events/5 GeV", 40, 0., 200.); 
    h1_mjj_2l2jsel         [pr]->Sumw2(true);
    h1_j1btag_2l2jsel      [pr] = new TH1F("h1_j1btag_2l2jsel_"+sProcess[pr]+"_"+sYear,";j1 DeepCSV; Events/0.04", 25, 0., 1.);
    h1_j1btag_2l2jsel      [pr]->Sumw2(true);
    h1_j2btag_2l2jsel      [pr] = new TH1F("h1_j2btag_2l2jsel_"+sProcess[pr]+"_"+sYear,";j2 DeepCSV; Events/0.04", 25, 0., 1.);
    h1_j2btag_2l2jsel      [pr]->Sumw2(true);
    h1_j1pT_20-200_2l2jsel [pr] = new TH1F("h1_j1pT_20-200_2l2jsel_"+sProcess[pr]+"_"+sYear,";j1 pT(GeV); Events/5 GeV", 36, 20., 200.);
    h1_j1pT_20-200_2l2jsel [pr]->Sumw2(true);
    h1_j2pT_20-200_2l2jsel [pr] = new TH1F("h1_j2pT_20-200_2l2jsel_"+sProcess[pr]+"_"+sYear,";j2 pT(GeV); Events/5 GeV", 36, 20., 200.);
    h1_j2pT_20-200_2l2jsel [pr]->Sumw2(true);
    h1_j1pT_20-1200_2l2jsel[pr] = new TH1F("h1_j1pT_20-1200_2l2jsel_"+sProcess[pr]+"_"+sYear,";j1 pT(GeV); Events/10 GeV", 118, 20., 1200.);
    h1_j1pT_20-1200_2l2jsel[pr]->Sumw2(true);
    h1_j2pT_20-1200_2l2jsel[pr] = new TH1F("h1_j2pT_20-1200_2l2jsel_"+sProcess[pr]+"_"+sYear,";j2 pT(GeV); Events/10 GeV", 118, 20., 1200.);
    h1_j2pT_20-1200_2l2jsel[pr]->Sumw2(true);
  }

  


  int currentProcess;  


  //--- loop over all datasets
  for(int d=0; d<nDatasets; d++){

    currentProcess = -1;

    if(datasets[d]=="AllData") currentProcess = Data;
    if(datasets[d]=="TTZToLLNuNu_M10" || 
       datasets[d]=="TTWJetsToLNu") currentProcess = TTV;
    if(datasets[d]=="WWZ" ||
       datasets[d]=="WZZ" ||
       datasets[d]=="ZZZ") currentProcess = VVV;
    if(datasets[d]=="DYJetsToLL_M50") currentProcess = DY;
    if(datasets[d]=="TTTo2L2Nu") currentProcess = TTbar;
    

    // select input file
    TString inputFileName;
    if(currentProcess == Data) inputFileName = inDataPath + datasets[d] + "/ZZXAnalysis.root";
    else inputFileName = inFilePath + datasets[d] + "/ZZXAnalysis.root";
    cout<<"Opening file "<<inputFileName<<" ..."<<endl;
    inputFile[d] = TFile::Open(inputFileName);


    hCounters[d] = (TH1F*)inputFile[d]->Get("ZTree/Counters");
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(1);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d];
    inputTree[d] = (TTree*)inputFile[d]->Get("ZTree/candTree");

    cout<<"debug "<<datasets[d]<<" "<<currentProcess<<" "<<sProcess[currentProcess]<<endl;

    // set branch addresses
    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    if(currentProcess != Data) inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    if(currentProcess != Data) inputTree[d]->SetBranchAddress("xsec", &xsec);
    if(currentProcess != Data) inputTree[d]->SetBranchAddress("L1prefiringWeight", &L1prefiringWeight);

    inputTree[d]->SetBranchAddress("Zsel", &Zsel);
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("ZMass", &ZMass);  
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("JetPhi",  &JetPhi);
    inputTree[d]->SetBranchAddress("JetBTagger",  &JetBTagger);
    inputTree[d]->SetBranchAddress("JetHadronFlavour",  &JetHadronFlavour);


    // --------------------------------------------------------
    // --- first loop over input tree to get norm for BtagSF
    Long64_t entries1 = inputTree[d]->GetEntries();    
    cout<<"First loop over input files to get norm for BTagSF ..."<<endl;
    cout<<"Processing file: "<< datasets[d] << " (" << entries1 <<" entries) ..."<< endl;    
    for(Long64_t z=0; z<entries1; z++){

      inputTree[d]->GetEntry(z);
  
      if( Zsel < 0 ) continue;


      if(datasets[d] == "AllData"){
        sum_events[d] = 1.;
    	sum_BTagSF[d] = 1.;
      }
      else{

        // compute SF
        double * scaleFactors;
        scaleFactors = evalEventSF( int(JetPt->size()), JetHadronFlavour, JetEta, JetPt, JetBTagger, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down );

        // total counters for BTagSF norm --- 2l2jsel
        sum_events[d] += 1.; 
        sum_BTagSF[d] += scaleFactors[0]; 

      }// end else
    
    } // end first loop over entries 

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

  
      if( Zsel < 0 ) continue;

      
      // compute SF
      double * scaleFactors;
      if(currentProcess != Data){
        scaleFactors = evalEventSF( int(JetPt->size()), JetHadronFlavour, JetEta, JetPt, JetBTagger, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down );
      }

      // --- event weights
      Double_t eventWeight = 1.;
      if(currentProcess != Data) eventWeight = partialSampleWeight[d] *xsec *overallEventWeight *L1prefiringWeight *scaleFactors[0] *sum_events[d] /sum_BTagSF[d];

      
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


      
      // --- fill histos
      h1_j1btag_2l2jsel[currentProcess]->Fill(JetBTagger->at(dj1_), eventWeight);
      h1_j2btag_2l2jsel[currentProcess]->Fill(JetBTagger->at(dj2_), eventWeight);
      h1_j1pT_2l2jsel  [currentProcess]->Fill(JetPt->at(dj1_),      eventWeight); 
      h1_j2pT_2l2jsel  [currentProcess]->Fill(JetPt->at(dj2_),      eventWeight);
      h1_mjj_2l2jsel   [currentProcess]->Fill(bbMass,               eventWeight);
      h1_mll_2l2jsel   [currentProcess]->Fill(ZMass,                eventWeight);
    


    

    }//end loop over tree events

  }//end loop over datasets




  //---save 1D histos in a root file
  TString fout_1Dhistos_name = "f_histos_h1_2l2jsel_" + sYear +".root";
  TFile* fout_1Dhistos = new TFile(fout_1Dhistos_name, "recreate");
  fout_1Dhistos->cd();
  for(int pr=0; pr<nProcesses; pr++){
    h1_j1btag_2l2jsel[pr]->Write();
    h1_j2btag_2l2jsel[pr]->Write();
    h1_j1pT_2l2jsel  [pr]->Write();
    h1_j2pT_2l2jsel  [pr]->Write();
    h1_mjj_2l2jsel   [pr]->Write();
    h1_mll_2l2jsel   [pr]->Write();
  }
  fout_1Dhistos->Close();
    

}//end doHistos function




//*********************************
//*** doPlots_2l2jsel function ***
//*********************************
void doPlots_2l2jsel(){

 cout<<"do plots after 2l2j sel..."<<endl;

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


  TString outPath_2l2jselplots = "plots_2l2jsel_" + sYear;
  cout<<"creating output dir "<<outPath_2l2jselplots<<" ... "<<endl;
  gSystem->Exec(("mkdir -p "+string(outPath_2l2jselplots)).c_str()); // create output dir



  // retrieve histos from file
  TString inFileName = "f_histos_h1_2l2jsel_" + sYear + ".root";     
  cout<<"Retrieving Data and MC histograms from file "<<inFileName<<" ..."<<endl;
  TFile* fInhistos = TFile::Open(inFileName);
  
  // --- take histos from file
  Int_t nPlots = 6;
  TString sPlots[] = {
    "mll_2l2jsel",
    "mjj_2l2jsel",
    "j1btag_2l2jsel",
    "j2btag_2l2jsel",
    "j1pT_20-200_2l2jsel",
    "j2pT_20-200_2l2jsel",
    "j1pT_20-1200_2l2jsel",
    "j2pT_20-1200_2l2jsel",
  };
  TH1F* h1_2l2jsel[nPlots][nProcesses];
  for(int pl=0; pl<nPlots; pl++){
    for(int pr=0; pr<nProcesses; pr++){
      h1_2l2jsel[pl][pr] = (TH1F*)fInhistos->Get("h1_"+sPlots[pl]+"_"+sProcess[pr]+"_"+sYear);
      cout<<h1_2l2jsel[pl][pr]->GetName()<<endl;
    }
  }   

  // --- define canvas, hstack and pads for BDT input plots
  TCanvas* c_2l2jsel      [nPlots];
  THStack* hs_2l2jsel     [nPlots];
  TPad*    pad1_2l2jsel   [nPlots];
  TLegend* leg_2l2jsel    [nPlots];
  TH1F*    hMCtot_2l2jsel [nPlots];
  TPad*    pad2_2l2jsel   [nPlots];
  TH1F*    rp_2l2jsel     [nPlots];
  TH1F*    hUncMC_2l2jsel [nPlots];

  
  //2l2jsel plots
  for(int pl=0; pl<nPlots; pl++){
    // canvas
    c_2l2jsel[pl] = new TCanvas("c_"+sPlots[pl]+"_"+sYear,"c_"+sPlots[pl]+"_"+sYear,800,800);
    // hstack
    hs_2l2jsel[pl] = new THStack("hs_"+sPlots[pl],"");
    // TTbar process
    h1_2l2jsel[pl][TTbar]->SetFillColor(kRed);
    h1_2l2jsel[pl][TTbar]->SetLineColor(kRed+4);
    hs_2l2jsel[pl]->Add(h1_2l2jsel[pl][TTbar]); //add to hs
    // TTV process
    h1_2l2jsel[pl][TTV]->SetFillColor(kBlue+3);
    h1_2l2jsel[pl][TTV]->SetLineColor(kBlue+3);
    hs_2l2jsel[pl]->Add(h1_2l2jsel[pl][TTV]); //add to hs
    // VVV process
    h1_2l2jsel[pl][VVV]->SetFillColor(kGreen-3);
    h1_2l2jsel[pl][VVV]->SetLineColor(kGreen-1);
    hs_2l2jsel[pl]->Add(h1_2l2jsel[pl][VVV]); //add to hs
    // DY process
    h1_2l2jsel[pl][DY]->SetFillColor(kCyan);
    h1_2l2jsel[pl][DY]->SetLineColor(kCyan+4);
    hs_2l2jsel[pl]->Add(h1_2l2jsel[pl][DY]); //add to hs

    // data
    h1_2l2jsel[pl][Data]->SetMarkerColor(kBlack);
    h1_2l2jsel[pl][Data]->SetLineColor(kBlack);
    h1_2l2jsel[pl][Data]->SetMarkerStyle(20);

    // --- upper plot pad
    pad1_2l2jsel[pl] = new TPad("pad1_"+sPlots[pl],"pad1_"+sPlots[pl], 0, 0.3, 1, 1.0);
    pad1_2l2jsel[pl]->Draw();
    pad1_2l2jsel[pl]->cd();

    hs_2l2jsel[pl]->SetMaximum(10e07);
    hs_2l2jsel[pl]->SetMinimum(10e-01);
      
    hs_2l2jsel[pl]->Draw("histo");
    h1_2l2jsel[pl][Data]->Draw("samepe");

    hs_2l2jsel[pl]->GetXaxis()->SetLabelFont(43);
    hs_2l2jsel[pl]->GetXaxis()->SetLabelSize(15);
    hs_2l2jsel[pl]->GetXaxis()->SetTitle(h1_2l2jsel[pl][DY]->GetXaxis()->GetTitle());
    hs_2l2jsel[pl]->GetYaxis()->SetTitleSize(20);
    hs_2l2jsel[pl]->GetYaxis()->SetTitleFont(43);
    hs_2l2jsel[pl]->GetYaxis()->SetTitleOffset(1.4);
    hs_2l2jsel[pl]->GetYaxis()->SetLabelFont(43);
    hs_2l2jsel[pl]->GetYaxis()->SetLabelSize(15);
    hs_2l2jsel[pl]->GetYaxis()->SetTitle(h1_2l2jsel[pl][DY]->GetYaxis()->GetTitle());

    // --- legend
    leg_2l2jsel[pl] = new TLegend(0.78,0.65,0.94,0.87);
    leg_2l2jsel[pl]->AddEntry(h1_2l2jsel[pl][Data], "Data",       "lp");
    leg_2l2jsel[pl]->AddEntry(h1_2l2jsel[pl][DY],   "DY",         "f");
    leg_2l2jsel[pl]->AddEntry(h1_2l2jsel[pl][VVV],  "VVV; V=Z,W", "f");
    leg_2l2jsel[pl]->AddEntry(h1_2l2jsel[pl][TTV],  "TTV; V=Z,W", "f");
    leg_2l2jsel[pl]->AddEntry(h1_2l2jsel[pl][TTbar],"tt",         "f");
    leg_2l2jsel[pl]->SetFillColor(kWhite);
    leg_2l2jsel[pl]->SetLineColor(kBlack);
    leg_2l2jsel[pl]->SetTextFont(43);
    leg_2l2jsel[pl]->Draw();

    c_2l2jsel[pl]->Update();

    pad1_2l2jsel[pl]->SetLogy();

    c_2l2jsel[pl]->Update();

    // --- tot hist for all MC
    hMCtot_2l2jsel[pl] = (TH1F*)h1_2l2jsel[pl][TTbar]->Clone("hMCtot_2l2jsel_"+sPlots[pl]);
    hMCtot_2l2jsel[pl]->Add(h1_2l2jsel[pl][TTV]);
    hMCtot_2l2jsel[pl]->Add(h1_2l2jsel[pl][VVV]);
    hMCtot_2l2jsel[pl]->Add(h1_2l2jsel[pl][DY]);

    // --- lower pad plot
    c_2l2jsel[pl]->cd();
    pad2_2l2jsel[pl] = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
    pad2_2l2jsel[pl]->SetGridy();
    pad2_2l2jsel[pl]->Draw();
    pad2_2l2jsel[pl]->cd();

    // --- define ratio plot
    rp_2l2jsel[pl] = (TH1F*)h1_2l2jsel[pl][Data]->Clone("rp_2l2jsel_"+sPlots[pl]);
    rp_2l2jsel[pl]->SetLineColor(kBlack);
    rp_2l2jsel[pl]->SetMinimum(0.);
    rp_2l2jsel[pl]->SetMaximum(2.);
    rp_2l2jsel[pl]->SetStats(0);
    rp_2l2jsel[pl]->Divide(hMCtot_2l2jsel[pl]); //divide histo rp/MC
    rp_2l2jsel[pl]->SetMarkerStyle(20);
    rp_2l2jsel[pl]->SetMarkerColor(kBlack);
    rp_2l2jsel[pl]->SetTitle("");

    rp_2l2jsel[pl]->SetYTitle("Data/MC");
    rp_2l2jsel[pl]->GetYaxis()->SetNdivisions(505);
    rp_2l2jsel[pl]->GetYaxis()->SetTitleSize(20);
    rp_2l2jsel[pl]->GetYaxis()->SetTitleFont(43);
    rp_2l2jsel[pl]->GetYaxis()->SetTitleOffset(1.4);
    rp_2l2jsel[pl]->GetYaxis()->SetLabelFont(43);
    rp_2l2jsel[pl]->GetYaxis()->SetLabelSize(15);

    rp_2l2jsel[pl]->GetXaxis()->SetTitleSize(20);
    rp_2l2jsel[pl]->GetXaxis()->SetTitleFont(43);
    rp_2l2jsel[pl]->GetXaxis()->SetTitleOffset(4.);
    rp_2l2jsel[pl]->GetXaxis()->SetLabelFont(43);
    rp_2l2jsel[pl]->GetXaxis()->SetLabelSize(15);

    // --- define mc shadow unc plot
    hUncMC_2l2jsel[pl] = (TH1F*)hMCtot_2l2jsel[pl]->Clone("hUncMC_2l2jsel_"+sPlots[pl]);
    for(int xbin=1; xbin < hUncMC_2l2jsel[pl]->GetXaxis()->GetNbins() + 1; xbin++){
      float err = 0.;
      if(hUncMC_2l2jsel[pl]->GetBinContent(xbin) == 0) continue;
      err = hUncMC_2l2jsel[pl]->GetBinError(xbin) / hUncMC_2l2jsel[pl]->GetBinContent(xbin);
      hUncMC_2l2jsel[pl]->SetBinContent(xbin, 1.);
      hUncMC_2l2jsel[pl]->SetBinError(xbin, err);
    }
    hUncMC_2l2jsel[pl]->SetLineColor(1);
    hUncMC_2l2jsel[pl]->SetFillStyle(3005);
    hUncMC_2l2jsel[pl]->SetFillColor(kGray+3);
    hUncMC_2l2jsel[pl]->SetMarkerColor(1);
    hUncMC_2l2jsel[pl]->SetMarkerStyle(1);
    hUncMC_2l2jsel[pl]->SetTitle("");
    hUncMC_2l2jsel[pl]->SetStats(0);

    // ---draw
    rp_2l2jsel[pl]->Draw("ep");
    hUncMC_2l2jsel[pl]->Draw("e2 same");

    c_2l2jsel[pl]->Update();

    // --- draw CMS and lumi text
    writeExtraText = true;
    extraText      = "Preliminary";
    lumi_sqrtS     = lumiText + " (13 TeV)";
    cmsTextSize    = 0.6;
    lumiTextSize   = 0.46;
    extraOverCmsTextSize = 0.75;
    relPosX = 0.12;
    CMS_lumi(pad1_2l2jsel[pl], 0, 0);


    c_2l2jsel[pl]->SaveAs(outPath_2l2jselplots + "/" + c_2l2jsel[pl]->GetName() + ".png");
    c_2l2jsel[pl]->SaveAs(outPath_2l2jselplots + "/" + c_2l2jsel[pl]->GetName() + ".pdf");
      

  }



} // end function doPlots_2l2jsel




//*********************
//*** main function ***
//*********************
void analysis_4lbb_2bjet_2l2jselPlots()
{


  if(REDOHISTOS) doHistos();

  doPlots_2l2jsel();

}
