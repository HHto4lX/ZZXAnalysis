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

// includes for btag SF
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#include "evalEventSF.C"


using namespace std;



void reweightingBTAG(){

  ///////////// SET UP B-TAG CALIBRATION ///////////////
    
  // set up calibration + reader
  cout << "Loading the .csv file..." << endl;
    
  std::string inputCSVfile        = "../../data/BTagging/DeepCSV_102XSF_V1.csv";  
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




  // variables for tree 

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
  //  inputTree->SetBranchAddress("xsec", &xsec)
  xsec = 0.00001017; //pb Angela

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
    vector<float>          JetpT_; 
    vector<float>          Jeteta_; 
    vector<float>          Jethadronflavour_; 
    vector<float>          JetBtagger_; 

    int d1_ = -999; // position of higest btagger jet
    int d2_ = -999; // position of highest pT jet


    for (UInt_t j = 0; j < JetPt->size(); j++)
    {
      if ( (fabs ( JetEta->at(j) ) > 2.4) || (JetPt->at(j) < 20 ) ) continue; // pt cut 20GeV from ntuplizer 
	  
      TLorentzVector temp;
      temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
      JetVec_          .push_back(temp);
      JetpT_           .push_back(JetPt->at(j));
      Jeteta_          .push_back(JetEta->at(j));
      Jethadronflavour_.push_back(JetHadronFlavour->at(j));
      JetBtagger_      .push_back(JetBTagger->at(j));
    }


    // at least 2 jets in the acceptance
    if (JetVec_.size() < 2) continue;


    cout<<"ciao"<<endl;


    double * scaleFactors;
    scaleFactors = evalEventSF( int(JetVec_.size()), Jethadronflavour_, Jeteta_, JetpT_, JetBtagger_, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down );

    cout<<"SF "<<scaleFactors[0]<<endl;

    // fill eventweight
    Double_t eventWeight = partialSampleWeight * xsec * overallEventWeight ;  //kfactor e l1prefiring non ci sono, tanto uso solo HH sample


  }// end loop over tree entries

  




}
