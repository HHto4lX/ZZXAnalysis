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
enum Process {Data=0, HH=1, ggH=2, VBF=3, VH=4, ttH=5, bbH=6, qqZZ=7, ggZZ=8, TTZ=9, TTW=10, VVV=11, HWW=12, ZXbkg=13}; 
const int nProcesses = 14;
TString sProcess[nProcesses] = {"Data", "HH", "ggH", "VBF", "VH", "ttH", "bbH", "qqZZ", "ggZZ", "TTZ", "TTW", "VVV", "HWW", "ZXbkg"};
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
    lumi       = 35.9; //fb-1 2016
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
  int nDatasets = 28;
  TString datasets[] = {
    "AllData", 
    "HH4lbb_Angela",
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "bbH125",
    "ttH125",
    "ZZTo4lamcatnlo",
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
    "HToWW125_ggH",
    "HToWW125_VBFH",
    "HToWW125_HWplusJ",
    "HToWW125_HWminusJ",
    "HToWW125_HZJ",
    "HToWW125_bbH",
    "ZXbkg",
  };




  // input file branches
  TFile* inputFile[nDatasets];
  TTree* inputTree[nDatasets];
  TH1F* hCounters[nDatasets];
  Long64_t NGenEvt[nDatasets];
  Float_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];
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

  // define yields histos and nEvents histos (no weight)
  TH1F* hYields_4lsel  [nProcesses][nFinalStates+1];
  TH1F* hYields_4ljjsel[nProcesses][nFinalStates+1];
  TH1F* hEvents_4lsel  [nProcesses][nFinalStates+1];
  TH1F* hEvents_4ljjsel[nProcesses][nFinalStates+1];
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
    }
  }

  // define 1D histos
  TH1F* h1_m4l_4lsel    [nProcesses][nFinalStates+1]; 
  TH1F* h1_m4l_4ljjsel  [nProcesses][nFinalStates+1]; 
  TH1F* h1_mbb_4ljjsel  [nProcesses][nFinalStates+1];
  TH1F* h1_j1Eta_4ljjsel[nProcesses][nFinalStates+1];
  TH1F* h1_j2Eta_4ljjsel[nProcesses][nFinalStates+1];
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      h1_m4l_4lsel    [pr][fs] = new TH1F("h1_m4l_4lsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events/2 GeV", 65, 70., 200.);
      h1_m4l_4lsel    [pr][fs]->Sumw2(true);
      h1_m4l_4ljjsel  [pr][fs] = new TH1F("h1_m4l_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events/2 GeV", 65, 70., 200.);
      h1_m4l_4ljjsel  [pr][fs]->Sumw2(true);
      h1_mbb_4ljjsel  [pr][fs] = new TH1F("h1_mbb_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{jj} (GeV); Events/2 GeV", 100, 0., 200.); 
      h1_mbb_4ljjsel  [pr][fs]->Sumw2(true);
      h1_j1Eta_4ljjsel[pr][fs] = new TH1F("h1_j1Eta_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j1 #eta ; Events",56, -2.8, 2.8);
      h1_j1Eta_4ljjsel[pr][fs]->Sumw2(true);
      h1_j2Eta_4ljjsel[pr][fs] = new TH1F("h1_j2Eta_4ljjsel_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";j2 #eta ; Events",56, -2.8, 2.8);
      h1_j2Eta_4ljjsel[pr][fs]->Sumw2(true);      
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
    if(datasets[d]=="ggH125") currentProcess = ggH;
    if(datasets[d]=="VBFH125") currentProcess = VBF;
    if(datasets[d]=="WplusH125" ||
       datasets[d]=="WminusH125" || 
       datasets[d]=="ZH125") currentProcess = VH; 
    if(datasets[d]=="ttH125") currentProcess = ttH;
    if(datasets[d]=="bbH125") currentProcess = bbH;
    if(datasets[d]=="ZZTo4lamcatnlo") currentProcess = qqZZ;
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
    if(datasets[d]=="HToWW125_ggH" ||
       datasets[d]=="HToWW125_VBFH" ||
       datasets[d]=="HToWW125_HWplusJ" ||
       datasets[d]=="HToWW125_HWminusJ" ||
       datasets[d]=="HToWW125_HZJ" ||
       datasets[d]=="HToWW125_bbH") currentProcess = HWW;
    if(datasets[d]=="ZXbkg") currentProcess = ZXbkg;
    

    // select input file
    TString inputFileName;
    if(currentProcess == Data) inputFileName = inDataPath + datasets[d] + "/ZZXAnalysis.root";
    else inputFileName = inFilePath + datasets[d] + "/ZZXAnalysis.root";
    cout<<"Opening file "<<inputFileName<<" ..."<<endl;
    inputFile[d] = TFile::Open(inputFileName);

    if(currentProcess == ZXbkg){
      hCounters[d] = 0;
      NGenEvt[d] = 0;
      gen_sumWeights[d] = 0.;
      partialSampleWeight[d] = 0;
      inputTree[d] = (TTree*)inputFile[d]->Get("candTree");
    }
    else{
      hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
      NGenEvt[d] = (Long64_t)hCounters[d]->GetBinContent(1);
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
    inputTree[d]->SetBranchAddress("JetIsBtagged",  &JetIsBtagged);  
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

    // loop over input tree
    Long64_t entries = inputTree[d]->GetEntries();    
    cout<<"Processing file: "<< datasets[d] << " (" << entries <<" entries) ..."<< endl;    

    for(Long64_t z=0; z<entries; z++){

      inputTree[d]->GetEntry(z);

  
      if( !(ZZsel>=0) ) continue;

      // 4l selection
      if(LepEta->size() != 4){
	cout << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
      }

      
      // --- fill k factors and event weights
      Float_t kfactor = 1.;
      if(currentProcess == qqZZ) { kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; } // qqZZ sample                      
      else if(currentProcess == ggZZ) { kfactor = KFactor_QCD_ggZZ_Nominal; } //ggZZ samples 

      Double_t eventWeight = 1.;
      if(currentProcess != Data && currentProcess != ZXbkg) eventWeight = partialSampleWeight[d] * xsec * kfactor * overallEventWeight;
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

      // --- fill yields after 4l sel
      hYields_4lsel[currentProcess][currentFinalState]->Fill(0.5, eventWeight);

      // --- fill number of non weighted events selected after 4l sel
      hEvents_4lsel[currentProcess][currentFinalState]->Fill(0.5, 1.);




      // mass cut: signal region
       if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV

 


      //JETSELECTION--------------------------------------------------  
      // at least 2 jets in the acceptance
      if (JetPt->size() < 2) continue;   
          
      cout<<"jetptsize: "<<JetPt->size()<<" jetbtaggersize: "<<JetBTagger->size()<<endl;
      // index of jet with max btagger value (j1)
      int dj1_ = distance( JetBTagger->begin(), max_element(JetBTagger->begin(), JetBTagger->end()));
      
      // index of jet with 2nd highest btagger value (j2)  
      float max2JetBtag = -9999.;
      int dj2_ = -9999;
      for(UInt_t i=0; i<JetBTagger->size(); i++){
        cout<<JetBTagger->at(i)<<" ";
  
        if(i == dj1_) continue;
        float temp = JetBTagger->at(i);
        if (temp > max2JetBtag){
          max2JetBtag = temp;
          dj2_ = i;
        }
      }
      cout<<" indici (1,2) "<<dj1_ <<" "<<dj2_<<" valori btag(1); btag(2) "<<JetBTagger->at(dj1_)<<" "<<JetBTagger->at(dj2_)<<endl;
      if(JetBTagger->at(dj1_)>=JetBTagger->at(dj2_)) cout<<"FUNZIONA"<<endl;
      else cout<<"NON FUNZIONA"<<endl;


      
      // build 2 jets tlorentzvectors
      TLorentzVector tlzvec_j1_;
      tlzvec_j1_.SetPtEtaPhiM(JetPt->at(dj1_), JetEta->at(dj1_), JetPhi->at(dj1_), JetMass->at(dj1_));
      TLorentzVector tlzvec_j2_;
      tlzvec_j2_.SetPtEtaPhiM(JetPt->at(dj2_), JetEta->at(dj2_), JetPhi->at(dj2_), JetMass->at(dj2_));
      
      
      // build H->bb tlorentzvector
      TLorentzVector Hbb_tlzvec = tlzvec_j1_ + tlzvec_j2_;
      Float_t bbMass = Hbb_tlzvec.M();


      // --- fill yields after 4ljj sel
      hYields_4ljjsel[currentProcess][currentFinalState]->Fill(0.5, eventWeight);
      // --- fill number of non weighted events selected after 4ljj sel
      hEvents_4ljjsel[currentProcess][currentFinalState]->Fill(0.5, 1.);
      

      // --- fill histos after 4ljj sel
      h1_m4l_4ljjsel  [currentProcess][currentFinalState]->Fill(ZZMass, eventWeight);
      h1_mbb_4ljjsel  [currentProcess][currentFinalState]->Fill(bbMass, eventWeight);
      h1_j1Eta_4ljjsel[currentProcess][currentFinalState]->Fill(JetEta->at(dj1_), eventWeight);
      h1_j2Eta_4ljjsel[currentProcess][currentFinalState]->Fill(JetEta->at(dj2_), eventWeight);
      //2D histo
      h2_m4lvsmbb_4ljjsel[currentProcess][currentFinalState]->Fill(ZZMass, bbMass, eventWeight);



      // // build H-H DeltaR
      // float DeltaPhi = ZZPhi - Hbb_Vec.Phi();
      // if( fabs(DeltaPhi) > acos(-1) ) { DeltaPhi = (2*acos(-1)) - fabs(DeltaPhi); }
      // float DeltaEta = ZZEta - Hbb_Vec.Eta();
      // float DeltaR = sqrt( DeltaPhi*DeltaPhi + DeltaEta*DeltaEta );
  
      // // save DeltaR
      // f_deltarsignal = DeltaR;
    


      // // build HZZ tlorentzvector
      // TLorentzVector HZZ_Vec;
      // HZZ_Vec.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, ZZMass);
    
      // // build HH tlorentzvector
      // TLorentzVector HH_Vec = HZZ_Vec + Hbb_Vec;


      // // save bb and HH masses
      // f_bbmass = Hbb_Vec.M();
      // f_HHmass = HH_Vec.M();

      // // save jet jet inv mass
      // f_massjetjet = Hbb_Vec.M();

      //JETSELECTION---------------------------------------------------






    }//end loop over tree events

  }//end loop over datasets


  //---fill inclusive yields and histos
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates; fs++){
      hYields_4lsel   [pr][nFinalStates]->Add(hYields_4lsel   [pr][fs]);
      hYields_4ljjsel [pr][nFinalStates]->Add(hYields_4ljjsel [pr][fs]);
      hEvents_4lsel   [pr][nFinalStates]->Add(hEvents_4lsel   [pr][fs]);
      hEvents_4ljjsel [pr][nFinalStates]->Add(hEvents_4ljjsel [pr][fs]);
      h1_m4l_4lsel    [pr][nFinalStates]->Add(h1_m4l_4lsel    [pr][fs]);
      h1_m4l_4ljjsel  [pr][nFinalStates]->Add(h1_m4l_4ljjsel  [pr][fs]);
      h1_mbb_4ljjsel  [pr][nFinalStates]->Add(h1_mbb_4ljjsel  [pr][fs]);
      h1_j1Eta_4ljjsel[pr][nFinalStates]->Add(h1_j1Eta_4ljjsel[pr][fs]);
      h1_j2Eta_4ljjsel[pr][nFinalStates]->Add(h1_j2Eta_4ljjsel[pr][fs]);
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
    }
  }
  fout_yields->Close();

  //---save histos in a root file
  TString fout_1Dhistos_name = "f_histos_h1_" + sYear +".root";
  TFile* fout_1Dhistos = new TFile(fout_1Dhistos_name, "recreate");
  fout_1Dhistos->cd();
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      h1_m4l_4lsel    [pr][fs]->Write();
      h1_m4l_4ljjsel  [pr][fs]->Write();
      h1_mbb_4ljjsel  [pr][fs]->Write();
      h1_j1Eta_4ljjsel[pr][fs]->Write();
      h1_j2Eta_4ljjsel[pr][fs]->Write();      
    }
  }
  fout_1Dhistos->Close();
    
  //---save 2D histos in a root file
  TString fout_2Dhistos_name = "f_histos_h2_" + sYear + ".root";
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
    }
  }


  // **********************************************
  // *** PRINT NUMBER OF EVENTS SEL (no weight) ***
  // **********************************************
  // --- print number of event selected after 4l sel for sync (no weight)
  cout<<"print number of event selected after 4l sel for sync (no weight) ... "<<endl; 
  ofstream f_events4lsel_sync;
  TString f_events4lsel_sync_name = "nEvents_4lsel_perSync_"+ sYear + ".txt";
  f_events4lsel_sync.open(f_events4lsel_sync_name);
  f_events4lsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF+H->WW) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_events4lsel_sync<<"|"<<sFinalState[fs]<<" |"<<nEvent_4lsel[HH][fs]<<" |"<<nEvent_4lsel[TTZ][fs]<<" |"<<nEvent_4lsel[ttH][fs]<<" |"<<nEvent_4lsel[qqZZ][fs]+nEvent_4lsel[ggZZ][fs]<<" |"<<nEvent_4lsel[ggH][fs]+nEvent_4lsel[VBF][fs]+nEvent_4lsel[HWW][fs]<<" |"<<nEvent_4lsel[VVV][fs]+nEvent_4lsel[VH][fs]+nEvent_4lsel[TTW][fs]<<" |"<<nEvent_4lsel[ZXbkg][fs]<<" |"<<endl;
  }
  f_events4lsel_sync.close();

  // --- print number of event selected after 4ljj sel for sync (no weight)
  cout<<"print number of event selected after 4ljj sel for sync (no weight) ... "<<endl; 
  ofstream f_events4ljjsel_sync;
  TString f_events4ljjsel_sync_name = "nEvents_4ljjsel_perSync_"+ sYear + ".txt";
  f_events4ljjsel_sync.open(f_events4ljjsel_sync_name);
  f_events4ljjsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF+H->WW) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_events4ljjsel_sync<<"|"<<sFinalState[fs]<<" |"<<nEvent_4ljjsel[HH][fs]<<" |"<<nEvent_4ljjsel[TTZ][fs]<<" |"<<nEvent_4ljjsel[ttH][fs]<<" |"<<nEvent_4ljjsel[qqZZ][fs]+nEvent_4ljjsel[ggZZ][fs]<<" |"<<nEvent_4ljjsel[ggH][fs]+nEvent_4ljjsel[VBF][fs]+nEvent_4ljjsel[HWW][fs]<<" |"<<nEvent_4ljjsel[VVV][fs]+nEvent_4ljjsel[VH][fs]+nEvent_4ljjsel[TTW][fs]<<" |"<<nEvent_4ljjsel[ZXbkg][fs]<<" |"<<endl;
  }
  f_events4ljjsel_sync.close();



  // ********************
  // *** PRINT YIELDS ***
  // ********************
  // --- print yields after 4l sel for sync 
  cout<<"print yields after 4l sel for sync ... "<<endl;
  ofstream f_yields4lsel_sync;
  TString f_yields4lsel_sync_name = "yields4lsel_perSync_"+ sYear + ".txt";
  f_yields4lsel_sync.open(f_yields4lsel_sync_name);
  f_yields4lsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF+H->WW) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4lsel_sync<<"|"<<sFinalState[fs]<<" |"<<yield_4lsel[HH][fs]<<" |"<<yield_4lsel[TTZ][fs]<<" |"<<yield_4lsel[ttH][fs]<<" |"<<yield_4lsel[qqZZ][fs]+yield_4lsel[ggZZ][fs]<<" |"<<yield_4lsel[ggH][fs]+yield_4lsel[VBF][fs]+yield_4lsel[HWW][fs]<<" |"<<yield_4lsel[VVV][fs]+yield_4lsel[VH][fs]+yield_4lsel[TTW][fs]<<" |"<<yield_4lsel[ZXbkg][fs]<<" |"<<endl;
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
  f_yields4ljjsel_sync<<"|Final state |signal HH |ttZ |ttH |ZZ(=qqZZ+ggZZ) |Higgs+VBF(=ggH+VBF+H->WW) |others(=VVV+VH+TTW) |Z+X |"<<endl;
  for(int fs=0; fs<nFinalStates+1; fs++){
    f_yields4ljjsel_sync<<"|"<<sFinalState[fs]<<" |"<<yield_4ljjsel[HH][fs]<<" |"<<yield_4ljjsel[TTZ][fs]<<" |"<<yield_4ljjsel[ttH][fs]<<" |"<<yield_4ljjsel[qqZZ][fs]+yield_4ljjsel[ggZZ][fs]<<" |"<<yield_4ljjsel[ggH][fs]+yield_4ljjsel[VBF][fs]+yield_4ljjsel[HWW][fs]<<" |"<<yield_4ljjsel[VVV][fs]+yield_4ljjsel[VH][fs]+yield_4ljjsel[TTW][fs]<<" |"<<yield_4ljjsel[ZXbkg][fs]<<" |"<<endl;
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



//*********************
//*** main function ***
//*********************
void analysis_4lbb_2bjet()
{

  // string outputFilePath = "out";
  // //string outputFilePath = "200212_mvaNtuples_1bjet1Pt_4l";
  // gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir

  if(REDOHISTOS) doHistos();

  printYields_forSync();

  printYields_forCards();

}
