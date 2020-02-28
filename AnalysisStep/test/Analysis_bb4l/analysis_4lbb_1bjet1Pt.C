// ***********************
// prepare ntuples for mva
// run with:
//
// root -l -b -q analysis_4lbb_1bjet1Pt.C++
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



//******************
//int year = 2016;
//int year = 2017;
int year = 2018;
//******************



//*************************************************************************************
// FINAL STATES
enum FinalState {fs_4mu=0, fs_4e=1, fs_2mu2e=3};  // 4mu, 4e, 2e2mu
const int nFinalStates = 3;
TString sFinalState[nFinalStates+1] = {"4mu", "4e","2e2mu","4l"};
//*************************************************************************************
// PROCESSES
enum Process {Data=0, HH=1, ggH=2, VBF=3, VH=4, ttH=5, bbH=6, qqZZ=7, ggZZ=8, TTZ=9, TTW=10, VVV=11, ZX=12, HWW=13}; 
const int nProcesses = 14;
TString sProcess[nProcesses] = {"HH", "ggH", "VBF", "VH", "ttH", "bbH", "qqZZ", "ggZZ", "TTZ", "TTW", "VVV", "ZX", "HWW"};
//*************************************************************************************




void doHistos()
{

  TH1::SetDefaultSumw2(true);

  //---lumi and input path
  float lumi = 0.;
  TString inFilePath;
  TString sYear;
  if(year==2016){
    lumi       = 35.9; //fb-1 2016
    sYear      = "2016";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2016/";
  }
  else if(year==2017){
    lumi       = 41.5; //fb-1 2017
    sYear      = "2017";
    inFilePath = "/eos/user/a/acappati/samples_HH4lbb/samples_2017/";
  }
  else if(year==2018){
    lumi       = 59.7; //fb-1 2018
    sYear      = "2018";
    inFilePath = "/eos/user/a/acappati/samples_4lX/20200205_bestKD_samples2018/";
    inDataPath = "/eos/user/a/acappati/samples_4lX/allsamples/";
  }
  else{ 
    cout<<"wrong year selected!"<<endl;
    break;
  }



  //---datasets
  int nDatasets = 22;
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
    "TTZJets_M10_MLMext1",
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


  // define 1D histos
  TH1F* h1_m4l_4lselOnly[nFinalStates+1][nProcesses]; 
  //  TH1F* h1_m4l_4ljjsel[nFinalStates+1][nProcesses]; 
  //  TH1F* h1_mbb_4ljjsel[nFinalStates+1][nProcesses];
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int pr=0; pr<nProcesses; pr++){
      h1_m4l_4lselOnly[fs][pr] = new TH1F("h1_m4l_4lselOnly_"+sProcess[pr]+"_"+sFinalState[fs]+"_"+sYear,";m_{4l} (GeV); Events / 2 GeV", 65, 70., 200.);
      h1_m4l_4lselOnly[fs][pr]->Sumw2(true);
      //  h1_m4l_4ljjsel[fs][pr]; 
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
       datasets[d]=="ggTo2mu2tau_Contin_MCFM70") currentProcess = ggZZ;
    if(datasets[d]=="TTZToLLNuNu_M10") currentProcess = TTZ; 
    if(datasets[d]=="TTWJetsToLNu") currentProcess = TTW;
    

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



    //save only events for 1 final state at the time
    //    if(currentFinalState != fs_4mu)   continue;  // save 4mu only
    //    if(currentFinalState != fs_4e)    continue;  // save 4e only
    //    if(currentFinalState != fs_2e2mu) continue;  // save 2e2mu only
    //    cout<<currentFinalState<<endl;


    // mass cut: signal region
     if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV
    //    if(ZZMass < 118 || ZZMass > 130) continue;

 


    //    if(JETSELECTION){

    //JETSELECTION---------------------------------------------------
      // jet quantities
      vector<TLorentzVector> JetVec; // TLorentz vector with all Jets per Event
      vector<TLorentzVector> JetPair; // TLorentz vector with all Jets Pairs
      vector<float> JetBinfo; // vector with b-tag info per each Jet of the Event
      
  
      for (UInt_t j = 0; j < JetPt->size(); j++)
        {
	  if ( fabs( JetEta->at(j) ) > 2.4 ) continue;
          if (JetPt->at(j) < 20 )  continue; // pt cut 20GeV from ntuplizer 
  	  
  	TLorentzVector temp;
  	temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
  	JetVec.push_back(temp);
  	JetBinfo.push_back(JetBTagger->at(j));
        }
  
  
  
      // at least 2 jets in the acceptance
      if (JetVec.size() < 2) continue;   
          
      
  
      // get and save vector with max btagger value
      f_bdiscjet1signal = *max_element(JetBinfo.begin(), JetBinfo.end());
  
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
      f_bdiscjet2signal = JetBinfo.at(d2_maxJetPt);

      // save jet pT
      f_ptjet1 = JetVec.at(d1_maxbtag).Pt();
      f_ptjet2 = JetVec.at(d2_maxJetPt).Pt();
  
  
      // build H->bb tlorentzvector
      TLorentzVector Hbb_Vec = JetVec.at(d1_maxbtag) + JetVec.at(d2_maxJetPt);
  
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

      //    } // end if JETSELECTION
    //JETSELECTION---------------------------------------------------





    // fill k factors and event weights
    Float_t kfactor = 1.;
    // qqZZ sample
    if(inFile.Contains("ZZTo4l")) { kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; }
    //ggZZ samples                       
    else if(inFile.Contains("ggTo")) { kfactor = KFactor_QCD_ggZZ_Nominal; }

    //    if (isHH) xsec = 0.00001017; // fb Angela

    Double_t eventWeight = 1.;
    if(!isDATA && !isZX) eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight;
    if(isZX) eventWeight = weight; //ZX weight

    //    cout<<" -- xsec "<<xsec<<endl;
    //    break;

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

  }//end loop over datasets
}


void prepareNtupleMVA_1bjet1Pt()
{

  float lumi = 59.74; //fb-1 2018
  //float lumi = 41.5; //fb-1 2017

  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/20200205_bestKD_samples2018/";
  //  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/20200209_samples2017/";
  //  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/allsamples/";
  TString inputFileName[] = {// "HH4lbb_Angela",
                             // "ggH125",
                             // "VBFH125",
                             // "WplusH125",
                             // "WminusH125",
                             // "ZH125",
                             // "bbH125",
                             // "ttH125",
                             // "ZZTo4lamcatnlo",
                             // "TTZJets_M10_MLMext1",
                             // "TTZToLL_M1to1O_MLM",
                             // "TTWJetsToLNu",
                             // "ggTo4e_Contin_MCFM701",
                             // "ggTo4mu_Contin_MCFM701",
                             // "ggTo4tau_Contin_MCFM701",
                             // "ggTo2e2mu_Contin_MCFM701",
                             // "ggTo2e2tau_Contin_MCFM701",
                             // "ggTo2mu2tau_Contin_MCFM701",
                             // "WWZ",
                             // "WZZ",
                             // "ZZZ",
                             "ZX",
			     //"AllData", 
                             };


  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<< "number of input files: " << nInputFiles<<endl;


  string outputFilePath = "out";
  //string outputFilePath = "200212_mvaNtuples_1bjet1Pt_4l";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir


  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doNtuplesForMVA(inputFilePath + inputFileName[i] + "/ZZXAnalysis.root", outputFilePath + "/reduced_" + inputFileName[i] + ".root", lumi);
  }


}
