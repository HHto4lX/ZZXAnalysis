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

#define MERGE2E2MU 1
#define DOBLINDHISTO 1
#define VERBOSE 0

double DELTAR = .4;

// Final State scan

enum FinalState {fs_4mu=0, fs_4e=1, fs_2e2mu=2, fs_2mu2e=3};  // 4mu, 4e, 2e2mu, 2mu2e
const int nFinalState = 4;
string FinalState[nFinalState+1] = {"4mu", "4e","2e2mu","2mu2e", "4L"};
string FinalStateLabel[nFinalState+1] = {"4#mu", "4e", "2e2#mu", "2#mu2e", "4L"};
int FinalStateColor[nFinalState+1] = {2, 4, 8, 3, 1};
int FinalStateMarkerStyle[nFinalState+1] = {20, 22, 21, 33, 29};

struct Histo1D 
{
  string Name; 
  string XLabel; 
  string YLabel; 
  string CutLabel;
  double Min;
  double Max;
  int nBin;
  bool isLogx;
  bool isLogy;
};

const int nHisto = 41;

Histo1D myHisto1D[nHisto] = {
  {"M4L_CR4Lonly", "m_{4l} (GeV)", "Events / 5 GeV", "", 70, 1070, 200, 1, 0},  
  {"MZ1_CR4Lonly", "m_{Z1} (GeV)", "Events / 2 GeV", "", 40, 120, 40, 0, 0},  
  {"MZ2_CR4Lonly", "m_{Z2} (GeV)", "Events / 2 GeV", "", 12, 120, 54, 0, 0}, 
  {"pt4L_CR4Lonly", "p_{T}^{4l}", "Events / 7 GeV", "", 0, 700, 100, 0, 0}, 
  {"eta4L_CR4Lonly", "#eta^{4l}", "Events", "", -4, 4, 40, 0, 0}, 
  {"phi4L_CR4Lonly", "#phi^{4l}", "Events", "", -4, 4, 40, 0, 0},
  {"NJets","# jets", "Events", "", 0, 20, 20, 0, 0}, 
  {"NJetsnoB","# jets not b-tagged", "Events", "", 0, 20, 20, 0, 0}, 
  {"NBJets","# b-jets", "Events", "", 0, 20, 20, 0, 0}, 
  {"M4L", "m_{4l} (GeV)", "Events / 10 GeV", "", 70, 1070, 100, 1, 0},  
  {"MZ1", "m_{Z1} (GeV)", "Events / 2 GeV", "", 40, 120, 40, 0, 0},  
  {"MZ2", "m_{Z2} (GeV)", "Events / 2 GeV", "", 12, 120, 54, 0, 0}, 
  {"pt4L", "p_{T}^{4l}", "Events / 28 GeV", "", 0, 700, 25, 0, 0}, 
  {"eta4L", "#eta^{4l}", "Events", "", -4, 4, 20, 0, 0}, 
  {"phi4L", "#phi^{4l}", "Events", "", -4, 4, 20, 0, 0}, 
  {"etabb_m", "#eta^{bb}", "Events", "", -4, 4, 40, 0, 0},
  {"phibb_m", "#phi^{bb}", "Events", "", -4, 4, 40, 0, 0},
  {"etabb_pt", "#eta^{bb}", "Events", "", -4, 4, 40, 0, 0},
  {"phibb_pt", "#phi^{bb}", "Events", "", -4, 4, 40, 0, 0},
  {"Deltaetabb_m", "#eta^{bb}", "Events", "", -8, 8, 64, 0, 0},
  {"Deltaphibb_m", "#phi^{bb}", "Events", "", -8, 8, 64, 0, 0},
  {"Deltaetabb_pt", "#eta^{bb}", "Events", "", -8, 8, 64, 0, 0},
  {"Deltaphibb_pt", "#phi^{bb}", "Events", "", -8, 8, 64, 0, 0},
  {"methodPTjet_M", "m_{pairJ}", "Events", "", 0, 500, 100, 0, 0},
  {"methodPTjet_PT", "p_{T}^{pairJ}", "Events", "", 0, 1000, 200, 0, 0},
  {"methodPTjet_Binfo", "jet b-tag", "Events", "", -1, 4, 5, 0, 0},
  {"methodM_M", "m_{pairJ}", "Events", "", 0, 500, 100, 0, 0},  
  {"methodM_PT", "p_{T}^{pairJ}", "Events", "", 0, 1000, 200, 0, 0},
  {"methodM_Binfo", "jet b-tag", "Events", "", -1, 4, 5, 0, 0},
  {"methodBM_M", "m_{pairJ}", "Events", "", 0, 500, 100, 0, 0},  
  {"methodBM_PT", "p_{T}^{pairJ}", "Events", "", 0, 1000, 200, 0, 0},
  {"methodBPT_M", "m^{jj}", "Events / 20 GeV", "", 0, 500, 25, 0, 0},  
  {"methodBPT_PT", "p_{T}^{jj}", "Events / 28 GeV", "", 0, 700, 25, 0, 0},
  {"GENRECOeff_methodPTjet", "GEN-RECO efficiency ptjet method", "Events", "", 0, 3, 3, 0, 0},
  {"GENRECOeff_methodM", "GEN-RECO efficiency madd method", "Events", "", 0, 3, 3, 0, 0},
  {"methodM_Pruned1", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodM_Pruned2", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodM_PrunedTOT", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodPTjet_Pruned1", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodPTjet_Pruned2", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodPTjet_PrunedTOT", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
};



void doHisto(TString inputFileMC, TString outputFile, double lumi=1)
{

  bool GEN = false;
  bool ISDATA = false;
  if ( inputFileMC.Contains("AllData")) ISDATA = true ;

  std::cout << "ISDATA: " << ISDATA << endl; 
  std::cout << "DOBLINDHISTO: " << DOBLINDHISTO << endl;

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

  Float_t KFactor_QCD_ggZZ_Nominal;
  Float_t KFactor_EW_qqZZ;
  Float_t KFactor_QCD_qqZZ_dPhi;
  Float_t KFactor_QCD_qqZZ_M;
  Float_t KFactor_QCD_qqZZ_Pt;

  Short_t ZZsel;
  vector<Float_t> *LepEta = 0;
  Float_t ZZMass;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Float_t ZZPt;
  Float_t ZZEta;
  Float_t ZZPhi;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t MPairJ = 0;
  Float_t PtPairJ = 0;
  
  Float_t Deltabb_eta_m = 0;
  Float_t Deltabb_phi_m = 0; 
  Float_t Deltabb_eta_pt = 0;
  Float_t Deltabb_phi_pt = 0;
 
  vector<Float_t> *JetIsBTagged = 0;
  vector<Float_t> *JetPt     = 0;
  vector<Float_t> *JetEta    = 0;
  vector<Float_t> *JetPhi    = 0;
  vector<Float_t> *JetMass   = 0;
  vector<Float_t> *GENjetParentID   = 0;
  vector<Float_t> *prunedGenPartID   = 0;
  vector<Float_t> *prunedGenPartPt   = 0;
  vector<Float_t> *prunedGenPartEta   = 0;
  vector<Float_t> *prunedGenPartPhi   = 0;
  vector<Float_t> *prunedGenPartMass   = 0;
  vector<Float_t> *prunedGenMotherID   = 0;

  Float_t BESTpair_methodPTjet_M = 0;
  Float_t BESTpair_methodPTjet_PT = 0;
  Float_t BESTpair_methodPTjet_ETA = 0;
  Float_t BESTpair_methodPTjet_PHI = 0;
  Int_t BESTpair_methodPTjet_Binfo = 0;
  Int_t BESTpair_methodPTjet_index1 = 0;
  Int_t BESTpair_methodPTjet_index2 = 0;

  Float_t BESTpair_methodM_M = 0;
  Float_t BESTpair_methodM_PT = 0;
  Float_t BESTpair_methodM_ETA = 0;
  Float_t BESTpair_methodM_PHI = 0;
  Int_t BESTpair_methodM_Binfo = 0;
  Int_t BESTpair_methodM_index1 = 0;
  Int_t BESTpair_methodM_index2 = 0;

  Float_t BESTpair_methodBM_M = 0;
  Float_t BESTpair_methodBM_PT = 0;
  Float_t BESTpair_methodBM_ETA = 0;
  Float_t BESTpair_methodBM_PHI = 0;
  Int_t BESTpair_methodBM_Binfo = 0;
  Int_t BESTpair_methodBM_index1 = 0;
  Int_t BESTpair_methodBM_index2 = 0;

  Float_t BESTpair_methodBPT_M = 0;
  Float_t BESTpair_methodBPT_PT = 0;
  Float_t BESTpair_methodBPT_ETA = 0;
  Float_t BESTpair_methodBPT_PHI = 0;
  Int_t BESTpair_methodBPT_Binfo = 0;
  Int_t BESTpair_methodBPT_index1 = 0;
  Int_t BESTpair_methodBPT_index2 = 0;

  Float_t BESTpair_methodM_M_0H = 0;
  Float_t BESTpair_methodM_M_1H_eff1 = 0;
  Float_t BESTpair_methodM_M_1H_eff0 = 0;
  Float_t BESTpair_methodM_M_2H_eff2 = 0;
  Float_t BESTpair_methodM_M_2H_eff1 = 0;
  Float_t BESTpair_methodM_M_2H_eff0 = 0;

  Float_t BESTpair_methodPTjet_M_0H = 0;
  Float_t BESTpair_methodPTjet_M_1H_eff1 = 0;
  Float_t BESTpair_methodPTjet_M_1H_eff0 = 0;
  Float_t BESTpair_methodPTjet_M_2H_eff2 = 0;
  Float_t BESTpair_methodPTjet_M_2H_eff1 = 0;
  Float_t BESTpair_methodPTjet_M_2H_eff0 = 0;

  Int_t methodM_eff_pruned1 = 0;
  Int_t methodM_eff_pruned2 = 0;
  Int_t methodM_eff_prunedTOT = 0;
  Int_t methodPTjet_eff_pruned1 = 0;
  Int_t methodPTjet_eff_pruned2 = 0;
  Int_t methodPTjet_eff_prunedTOT = 0;

  Float_t efficiency_mass = 0;
  Float_t efficiency_pt = 0;

  Int_t count_0H = 0;
  Float_t yield_0H = 0;
  Int_t count_1H = 0;
  Float_t yield_1H = 0;
  Int_t count_2H = 0;
  Float_t yield_2H = 0;

  double yield = 0;

  TH1F* h1[nHisto][nFinalState+1];
  
  for (int fs = 0; fs < (nFinalState+1); fs++ )
    {
      for (int h = 0; h < nHisto; h++)
	{
	  h1[h][fs] = new TH1F( 
			       Form("h1_%s_%s", myHisto1D[h].Name.c_str(), FinalState[fs].c_str()), 
			       Form(";%s;%s", myHisto1D[h].XLabel.c_str(), myHisto1D[h].YLabel.c_str()),
			       myHisto1D[h].nBin, myHisto1D[h].Min, myHisto1D[h].Max );

	  h1[h][fs]->Sumw2(true);
	}       
    }
  
  int currentFinalState;

  inputFile =  TFile::Open( inputFileMC );
  
  hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
  NGenEvt = (Float_t)hCounters->GetBinContent(1);
  gen_sumWeights = (Float_t)hCounters->GetBinContent(40);
  partialSampleWeight = lumi * 1000 / gen_sumWeights;
  inputTree = (TTree*)inputFile->Get("ZZTree/candTree");
  
  // set address of all branches
  if (!ISDATA) inputTree->SetBranchAddress("genBR", &genBR);
  inputTree->SetBranchAddress("RunNumber", &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber", &nLumi);
  if (!ISDATA) inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
  if (!ISDATA) inputTree->SetBranchAddress("xsec", &xsec);
  inputTree->SetBranchAddress("ZZsel", &ZZsel);
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
  inputTree->SetBranchAddress("JetIsBtagged",  &JetIsBTagged);
  if (GEN) inputTree->SetBranchAddress("GENjetParentID",  &GENjetParentID);
  if (GEN) inputTree->SetBranchAddress("prunedGenPartEta", &prunedGenPartEta );
  if (GEN) inputTree->SetBranchAddress("prunedGenPartPhi", &prunedGenPartPhi );
  if (GEN) inputTree->SetBranchAddress("prunedGenPartPt", &prunedGenPartPt );
  if (GEN) inputTree->SetBranchAddress("prunedGenPartMass", &prunedGenPartMass );
  if (GEN) inputTree->SetBranchAddress("prunedGenPartID", &prunedGenPartID );
  if (GEN) inputTree->SetBranchAddress("prunedGenMotherID", &prunedGenMotherID );
  if(inputFileMC.Contains("ggTo"))    //ggZZ samples                                                                                                                                                        
    {
      inputTree->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
    }
  if(inputFileMC.Contains("ZZTo4l"))   //qqZZ samples                                                                                                                                                       
    {
      inputTree->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
      inputTree->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi);
      inputTree->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
      inputTree->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt);
    }

  int entries = inputTree->GetEntries();
  std::cout<<"Processing file: "<< inputFileMC << "\nNumber of entries: " << entries << endl;
  
  //  entries = 10;

  for (Long64_t entry = 0; entry < entries; entry++)
    {  
      int nbjet = 0;
      int njet = 0;
      int njet_noB = 0;
      int counterH = 0;
      bool event_2H = false;
      bool event_1H = false;
      bool event_0H = false;
 
      inputTree->GetEntry(entry);
      
      if(LepEta->size() != 4)
	{
	  cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	  continue;
	}
      if( !(ZZsel >= 90) ) continue;

      Float_t kfactor = 1.;
      // qqZZ sample                                                                                                                                                                                        
      if(inputFileMC.Contains("ZZTo4l")) { kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; }
      //ggZZ samples                                                                                                                                                                                        
      else if(inputFileMC.Contains("ggTo")) { kfactor = KFactor_QCD_ggZZ_Nominal; }

      Double_t eventWeight = 0.;
      if (!ISDATA) eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight ;
      else eventWeight = 1.;

      if (VERBOSE) 
	{std::cout << "--------- event: " << entry << " weight: " << eventWeight << endl;
	  /*
	    std::cout << "-------------------\nlumi " << lumi << endl;
	    std::cout << "xsec " << xsec << endl;
	    std::cout << "genBR " << genBR << endl;
	    std::cout << "genSumwweight " << gen_sumWeights << endl;
	    std::cout << "overallEventwweight " << overallEventWeight << endl;
	    std::cout << "Bin 0 " << NGenEvt << endl;
	  */
	}

      currentFinalState = -1;
      if (Z1Flav == -121)
	{
	  if (Z2Flav == -121)
	    currentFinalState = fs_4e;
	  else if ( Z2Flav == -169)
	    currentFinalState = fs_2e2mu;
	  else
	    cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z2Flav = " << Z2Flav << endl;
	}
      else if (Z1Flav == -169)
	{
	  if (Z2Flav == -121)
	    currentFinalState = fs_2mu2e;
	  else if (Z2Flav == -169)
	    currentFinalState = fs_4mu;
	  else
	    cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z2Flav = " << Z2Flav << endl;
	}
      else
	{
	  cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; Z1Flav = " << Z1Flav << endl;
	}
      if(MERGE2E2MU && ( currentFinalState == fs_2mu2e )) currentFinalState = fs_2e2mu;

      // Mass Cut      
      //      if (ZZMass < 115 || ZZMass > 135) continue;

      //      std::cout << "ZZ Mass: " << ZZMass << endl;


      Float_t histo[nHisto];

      bool SHOW =  ( (ZZMass > 135) || (ZZMass < 115) ) || (!DOBLINDHISTO) ;

      for(int v = 0; v < nHisto; v++)
	{
	  string histoString = myHisto1D[v].Name.c_str();
	  if      (histoString == "M4L_CR4Lonly" && ( !ISDATA || SHOW ) ) histo[v] = ZZMass;
	  else if (histoString == "MZ1_CR4Lonly") histo[v] = Z1Mass;
	  else if (histoString == "MZ2_CR4Lonly") histo[v] = Z2Mass;
	  else if (histoString == "pt4L_CR4Lonly") histo[v] = ZZPt;
	  else if (histoString == "eta4L_CR4Lonly") histo[v] = ZZEta;
	  else if (histoString == "phi4L_CR4Lonly") histo[v] = ZZPhi;
	  else continue;
	  	
	  h1[v][currentFinalState]->Fill(histo[v], eventWeight);
	}



      // H(bb) selection
      
      vector<TLorentzVector> JetVec; // TLorentz vector with all Jets per Event
      vector<TLorentzVector> JetPair; // TLorentz vector with all Jets Pairs
      vector<int> JetBinfo; // vector with b-tag info per each Jet of the Event
      vector<int> JetPairBinfo; // vector with b-tag info per each Jet Pair of the Event
      vector<int> JetPair_index1; // vector with the index of each first jet of the pair 
      vector<int> JetPair_index2; // vector with the index of each second jet of the pair
      vector<int> index_GEN;
     
      int btag;      
       
      for (UInt_t j = 0; j < JetPt->size(); j++)
	{
	  btag = 0;
  	  if ( fabs ( JetEta->at(j) ) > 2.4 ) continue; // pt cut 20GeV from ntuplizer
	  
	  njet++;
	  if (JetIsBTagged->at(j) >0) 
	    {
	      btag = 1;
	      nbjet ++;
	    }
	  else njet_noB++;
	  if (VERBOSE)
	    {
	      std::cout << "Jet is: " << JetIsBTagged->at(j) << endl;
	      std::cout << "btag  : " << btag << endl;
	    }
	  TLorentzVector temp;
	  temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
	  JetVec.push_back(temp);
	  JetBinfo.push_back(btag);
	}

      if (VERBOSE)
	{
	  std::cout << "Jet Vec size: " << JetVec.size() << endl;
	  std::cout << "B info  size: " << JetBinfo.size() << endl;
	  std::cout << "nbjet       : " << nbjet << endl;
	  std::cout << "b info : ";
	  for (int i = 0; i < JetBinfo.size(); i++)
	    {
	      std::cout << JetBinfo[i] << "  ";
	    }      
	  std::cout << endl;
	}
      
      // remove events with <2 jet OR no btagged jets
      
      if (JetVec.size() < 2) continue;      
      if (nbjet < 1) continue;
      if (VERBOSE) std::cout << "passed: >= 2 jet, >= 1 b-jet" << endl;
      
      // build pairs of jets
      
      for (UInt_t i = 0; i < JetVec.size() -1; i++)
	{
	  TLorentzVector jet0 = JetVec[i];
	  for (UInt_t j = i+1; j < JetVec.size(); j++)
	    {
	      TLorentzVector jet1 = JetVec[j];
	      JetPair.push_back(jet0 + jet1);
	      int temp_int = JetBinfo[i] + JetBinfo[j];
	      JetPairBinfo.push_back(temp_int); 
	      JetPair_index1.push_back(i);
	      JetPair_index2.push_back(j);
	    }	  
	}
      
      if (VERBOSE)
	{
	  std::cout << "Pair Binfo  size: " << JetPairBinfo.size() << endl;
	  std::cout << "Pair index1 size: " << JetPair_index1.size() << endl;
	  std::cout << "Pair index2 size: " << JetPair_index2.size() << endl;
	  std::cout << "b-info : ";
	  for (int i = 0; i < JetPairBinfo.size(); i++)
	    {
	      std::cout << JetPairBinfo[i] << " ";
	    }
	  std::cout << "\nindex1 : " ;
	  for (int i = 0; i < JetPair_index1.size(); i++)
	    {
	      std::cout << JetPair_index1[i] << " ";
	    }
	  std::cout << "\nindex2 : " ;
	  for (int i = 0; i < JetPair_index1.size(); i++)
	    {
	      std::cout << JetPair_index2[i] << " ";
	    }
	  std::cout << endl;
	}
      // remove events with less than 2 pairs of jets

      if (JetPair.size() < 1) continue;      
      if (VERBOSE) std::cout << "passed: >= 1 jetpair" << endl;

      double deltaM = 100000.;
      double deltaM_B = 100000.;
      double sumPT = 0.;

      BESTpair_methodPTjet_M = 0;
      BESTpair_methodPTjet_PT = 0;
      BESTpair_methodPTjet_ETA = 0;
      BESTpair_methodPTjet_PHI = 0;
      BESTpair_methodPTjet_Binfo = -1;
      BESTpair_methodPTjet_index1 = -1;
      BESTpair_methodPTjet_index2 = -1;
      
      BESTpair_methodM_M = 0;
      BESTpair_methodM_PT = 0;
      BESTpair_methodM_ETA = 0;
      BESTpair_methodM_PHI = 0;
      BESTpair_methodM_Binfo = -1;
      BESTpair_methodM_index1 = -1;
      BESTpair_methodM_index2 = -1;

      BESTpair_methodBM_M = 0;
      BESTpair_methodBM_PT = 0;
      BESTpair_methodBM_ETA = 0;
      BESTpair_methodBM_PHI = 0;
      BESTpair_methodBM_Binfo = -1;
      BESTpair_methodBM_index1 = -1;
      BESTpair_methodBM_index2 = -1;

      BESTpair_methodBPT_M = 0;
      BESTpair_methodBPT_PT = 0;
      BESTpair_methodBPT_ETA = 0;
      BESTpair_methodBPT_PHI = 0;
      BESTpair_methodBPT_Binfo = -1;
      BESTpair_methodBPT_index1 = -1;
      BESTpair_methodBPT_index2 = -1;

      methodM_eff_pruned1 = 0;
      methodM_eff_pruned2 = 0;
      methodM_eff_prunedTOT = 0;
      methodPTjet_eff_pruned1 = 0;
      methodPTjet_eff_pruned2 = 0;
      methodPTjet_eff_prunedTOT = 0;

      Deltabb_eta_m = 0;
      Deltabb_phi_m = 0;
      Deltabb_eta_pt = 0;
      Deltabb_phi_pt = 0;

      int size_binfo_2 = 0;
      int size_binfo_1 = 0;

      // loop on pairs to find best pairs

      for (UInt_t s = 0; s<JetPair.size(); s++)
	{
	  if (JetPairBinfo[s] > 1) size_binfo_2++;
	  if (JetPairBinfo[s] == 1) size_binfo_1++;
	}    

      if (VERBOSE)
	{
	  std::cout << "# 2 b-tag : " << size_binfo_2 << endl;
	  std::cout << "# 1 b-tag : " << size_binfo_1 << endl;
	}
      for (UInt_t p = 0; p<JetPair.size(); p++)
	{

	  if (VERBOSE) std::cout << "loop i : " << p << endl;

	  if ( JetPairBinfo[p] > 0 )
	    {
	      if (BESTpair_methodPTjet_M == 0)
		{
		  if (VERBOSE) std::cout << "best PTjet" << endl;
		  BESTpair_methodPTjet_M = JetPair[p].M();
		  BESTpair_methodPTjet_ETA = JetPair[p].Eta();
		  BESTpair_methodPTjet_PHI = JetPair[p].Phi();
		  BESTpair_methodPTjet_PT = JetPair[p].Pt();
		  BESTpair_methodPTjet_Binfo = JetPairBinfo[p];
		  
		  BESTpair_methodPTjet_index1 = JetPair_index1[p];
		  BESTpair_methodPTjet_index2 = JetPair_index2[p];
		}  
	    }
	  
	  if (size_binfo_2 == 1 && JetPairBinfo[p] == 2) 
	    {

	      if (VERBOSE) std::cout << "#2b OK + btag2" << endl;

	      BESTpair_methodBM_M = JetPair[p].M();
	      BESTpair_methodBM_ETA = JetPair[p].Eta();
	      BESTpair_methodBM_PHI = JetPair[p].Phi();
	      BESTpair_methodBM_PT = JetPair[p].Pt();
	      BESTpair_methodBM_Binfo = JetPairBinfo[p];
	      BESTpair_methodBM_index1 = JetPair_index1[p];
	      BESTpair_methodBM_index2 = JetPair_index2[p];

	      BESTpair_methodBPT_M = JetPair[p].M();
	      BESTpair_methodBPT_ETA = JetPair[p].Eta();
	      BESTpair_methodBPT_PHI = JetPair[p].Phi();
	      BESTpair_methodBPT_PT = JetPair[p].Pt();
	      BESTpair_methodBPT_Binfo = JetPairBinfo[p];
	      BESTpair_methodBPT_index1 = JetPair_index1[p];
	      BESTpair_methodBPT_index2 = JetPair_index2[p];
	    }
	  
	  if (size_binfo_2 > 1 && JetPairBinfo[p] == 2 && (fabs (JetPair[p].M() - 125.) < deltaM_B) )
	    {

	      if (VERBOSE) std::cout << "MORE THAN 2 B-TAg JET ************************\n#2b PLUS  + 2Btag + mass" << endl;

	      deltaM_B = fabs(JetPair[p].M() - 125.);

              BESTpair_methodBM_M = JetPair[p].M();
              BESTpair_methodBM_ETA = JetPair[p].Eta();
              BESTpair_methodBM_PHI = JetPair[p].Phi();
              BESTpair_methodBM_PT = JetPair[p].Pt();
              BESTpair_methodBM_Binfo = JetPairBinfo[p];

              BESTpair_methodBM_index1 = JetPair_index1[p];
              BESTpair_methodBM_index2 = JetPair_index2[p];
	    }

	  if (size_binfo_2 > 1 && JetPairBinfo[p] == 2 && JetPair[p].Pt() > sumPT )
	    {

	      if (VERBOSE) std::cout << "MORE THAN 2 B-TAg JET ************************\n#2b PLUS  + 2Btag + pt" << endl;

	      sumPT = JetPair[p].Pt();

              BESTpair_methodBPT_M = JetPair[p].M();
              BESTpair_methodBPT_ETA = JetPair[p].Eta();
              BESTpair_methodBPT_PHI = JetPair[p].Phi();
              BESTpair_methodBPT_PT = JetPair[p].Pt();
              BESTpair_methodBPT_Binfo = JetPairBinfo[p];

              BESTpair_methodBPT_index1 = JetPair_index1[p];
              BESTpair_methodBPT_index2 = JetPair_index2[p];
	    }


	  if (size_binfo_2 == 0 && JetPairBinfo[p] > 0 && (fabs (JetPair[p].M() - 125.) < deltaM_B) )
            {

	      if (VERBOSE) std::cout << "#1b OK + btag1 + mass" << endl;
	      
              deltaM_B = fabs(JetPair[p].M() - 125.);

              BESTpair_methodBM_M = JetPair[p].M();
              BESTpair_methodBM_ETA = JetPair[p].Eta();
              BESTpair_methodBM_PHI = JetPair[p].Phi();
              BESTpair_methodBM_PT = JetPair[p].Pt();
              BESTpair_methodBM_Binfo = JetPairBinfo[p];
              BESTpair_methodBM_index1 = JetPair_index1[p];
              BESTpair_methodBM_index2 = JetPair_index2[p];
            }


	  if (size_binfo_2 == 0 && JetPairBinfo[p] > 0 && JetPair[p].Pt() > sumPT )
            {

	      if (VERBOSE) std::cout << "#1b OK + btag1 + pt" << endl;
	      
              sumPT = JetPair[p].Pt();

              BESTpair_methodBPT_M = JetPair[p].M();
              BESTpair_methodBPT_ETA = JetPair[p].Eta();
              BESTpair_methodBPT_PHI = JetPair[p].Phi();
              BESTpair_methodBPT_PT = JetPair[p].Pt();
              BESTpair_methodBPT_Binfo = JetPairBinfo[p];
              BESTpair_methodBPT_index1 = JetPair_index1[p];
              BESTpair_methodBPT_index2 = JetPair_index2[p];
            }
	  
	  
	  if ( (fabs (JetPair[p].M() - 125.) < deltaM) && JetPairBinfo[p] > 0 ) 
	    {

	      if (VERBOSE) std::cout << "best mass" << endl;

	      deltaM = fabs(JetPair[p].M() - 125.);

	      BESTpair_methodM_M = JetPair[p].M();
	      BESTpair_methodM_ETA = JetPair[p].Eta();
	      BESTpair_methodM_PHI = JetPair[p].Phi();
	      BESTpair_methodM_PT = JetPair[p].Pt();
	      BESTpair_methodM_Binfo = JetPairBinfo[p];
	      BESTpair_methodM_index1 = JetPair_index1[p];
	      BESTpair_methodM_index2 = JetPair_index2[p];
	    }

	}

      if (GEN)
	{
	  for (UInt_t pr = 0; pr < prunedGenPartEta->size(); pr ++)
	    {
	      if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25)
		{
		  double M_deltaPhi1 = JetVec[BESTpair_methodBM_index1].Phi() - prunedGenPartPhi->at(pr);
		  if (fabs (M_deltaPhi1) > 3.14) M_deltaPhi1 = 2 * 3.14 - fabs(M_deltaPhi1);
		  double M_deltaEta1 = JetVec[BESTpair_methodBM_index1].Eta() - prunedGenPartEta->at(pr);
		  
		  double M_DeltaR1 = sqrt(M_deltaEta1 * M_deltaEta1 + M_deltaPhi1 * M_deltaPhi1);
		  
		  if (M_DeltaR1 < DELTAR) methodM_eff_pruned1++;
		  
		  double M_deltaPhi2 = JetVec[BESTpair_methodBM_index2].Phi() - prunedGenPartPhi->at(pr);
		  if (fabs (M_deltaPhi2) > 3.14) M_deltaPhi2 = 2 * 3.14 - fabs(M_deltaPhi2);
		  double M_deltaEta2 = JetVec[BESTpair_methodBM_index2].Eta() - prunedGenPartEta->at(pr);
		  
		  double M_DeltaR2 = sqrt(M_deltaEta2 * M_deltaEta2 + M_deltaPhi2 * M_deltaPhi2);
		  
		  if (M_DeltaR2 < DELTAR) methodM_eff_pruned2++;
		  
		  
		  double PTjet_deltaPhi1 = JetVec[BESTpair_methodBPT_index1].Phi() - prunedGenPartPhi->at(pr);
		  if (fabs (PTjet_deltaPhi1) > 3.14) PTjet_deltaPhi1 = 2 * 3.14 - fabs(PTjet_deltaPhi1);
		  double PTjet_deltaEta1 = JetVec[BESTpair_methodBPT_index1].Eta() - prunedGenPartEta->at(pr);
		  
		  double PTjet_DeltaR1 = sqrt(PTjet_deltaEta1 * PTjet_deltaEta1 + PTjet_deltaPhi1 * PTjet_deltaPhi1);
		  
		  if (PTjet_DeltaR1 < DELTAR) methodPTjet_eff_pruned1++;
		  
		  double PTjet_deltaPhi2 = JetVec[BESTpair_methodBPT_index2].Phi() - prunedGenPartPhi->at(pr);
		  if (fabs (PTjet_deltaPhi2) > 3.14) PTjet_deltaPhi2 = 2 * 3.14 - fabs(PTjet_deltaPhi2);
		  double PTjet_deltaEta2 = JetVec[BESTpair_methodBPT_index2].Eta() - prunedGenPartEta->at(pr);
		  
		  double PTjet_DeltaR2 = sqrt(PTjet_deltaEta2 * PTjet_deltaEta2 + PTjet_deltaPhi2 * PTjet_deltaPhi2);
		  
		  if (PTjet_DeltaR2 < DELTAR) methodPTjet_eff_pruned2++;
		  
		}
	    }
	  
	  methodM_eff_prunedTOT = methodM_eff_pruned1 + methodM_eff_pruned2;
	  methodPTjet_eff_prunedTOT = methodPTjet_eff_pruned1 + methodPTjet_eff_pruned2; 
	  
	}
      
      // check efficiency for methodPT
     
      Deltabb_eta_m = ZZEta + BESTpair_methodM_ETA;    
      Deltabb_phi_m = ZZPhi - BESTpair_methodM_PHI;
      Deltabb_eta_pt = ZZEta + BESTpair_methodPTjet_ETA;
      Deltabb_phi_pt = ZZPhi - BESTpair_methodPTjet_PHI;
    


      for(int v = 0; v < nHisto; v++)
	{
	  string histoString = myHisto1D[v].Name.c_str();
	  if      (histoString == "M4L" && ( !ISDATA || SHOW ) ) histo[v] = ZZMass;
	  else if (histoString == "MZ1") histo[v] = Z1Mass;
	  else if (histoString == "MZ2") histo[v] = Z2Mass;
	  else if (histoString == "pt4L") histo[v] = ZZPt;
	  else if (histoString == "eta4L") histo[v] = ZZEta;
	  else if (histoString == "phi4L") histo[v] = ZZPhi;
	  else if (histoString == "etabb_m") histo[v] = BESTpair_methodM_ETA;
	  else if (histoString == "phibb_m") histo[v] = BESTpair_methodM_PHI;
	  else if (histoString == "etabb_pt") histo[v] = BESTpair_methodPTjet_ETA;
	  else if (histoString == "phibb_pt") histo[v] = BESTpair_methodPTjet_PHI;
	  else if (histoString == "Deltaetabb_m") histo[v] = Deltabb_eta_m;
	  else if (histoString == "Deltaphibb_m") histo[v] = Deltabb_phi_m;
 	  else if (histoString == "Deltaetabb_pt") histo[v] = Deltabb_eta_pt;
	  else if (histoString == "Deltaphibb_pt") histo[v] = Deltabb_phi_pt;
	  else if (histoString == "NJets") histo[v] = njet;
	  else if (histoString == "NJets_noB") histo[v] = njet_noB;
	  else if (histoString == "NBJets") histo[v] = nbjet;
	  else if (histoString == "methodPTjet_M") histo[v] = BESTpair_methodPTjet_M; 
	  else if (histoString == "methodPTjet_PT") histo[v] = BESTpair_methodPTjet_PT;
	  else if (histoString == "methodPTjet_ETA") histo[v] = BESTpair_methodPTjet_ETA;
	  else if (histoString == "methodPTjet_PHI") histo[v] = BESTpair_methodPTjet_PHI;
	  else if (histoString == "methodPTjet_Binfo") histo[v] = BESTpair_methodPTjet_Binfo;
 	  else if (histoString == "methodM_M") histo[v] = BESTpair_methodM_M;
	  else if (histoString == "methodM_PT") histo[v] = BESTpair_methodM_PT;
 	  else if (histoString == "methodM_ETA") histo[v] = BESTpair_methodM_ETA;
 	  else if (histoString == "methodM_PHI") histo[v] = BESTpair_methodM_PHI;
	  else if (histoString == "methodM_Binfo") histo[v] = BESTpair_methodM_Binfo;
	  else if (histoString == "methodM_Pruned1" && GEN ) histo[v] = methodM_eff_pruned1;
	  else if (histoString == "methodM_Pruned2" && GEN ) histo[v] = methodM_eff_pruned2;
	  else if (histoString == "methodM_PrunedTOT" && GEN ) histo[v] = methodM_eff_prunedTOT;
	  else if (histoString == "methodPTjet_Pruned1" && GEN ) histo[v] = methodPTjet_eff_pruned1;
	  else if (histoString == "methodPTjet_Pruned2" && GEN ) histo[v] = methodPTjet_eff_pruned2;
	  else if (histoString == "methodPTjet_PrunedTOT" && GEN ) histo[v] = methodPTjet_eff_prunedTOT;
	  else if (histoString == "methodBM_M") histo[v] = BESTpair_methodBM_M;
          else if (histoString == "methodBM_PT") histo[v] = BESTpair_methodBM_PT;
          else if (histoString == "methodBM_ETA") histo[v] = BESTpair_methodBM_ETA;
          else if (histoString == "methodBM_PHI") histo[v] = BESTpair_methodBM_PHI;
          else if (histoString == "methodBM_Binfo") histo[v] = BESTpair_methodBM_Binfo;
	  else if (histoString == "methodBPT_M" && SHOW ) histo[v] = BESTpair_methodBPT_M;
          else if (histoString == "methodBPT_PT") histo[v] = BESTpair_methodBPT_PT;
          else if (histoString == "methodBPT_ETA") histo[v] = BESTpair_methodBPT_ETA;
          else if (histoString == "methodBPT_PHI") histo[v] = BESTpair_methodBPT_PHI;
          else if (histoString == "methodBPT_Binfo") histo[v] = BESTpair_methodBPT_Binfo;
	  else continue;
	  	
	  h1[v][currentFinalState]->Fill(histo[v], eventWeight);
	}
      
      JetVec.clear();
      JetPair.clear();
     
    } //end loop entries
  
  // fill inclusive histo (4L)
  
  for(int n = 0; n < nHisto; n++)
    {
      for(int f = 0; f < nFinalState; f++)
	{
	  h1[n][nFinalState] -> Add ( h1[n][f] );
	}
    } 
  
  
  std::cout << "Output file "<< outputFile << endl;
  TFile* fOut = new TFile(outputFile,"recreate");  
  fOut->cd();
  
  for (int nfs = 0; nfs < nFinalState+1; nfs ++)
    {
      for (int v = 0; v < nHisto; v ++)
	{
	  h1[v][nfs]->Write();
	}
    }
  
  std::cout << "YIELD " << h1[11][4]->Integral(0,500) << endl; 
    
} //end doHisto


void ZZbb_analysis() 
{

  double lumi = 140; // full Run2 Lumi
  if (DOBLINDHISTO) lumi = 59.74;
  std::cout << "Lumi: " << lumi << endl;

  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/190829/";
  //TString inputFilePath = "/eos/user/a/acappati/samples_4lX/190924/";
  //TString inputFilePath = "/eos/user/a/acappati/samples_4lX/190626/";
  TString inputFileName[] = {"HH4lbb",
			     "ggH125",
			     "VBFH125",
			     "WplusH125",
			     "WminusH125",
			     "ZH125",
			     "bbH125",
			     "ttH125",
                             //"ggTo4e_Contin_MCFM701",
			     //"ggTo4mu_Contin_MCFM701",
			     //"ggTo4tau_Contin_MCFM701",
			     //"ggTo2e2mu_Contin_MCFM701",
			     //"ggTo2e2tau_Contin_MCFM701",
			     //"ggTo2mu2tau_Contin_MCFM701",
			     "ZZTo4lext1",
			     "TTZJets_M10_MLMext1",
			     "TTZToLL_M1to1O_MLM",
			     "TTWJetsToLNu",
			     //"DY2JetsToLL_M50",
			     "AllData"
                            };

  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<< "number of input files: " << nInputFiles<<endl;

  string outputFilePath = "histos_4lbb";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir
  

  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doHisto(inputFilePath + inputFileName[i] + "/ZZXAnalysis.root", outputFilePath + "/histos_" + inputFileName[i] + ".root", lumi);
  }
  
}
