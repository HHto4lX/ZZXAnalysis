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

double DELTAR = 1.2;

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

const int nHisto = 43;

Histo1D myHisto1D[nHisto] = {
  {"NJets","# jets", "Events", "", 0, 20, 20, 0, 0},
  {"NJetsnoB","# jets not b-tagged", "Events", "", 0, 20, 20, 0, 0},
  {"NBJets","# b-jets", "Events", "", 0, 20, 20, 0, 0},
  {"M4L", "m_{4L}", "Events", "", 0, 500, 100, 0, 0},
  {"MZ1", "m_{Z1}", "Events", "", 0, 200, 100, 0, 0},
  {"MZ2", "m_{Z2}", "Events", "", 0, 200, 100, 0, 0},
  {"pt4L", "p_{T}^{4L}", "Events", "", 0, 1000, 200, 0, 0},
  {"eta4L", "#eta^{4L}", "Events", "", -4, 4, 40, 0, 0},
  {"phi4L", "#phi^{4L}", "Events", "", -4, 4, 40, 0, 0},
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
  {"GENRECOeff_methodPTjet", "GEN-RECO efficiency ptjet method", "Events", "", 0, 3, 3, 0, 0},
  {"GENRECOeff_methodM", "GEN-RECO efficiency madd method", "Events", "", 0, 3, 3, 0, 0},
  {"methodM_Pruned1", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodM_Pruned2", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodM_PrunedTOT", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodPTjet_Pruned1", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodPTjet_Pruned2", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodPTjet_PrunedTOT", "Pruned gen efficiency", "Events", "", 0, 5, 5, 0, 0},
  {"methodM_M_0H","methodM_M_0H","Events","", 0, 500, 100, 0, 0},
  {"methodM_M_1H_eff1","methodM_M_1H_eff1","Events","", 0, 500, 100, 0, 0},
  {"methodM_M_1H_eff0","methodM_M_1H_eff0","Events","", 0, 500, 100, 0, 0},
  {"methodM_M_2H_eff2","methodM_M_2H_eff2","Events","", 0, 500, 100, 0, 0},
  {"methodM_M_2H_eff1","methodM_M_2H_eff1","Events","", 0, 500, 100, 0, 0},
  {"methodM_M_2H_eff0","methodM_M_2H_eff0","Events","", 0, 500, 100, 0, 0},
  {"methodPTjet_M_0H","methodPTjet_M_0H","Events","", 0, 500, 100, 0, 0},
  {"methodPTjet_M_1H_eff1","methodPTjet_M_1H_eff1","Events","", 0, 500, 100, 0, 0},
  {"methodPTjet_M_1H_eff0","methodPTjet_M_1H_eff0","Events","", 0, 500, 100, 0, 0},
  {"methodPTjet_M_2H_eff2","methodPTjet_M_2H_eff2","Events","", 0, 500, 100, 0, 0},
  {"methodPTjet_M_2H_eff1","methodPTjet_M_2H_eff1","Events","", 0, 500, 100, 0, 0},
  {"methodPTjet_M_2H_eff0","methodPTjet_M_2H_eff0","Events","", 0, 500, 100, 0, 0}
};
  
void doHisto(const std::string inputFileMC, const std::string outputFile, double lumi=1)
{
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

  inputFile =  TFile::Open( inputFileMC.c_str() );
  
  hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
  NGenEvt = (Float_t)hCounters->GetBinContent(1);
  gen_sumWeights = (Float_t)hCounters->GetBinContent(40);
  partialSampleWeight = lumi * 1000 / gen_sumWeights;
  inputTree = (TTree*)inputFile->Get("ZZTree/candTree");
  
  // set address of all branches
  inputTree->SetBranchAddress("genBR", &genBR);
  inputTree->SetBranchAddress("RunNumber", &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber", &nLumi);
  inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
  inputTree->SetBranchAddress("xsec", &xsec);
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
  inputTree->SetBranchAddress("GENjetParentID",  &GENjetParentID);
  inputTree->SetBranchAddress("prunedGenPartEta", &prunedGenPartEta );
  inputTree->SetBranchAddress("prunedGenPartPhi", &prunedGenPartPhi );
  inputTree->SetBranchAddress("prunedGenPartPt", &prunedGenPartPt );
  inputTree->SetBranchAddress("prunedGenPartMass", &prunedGenPartMass );
  inputTree->SetBranchAddress("prunedGenPartID", &prunedGenPartID );
  inputTree->SetBranchAddress("prunedGenMotherID", &prunedGenMotherID );

  int entries = inputTree->GetEntries();
  std::cout<<"Processing file: "<< inputFileMC.c_str() << "\nNumber of entries: " << entries << endl;
  
  //  entries = 100;

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
      
      Double_t eventWeight = partialSampleWeight * xsec * genBR * kfactor * overallEventWeight ;
    
      //Double_t eventWeight = 1;
      std::cout << "--------- event: " << entry << " weight: " << eventWeight << endl;

      /*
      std::cout << "-------------------\nlumi " << lumi << endl;
      std::cout << "xsec " << xsec << endl;
      std::cout << "genBR " << genBR << endl;
      std::cout << "genSumwweight " << gen_sumWeights << endl;
      std::cout << "overallEventwweight " << overallEventWeight << endl;
      std::cout << "Bin 0 " << NGenEvt << endl;
      */

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
      

      // H(bb) selection
      
      vector<TLorentzVector> JetVec;
      vector<TLorentzVector> JetPair;
      vector<int> JetBinfo;
      vector<int> JetPairBinfo;
      vector<int> JetPair_index1;
      vector<int> JetPair_index2;
      vector<int> index_GEN;
      
      int btag = 0;      
      
      efficiency_mass = 0;
      efficiency_pt = 0;

      for(UInt_t i=0; i<GENjetParentID->size(); i++)
	{
	  if (TMath::Abs(JetEta->at(i) > 2.4)) continue; 
	
	  if ( GENjetParentID->at(i)==25 ) 
	    {
	      index_GEN.push_back(i);
	      counterH++;
	    }
	}


      if (counterH == 2) event_2H = true;
      else if (counterH == 1) event_1H = true;
      else if (counterH == 0) event_0H = true;
 
      for (int j = 0; j < JetPt->size(); j++)
	{
	  if ( fabs ( JetEta->at(j) ) > 2.4 ) continue;
	  
	  njet++;
	  if (JetIsBTagged->at(j) >0) 
	    {
	      btag = 1;
	      nbjet ++;
	    }
	  else njet_noB++;
	  
	  TLorentzVector temp;
	  temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
	  JetVec.push_back(temp);
	  JetBinfo.push_back(btag);
	}
      
      // remove events with <2 jet OR no btagged jets
      
      if (JetVec.size() < 2) continue;      
      if (nbjet < 1) continue;

      // build pairs of jets
      
      for (int i = 0; i < JetVec.size() -1; i++)
	{
	  TLorentzVector jet0 = JetVec[i];
	  for (int j = i+1; j < JetVec.size(); j++)
	    {
	      TLorentzVector jet1 = JetVec[j];
	      
	      JetPair.push_back(jet0 + jet1);
	      JetPairBinfo.push_back(JetBinfo[i] + JetBinfo[j]); 
	      JetPair_index1.push_back(i);
	      JetPair_index2.push_back(j);
	    }	  
	}
      
      // remove events with less than 2 pairs of jets

      if (JetPair.size() < 1) continue;      

      double deltaM = 100000.;
      
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

      BESTpair_methodM_M_0H = -1.;
      BESTpair_methodM_M_1H_eff1 = -1.;
      BESTpair_methodM_M_1H_eff0 = -1.;
      BESTpair_methodM_M_2H_eff2 = -1.;
      BESTpair_methodM_M_2H_eff1 = -1.;
      BESTpair_methodM_M_2H_eff0 = -1.;

      BESTpair_methodPTjet_M_0H = -1.;
      BESTpair_methodPTjet_M_1H_eff1 = -1.;
      BESTpair_methodPTjet_M_1H_eff0 = -1.;
      BESTpair_methodPTjet_M_2H_eff2 = -1.;
      BESTpair_methodPTjet_M_2H_eff1 = -1.;
      BESTpair_methodPTjet_M_2H_eff0 = -1.;

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

      // loop on pairs to find best pairs
    

      for (int p = 0; p<JetPair.size(); p++)
	{
	  if ( JetPairBinfo[p] > 0 )
	    {
	      if (BESTpair_methodPTjet_M == 0)
		{
		  BESTpair_methodPTjet_M = JetPair[p].M();
		  BESTpair_methodPTjet_ETA = JetPair[p].Eta();
		  BESTpair_methodPTjet_PHI = JetPair[p].Phi();
		  BESTpair_methodPTjet_PT = JetPair[p].Pt();
		  BESTpair_methodPTjet_Binfo = JetPairBinfo[p];
		  
		  BESTpair_methodPTjet_index1 = JetPair_index1[p];
		  BESTpair_methodPTjet_index2 = JetPair_index2[p];
		}  
	    }
	  
	  if ( (fabs (JetPair[p].M() - 125.) < deltaM) && JetPairBinfo[p] > 0 ) 
	    {
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

      
      for (int pr = 0; pr < prunedGenPartEta->size(); pr ++)
	{
	  if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25)
	    {
	      double M_deltaPhi1 = JetVec[BESTpair_methodM_index1].Phi() - prunedGenPartPhi->at(pr);
	      if (fabs (M_deltaPhi1) > 3.14) M_deltaPhi1 = 2 * 3.14 - fabs(M_deltaPhi1);
	      double M_deltaEta1 = JetVec[BESTpair_methodM_index1].Eta() - prunedGenPartEta->at(pr);
	      
	      double M_DeltaR1 = sqrt(M_deltaEta1 * M_deltaEta1 + M_deltaPhi1 * M_deltaPhi1);

	      if (M_DeltaR1 < DELTAR) methodM_eff_pruned1++;

	      double M_deltaPhi2 = JetVec[BESTpair_methodM_index2].Phi() - prunedGenPartPhi->at(pr);
	      if (fabs (M_deltaPhi2) > 3.14) M_deltaPhi2 = 2 * 3.14 - fabs(M_deltaPhi2);
	      double M_deltaEta2 = JetVec[BESTpair_methodM_index2].Eta() - prunedGenPartEta->at(pr);
	      
	      double M_DeltaR2 = sqrt(M_deltaEta2 * M_deltaEta2 + M_deltaPhi2 * M_deltaPhi2);

	      if (M_DeltaR2 < DELTAR) methodM_eff_pruned2++;


	      double PTjet_deltaPhi1 = JetVec[BESTpair_methodPTjet_index1].Phi() - prunedGenPartPhi->at(pr);
	      if (fabs (PTjet_deltaPhi1) > 3.14) PTjet_deltaPhi1 = 2 * 3.14 - fabs(PTjet_deltaPhi1);
	      double PTjet_deltaEta1 = JetVec[BESTpair_methodPTjet_index1].Eta() - prunedGenPartEta->at(pr);
	      
	      double PTjet_DeltaR1 = sqrt(PTjet_deltaEta1 * PTjet_deltaEta1 + PTjet_deltaPhi1 * PTjet_deltaPhi1);

	      if (PTjet_DeltaR1 < DELTAR) methodPTjet_eff_pruned1++;

	      double PTjet_deltaPhi2 = JetVec[BESTpair_methodPTjet_index2].Phi() - prunedGenPartPhi->at(pr);
	      if (fabs (PTjet_deltaPhi2) > 3.14) PTjet_deltaPhi2 = 2 * 3.14 - fabs(PTjet_deltaPhi2);
	      double PTjet_deltaEta2 = JetVec[BESTpair_methodPTjet_index2].Eta() - prunedGenPartEta->at(pr);
	      
	      double PTjet_DeltaR2 = sqrt(PTjet_deltaEta2 * PTjet_deltaEta2 + PTjet_deltaPhi2 * PTjet_deltaPhi2);

	      if (PTjet_DeltaR2 < DELTAR) methodPTjet_eff_pruned2++;

	    }
	}

      methodM_eff_prunedTOT = methodM_eff_pruned1 + methodM_eff_pruned2;
      methodPTjet_eff_prunedTOT = methodPTjet_eff_pruned1 + methodPTjet_eff_pruned2; 

      // check efficiency for methodPT

      if (event_0H == true)
	{
	  BESTpair_methodM_M_0H = BESTpair_methodM_M;
	  BESTpair_methodPTjet_M_0H = BESTpair_methodPTjet_M;
	  yield_0H += eventWeight;
	  count_0H++;
	}
      if (event_1H == true)
	{
	  yield_1H += eventWeight;
	  count_1H++;
	}
      if (event_2H == true)
	{
	  yield_2H += eventWeight;
	  count_2H++;
	}

      if ((event_1H == true) || (event_2H == true))
	{
	  for (int i = 0; i < index_GEN.size(); i++)
	    {
	      if ((BESTpair_methodPTjet_index1 == index_GEN[i]) || (BESTpair_methodPTjet_index2 == index_GEN[i])) efficiency_pt++;
	      if ((BESTpair_methodM_index1 == index_GEN[i]) || (BESTpair_methodM_index2 == index_GEN[i])) efficiency_mass++;	  
	    }
	}

      if ((event_1H == true) && (efficiency_mass == 1))
	{
	  BESTpair_methodM_M_1H_eff1 = BESTpair_methodM_M;
	}
      else if ((event_1H == true) && (efficiency_mass == 0))
	{
	  BESTpair_methodM_M_1H_eff0 = BESTpair_methodM_M;
	}
      else if ((event_2H == true) && (efficiency_mass == 2))
	{
	  BESTpair_methodM_M_2H_eff2 = BESTpair_methodM_M;
	}
      else if ((event_2H == true) && (efficiency_mass == 1))
	{
	  BESTpair_methodM_M_2H_eff1 = BESTpair_methodM_M;
	}
      else if ((event_2H == true) && (efficiency_mass == 0))
	{
	  BESTpair_methodM_M_2H_eff0 = BESTpair_methodM_M;
	}



      if ((event_1H == true) && (efficiency_pt == 1))
	{
	  BESTpair_methodPTjet_M_1H_eff1 = BESTpair_methodPTjet_M;
	}
      else if ((event_1H == true) && (efficiency_pt == 0))
	{
	  BESTpair_methodPTjet_M_1H_eff0 = BESTpair_methodPTjet_M;
	}
      else if ((event_2H == true) && (efficiency_pt == 2))
	{
	  BESTpair_methodPTjet_M_2H_eff2 = BESTpair_methodPTjet_M;
	}
      else if ((event_2H == true) && (efficiency_pt == 1))
	{
	  BESTpair_methodPTjet_M_2H_eff1 = BESTpair_methodPTjet_M;
	}
      else if ((event_2H == true) && (efficiency_pt == 0))
	{
	  BESTpair_methodPTjet_M_2H_eff0 = BESTpair_methodPTjet_M;
	}
     
      Deltabb_eta_m = ZZEta + BESTpair_methodM_ETA;    
      Deltabb_phi_m = ZZPhi - BESTpair_methodM_PHI;
      Deltabb_eta_pt = ZZEta + BESTpair_methodPTjet_ETA;
      Deltabb_phi_pt = ZZPhi - BESTpair_methodPTjet_PHI;

      Float_t histo[nHisto];
      
      for(int v = 0; v < nHisto; v++)
	{
	  string histoString = myHisto1D[v].Name.c_str();
	  if      (histoString == "M4L") histo[v] = ZZMass;
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
	  else if (histoString == "GENRECOeff_methodM") histo[v] = efficiency_mass;
	  else if (histoString == "GENRECOeff_methodPTjet") histo[v] = efficiency_pt;
	  else if (histoString == "methodM_Pruned1") histo[v] = methodM_eff_pruned1;
	  else if (histoString == "methodM_Pruned2") histo[v] = methodM_eff_pruned2;
	  else if (histoString == "methodM_PrunedTOT") histo[v] = methodM_eff_prunedTOT;
	  else if (histoString == "methodPTjet_Pruned1") histo[v] = methodPTjet_eff_pruned1;
	  else if (histoString == "methodPTjet_Pruned2") histo[v] = methodPTjet_eff_pruned2;
	  else if (histoString == "methodPTjet_PrunedTOT") histo[v] = methodPTjet_eff_prunedTOT;
	  else if (histoString == "methodM_M_0H") histo[v] = BESTpair_methodM_M_0H;
	  else if (histoString == "methodM_M_1H_eff1") histo[v] = BESTpair_methodM_M_1H_eff1;
	  else if (histoString == "methodM_M_1H_eff0") histo[v] = BESTpair_methodM_M_1H_eff0;
	  else if (histoString == "methodM_M_2H_eff2") histo[v] = BESTpair_methodM_M_2H_eff2;
	  else if (histoString == "methodM_M_2H_eff1") histo[v] = BESTpair_methodM_M_2H_eff1;
	  else if (histoString == "methodM_M_2H_eff0") histo[v] = BESTpair_methodM_M_2H_eff0;
	  else if (histoString == "methodPTjet_M_0H") histo[v] = BESTpair_methodPTjet_M_0H;
	  else if (histoString == "methodPTjet_M_1H_eff1") histo[v] = BESTpair_methodPTjet_M_1H_eff1;
	  else if (histoString == "methodPTjet_M_1H_eff0") histo[v] = BESTpair_methodPTjet_M_1H_eff0;
	  else if (histoString == "methodPTjet_M_2H_eff2") histo[v] = BESTpair_methodPTjet_M_2H_eff2;
	  else if (histoString == "methodPTjet_M_2H_eff1") histo[v] = BESTpair_methodPTjet_M_2H_eff1;
	  else if (histoString == "methodPTjet_M_2H_eff0") histo[v] = BESTpair_methodPTjet_M_2H_eff0;
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
  TFile* fOut = new TFile(outputFile.c_str(),"recreate");  
  fOut->cd();
  
  for (int nfs = 0; nfs < nFinalState+1; nfs ++)
    {
      for (int v = 0; v < nHisto; v ++)
	{
	  h1[v][nfs]->Write();
	}
    }


  std::cout << "yield 0H " << yield_0H << endl;
  std::cout << "count 0H " << count_0H << endl;
  std::cout << "yield 1H " << yield_1H << endl;
  std::cout << "count 1H " << count_1H << endl;
  std::cout << "yield 2H " << yield_2H << endl;
  std::cout << "count 2H " << count_2H << endl;
  
  std::cout << "YIELD " << h1[11][4]->Integral(0,500) << endl; 
    
} //end doHisto


/*
void doPlot(const std::string outputFile)
{
  

}
*/


void ZZbb_analysis(const std::string inputFileMC, const std::string outputFile, double lumi=140) 
{
  doHisto(inputFileMC, outputFile, lumi);

  
  //  doPlot(ouputFile);
}
