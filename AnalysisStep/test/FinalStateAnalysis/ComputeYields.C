// *********************************************
//
// run with: root -l -b -q ComputeYields.C++
//
// *********************************************

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

#include "ZZXAnalysis/AnalysisStep/src/Category.cc"

using namespace std;

#define MERGE2E2MU 1

double DELTAR = .4;

// Final State scan
enum FinalState {fs_4mu=0, fs_4e=1, fs_2e2mu=2, fs_2mu2e=3};  // 4mu, 4e, 2e2mu, 2mu2e
const int nFinalState = 4;
string FinalState[nFinalState+1] = {"4mu", "4e","2e2mu","2mu2e", "4L"};
string FinalStateLabel[nFinalState+1] = {"4#mu", "4e", "2e2#mu", "2#mu2e", "4L"};
int FinalStateColor[nFinalState+1] = {2, 4, 8, 3, 1};
int FinalStateMarkerStyle[nFinalState+1] = {20, 22, 21, 33, 29};


// Categories
const int nCat = 5;
string sCategory[nCat] = {
  "HHUntagged",
  "HH4lbbTagged",
  "HH4lWWTagged",
  "HH4lgammagammaTagged",
  "HH4ltautauTagged"
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



  TH1F* h1[nFinalState+1][nCat];
  
  for (int fs = 0; fs < (nFinalState+1); fs++ )
    {
      for (int c = 0; c < nCat; c++)
	{
	  h1[fs][c] = new TH1F(Form("h1_%s_%s", FinalState[fs].c_str(), sCategory[c].c_str()),"",1,0,1);
	  h1[fs][c]->Sumw2(true);
	}       
    }
  
  int currentFinalState;
  int currentCategory;

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
  
  // --- loop over tree entries
  for (Long64_t entry = 0; entry < entries; entry++)
    {  
 
      inputTree->GetEntry(entry);
      

      if(LepEta->size() != 4)
	{
	  cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	  continue;
	}
      if( !(ZZsel >= 90) ) continue;


      Float_t kfactor = 1.;
      

      Double_t eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight ;

 
      // 4l final state check    
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


      // --- Mass Cut      
      if (ZZMass < 115 || ZZMass > 135) continue;


      // --- choose category
      currentCategory = categoryHH(JetPt,
                                   JetEta,
                                   JetIsBTagged
      				  );


      // --- fill yields histos
      h1[currentFinalState][currentCategory]->Fill(0.5, eventWeight);

      
    } //end loop entries
  

  // fill inclusive histo (4L)  
  for(int f = 0; f < nFinalState; f++)
  {
    for(int c = 0; c < nCat; c++)
    {
      h1[nFinalState][c]->Add( h1[f][c] );
    }
  }
 
  
  
  std::cout << "Output file "<< outputFile << endl;
  TFile* fOut = new TFile(outputFile.c_str(),"recreate");  
  fOut->cd();
  
  for (int nfs = 0; nfs < nFinalState+1; nfs ++)
    {
      for (int c = 0; c < nCat; c++)
	{
	  h1[nfs][c]->Write();
	}
    }
  
 
  // print yields per category
  for(int c = 0; c < nCat; c++)
  {
    std::cout << "YIELD "<< sCategory[c] << " " << h1[4][c]->GetBinContent(1) << endl; 
  }
    
} //end doHisto


void ComputeYields() 
{

  double lumi = 140; // full Run2 Lumi

  string inputFilePath = "/eos/user/a/acappati/samples_4lX/190829/";
  string inputFileName[] = {"HH4lbb",
                            "ggH125",
			    "VBFH125",
			    "WplusH125",
			    "WminusH125",
			    "ZH125",
			    "bbH125",
			    "ttH125",
			    // "ggTo4e_Contin_MCFM701",
			    // "ggTo4mu_Contin_MCFM701",
			    // "ggTo4tau_Contin_MCFM701",
			    // "ggTo2e2mu_Contin_MCFM701",
			    // "ggTo2e2tau_Contin_MCFM701",
			    // "ggTo2mu2tau_Contin_MCFM701",
			    "ZZTo4lext1",
			    "TTZJets_M10_MLMext1",
			    "TTZToLL_M1to1O_MLM",
			    "TTWJetsToLNu",
                            //"DYJetsToLL_M50",
                            //"TTTo2L2Nu",
                            };

  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<<nInputFiles<<endl;

  string outputFilePath = "Yields_histos";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir
  

  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doHisto(Form("%s%s%s",inputFilePath.c_str(),inputFileName[i].c_str(),"/ZZXAnalysis.root"), Form("%s%s%s%s",outputFilePath.c_str(), "/histos_", (inputFileName[i]).c_str(), ".root"), lumi);
  }
  
}
