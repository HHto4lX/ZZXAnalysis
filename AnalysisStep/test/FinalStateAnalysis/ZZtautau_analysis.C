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

#define DEBUG 0
#define VERBOSE2 1
#define MERGE2E2MU 1

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

const int nHisto = 7;

Histo1D myHisto1D[nHisto] = {
  {"M4L", "m_{4L}", "Events", "", 0, 500, 100, 0, 0},
  {"MZ1", "m_{Z1}", "Events", "", 0, 200, 100, 0, 0},
  {"MZ2", "m_{Z2}", "Events", "", 0, 200, 100, 0, 0},
  {"pt4L", "p_{T}^{4L}", "Events", "", 0, 1000, 200, 0, 0},
  {"eta4L", "#eta^{4L}", "Events", "", -4, 4, 80, 0, 0},
  {"yield", "", "Events", "", 0, 10, 10, 0, 0},
  {"Mtau", "m_{#tautau}", "Events", "", 0, 200, 100, 0, 0},
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
  Short_t Z1Flav;
  Short_t Z2Flav;
  
  vector<Float_t> *tauPt = 0;
  vector<Float_t> *tauEta = 0;

  Int_t ntau = 0;
  
  Int_t YIELD = 0;
  
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
  inputTree->SetBranchAddress("ZZEta", &ZZEta);
  inputTree->SetBranchAddress("tauPt", &tauPt);
  inputTree->SetBranchAddress("tauEta", &tauEta);

  int entries = inputTree->GetEntries();
  std::cout<<"Processing file: "<< inputFileMC.c_str() << "\nNumber of entries: " << entries << endl;
  
  //  entries = 100;

  for (Long64_t entry = 0; entry < entries; entry++)
    {  
      if (DEBUG && (entry > 100)) break;
      
      inputTree->GetEntry(entry);
      
      if(LepEta->size() != 4)
	{
	  cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	  continue;
	}
      if( !(ZZsel >= 90) ) continue;

      Float_t kfactor = 1.;
      
      Double_t eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight ; // add genBR for SIGNAL
    
      //Double_t eventWeight = 1;
      std::cout << "-----------------------------------------------------------------------\nevent: " << entry << " weight: " << eventWeight << endl;
      
      std::cout << "lumi " << lumi << endl;
      std::cout << "xsec " << xsec << endl;
      std::cout << "genBR " << genBR << endl;
      std::cout << "genSumwweight " << gen_sumWeights << endl;
      std::cout << "overallEventwweight " << overallEventWeight << endl;
      std::cout << "Bin 0 " << NGenEvt << endl;
      
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
      

      ntau = 0;
      YIELD = 0;

      std::cout <<"\n tau vec size " << tauPt->size() << endl;
      
      if (tauPt->size() > 1)
	{
	  for (int j = 0; j < tauPt->size(); j++)
	    {
	      if (tauPt->at(j) > 20)
		ntau++;
	    }	  
	}
      std::cout << "N tau in acceptance " << ntau << endl;

      if (ntau > 1) YIELD++;


      std::cout << "YIELD " << YIELD << endl;
    
            
      Float_t histo[nHisto];
      
      for(int v = 0; v < nHisto; v++)
	{
	  string histoString = myHisto1D[v].Name.c_str();
	  if      (histoString == "M4L") histo[v] = ZZMass;
	  else if (histoString == "MZ1") histo[v] = Z1Mass;
	  else if (histoString == "MZ2") histo[v] = Z2Mass;
	  else if (histoString == "pt4L") histo[v] = ZZPt;
	  else if (histoString == "eta4L") histo[v] = ZZEta;
	  else if (histoString == "eta4L") histo[v] = ZZEta;
	  else if (histoString == "yield") histo[v] = YIELD;
	  else continue;
		 
	  h1[v][currentFinalState]->Fill(histo[v], eventWeight);
	}
      
      
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
  
    
} //end doHisto


/*
void doPlot(const std::string outputFile)
{
  

}
*/


void ZZtautau_analysis() 
{
  double lumi = 140; // full Run2 Lumi                                                                                                                                            

  string inputFilePath = "/eos/user/a/acappati/samples_4lX/190626/";
  string inputFileName[] = {//"HH4Ltautau",
			    "ggH125",
			    "VBFH125",
			    "WplusH125",
			    "WminusH125",
			    "ZH125",
			    "bbH125",
			    "ttH125",
			    "ggTo4e_Contin_MCFM701",
			    "ggTo4mu_Contin_MCFM701",
			    "ggTo4tau_Contin_MCFM701",
			    "ggTo2e2mu_Contin_MCFM701",
			    "ggTo2e2tau_Contin_MCFM701",
			    "ggTo2mu2tau_Contin_MCFM701",
			    "ZZTo4lext1",
			    "TTZJets_M10_MLMext1",
			    "TTZToLL_M1to1O_MLM",
			    "TTWJetsToLNu",
			    //"DYJetsToLL_M50",                                                                                                                                  
			    //"TTTo2L2Nu",                                                                                                                                        
  };

  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<<nInputFiles<<endl;

  string outputFilePath = "histos_4ltautau";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir                                                                                                        

  //call function                                                                                                                                                                  
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doHisto(Form("%s%s%s",inputFilePath.c_str(),inputFileName[i].c_str(),"/ZZXAnalysis.root"), Form("%s%s%s%s",outputFilePath.c_str(), "/histos_", (inputFileName[i]).c_str(), ".root"), lumi);
  }





  //  doPlot(ouputFile);
}
