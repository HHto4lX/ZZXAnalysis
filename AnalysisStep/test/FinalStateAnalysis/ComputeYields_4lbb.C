// *********************************************
//
// run with: root -l -b -q ComputeYields_4lbb.C++
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
#include "map"


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

  
void doHisto(TString inputFileMC, TString outputFile, double lumi=1)
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
  Float_t weight;

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
  
  // variables for categories
  float jetPt[99];        
  float jetEta[99];       
  float jetIsBTagged[99];

  TH1F* h1[nFinalState+1];

  TH1F* h_Event_4L[nFinalState+1];
  TH1F* h_Event_4LMassCut[nFinalState+1];
  TH1F* h_Event_4LCR[nFinalState+1];
  TH1F* h_Event_4LMassCut2Jets[nFinalState+1];
  TH1F* h_Event_4LMassCut2Jets1BTag[nFinalState+1];
  
  for (int fs = 0; fs < (nFinalState+1); fs++ )
    {
      h1[fs] = new TH1F(Form("h1_%s", FinalState[fs].c_str()),"",1,0,1);
      h1[fs]->Sumw2(true); 
      h_Event_4L[fs] = new TH1F(Form("h_Events_4L_%s", FinalState[fs].c_str()),"",1,0,1);
      h_Event_4L[fs]->Sumw2(true);
      h_Event_4LCR[fs] = new TH1F(Form("h_Event_4LCR_%s", FinalState[fs].c_str()),"",1,0,1);
      h_Event_4LCR[fs]->Sumw2(true);
      h_Event_4LMassCut[fs] = new TH1F(Form("h_Event_4LMassCut_%s", FinalState[fs].c_str()),"",1,0,1);
      h_Event_4LMassCut[fs]->Sumw2(true);
      h_Event_4LMassCut2Jets[fs] = new TH1F(Form("h_Event_4LMassCut2Jets_%s", FinalState[fs].c_str()),"",1,0,1);
      h_Event_4LMassCut2Jets[fs]->Sumw2(true);
      h_Event_4LMassCut2Jets1BTag[fs] = new TH1F(Form("h_Event_4LMassCut2Jets1BTag_%s", FinalState[fs].c_str()),"",1,0,1);
      h_Event_4LMassCut2Jets1BTag[fs]->Sumw2(true);
    }
  
  int currentFinalState;

  inputFile =  TFile::Open( inputFileMC );
  
  if( !(inputFileMC.Contains("Z+X")) )
    { 
      hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
      NGenEvt = (Float_t)hCounters->GetBinContent(1);
      gen_sumWeights = (Float_t)hCounters->GetBinContent(40);
      partialSampleWeight = lumi * 1000 / gen_sumWeights;
      inputTree = (TTree*)inputFile->Get("ZZTree/candTree");
    }
  
  else if (inputFileMC.Contains("Z+X"))
    {
      inputTree = (TTree*)inputFile->Get("candTree");
    }

  // set address of all branches
  if( !(inputFileMC.Contains("Z+X")) )
    {
      //      std::cout << "NOT Z+X" << endl;
      inputTree->SetBranchAddress("genBR", &genBR);
      inputTree->SetBranchAddress("RunNumber", &nRun);
      inputTree->SetBranchAddress("EventNumber", &nEvent);
      inputTree->SetBranchAddress("LumiNumber", &nLumi);
      inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
      inputTree->SetBranchAddress("xsec", &xsec);
    }
  if( inputFileMC.Contains("Z+X") )
    {
      inputTree->SetBranchAddress("weight", &weight);
    }
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
  
  std::cout << "entries: " << entries << endl;
  
  // --- loop over tree entries
  for (Long64_t entry = 0; entry < entries; entry++)
    { 
      inputTree->GetEntry(entry);
      
      if(LepEta->size() != 4)
	{
	  cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	  continue;
	}
      if( !(ZZsel >= 0) ) continue;
      
      Float_t kfactor = 1.;
      if(inputFileMC.Contains("ggTo"))       { kfactor = KFactor_QCD_ggZZ_Nominal; }   //ggZZ samples
      else if(inputFileMC.Contains("ZZTo4l")){ kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; }               //qqZZ samples
      
      Double_t eventWeight = 1.;
      
      if( !(inputFileMC.Contains("Z+X")) ) 
      	{ 
      	  eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight ;
      	}
      if( inputFileMC.Contains("Z+X") ) 
      	{
      	  eventWeight = weight ;
      	}
      
      
      // 4l final state check    
      currentFinalState = -1;
      if( !(inputFileMC.Contains("Z+X")) )
	{
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
	}

      else  if ( inputFileMC.Contains("Z+X") )
	{
          if (Z1Flav == -121)
            {
              if (Z2Flav == +121)
                currentFinalState = fs_4e;
              else if ( Z2Flav == +169)
                currentFinalState = fs_2e2mu;
              else
                cerr << "error in event " <<  endl;
            }
          else if (Z1Flav == -169)
            {
              if (Z2Flav == +121)
                currentFinalState = fs_2mu2e;
              else if (Z2Flav == +169)
                currentFinalState = fs_4mu;
              else
		cerr << "error in event " <<  endl;
            }
          else
            {
	      cerr << "error in event " <<  endl;
            }
          if(MERGE2E2MU && ( currentFinalState == fs_2mu2e )) currentFinalState = fs_2e2mu;
        }
      

      // FILL CUT #1      
      h_Event_4L[currentFinalState]->Fill(0.5, eventWeight);


      // --- Mass Cut      
      if (ZZMass < 115 || ZZMass > 135) 
	{
	  h_Event_4LCR[currentFinalState]->Fill(0.5, eventWeight);
	  continue;
	}

      // FILL CUT #2
      h_Event_4LMassCut[currentFinalState]->Fill(0.5, eventWeight);
      
      short nJets20 = 0.;
      for(UInt_t i = 0; i < JetPt->size(); i++)
	{
	  jetPt[i]        = JetPt->at(i);
	  jetEta[i]       = JetEta->at(i);
	  jetIsBTagged[i] = JetIsBTagged->at(i);
	  nJets20++;
	}

      // check on Jets Variables and btagging                                                                                                                                                                  
      Int_t nBtaggedJets = 0;
      Int_t nGoodJets    = 0;

      for(int i = 0; i < nJets20; i++)
	{
	  if( abs(jetEta[i]) > 2.4 ) continue;
	  if( jetPt[i] < 30. ) continue;

	  nGoodJets++;
	  if(jetIsBTagged[i] > 0) nBtaggedJets++; 
	}
      
      
      // FILL CUT #3 
      if( nGoodJets > 1 ) h_Event_4LMassCut2Jets[currentFinalState]->Fill(0.5, eventWeight);

      // FILL CUT #4
      if( nGoodJets > 1 && nBtaggedJets > 0 ) h_Event_4LMassCut2Jets1BTag[currentFinalState]->Fill(0.5, eventWeight);
      

      // --- fill yields histos      
      if (nGoodJets > 1 && nBtaggedJets > 0)
      h1[currentFinalState]->Fill(0.5, eventWeight);

      
    } //end loop entries
  

  // fill inclusive histo (4L)  
  for(int f = 0; f < nFinalState; f++) 
    {
      h1[nFinalState]->Add( h1[f] );
      h_Event_4L[nFinalState]->Add ( h_Event_4L[f] );
      h_Event_4LCR[nFinalState]->Add ( h_Event_4LCR[f] );
      h_Event_4LMassCut[nFinalState]->Add ( h_Event_4LMassCut[f] ) ;
      h_Event_4LMassCut2Jets[nFinalState]->Add ( h_Event_4LMassCut2Jets[f] ) ;
      h_Event_4LMassCut2Jets1BTag[nFinalState]->Add ( h_Event_4LMassCut2Jets1BTag[f] ) ;
    }
  

  std::cout << "Output file "<< outputFile << endl;
  TFile* fOut = new TFile(outputFile ,"recreate");  
  fOut->cd();
  
  for (int nfs = 0; nfs < nFinalState+1; nfs ++)
    {
      h1[nfs]->Write();
    }

  std::cout << " FINAL STATE : " << nFinalState << endl;  
  for (int nfs = 0; nfs < nFinalState+1; nfs ++)
    {
      std::cout << "\n\nYIELDS for final state: " << nfs
		<< "\nnEvents 4L: " << h_Event_4L[nfs]->GetBinContent(1)
		<< "\nnEvents 4L CR: " << h_Event_4LCR[nfs]->GetBinContent(1)
		<< "\nnEvents 4L Mass Cut: " << h_Event_4LMassCut[nfs]->GetBinContent(1)
		<< "\nnEvents 4L Mass Cut 2 Jets: " << h_Event_4LMassCut2Jets[nfs]->GetBinContent(1) 
		<< "\nnEvents 4L Mass Cut 2 Jets 1 BTag: " << h_Event_4LMassCut2Jets1BTag[nfs]->GetBinContent(1) << endl;
    }
  
  std::cout << "---------" 
	    << "\nEventi old histo 4L: " << h1[4]->GetBinContent(1) << endl;
 
  
  fOut->Close(); 
  
} //end doHisto


/*
void doCount(TString inputfile, TString inputSample, std::vector<double> > & map_chan, std::map<TString,std::vector<TString> > & map_names)
{
  TFile * inputFile =  TFile::Open( inputfile );
  
  enum FinalState {fs_4mu=0, fs_4e=1, fs_2e2mu=2, fs_4L=3};  // 4mu, 4e, 2e2mu, 2mu2e
  const int nFinalState = 4;
  string FinalState[nFinalState] = {"4mu", "4e","2e2mu", "4L"};

  for (int nF = 0; nF < nFinalState; nF++ )
    {
      std::string histo_name = "h1_" + FinalState[nF] + "_HH4lbbTagged";
      TH1F *h_temp = (TH1F*)inputFile->Get( histo_name.c_str() );
      
      double counter = h_temp->GetBinContent(1);
      
      map_chan[FinalState[nF]].push_back(counter);
      map_names[FinalState[nF]].push_back(inputSample); 
    }
  
  inputFile->Close();
  
} //end doCount
*/


void ComputeYields_4lbb() 
{
  
  //double lumi = 140; // full Run2 Lumi
  //  double lumi = 35.92; // 2016 data
  double lumi = 59.74; // 2018 data 
  
  //TString inputFilePath = "/eos/user/a/acappati/samples_4lX/190626/";

  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/allsamples/";
  //  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/190924/";
  TString inputFileName[] = {"HH4lbb",
                            "ggH125",
			    "VBFH125",
			    "ZH125",
			    "bbH125",
 			    "ttH125",
                             // "ZZTo4lext1",
			     // "TTWJetsToLNu",
                             // "Z+X",
			     // "WWZ",
                             // "WZZ",
                             // "ZZZ",
			     // "WplusH125",
 			     // "WminusH125",
			     // "TTZJets_M10_MLMext1",
			     // "TTZToLL_M1to1O_MLM",
			     // "ggTo4e_Contin_MCFM701",
			     // "ggTo4mu_Contin_MCFM701",
			     // "ggTo4tau_Contin_MCFM701",
			     // "ggTo2e2mu_Contin_MCFM701",
			     // "ggTo2e2tau_Contin_MCFM701",
			     // "ggTo2mu2tau_Contin_MCFM701",
                            };


  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<<nInputFiles<<endl;

  string outputFilePath = "Yields_histos";
  gSystem->Exec(("mkdir -p "+ outputFilePath).c_str()); // create output dir

  /*  
  std::ofstream outFile_4mubb("Yields_4mubb.txt");
  std::ofstream outFile_4ebb("Yields_4ebb.txt");
  std::ofstream outFile_2e2mubb("Yields_2e2mubb.txt");
  std::ofstream outFile_4Lbb("Yields_4Lbb.txt");


  std::map<TString, std::vector<double> > map_chan;
  std::map<TString, std::vector<TString> > map_names;
  */
  
  //call function
  for(UInt_t i=0; i<nInputFiles; i++)
    {
      cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
  
      doHisto(inputFilePath + inputFileName[i] + "/ZZXAnalysis.root" , outputFilePath + "/histos_" + inputFileName[i] + ".root", lumi);

      //    doCount(outputFilePath + "/histos_" + inputFileName[i] + ".root", inputFileName[i], map_chan, map_names);
    }

  /*
  int sizes = map_names["4mu"].size();
  for (int i=0; i< sizes; ++i )
    {
      outFile_4mubb   << map_names["4mu"][i]   << '\t';
      outFile_4ebb    << map_names["4e"][i]   << '\t';
      outFile_2e2mubb << map_names["2e2mu"][i]   << '\t';
      outFile_4Lbb    << map_names["4L"][i]   << '\t';
    } 
  
  for (int i=0; i< sizes; ++i ) 
    {
      outFile_4mubb   << map_chan["4mu"][i]   << '\t'; 
      outFile_4ebb    << map_chan["4e"][i]    << '\t';
      outFile_2e2mubb << map_chan["2e2mu"][i] << '\t';
      outFile_4Lbb    << map_chan["4L"][i]    << '\t';
    }  
  */

}
