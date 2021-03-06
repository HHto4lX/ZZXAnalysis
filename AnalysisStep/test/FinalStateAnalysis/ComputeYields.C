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
#include "map"

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
 
  vector<Float_t> *JetIsBTaggedWithSF = 0;
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

  Short_t nExtraLep = 0;;
  vector<Float_t> *ExtraLepPt = 0;
  vector<Float_t> *ExtraLepEta = 0;
  vector<Int_t> *ExtraLepLepId = 0;

  Short_t nPhotons = 0;

  Short_t nTaus  =  0;

  // variables for categories
  float jetPt[99];        
  float jetEta[99];       
  float jetIsBTaggedWithSF[99];
  float extraLepPt[99];   
  float extraLepEta[99];  
  int   extraLepLepId[99];  


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

  inputFile =  TFile::Open( inputFileMC );
  
  if( !(inputFileMC.Contains("Z+X")) )
    { 
      //      std::cout << "NOT Z+X" << endl;
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
  inputTree->SetBranchAddress("JetIsBtaggedWithSF",  &JetIsBTaggedWithSF);
  if( !(inputFileMC.Contains("Z+X")) )
    {
      //      std::cout << "NOT Z+X" << endl;
      inputTree->SetBranchAddress("nExtraLep", &nExtraLep);
      inputTree->SetBranchAddress("ExtraLepPt", &ExtraLepPt);
      inputTree->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
      inputTree->SetBranchAddress("ExtraLepLepId", &ExtraLepLepId);
      inputTree->SetBranchAddress("nPhotons", &nPhotons);
      inputTree->SetBranchAddress("nTaus", &nTaus);
    }
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
      
      // --- Mass Cut      
      //      if (ZZMass < 115 || ZZMass > 135) continue;


      // --- find category
      short nJets20 = 0.;
      for(UInt_t i = 0; i < JetPt->size(); i++){
        jetPt[i]        = JetPt->at(i);
        jetEta[i]       = JetEta->at(i);
        jetIsBTaggedWithSF[i] = JetIsBTaggedWithSF->at(i);
        nJets20++;
      }
      if( !(inputFileMC.Contains("Z+X")) )
	{
	  //	  std::cout << "NOT Z+X" << endl;
	  for(int i = 0; i < nExtraLep; i++){
	    extraLepPt[i]    = ExtraLepPt->at(i);
	    extraLepEta[i]   = ExtraLepEta->at(i);
	    extraLepLepId[i] = ExtraLepLepId->at(i);
	  }
	}
      currentCategory = categoryHH(nJets20,
                                   jetPt,
                                   jetEta,
                                   jetIsBTaggedWithSF,
                                   nExtraLep,
                                   extraLepPt,
                                   extraLepEta,
                                   extraLepLepId
				   //        nPhotons
				   //         nTaus
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
  TFile* fOut = new TFile(outputFile ,"recreate");  
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

  fOut->Close(); 
    
} //end doHisto


void  doCount(TString inputfile, TString inputSample, std::map<TString,std::vector<double> > & map_chan, std::map<TString,std::vector<TString> > & map_names, double counter_VVV[4], double counter_WH[4], double counter_TTZ[4], double counter_ggZZ[4])
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
      
      if ( inputfile.Contains("ZZZ") || inputfile.Contains("WWZ") || inputfile.Contains("WZZ") ) counter_VVV[nF] += counter;
      else if ( inputfile.Contains("TTZToLL") || inputfile.Contains("TTZJets") ) counter_TTZ[nF] += counter;
      else if ( inputfile.Contains("WplusH") || inputfile.Contains("WminusH") ) counter_WH[nF] += counter;
      else if ( inputfile.Contains("ggTo") ) counter_ggZZ[nF] += counter;

      else 
	{
	  map_chan[FinalState[nF]].push_back(counter);
	  map_names[FinalState[nF]].push_back(inputSample);
	}
    }
  
  
  inputFile->Close();
  
} //end doCount

void ComputeYields() 
{

  //double lumi = 140; // full Run2 Lumi
  //  double lumi = 35.92; // 2016 data
  double lumi = 59.74; // 2018 data 

  //TString inputFilePath = "/eos/user/a/acappati/samples_4lX/190626/";

  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/allsamples/";
  //  TString inputFilePath = "/eos/user/a/acappati/samples_4lX/190924/";
  TString inputFileName[] = {"HH4lbb",
                            // "HH4lww",
                            // "HH4lgammagamma",
                            // "HH4ltautau",
                            "ggH125",
			    "VBFH125",
			    "ZH125",
			    "bbH125",
			    "ttH125",
                             "ZZTo4lext1",
			     "TTWJetsToLNu",
                             "Z+X",
			     "WWZ",
                             "WZZ",
                             "ZZZ",
			     "WplusH125",
 			     "WminusH125",
			     "TTZJets_M10_MLMext1",
			     "TTZToLL_M1to1O_MLM",
			     "ggTo4e_Contin_MCFM701",
			     "ggTo4mu_Contin_MCFM701",
			     "ggTo4tau_Contin_MCFM701",
			     "ggTo2e2mu_Contin_MCFM701",
			     "ggTo2e2tau_Contin_MCFM701",
			     "ggTo2mu2tau_Contin_MCFM701",
                            };


  size_t nInputFiles = sizeof(inputFileName)/sizeof(inputFileName[0]);
  cout<<nInputFiles<<endl;

  string outputFilePath = "Yields_histos";
  gSystem->Exec(("mkdir -p "+ outputFilePath).c_str()); // create output dir
  
  std::ofstream outFile_4mubb("Yields_forCombine_4mubb.txt");
  std::ofstream outFile_4ebb("Yields_forCombine_4ebb.txt");
  std::ofstream outFile_2e2mubb("Yields_forCombine_2e2mubb.txt");
  std::ofstream outFile_4Lbb("Yields_forCombine_4Lbb.txt");

  double counter_VVV[4] = {0,0,0,0};
  double counter_WH[4] = {0,0,0,0};  
  double counter_TTZ[4] = {0,0,0,0};
  double counter_ggZZ[4] = {0,0,0,0};

  // outFile_4mubb << "HH4Lbb\tggH\tVBFH\tZH\tbbH\tttH\tqqZ\tttW\tWH\tttZ\tggZZ\n"; 
  // outFile_4ebb << "HH4Lbb\tggH\tVBFH\tZH\tbbH\tttH\tqqZ\tttW\tWH\tttZ\tggZZ\n"; 
  // outFile_2e2mubb << "HH4Lbb\tggH\tVBFH\tZH\tbbH\tttH\tqqZ\tttW\tWH\tttZ\tggZZ\n"; 
  // outFile_4Lbb << "HH4Lbb\tggH\tVBFH\tZH\tbbH\tttH\tqqZ\tttW\tWH\tttZ\tggZZ\n";


  std::map<TString, std::vector<double> > map_chan;
  std::map<TString, std::vector<TString> > map_names;
   
   //call function
  for(UInt_t i=0; i<nInputFiles; i++)
    {
      cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    
      doHisto(inputFilePath + inputFileName[i] + "/ZZXAnalysis.root" , outputFilePath + "/histos_" + inputFileName[i] + ".root", lumi);
      //doCount(outputFilePath + "/histos_" + inputFileName[i] + ".root", outFile_4mubb, outFile_4ebb, outFile_2e2mubb, outFile_4Lbb, counter_WH, counter_TTZ, counter_ggZZ);
      doCount(outputFilePath + "/histos_" + inputFileName[i] + ".root", inputFileName[i], map_chan, map_names, counter_VVV, counter_WH, counter_TTZ, counter_ggZZ);
    }

  int sizes = map_names["4mu"].size();
  for (int i=0; i< sizes; ++i )
    {
      outFile_4mubb   << map_names["4mu"][i]   << '\t';
      outFile_4ebb    << map_names["4e"][i]   << '\t';
      outFile_2e2mubb << map_names["2e2mu"][i]   << '\t';
      outFile_4Lbb    << map_names["4L"][i]   << '\t';
    } 
  outFile_4mubb   <<  "VVV\tWH\tTTZ\tggZZ\n";
  outFile_4ebb    <<  "VVV\tWH\tTTZ\tggZZ\n";
  outFile_2e2mubb <<  "VVV\tWH\tTTZ\tggZZ\n";
  outFile_4Lbb    <<  "VVV\tWH\tTTZ\tggZZ\n";

  for (int i=0; i< sizes; ++i ) 
    {
      outFile_4mubb   << map_chan["4mu"][i]   << '\t'; 
      outFile_4ebb    << map_chan["4e"][i]    << '\t';
      outFile_2e2mubb << map_chan["2e2mu"][i] << '\t';
      outFile_4Lbb    << map_chan["4L"][i]    << '\t';
    }  
  for (int nF = 0; nF < nFinalState; nF++ )
    {
      if (nF == 0)  outFile_4mubb << counter_VVV[nF] << "\t" << counter_WH[nF] << "\t" << counter_TTZ[nF] << "\t" << counter_ggZZ[nF] << "\n";
      else if (nF == 1) outFile_4ebb << counter_VVV[nF] << "\t" << counter_WH[nF] << "\t" << counter_TTZ[nF] << "\t" << counter_ggZZ[nF] << "\n";
      else if (nF == 2) outFile_2e2mubb << counter_VVV[nF] << "\t" << counter_WH[nF] << "\t" << counter_TTZ[nF] << "\t" << counter_ggZZ[nF] << "\n";
      else if (nF == 3) outFile_4Lbb << counter_VVV[nF] << "\t" << counter_WH[nF] << "\t" << counter_TTZ[nF] << "\t" << counter_ggZZ[nF] << "\n";
    }
}
