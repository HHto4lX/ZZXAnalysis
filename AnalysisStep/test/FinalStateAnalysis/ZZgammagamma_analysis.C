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
#include <fstream>

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

const int nHisto = 11;

Histo1D myHisto1D[nHisto] = {
  {"NPhotons","# photons", "Events", "", 0, 20, 20, 0, 0},
  {"NPhotonsID","# photons with ID", "Events", "", 0, 20, 20, 0, 0},
  {"M4L", "m_{4L}", "Events", "", 0, 500, 100, 0, 0},
  {"MZ1", "m_{Z1}", "Events", "", 0, 200, 100, 0, 0},
  {"MZ2", "m_{Z2}", "Events", "", 0, 200, 100, 0, 0},
  {"pt4L", "p_{T}^{4L}", "Events", "", 0, 1000, 200, 0, 0},
  {"eta4L", "#eta^{4L}", "Events", "", -8, 8, 80, 0, 0},
  {"Mpairmgamma", "m_{pair#gamma}", "Events", "", 0, 500, 250, 0, 0},
  {"Mpairptgamma", "p_{T}^{pair#gamma}", "Events", "", 0, 1000, 200, 0, 0},
  {"NGenPhotons","# gen photons", "Events", "", 0, 20, 20, 0, 0},
  {"NGenPhotonsInAcc","# gen photons", "Events", "", 0, 20, 20, 0, 0},
};
  
float deltaR(float  p1, float p2, float e1, float e2)
{
 float dp = std::abs(p1 - p2);
 if (dp > M_PI) dp -= 2 * M_PI;
 return (e1 - e2) * (e1 - e2) + dp * dp;
}

void doHisto(const std::string inputFileMC, const std::string outputFile, double lumi=1)
{

  ofstream myfile;
  myfile.open ("HH4lgg.txt");

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
  vector<Float_t> *LepPhi = 0;
  vector<Float_t> *LepPt = 0;
  vector<Int_t> *LepLepId = 0;
  Float_t ZZMass;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Float_t ZZPt;
  Float_t ZZEta;
  Short_t Z1Flav;
  Short_t Z2Flav;
  
  vector<Float_t> *photonIsID = 0;
  vector<Float_t> *photonPt = 0;
  vector<Float_t> *photonEta = 0;
  vector<Float_t> *photonPhi= 0;
  vector<Float_t> *photonMass = 0;
  vector<Float_t> *photonPassElectronVeto = 0;

  vector<Float_t> *prunedGenPhoPt = 0;
  vector<Float_t> *prunedGenPhoEta = 0; 
  vector<Float_t> *prunedGenPhoPhi = 0; 
  vector<Float_t> *prunedGenPhoMass = 0;
  vector<Short_t> *prunedGenPhoID = 0;  
  vector<Short_t> *prunedGenPhoMotherID = 0;             
               


  Float_t MASSPairMass;
  Float_t MASSPairPt;

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
  inputTree->SetBranchAddress("LepPhi", &LepPhi);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepLepId", &LepLepId);
  inputTree->SetBranchAddress("ZZMass", &ZZMass);  
  inputTree->SetBranchAddress("Z1Flav", &Z1Flav);
  inputTree->SetBranchAddress("Z2Flav", &Z2Flav);
  inputTree->SetBranchAddress("Z1Mass", &Z1Mass);
  inputTree->SetBranchAddress("Z2Mass", &Z2Mass);
  inputTree->SetBranchAddress("ZZPt", &ZZPt);
  inputTree->SetBranchAddress("ZZEta", &ZZEta);
  inputTree->SetBranchAddress("photonPt", &photonPt);
  inputTree->SetBranchAddress("photonEta", &photonEta);
  inputTree->SetBranchAddress("photonPhi", &photonPhi);
  inputTree->SetBranchAddress("photonMass", &photonMass);
  inputTree->SetBranchAddress("photonPassElectronVeto", &photonPassElectronVeto);
  inputTree->SetBranchAddress("photonPhi",  &photonPhi);
  inputTree->SetBranchAddress("photonIsID",  &photonIsID);
  inputTree->SetBranchAddress("prunedGenPhoPt", &prunedGenPhoPt);
  inputTree->SetBranchAddress("prunedGenPhoEta", &prunedGenPhoEta);
  inputTree->SetBranchAddress("prunedGenPhoPhi", &prunedGenPhoPhi);
  inputTree->SetBranchAddress("prunedGenPhoMass", &prunedGenPhoMass);
  inputTree->SetBranchAddress("prunedGenPhoID", &prunedGenPhoID);
  inputTree->SetBranchAddress("prunedGenPhoMotherID", &prunedGenPhoMotherID);
 
  //inputTree->Draw("prunedGenPhoPt");
  int entries = inputTree->GetEntries();
 
  myfile<<"Processing file: "<< inputFileMC.c_str() << "\nNumber of entries: " << entries << endl;
  
  //  entries = 100;

  double yield = 0;
  for (Long64_t entry = 0; entry < entries; entry++)
    {  
      int nphoton = 0;
      int nphoton_ID = 0;
      int ngenphoton = 0;
      int ngenphoton_inacc = 0;

      if (DEBUG && (entry > 100)) break;
      
      inputTree->GetEntry(entry);
      
      if(LepEta->size() != 4)
	{
	  cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	  continue;
	}
      if( !(ZZsel >= 90) ) continue;

      Float_t kfactor = 1.;
      
      Double_t eventWeight = partialSampleWeight * xsec * kfactor * overallEventWeight ;
    
      //Double_t eventWeight = 1;
      myfile << "-----------------------------------------------------------------------" << endl;
      // myfile << "lumi " << lumi << endl;
      // myfile << "xsec " << xsec << endl;
      // myfile << "genBR " << genBR << endl;
      // myfile << "genSumwweight " << gen_sumWeights << endl;
      // myfile << "overallEventwweight " << overallEventWeight << endl;
      // myfile << "Bin 0 " << NGenEvt << endl;
      
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
      
      myfile <<FinalState[currentFinalState].c_str() <<" final state "<< endl;
    

      // H(gg) selection
      vector<TLorentzVector> PhotonVec;
      vector<TLorentzVector> PhotonPair;
      
      int photonID = 0;      
      float dR; 
      float dRmin;
      
      //count how many gen photons are there
      for (UInt_t k = 0; k < prunedGenPhoPt->size(); k++)
	{
	  ngenphoton++;
	  //Acceptance cuts at gen level
	  if(prunedGenPhoPt->at(k)>10 && fabs(prunedGenPhoEta->at(k))<2.5 && prunedGenPhoMotherID->at(k) == 25) ngenphoton_inacc++;
	 
	}
      

      myfile <<"nphoton size " << photonEta->size() << endl;

      for (UInt_t j = 0; j < photonPt->size(); j++)
	{
	  nphoton++;
	  //myfile << "Photon ID = " << photonIsID->at(j) << " pt = " << photonPt->at(j) << " eta = " << photonEta->at(j)<< " phi = " << photonPhi->at(j)<<endl;
	  dR = 0;
	  dRmin = 1000000;

	  //Look for the smallest deltaR between the photon and the leptons
	  for (UInt_t k = 0; k < LepPt->size(); k++)
	    {
	      dR = deltaR(LepPhi->at(k),photonPhi->at(j),LepEta->at(k),photonEta->at(j)); 
	      if(dR < dRmin) dRmin = dR;
	      //myfile << "dR = " << dR << endl; 
	    }
	   myfile << "dRmin = " << dRmin << endl; 
	      
	   //cuts for efficiency, at reco levels
	   //if ( dRmin > 0.01 && photonIsID->ats(j) >0 && fabs(photonEta->at(j))<2.5 && photonPt->at(j)>10)
	   if (photonPassElectronVeto->at(j)>0.5 && photonIsID->at(j) >0 && fabs(photonEta->at(j))<2.5 && photonPt->at(j)>10)
	     
	     {
	       photonID = 1;
	       nphoton_ID ++;
	       TLorentzVector temp;
	       temp.SetPtEtaPhiM(photonPt->at(j), photonEta->at(j), photonPhi->at(j), photonMass->at(j));
	       PhotonVec.push_back(temp);
	     }
	}
      myfile << PhotonVec.size()<<" photon passing ID" << endl;
      // myfile << "photon n \t" << nphoton << endl;
      // myfile << "ID n \t\t" << nphoton_ID << endl;
      // myfile << "vec size \t" << PhotonVec.size() << endl;
      

      
      // create gamma gamma pairs
      double deltaM = 100000.;
      MASSPairMass = 0;
      MASSPairPt = 0;
      
      //if (PhotonVec.size() < 2) continue;

      //Needed for computing the acceptance
      if (PhotonVec.size() < 2){ 
	MASSPairMass = -1000;
	MASSPairPt = -1000;
	  }      

      else{
	for (UInt_t i = 0; i < PhotonVec.size() -1; i++)
	  {
	    TLorentzVector ph0 = PhotonVec[i];
	    for (UInt_t j = i+1; j < PhotonVec.size(); j++)
	      {
		TLorentzVector ph1 = PhotonVec[j];
		PhotonPair.push_back(ph0 + ph1);
	      }
	  }
	
	myfile << "pair n \t\t" << PhotonPair.size() << endl;     
     
	
	if (PhotonPair.size() < 1) continue;
	
	for (UInt_t i = 0; i< PhotonPair.size(); i++)
	  {
	    myfile << i <<  "  mass: " << PhotonPair[i].M() << " Delta M: " << deltaM << endl; 
	    if ( fabs (PhotonPair[i].M() - 125.) < deltaM ) 
	      {
		deltaM = fabs(PhotonPair[i].M() - 125.);
		MASSPairMass = PhotonPair[i].M();
		MASSPairPt = PhotonPair[i].Pt();
	      }
	  }
      }
      myfile << "AFTER delta M " << deltaM << endl;
      myfile << "MASS " << MASSPairMass << endl;
      myfile << "Z1 Mass " << Z1Mass << " Z2Mass "<< Z2Mass << endl;

      
      Float_t histo[nHisto];
      
      yield += eventWeight;

      for(int v = 0; v < nHisto; v++)
	{
	  string histoString = myHisto1D[v].Name.c_str();
	  if      (histoString == "M4L") histo[v] = ZZMass;
	  else if (histoString == "MZ1") histo[v] = Z1Mass;
	  else if (histoString == "MZ2") histo[v] = Z2Mass;
	  else if (histoString == "pt4L") histo[v] = ZZPt;
	  else if (histoString == "eta4L") histo[v] = ZZEta;
	  else if (histoString == "NPhotons") histo[v] = nphoton;
	  else if (histoString == "NPhotonsID") histo[v] = nphoton_ID;
	  else if (histoString == "Mpairmgamma") histo[v] = MASSPairMass;
	  else if (histoString == "Mpairptgamma") histo[v] = MASSPairPt;
	  else if (histoString == "NGenPhotons") histo[v] = ngenphoton;
	  else if (histoString == "NGenPhotonsInAcc") histo[v] = ngenphoton_inacc;
	  else continue;
		 
	  h1[v][currentFinalState]->Fill(histo[v], eventWeight);
	  
	}
      
      PhotonVec.clear();
      PhotonPair.clear();
      
    } //end loop entries

  myfile << "TOT YIELD " << yield << endl;
  
  
  // fill inclusive histo (4L)
  
  for(int n = 0; n < nHisto; n++)
    {
      for(int f = 0; f < nFinalState; f++)
	{
	  h1[n][nFinalState] -> Add ( h1[n][f] );
	}
    } 
  
  
  myfile << "Output file "<< outputFile << endl;
  myfile.close();
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

void doAccEff(const std::string inputFileMC)
{

 TFile* inputFile;
 inputFile =  TFile::Open( inputFileMC.c_str() );
 
 TH1F* h_4mu = (TH1F*)inputFile->Get("h1_NGenPhotons_4mu");
 TH1F* h_4e = (TH1F*)inputFile->Get("h1_NGenPhotons_4e");
 TH1F* h_2e2mu = (TH1F*)inputFile->Get("h1_NGenPhotons_2e2mu");
 
 TH1F* h_4mu_inacc = (TH1F*)inputFile->Get("h1_NGenPhotonsInAcc_4mu");
 TH1F* h_4e_inacc = (TH1F*)inputFile->Get("h1_NGenPhotonsInAcc_4e");
 TH1F* h_2e2mu_inacc = (TH1F*)inputFile->Get("h1_NGenPhotonsInAcc_2e2mu");
 
 TH1F* h_4mu_reco = (TH1F*)inputFile->Get("h1_NPhotonsID_4mu");
 TH1F* h_4e_reco = (TH1F*)inputFile->Get("h1_NPhotonsID_4e");
 TH1F* h_2e2mu_reco = (TH1F*)inputFile->Get("h1_NPhotonsID_2e2mu");

 Double_t e[3];
 Double_t e_inacc[3];
 Double_t e_reco[3];
 
 Double_t integral[3];
 Double_t integral_inacc[3];
 Double_t integral_reco[3];

  
 integral[0] = h_4mu->IntegralAndError(0, -1, e[0], "");
 integral[1] = h_4e->IntegralAndError(0, -1, e[1], "");
 integral[2] = h_2e2mu->IntegralAndError(0, -1, e[2], "");
 
 integral_inacc[0] = h_4mu_inacc->IntegralAndError(3, 100, e_inacc[0], "");
 integral_inacc[1] = h_4e_inacc->IntegralAndError(3, 100, e_inacc[1], "");
 integral_inacc[2] = h_2e2mu_inacc->IntegralAndError(3, 100, e_inacc[2], "");
 
 integral_reco[0] = h_4mu_reco->IntegralAndError(3, 100, e_reco[0], "");
 integral_reco[1] = h_4e_reco->IntegralAndError(3, 100, e_reco[1], "");
 integral_reco[2] = h_2e2mu_reco->IntegralAndError(3, 100, e_reco[2], "");

 //for(int i =0; i<20; i++){
 //cout << "4mu bin " << i << " " << h_4mu_reco->GetBinContent(i) << " " << h_4mu_inacc->GetBinContent(i) << " " << h_4mu->GetXaxis()->FindBin(i) << endl;
   //cout << "4e bin " << i << " " << h_4e_reco->GetBinContent(i) << " " << h_4e_inacc->GetBinContent(i) << " " << h_4e->GetXaxis()->FindBin(i) << endl;
   //cout << "2e2mu bin " << i << " " << h_2e2mu_reco->GetBinContent(i) << " " << h_2e2mu_inacc->GetBinContent(i) << endl;

 //}

 cout << "4mu " <<integral[0] << "+-" << e[0] << " " << integral_inacc[0] << "+-" << e_inacc[0] << " " << integral_reco[0] << "+-" << e_reco[0] <<" ==> acc = " << integral_inacc[0]/integral[0] << " eff = " << integral_reco[0]/integral_inacc[0] << endl;  
 cout << "4e " << integral[1] << "+-" << e[1] << " " << integral_inacc[1] << "+-" << e_inacc[1] << " " << integral_reco[1] << "+-" << e_reco[1] <<" ==> acc = " << integral_inacc[1]/integral[1] << " eff = " << integral_reco[1]/integral_inacc[1] << endl;    
 cout << "2e2mu " << integral[2] << "+-" << e[2] << " " << integral_inacc[2] << "+-" << e_inacc[2] << " " << integral_reco[2] << "+-" << e_reco[2] <<" ==> acc = " << integral_inacc[2]/integral[2] << " eff = " << integral_reco[2]/integral_inacc[2] << endl; 
 cout << "number of events with more than 2 (reco) photons passng the analysis selection" << endl;
 cout << "4mu " << h_4mu_reco->IntegralAndError(4, 100, e_reco[0], "")/h_4mu_reco->IntegralAndError(0, -1, e_reco[0], "")<<endl;
 cout << "4e " << h_4e_reco->IntegralAndError(4, 100, e_reco[1], "")/h_4e_reco->IntegralAndError(0, -1, e_reco[0], "") << endl;
   cout << "2e2mu " << h_2e2mu_reco->IntegralAndError(4, 100, e_reco[2], "")/h_2e2mu_reco->IntegralAndError(0, -1, e_reco[0], "") << endl;

}
void doComparisonPlots(const std::string inputFileMC1,const std::string inputFileMC2,const std::string inputFileMC3,const std::string inputFileMC4, const std::string variable, Double_t xmin = 0, Double_t xmax = -1)
{
  TFile* inputFile1;
  TFile* inputFile2;
  TFile* inputFile3;
  TFile* inputFile4;
      
  inputFile1 =  TFile::Open( inputFileMC1.c_str() );
  inputFile2 =  TFile::Open( inputFileMC2.c_str() );
  inputFile3 =  TFile::Open( inputFileMC3.c_str() );
  inputFile4 =  TFile::Open( inputFileMC4.c_str() );
 
  TH1F* h1 = (TH1F*)inputFile1->Get(TString("h1_")+TString(variable));
  TH1F* h2= (TH1F*)inputFile2->Get(TString("h1_")+TString(variable));
  TH1F* h3= (TH1F*)inputFile3->Get(TString("h1_")+TString(variable));
  TH1F* h4 = (TH1F*)inputFile4->Get(TString("h1_")+TString(variable));
 
  Float_t max1 = h1->GetMaximum();
  Float_t max2 = h2->GetMaximum();
  Float_t max3 = h3->GetMaximum();
  Float_t max4 = h4->GetMaximum();

  std::vector<float> max_values = {max1,max2,max3,max4}; 
  std::sort(max_values.begin(), max_values.end());

  TCanvas *c = new TCanvas("c", "c");
  //TH1*h = c->DrawFrame(xmin,ymin,xmax,ymax);
  //c->cd();
     
  TLegend *legend = new TLegend(.45,.25,.75,.45);
  legend->SetBorderSize(1);
  legend->SetTextSize(0.03 );
  h1->Draw("hist");
  h1->SetLineColor(kBlack);
  
  h2->Draw("hist SAME");
  h2->SetLineColor(kBlue);
  h3->Draw("hist SAME");
  h3->SetLineColor(kRed);
  h4->Draw("hist SAME");
  h4->SetLineColor(kMagenta);
  h1->Draw("hist SAME");
  h1->GetXaxis()->SetTitle(variable.c_str()); 
  h1->GetYaxis()->SetTitle(TString("Events")); 
  h1->GetYaxis()->SetTitleOffset(1.4);  
  h1->SetTitle("");

  Int_t bin_xmin = h1->GetXaxis()->FindBin(xmin);
  Int_t bin_xmax =  h1->GetXaxis()->FindBin(xmax);
  h1->GetXaxis()->SetRange(bin_xmin,bin_xmax);
  h1->SetMaximum(max_values[3]*1.05);

  cout << bin_xmin << " " << bin_xmax << " " << bin_xmax - bin_xmin << endl;
  //c->Modified();
  TString file1 = TString(inputFileMC1.c_str()).ReplaceAll(".root","");
  TString file2 = TString(inputFileMC2.c_str()).ReplaceAll(".root","");
  TString file3 = TString(inputFileMC3.c_str()).ReplaceAll(".root","");
  TString file4 = TString(inputFileMC4.c_str()).ReplaceAll(".root","");

  legend->AddEntry(h1,file1,"le");
  legend->AddEntry(h2,file2,"le");
  legend->AddEntry(h3,file3,"le");
  legend->AddEntry(h4,file4,"le");
  legend->Draw("SAME");

  //c->Print(TString(file1+TString("_")+TString(variable)+TString(".png"))); 
  //c->Print(TString(file1+TString("_")+TString(variable)+TString(".pdf"))); 
  //c->Print(TString(file1+TString("_")+TString(variable)+TString(".root"))); 

  // for(int k =0; k<50 ;k++){
  //   cout << k << " " << h1->GetBinContent(k)<< endl;
  // }



  Double_t bins[5]= { 0.0001, 0.001, 0.01, 0.1 ,1.};
  Int_t totBins = bin_xmax - bin_xmin;
  TH2F *h_2D = new TH2F("h2","",totBins,xmin,xmax,3,bins);
 
  Float_t entries;
  
  for (Int_t j = 1; j < 5; j++) {
    for (Int_t i = 1; i < totBins+1; i++) {
      
      entries = 0;
          
      if(j == 1) entries = h1->GetBinContent(bin_xmin+i);
      else if(j == 2) entries = h2->GetBinContent(bin_xmin+i);
      else if(j == 3) entries = h3->GetBinContent(bin_xmin+i);
      else if(j == 4) entries = h4->GetBinContent(bin_xmin+i);

      h_2D->SetBinContent(i,j,entries);
      //cout << "x = " << h_2D->GetXaxis()->GetBinCenter(i) <<" y = " << h_2D->GetYaxis()->GetBinCenter(j) << " entries = " << h_2D->GetBinContent(i,j) << " " << entries << endl;
    }
  }
  TCanvas *c2 = new TCanvas("c2", "c2");
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  h_2D->GetXaxis()->SetRange(bin_xmin,bin_xmax);
  h_2D->GetYaxis()->SetRange(0,5);
  h_2D->Draw("colz2"); 
  h_2D->GetXaxis()->SetTitle(variable.c_str()); 
  h_2D->GetYaxis()->SetTitle(TString("Lep-pho #DeltaR cut")); 

  c2->Print(TString(file1+TString("_")+TString(variable)+TString("_2D.png"))); 
  c2->Print(TString(file1+TString("_")+TString(variable)+TString("_2D.pdf"))); 
  c2->Print(TString(file1+TString("_")+TString(variable)+TString("_2D.root"))); 
}
/*
void doPlot(const std::string outputFile)
{
  

}
*/


void ZZgammagamma_analysis() 
{

  double lumi = 140; // full Run2 Lumi

  string inputFilePath = "/eos/user/a/acappati/samples_4lX/190626/";
  string inputFileName[] = {"ggH125",
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

  string outputFilePath = "histos_4lgammagamma";
  gSystem->Exec(("mkdir -p "+outputFilePath).c_str()); // create output dir
  

  //call function
  for(UInt_t i=0; i<nInputFiles; i++){
    cout<<"Processing sample "<<inputFileName[i]<<" ... "<<endl;
    doHisto(Form("%s%s%s",inputFilePath.c_str(),inputFileName[i].c_str(),"/ZZXAnalysis.root"), Form("%s%s%s%s",outputFilePath.c_str(), "/histos_", (inputFileName[i]).c_str(), ".root"), lumi);
  }
  
  //  doPlot(ouputFile);
}
