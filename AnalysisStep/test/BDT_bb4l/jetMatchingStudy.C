// ******************************************
// use: root -l -b -q jetMatchingStudy.C++
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

using namespace std;



void jetMatchingStudy(){

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

  //  vector<Float_t> *GENjetParentID   = 0;
  vector<Float_t> *prunedGenPartID   = 0;
  vector<Float_t> *prunedGenPartPt   = 0;
  vector<Float_t> *prunedGenPartEta   = 0;
  vector<Float_t> *prunedGenPartPhi   = 0;
  vector<Float_t> *prunedGenPartMass   = 0;
  vector<Float_t> *prunedGenMotherID   = 0;

  // counters
  int count2match_Method1 = 0;
  int count1match_Method1 = 0;
  int count0match_Method1 = 0;
  int countTOTmatch_Method1 = 0;

  int count2match_Method2 = 0;
  int count1match_Method2 = 0;
  int count0match_Method2 = 0;
  int countTOTmatch_Method2 = 0;

  int count2match_Method3 = 0;
  int count1match_Method3 = 0;
  int count0match_Method3 = 0;
  int countTOTmatch_Method3 = 0;

  // yields
  float yield_Method1_ = 0.;
  float yield_Method2_ = 0.;
  float yield_Method3_ = 0.;


  // histos
  TH1F* h_matches_Method1 = new TH1F("h_matches_Method1","reco-GEN jet matches;# reco jet with GEN match;#", 3,-0.5,2.5);
  TH1F* h_matches_Method2 = new TH1F("h_matches_Method2","reco-GEN jet matches;# reco jet with GEN match;#", 3,-0.5,2.5);
  TH1F* h_matches_Method3 = new TH1F("h_matches_Method3","reco-GEN jet matches;# reco jet with GEN match;#", 3,-0.5,2.5);

  TH1F* h_matches_Method1_eff = new TH1F("h_matches_Method1_eff","reco-GEN jet matches;# reco jet with GEN match;eff", 3,-0.5,2.5);
  TH1F* h_matches_Method2_eff = new TH1F("h_matches_Method2_eff","reco-GEN jet matches;# reco jet with GEN match;eff", 3,-0.5,2.5);
  TH1F* h_matches_Method3_eff = new TH1F("h_matches_Method3_eff","reco-GEN jet matches;# reco jet with GEN match;eff", 3,-0.5,2.5);

  TH2F* h2_matches_Method1onx_Method2ony = new TH2F("h2_matches_Method1onx_Method2ony",";# reco jet with GEN match (2high pT);# reco jet with GEN match (2high btag)", 3,-0.5,2.5,3,-0.5,2.5);
  TH2F* h2_matches_Method1onx_Method3ony = new TH2F("h2_matches_Method1onx_Method3ony",";# reco jet with GEN match (2high pT);# reco jet with GEN match (1high pT+1high btag)", 3,-0.5,2.5,3,-0.5,2.5);
  TH2F* h2_matches_Method3onx_Method2ony = new TH2F("h2_matches_Method3onx_Method2ony",";# reco jet with GEN match (1high pT+1high btag);# reco jet with GEN match (2high btag)", 3,-0.5,2.5,3,-0.5,2.5);


  TH1F* h_bbmass_Method1 = new TH1F("h_bbmass_Method1", ";bbMass;events", 50, 0., 500.); h_bbmass_Method1->Sumw2(true);
  TH1F* h_bbmass_Method2 = new TH1F("h_bbmass_Method2", ";bbMass;events", 50, 0., 500.); h_bbmass_Method2->Sumw2(true);
  TH1F* h_bbmass_Method3 = new TH1F("h_bbmass_Method3", ";bbMass;events", 50, 0., 500.); h_bbmass_Method3->Sumw2(true);
  TH1F* h_bbmass_GEN     = new TH1F("h_bbmass_GEN",     ";bbMass;events", 100, 0., 300.);
  
  TH1F* h_weights1 = new TH1F("h_weights1","",100,-0.0001,0.0008);
  TH1F* h_weights2 = new TH1F("h_weights2","",100,-0.0001,0.0008);
  TH1F* h_weights3 = new TH1F("h_weights3","",100,-0.0001,0.0008);
  
  TString inFile = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/HH4lbb_Angela/ZZXAnalysis.root";
  //  TString inFile = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/HH4lbb_Roberto_MadGraphLO/ZZXAnalysis.root";
  //  TString inFile = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/HH4lbb_Roberto_PowhegNLO/ZZXAnalysis.root";
  cout<<"reading file "<<inFile<<" ..."<<endl;
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
  inputTree->SetBranchAddress("xsec", &xsec);

  inputTree->SetBranchAddress("ZZsel", &ZZsel);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("ZZMass", &ZZMass); 
  
  inputTree->SetBranchAddress("JetPt", &JetPt);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetMass", &JetMass);
  inputTree->SetBranchAddress("JetPhi",  &JetPhi);
  inputTree->SetBranchAddress("JetBTagger",  &JetBTagger);

  //  inputTree->SetBranchAddress("GENjetParentID",  &GENjetParentID);
  inputTree->SetBranchAddress("prunedGenPartEta", &prunedGenPartEta );
  inputTree->SetBranchAddress("prunedGenPartPhi", &prunedGenPartPhi );
  inputTree->SetBranchAddress("prunedGenPartPt", &prunedGenPartPt );
  inputTree->SetBranchAddress("prunedGenPartMass", &prunedGenPartMass );
  inputTree->SetBranchAddress("prunedGenPartID", &prunedGenPartID );
  inputTree->SetBranchAddress("prunedGenMotherID", &prunedGenMotherID );




  Long64_t entries = inputTree->GetEntries();
  std::cout<<"Processing file: "<< inFile << "\nNumber of entries: " << entries << endl;

  for (Long64_t entry = 0; entry < entries; entry++)
  {  

    inputTree->GetEntry(entry);

    if(entry==10) cout<<"xsec  " <<xsec<<endl;

    // trigger check
    if( !(ZZsel >= 0) ) continue;


    // check if only 4 lepton
    if(LepEta->size() != 4){
	cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
    }


    // mass cut: signal region
    if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV


    // fill eventweight
    Double_t eventWeight = partialSampleWeight * xsec * overallEventWeight ;  //kfactor e l1prefiring non ci sono, tanto uso solo HH sample
    



    // Delta R for Gen jet matching
    double DELTAR = 0.4;    
 
    // -------------------------------------------
    vector<TLorentzVector> JetVec_; // TLorentz vector with all Jets per Event
    vector<float> JetpT_; // vector with the pT per each Jet of the Event
    vector<float> JetBinfo_; // vector with the Btagger per each Jet of the Event

    for (UInt_t j = 0; j < JetPt->size(); j++)
    {
      if ( (fabs ( JetEta->at(j) ) > 2.4) || (JetPt->at(j) < 20 ) ) continue; // pt cut 20GeV from ntuplizer 
	  
      TLorentzVector temp;
      temp.SetPtEtaPhiM(JetPt->at(j), JetEta->at(j), JetPhi->at(j), JetMass->at(j));
      JetVec_.push_back(temp);
      JetpT_.push_back(JetPt->at(j));
      JetBinfo_.push_back(JetBTagger->at(j));
    }



    // ------------------------------------------
    // at least 2 jets in the acceptance
    if (JetVec_.size() < 2) continue;




    // --- jet selection METHOD1: 2 high pT jets
    int d1_Method1 = -999; // position of higest pT jet
    int d2_Method1 = -999; // position of second-highest pT jet


    bool bool_match1_Method1 = false;
    bool bool_match2_Method1 = false;

    // count yield
    yield_Method1_ += eventWeight;

    // count number of events with at least 2 jets
    countTOTmatch_Method1++;

    // find the position of the highest pT jet
    d1_Method1 = distance( JetpT_.begin(), max_element(JetpT_.begin(), JetpT_.end()));

    // find the position of the second-highest pT jet
    float max_Method1 = -9999.;
    for(UInt_t i=0; i<JetpT_.size(); i++){
      if(i == d1_Method1) continue;
      float temp = JetpT_.at(i);
      if(temp > max_Method1){
        max_Method1 = temp;
        d2_Method1 = i;
      }
    }

    //      cout<<d1_Method1<<" "<<d2_Method1<<endl;      
   
    // find match between 1st reco jets and GEN jets
    for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
      if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
   
        float deltaPhi1 = JetVec_.at(d1_Method1).Phi() - prunedGenPartPhi->at(pr);
        if(fabs(deltaPhi1) > acos(-1)){ deltaPhi1 = (2*acos(-1)) - fabs(deltaPhi1); }
        float deltaEta1 = JetVec_.at(d1_Method1).Eta() - prunedGenPartEta->at(pr);

        float deltaR1 = sqrt(deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1);

        if(deltaR1 < DELTAR){ bool_match1_Method1 = true ;}
      }
    }

    // find match between 2nd reco jets and GEN jets
    for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
      if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
   
        float deltaPhi2 = JetVec_.at(d2_Method1).Phi() - prunedGenPartPhi->at(pr);
        if(fabs(deltaPhi2) > acos(-1)){ deltaPhi2 = (2*acos(-1)) - fabs(deltaPhi2); }
        float deltaEta2 = JetVec_.at(d2_Method1).Eta() - prunedGenPartEta->at(pr);

        float deltaR2 = sqrt(deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2);

        if(deltaR2 < 0.4){ bool_match2_Method1 = true ;}
      }
    }

    //      cout<<bool_match1_Method1<<" "<<bool_match2_Method1<<endl;
    //      cout<<count0match_Method1<<" "<<count1match_Method1<<" "<<count2match_Method1<<" "<<countTOTmatch_Method1<<endl;


    TLorentzVector Hbb_vec_Method1 = JetVec_.at(d1_Method1) + JetVec_.at(d2_Method1);
    h_bbmass_Method1->Fill(Hbb_vec_Method1.M(), eventWeight);
    h_weights1->Fill(eventWeight);

    // --- end jet selection METHOD1: 2 high pT jets
    // ----------------------------------------------




    // -------------------------------------------
    // --- jet selection METHOD2: 2 high btagger jets
    int d1_Method2 = -999; // position of higest btagger jet
    int d2_Method2 = -999; // position of second-highest btagger jet
    
    bool bool_match1_Method2 = false;
    bool bool_match2_Method2 = false;

    // count yield
    yield_Method2_ += eventWeight;

    // count number of events with at least 2 jets
    countTOTmatch_Method2++;

    // find the position of the highest btagger jet
    d1_Method2 = distance( JetBinfo_.begin(), max_element(JetBinfo_.begin(), JetBinfo_.end()));

    // find the position of the second-highest btagger jet
    float max_Method2 = -9999.;
    for(UInt_t i=0; i<JetBinfo_.size(); i++){
      if(i == d1_Method2) continue;
      float temp = JetBinfo_.at(i);
      if(temp > max_Method2){
        max_Method2 = temp;
        d2_Method2 = i;
      }
    }

    //      cout<<d1_Method2<<" "<<d2_Method2<<endl;      


    // find match between 1st reco jets and GEN jets
    for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
      if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
 
        float deltaPhi1 = JetVec_.at(d1_Method2).Phi() - prunedGenPartPhi->at(pr);
        if(fabs(deltaPhi1) > acos(-1)){ deltaPhi1 = (2*acos(-1)) - fabs(deltaPhi1); }
        float deltaEta1 = JetVec_.at(d1_Method2).Eta() - prunedGenPartEta->at(pr);

        float deltaR1 = sqrt(deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1);

        if(deltaR1 < DELTAR){ bool_match1_Method2 = true ;}
      }
    }

    // find match between 2nd reco jets and GEN jets
    for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
      if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
 
        float deltaPhi2 = JetVec_.at(d2_Method2).Phi() - prunedGenPartPhi->at(pr);
        if(fabs(deltaPhi2) > acos(-1)){ deltaPhi2 = (2*acos(-1)) - fabs(deltaPhi2); }
        float deltaEta2 = JetVec_.at(d2_Method2).Eta() - prunedGenPartEta->at(pr);

        float deltaR2 = sqrt(deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2);

        if(deltaR2 < 0.4){ bool_match2_Method2 = true ;}
      }
    }

    //      cout<<bool_match1_Method2<<" "<<bool_match2_Method2<<endl;
    //      cout<<count0match_Method2<<" "<<count1match_Method2<<" "<<count2match_Method2<<" "<<countTOTmatch_Method2<<endl;

    TLorentzVector Hbb_vec_Method2 = JetVec_.at(d1_Method2) + JetVec_.at(d2_Method2);
    h_bbmass_Method2->Fill(Hbb_vec_Method2.M(), eventWeight);
    h_weights2->Fill(eventWeight);
    //    } 
    // --- end jet selection METHOD2: 2 high btagger jets
    // ----------------------------------------------



    // -------------------------------------------
    // --- jet selection METHOD3: higher btagger jet + higher pT jet
    int d1_Method3 = -999; // position of higest btagger jet
    int d2_Method3 = -999; // position of second-highest btagger jet

    bool bool_match1_Method3 = false;
    bool bool_match2_Method3 = false;

    // count yield
    yield_Method3_ += eventWeight;

    // count number of events with at least 2 jets
    countTOTmatch_Method3++;

    // find the position of the highest btagger jet
    d1_Method3 = distance( JetBinfo_.begin(), max_element(JetBinfo_.begin(), JetBinfo_.end()));

    // find the position of the highest pT jet
    float max_Method3 = -9999.;
    for(UInt_t i=0; i<JetpT_.size(); i++){
      if(i == d1_Method3) continue;
      float temp = JetpT_.at(i);
      if(temp > max_Method3){
        max_Method3 = temp;
        d2_Method3 = i;
      }
    }

    //      cout<<d1_Method3<<" "<<d2_Method3<<endl;      
  

    // find match between 1st reco jets and GEN jets
    for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
      if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
 
        float deltaPhi1 = JetVec_.at(d1_Method3).Phi() - prunedGenPartPhi->at(pr);
        if(fabs(deltaPhi1) > acos(-1)){ deltaPhi1 = (2*acos(-1)) - fabs(deltaPhi1); }
        float deltaEta1 = JetVec_.at(d1_Method3).Eta() - prunedGenPartEta->at(pr);

        float deltaR1 = sqrt(deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1);

        if(deltaR1 < DELTAR){ bool_match1_Method3 = true ;}
      }
    }

    // find match between 2nd reco jets and GEN jets
    for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
      if (fabs(prunedGenPartID->at(pr)) == 5 && prunedGenMotherID->at(pr) == 25){     // check that GEN particle is b quark and mother of b quark is Higgs
 
        float deltaPhi2 = JetVec_.at(d2_Method3).Phi() - prunedGenPartPhi->at(pr);
        if(fabs(deltaPhi2) > acos(-1)){ deltaPhi2 = (2*acos(-1)) - fabs(deltaPhi2); }
        float deltaEta2 = JetVec_.at(d2_Method3).Eta() - prunedGenPartEta->at(pr);

        float deltaR2 = sqrt(deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2);

        if(deltaR2 < 0.4){ bool_match2_Method3 = true ;}
      }
    }

    //      cout<<bool_match1_Method3<<" "<<bool_match2_Method3<<endl;
    //      cout<<count0match_Method3<<" "<<count1match_Method3<<" "<<count2match_Method3<<" "<<countTOTmatch_Method3<<endl;

    TLorentzVector Hbb_vec_Method3 = JetVec_.at(d1_Method3) + JetVec_.at(d2_Method3);
    h_bbmass_Method3->Fill(Hbb_vec_Method3.M(), eventWeight);
    h_weights3->Fill(eventWeight);
    // --- end jet selection METHOD3: higher btagger jet + higher pT jet
    // ----------------------------------------------



    //CONDIZIONI PER FILL HISTOS
    Float_t x_Method1 = -999.;
    Float_t x_Method2 = -999.;
    Float_t x_Method3 = -999.;

    //method1
    if(bool_match1_Method1 && bool_match2_Method1){
      count2match_Method1++; 
      x_Method1 = 2.;
    }
    else if( (bool_match1_Method1 && !bool_match2_Method1) || (!bool_match1_Method1 && bool_match2_Method1) ){
      count1match_Method1++;
      x_Method1 = 1.;
    }
    else{
      count0match_Method1++;
      x_Method1 = 0.;
    }

    //method2
    if(bool_match1_Method2 && bool_match2_Method2){
      count2match_Method2++; 
      x_Method2 = 2.;
    }
    else if( (bool_match1_Method2 && !bool_match2_Method2) || (!bool_match1_Method2 && bool_match2_Method2) ){
      count1match_Method2++;
      x_Method2 = 1.;
    }
    else{
      count0match_Method2++;
      x_Method2 = 0.;
    }

    //method3
    if(bool_match1_Method3 && bool_match2_Method3){
      count2match_Method3++; 
      x_Method3 = 2.;
    }
    else if( (bool_match1_Method3 && !bool_match2_Method3) || (!bool_match1_Method3 && bool_match2_Method3) ){
      count1match_Method3++;
      x_Method3 = 1.;
    }
    else{
      count0match_Method3++;
      x_Method3 = 0.;
    }


    //FILL HISTOS
    h_matches_Method1->Fill(x_Method1);
    h_matches_Method2->Fill(x_Method2);
    h_matches_Method3->Fill(x_Method3);
    //FILL HISTOS 2D
    // histo with Method1 on x axis and Method2 on y axis
    h2_matches_Method1onx_Method2ony->Fill(x_Method1,x_Method2);
    // histo with Method1 on x axis and Method3 on y axis
    h2_matches_Method1onx_Method3ony->Fill(x_Method1,x_Method3);
    // histo with Method3 on x axis and Method2 on y axis
    h2_matches_Method3onx_Method2ony->Fill(x_Method3,x_Method2);



  }// end loop over tree entries

  




  // count matching efficiencies 

  // method 1
  cout<<(float)count0match_Method1/(float)countTOTmatch_Method1<<" "<<(float)count1match_Method1/(float)countTOTmatch_Method1<<" "<<(float)count2match_Method1/(float)countTOTmatch_Method1<<" "<<(float)countTOTmatch_Method1/(float)countTOTmatch_Method1<<endl;
  cout<<(float)count0match_Method1<<" "<<(float)count1match_Method1<<" "<<(float)count2match_Method1<<" "<<(float)countTOTmatch_Method1<<endl;


  h_matches_Method1_eff->SetBinContent(1, (float)count0match_Method1/(float)countTOTmatch_Method1);
  h_matches_Method1_eff->SetBinContent(2, (float)count1match_Method1/(float)countTOTmatch_Method1);
  h_matches_Method1_eff->SetBinContent(3, (float)count2match_Method1/(float)countTOTmatch_Method1);

  // method 2
  cout<<(float)count0match_Method2/(float)countTOTmatch_Method2<<" "<<(float)count1match_Method2/(float)countTOTmatch_Method2<<" "<<(float)count2match_Method2/(float)countTOTmatch_Method2<<" "<<(float)countTOTmatch_Method2/(float)countTOTmatch_Method2<<endl;
  cout<<(float)count0match_Method2<<" "<<(float)count1match_Method2<<" "<<(float)count2match_Method2<<" "<<(float)countTOTmatch_Method2<<endl;

  h_matches_Method2_eff->SetBinContent(1, (float)count0match_Method2/(float)countTOTmatch_Method2);
  h_matches_Method2_eff->SetBinContent(2, (float)count1match_Method2/(float)countTOTmatch_Method2);
  h_matches_Method2_eff->SetBinContent(3, (float)count2match_Method2/(float)countTOTmatch_Method2);

  // method 3
  cout<<(float)count0match_Method3/(float)countTOTmatch_Method3<<" "<<(float)count1match_Method3/(float)countTOTmatch_Method3<<" "<<(float)count2match_Method3/(float)countTOTmatch_Method3<<" "<<(float)countTOTmatch_Method3/(float)countTOTmatch_Method3<<endl;
  cout<<(float)count0match_Method3<<" "<<(float)count1match_Method3<<" "<<(float)count2match_Method3<<" "<<(float)countTOTmatch_Method3<<endl;

  h_matches_Method3_eff->SetBinContent(1, (float)count0match_Method3/(float)countTOTmatch_Method3);
  h_matches_Method3_eff->SetBinContent(2, (float)count1match_Method3/(float)countTOTmatch_Method3);
  h_matches_Method3_eff->SetBinContent(3, (float)count2match_Method3/(float)countTOTmatch_Method3);


  


  // draw histos
  gStyle->SetOptStat(0);

  // jet matching
  TCanvas* c_matches = new TCanvas();
  c_matches->cd();

  h_matches_Method2->SetLineColor(kRed);
  h_matches_Method2->Draw("hist");

  h_matches_Method1->SetLineColor(kBlue);
  h_matches_Method1->Draw("hist same");

  h_matches_Method3->SetLineColor(kGreen+2);
  h_matches_Method3->Draw("hist same");


  TLegend* leg = new TLegend(0.65,0.21,0.94,0.39);
  leg->AddEntry(h_matches_Method1, "2 highest pT", "l");
  leg->AddEntry(h_matches_Method2, "2 highest btagger", "l");
  leg->AddEntry(h_matches_Method3, "highest btagger +", "l");
  leg->AddEntry((TObject*)0,"highest pT","");
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kBlack);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->Draw();

  c_matches->Update();

  c_matches->SaveAs("jetmatches.png");



  // jet matching eff
  TCanvas* c_matches_eff = new TCanvas();
  c_matches_eff->cd();

  h_matches_Method2_eff->SetLineColor(kRed);
  h_matches_Method2_eff->Draw("hist");

  h_matches_Method1_eff->SetLineColor(kBlue);
  h_matches_Method1_eff->Draw("hist same");

  h_matches_Method3_eff->SetLineColor(kGreen+2);
  h_matches_Method3_eff->Draw("hist same");


  TLegend* leg_eff = new TLegend(0.65,0.21,0.94,0.39);
  leg_eff->AddEntry(h_matches_Method1_eff, "2 highest pT", "l");
  leg_eff->AddEntry(h_matches_Method2_eff, "2 highest btagger", "l");
  leg_eff->AddEntry(h_matches_Method3_eff, "highest btagger +", "l");
  leg_eff->AddEntry((TObject*)0,"highest pT","");
  leg_eff->SetFillColor(kWhite);
  leg_eff->SetLineColor(kBlack);
  leg_eff->SetTextFont(43);
  leg_eff->SetTextSize(20);
  leg_eff->Draw();

  c_matches_eff->Update();

  c_matches_eff->SaveAs("jetmatches_eff.png");





  //get yields
  cout<<"yields "<<yield_Method1_<<" "<<yield_Method2_<<" "<<yield_Method3_<<endl;
  cout<<"mbb integral "<<h_bbmass_Method1->Integral(0,-1)<<" "<<h_bbmass_Method2->Integral(0,-1)<<" "<<h_bbmass_Method3->Integral(0,-1)<<endl;
  cout<<"entries "<<h_bbmass_Method1->GetEntries()<<" "<<h_bbmass_Method2->GetEntries()<<" "<<h_bbmass_Method3->GetEntries()<<endl;

  // bbmass
  TCanvas* c_Massbb = new TCanvas();
  c_Massbb->cd();

  h_bbmass_Method2->SetLineColor(kRed);
  h_bbmass_Method2->Draw("hist");

  h_bbmass_Method3->SetLineColor(kGreen+2);
  h_bbmass_Method3->Draw("hist same");

  h_bbmass_Method1->SetLineColor(kBlue);
  h_bbmass_Method1->Draw("hist same");


  TLegend* leg2 = new TLegend(0.65,0.21,0.94,0.39);
  leg2->AddEntry(h_bbmass_Method1, "2 highest pT", "l");
  leg2->AddEntry(h_bbmass_Method2, "2 highest btagger", "l");
  leg2->AddEntry(h_bbmass_Method3, "highest btagger +", "l");
  leg2->AddEntry((TObject*)0,      "highest pT","");
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kBlack);
  leg2->SetTextFont(43);
  leg2->SetTextSize(20);
  leg2->Draw();

  c_Massbb->Update();

  c_Massbb->SaveAs("bbMass.png");





  //weights
  TCanvas* c_weights = new TCanvas();
  c_weights->cd();
  h_weights1->SetLineColor(kBlue);
  h_weights1->Draw("hist");
  h_weights2->SetLineColor(kRed);
  h_weights2->Draw("hist");
  h_weights3->SetLineColor(kGreen+3);
  h_weights3->Draw("hist");

  c_weights->SaveAs("weights.png");



  //2D HISTOS
  TCanvas* c_h2_jetmatching_1onx_2ony = new TCanvas();
  c_h2_jetmatching_1onx_2ony->cd();
  h2_matches_Method1onx_Method2ony->Draw("COLZ");
  c_h2_jetmatching_1onx_2ony->SaveAs("jetmatching_2D_1onx_2ony.png");

  TCanvas* c_h2_jetmatching_1onx_3ony = new TCanvas();
  c_h2_jetmatching_1onx_3ony->cd();
  h2_matches_Method1onx_Method3ony->Draw("COLZ");
  c_h2_jetmatching_1onx_3ony->SaveAs("jetmatching_2D_1onx_3ony.png");

  TCanvas* c_h2_jetmatching_3onx_2ony = new TCanvas();
  c_h2_jetmatching_3onx_2ony->cd();
  h2_matches_Method3onx_Method2ony->Draw("COLZ");
  c_h2_jetmatching_3onx_2ony->SaveAs("jetmatching_2D_3onx_2ony.png");



}
