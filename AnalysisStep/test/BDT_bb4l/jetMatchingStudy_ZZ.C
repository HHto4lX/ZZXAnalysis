// ******************************************
// use: root -l -b -q jetMatchingStudy_ZZ.C++
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

// ---------------------------------------
// --- sampleType
TString sampleType = "GGH";
//TString sampleType = "GGZZ";
// ---------------------------------------
// --- massRegion
TString massRegion = "SR";
//TString massRegion = "fullmass";
// ---------------------------------------


void jetMatchingStudy_ZZ(){

  float lumi = 59.7; //fb-1


  // -------------------------------------------
  // datasets
  TString datasets[] = {
    "ggH125",  // ggH
    // "ggTo4e_Contin_MCFM701",
    // "ggTo4mu_Contin_MCFM701",
    // "ggTo4tau_Contin_MCFM701",
    // "ggTo2e2mu_Contin_MCFM701",
    // "ggTo2e2tau_Contin_MCFM701",
    // "ggTo2mu2tau_Contin_MCFM701",
   };

  static size_t nDatasets = sizeof(datasets)/sizeof(datasets[0]);
  cout<< "number of input files: " << nDatasets<<endl;
  // -------------------------------------------    



  TFile* inputFile[nDatasets]; 
  TTree* inputTree[nDatasets]; 
  TH1F* hCounters[nDatasets]; 
  int NGenEvt[nDatasets];
  Float_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];

  Float_t genBR;
  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t overallEventWeight;
  Float_t xsec;

  Float_t KFactor_QCD_ggZZ_Nominal;

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
  int count0match_Method2 = 0;
  int count1match_Method2 = 0;
  int count2match_Method2 = 0;
  int countTOTmatch_Method2 = 0;


  // yields
  float yield_Method2_ = 0.;



  // histos
  TH1F* h_matches_Method2 = new TH1F("h_matches_Method2","reco-GEN jet matches;# reco jet with GEN match;eff", 3,-0.5,2.5);
  TH1F* h_bbmass_Method2 = new TH1F("h_bbmass_Method2", ";bbMass;events", 50, 0., 500.); h_bbmass_Method2->Sumw2(true);  
  TH1F* h_j1btag_Method2 = new TH1F("h_j1btag_Method2", ";j1 DeepCSV;events", 25, 0., 1.); h_j1btag_Method2->Sumw2(true);  
  TH1F* h_j2btag_Method2 = new TH1F("h_j2btag_Method2", ";j2 DeepCSV;events", 25, 0., 1.); h_j2btag_Method2->Sumw2(true);  
  TH1F* h_weights2 = new TH1F("h_weights2","",100,-0.0001,0.0008);
  



  for(int d=0; d<nDatasets; d++){  

    TString inputFileName = "/eos/user/a/acappati/samples_HH4lbb/samples_2018/" + datasets[d] + "/ZZXAnalysis.root";
    cout<<"Opening file "<<inputFileName<<" ..."<<endl;
    inputFile[d] = TFile::Open(inputFileName);

    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
    NGenEvt[d] = (Float_t)hCounters[d]->GetBinContent(1);
    gen_sumWeights[d] = (Float_t)hCounters[d]->GetBinContent(40);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d];
    inputTree[d] = (TTree*)inputFile[d]->Get("ZZTree/candTree");


    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("xsec", &xsec);

    inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass); 
    
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("JetPhi",  &JetPhi);
    inputTree[d]->SetBranchAddress("JetBTagger",  &JetBTagger);

    //  inputTree->SetBranchAddress("GENjetParentID",  &GENjetParentID);
    inputTree[d]->SetBranchAddress("prunedGenPartEta", &prunedGenPartEta );
    inputTree[d]->SetBranchAddress("prunedGenPartPhi", &prunedGenPartPhi );
    inputTree[d]->SetBranchAddress("prunedGenPartPt", &prunedGenPartPt );
    inputTree[d]->SetBranchAddress("prunedGenPartMass", &prunedGenPartMass );
    inputTree[d]->SetBranchAddress("prunedGenPartID", &prunedGenPartID );
    inputTree[d]->SetBranchAddress("prunedGenMotherID", &prunedGenMotherID );

    if(sampleType == "GGZZ"){
      inputTree[d]->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
    }




    Long64_t entries = inputTree[d]->GetEntries();
    cout<<"Processing file: "<< datasets[d] << " (" << entries <<" entries) ..."<< endl;

    for (Long64_t entry = 0; entry < entries; entry++)
    {   

      inputTree[d]->GetEntry(entry);

      if(entry==10) cout<<"xsec  " <<xsec<<endl;

      // trigger check
      if( !(ZZsel >= 0) ) continue;


      // check if only 4 lepton
      if(LepEta->size() != 4){
        cerr << "error in event " << nRun << ":" << nLumi << ":" << nEvent << "; stored leptons: "<< LepEta->size() << endl;
	continue;
      }



      if(massRegion == "SR"){
	// mass cut: signal region
	if(ZZMass < 115 || ZZMass > 135) continue; // 115 < ZZMass < 135 GeV
      }
      // else if(massRegion == "sidebands"){
      // 	// mass cut: side bands
      // 	if(ZZMass >= 115 && ZZMass <= 135) continue; // ZZMass < 115 or ZZMass > 135 GeV
      // }
      else if(massRegion == "fullmass"){
	// full mass range
      }
      else cerr<<"wrong mass region!"<<endl;



      // kfactors
      Float_t kfactor = 1.;
      if(sampleType == "GGZZ"){
        kfactor = KFactor_QCD_ggZZ_Nominal;
      }


      // fill eventweight
      Double_t eventWeight = partialSampleWeight[d] * xsec * kfactor * overallEventWeight;
    



      // Delta R for Gen jet matching
      double DELTAR = 0.4;    
 




      // -------------------------------------------
      // --- jet selection METHOD2: 2 high btagger jets
      vector<TLorentzVector> JetPair_Method2; // TLorentz vector with all Jets Pairs

      int d1_Method2 = -999; // position of higest btagger jet
      int d2_Method2 = -999; // position of second-highest btagger jet

      // at least 2 jets in the acceptance
      if (JetPt->size() >= 2){

        // count yield
        yield_Method2_ += eventWeight;

        // count number of events with at least 2 jets
        countTOTmatch_Method2++;

        // find the position of the highest btagger jet
        d1_Method2 = distance( JetBTagger->begin(), max_element(JetBTagger->begin(), JetBTagger->end()));

        // find the position of the second-highest btagger jet
        float max = -9999.;
        for(UInt_t i=0; i<JetBTagger->size(); i++){
          if(i == d1_Method2) continue;
          float temp = JetBTagger->at(i);
          if(temp > max){
            max = temp;
            d2_Method2 = i;
          }
        }


        // find match between 1st reco jets and GEN jets
        bool bool_match1_Method2 = false;
        for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
          if (fabs(prunedGenPartID->at(pr)) == 5 ){     // check that GEN particle is b quark and mother of b quark is Higgs
   
            float deltaPhi1 = JetPhi->at(d1_Method2) - prunedGenPartPhi->at(pr);
            if(fabs(deltaPhi1) > acos(-1)){ deltaPhi1 = (2*acos(-1)) - fabs(deltaPhi1); }
            float deltaEta1 = JetEta->at(d1_Method2) - prunedGenPartEta->at(pr);

            float deltaR1 = sqrt(deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1);

            if(deltaR1 < DELTAR){ bool_match1_Method2 = true ;}
          }
        }

        // find match between 2nd reco jets and GEN jets
        bool bool_match2_Method2 = false;
        for(UInt_t pr = 0; pr<prunedGenPartEta->size(); pr++){
          if (fabs(prunedGenPartID->at(pr)) == 5 ){     // check that GEN particle is b quark and mother of b quark is Higgs
   
            float deltaPhi2 = JetPhi->at(d2_Method2) - prunedGenPartPhi->at(pr);
            if(fabs(deltaPhi2) > acos(-1)){ deltaPhi2 = (2*acos(-1)) - fabs(deltaPhi2); }
            float deltaEta2 = JetEta->at(d2_Method2) - prunedGenPartEta->at(pr);

            float deltaR2 = sqrt(deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2);

            if(deltaR2 < 0.4){ bool_match2_Method2 = true ;}
          }
        }

      //      cout<<bool_match1_Method2<<" "<<bool_match2_Method2<<endl;

        if(bool_match1_Method2 && bool_match2_Method2) {count2match_Method2++;}
        else if( (bool_match1_Method2 && !bool_match2_Method2) || (!bool_match1_Method2 && bool_match2_Method2) ) {count1match_Method2++;}
        else {count0match_Method2++;}

        //      cout<<count0match_Method2<<" "<<count1match_Method2<<" "<<count2match_Method2<<" "<<countTOTmatch_Method2<<endl;

        // build 2 jets tlorentzvectors
        TLorentzVector tlzvec_j1_;
        tlzvec_j1_.SetPtEtaPhiM(JetPt->at(d1_Method2), JetEta->at(d1_Method2), JetPhi->at(d1_Method2), JetMass->at(d1_Method2));
        TLorentzVector tlzvec_j2_;
        tlzvec_j2_.SetPtEtaPhiM(JetPt->at(d2_Method2), JetEta->at(d2_Method2), JetPhi->at(d2_Method2), JetMass->at(d2_Method2));


        TLorentzVector Hbb_vec_Method2 = tlzvec_j1_ + tlzvec_j2_;
        h_bbmass_Method2->Fill(Hbb_vec_Method2.M(), eventWeight);
        h_weights2->Fill(eventWeight);
        h_j1btag_Method2->Fill(JetBTagger->at(d1_Method2), eventWeight);
        h_j2btag_Method2->Fill(JetBTagger->at(d2_Method2), eventWeight);
      } 
      // --- end jet selection METHOD2: 2 high btagger jets
      // ----------------------------------------------



    }// end loop over tree entries

  } //end loop over datasets

  




  // count matching efficiencies 


  // method 2
  cout<<(float)count0match_Method2/(float)countTOTmatch_Method2<<" "<<(float)count1match_Method2/(float)countTOTmatch_Method2<<" "<<(float)count2match_Method2/(float)countTOTmatch_Method2<<" "<<(float)countTOTmatch_Method2/(float)countTOTmatch_Method2<<endl;
  cout<<(float)count0match_Method2<<" "<<(float)count1match_Method2<<" "<<(float)count2match_Method2<<" "<<(float)countTOTmatch_Method2<<endl;

  h_matches_Method2->SetBinContent(1, (float)count0match_Method2/(float)countTOTmatch_Method2);
  h_matches_Method2->SetBinContent(2, (float)count1match_Method2/(float)countTOTmatch_Method2);
  h_matches_Method2->SetBinContent(3, (float)count2match_Method2/(float)countTOTmatch_Method2);

  


  // draw histos
  gStyle->SetOptStat(0);

  // jet matching
  TCanvas* c_matches = new TCanvas();
  c_matches->cd();

  h_matches_Method2->SetLineColor(kRed);
  h_matches_Method2->Draw("hist");



  // TLegend* leg = new TLegend(0.65,0.21,0.94,0.39);
  // leg->AddEntry(h_matches_Method2, "2 highest btagger", "l");
  // leg->SetFillColor(kWhite);
  // leg->SetLineColor(kBlack);
  // leg->SetTextFont(43);
  // leg->SetTextSize(20);
  // leg->Draw();

  c_matches->Update();

  c_matches->SaveAs("jetmatches_" + sampleType + "_" + massRegion + ".png");



  //get yields
  cout<<"yields "<<yield_Method2_<<endl;
  cout<<"mbb integral "<<h_bbmass_Method2->Integral(0,-1)<<endl;
  cout<<"entries "<<h_bbmass_Method2->GetEntries()<<endl;

  // bbmass
  TCanvas* c_Massbb = new TCanvas();
  c_Massbb->cd();

  h_bbmass_Method2->SetLineColor(kRed);
  h_bbmass_Method2->Draw("hist");


  // TLegend* leg2 = new TLegend(0.65,0.21,0.94,0.39);
  // leg2->AddEntry(h_bbmass_Method2, "2 highest btagger", "l");
  // leg2->SetFillColor(kWhite);
  // leg2->SetLineColor(kBlack);
  // leg2->SetTextFont(43);
  // leg2->SetTextSize(20);
  // leg2->Draw();

  c_Massbb->Update();

  c_Massbb->SaveAs("bbMass_" + sampleType + "_" + massRegion + ".png");

  // j1btagger
  TCanvas* c_j1btagger = new TCanvas();
  c_j1btagger->cd();
  h_j1btag_Method2->SetLineColor(kRed);
  h_j1btag_Method2->Draw("hist");
  c_j1btagger->SaveAs("btaggerj1_" + sampleType + "_" + massRegion + ".png");

  // j2btagger
  TCanvas* c_j2btagger = new TCanvas();
  c_j2btagger->cd();
  h_j2btag_Method2->SetLineColor(kRed);
  h_j2btag_Method2->Draw("hist");
  c_j2btagger->SaveAs("btaggerj2_" + sampleType + "_" + massRegion + ".png");

  //weights
  TCanvas* c_weights = new TCanvas();
  c_weights->cd();
  h_weights2->SetLineColor(kRed);
  h_weights2->Draw("hist");

  c_weights->SaveAs("weights_" + sampleType + "_" + massRegion + ".png");

}
