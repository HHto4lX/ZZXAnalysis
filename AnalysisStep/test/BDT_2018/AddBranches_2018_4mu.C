#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
//#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

//#include "/lustre/home/taliercio/SL7/CMSSW_9_3_2/src/NeuralNetwork_new/dataset_2016_2e2mu_massjetjet/weights/TMVAClassification_BDT_Grad_400_2.class.C"
//#include "/lustre/home/taliercio/SL7/CMSSW_9_3_2/src/NeuralNetwork_new/dataset_2016_2e2mu_massjetjet/weights/TMVAClassification_MLP_500_N-7.class.C"
//#include "/lustre/home/taliercio/CMSSW_9_3_2/src/dataset_ggH8_utilizzatobjet/weights/TMVAClassification_BDT.class.C"

//adding the BDTs
//#include "/lustre/home/taliercio/SL7/CMSSW_9_3_2/src/NeuralNetwork_new/dataset_new/weights/TMVAClassification_BDT_Grad_200_2.class.C"
//#include "/lustre/home/taliercio/SL7/CMSSW_9_3_2/src/NeuralNetwork_new/dataset_study_4mu_2016_4_new/weights/TMVAClassification_BDT_Grad_400_2_4.class.C"
//#include "/lustre/home/taliercio/SL7/CMSSW_9_3_2/src/NeuralNetwork_new/dataset_study_4mu_2016_10_new/weights/TMVAClassification_BDT_Grad_400_2_10.class.C"
//#include "/lustre/home/taliercio/SL7/CMSSW_9_3_2/src/NeuralNetwork_new/dataset_study_4mu_2016_11_new/weights/TMVAClassification_BDT_Grad_400_2_11.class.C"
#include "/4mu/TMVAClassification_BDT_Grad_400_2.class.C"


void AddBranches_2018(){
using namespace std;
//int main (int argc, char ** argv){
  //signal

  //TFile *signal_4mu = new TFile("/lustre/home/taliercio/CMSSW_9_3_2/src/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root","update");

  
  std::ifstream infile_data, infile_mc;
  infile_data.open("bkg_2018_4mu.txt");//bari_data_2e2mu_good.txt
  infile_mc.open("bkg_2018_4mu.txt");

  std::string inputfilename_data, inputfilename_mc;
  std::vector < TFile* > data, mc;
  /*
  while(std::getline(infile_mc,inputfilename_mc)){
    mc.push_back(new TFile(static_cast<TString>(inputfilename_mc)));
  }
  infile_mc.close();
  */

  while(std::getline(infile_data,inputfilename_data)){
    data.push_back(new TFile(static_cast<TString>(inputfilename_data), "update"));
  }
  infile_data.close();
  /*
  while(std::getline(infile_mc,inputfilename_mc)){
    mc.push_back(new TFile(static_cast<TString>(inputfilename_mc)));
  }
  infile_mc.close();
  */

  cout << "fuori dal for " << "\n";
  for(auto flag : data){
  flag->cd();  
  cout << "name " << flag->GetName() << "\n";
  TTree *treesignal_4mu = (TTree*)flag->Get("reducedTree");
  float f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_lept4_ptsignal_4mu, f_deltaphisignal_4mu, f_deltarsignal_4mu, f_METsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu, f_etajet1signal_4mu, f_etajet2signal_4mu, f_ptjet1signal_4mu, f_ptjet2signal_4mu, f_massjetjetsignal_4mu, jet1_pt,
    //signal_4mu variable
    weightsignal_4mu, masssignal_4mu, jetsignal_4mu, pt3signal_4mu;
  
  
  int Ngoodsignal_4mu, bjetsignal_4mu, Nbjetssignal_4mu;
  int cont_4mu =0, cont_4mu_all=0;
  float BDT_shape_fill, BDT_shape_fill_4_conf, BDT_shape_fill_10_conf, BDT_shape_fill_11_conf;
  treesignal_4mu->SetBranchAddress("f_lept1_pt",&f_lept1_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_lept2_pt",&f_lept2_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_lept3_pt",&f_lept3_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_lept4_pt",&f_lept4_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_massjetjet",&f_massjetjetsignal_4mu);
  //treesignal_4mu->SetBranchAddress("f_deltaphi_norm",&f_deltaphisignal_4mu);
  treesignal_4mu->SetBranchAddress("f_bdiscjet1",&f_bdiscjet1signal_4mu);
  treesignal_4mu->SetBranchAddress("f_bdiscjet2",&f_bdiscjet2signal_4mu);
  //treesignal_4mu->SetBranchAddress("f_etajet1",&f_etajet1signal_4mu);
  //treesignal_4mu->SetBranchAddress("f_etajet2",&f_etajet2signal_4mu);
  //treesignal_4mu->SetBranchAddress("f_Nbjets",&Nbjetssignal_4mu);
  treesignal_4mu->SetBranchAddress("f_deltar_norm",&f_deltarsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_MET_norm",&f_METsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_ptjet1",&f_ptjet1signal_4mu);
  treesignal_4mu->SetBranchAddress("f_ptjet2",&f_ptjet2signal_4mu);
  //Set other branches for signal_4mu
  //treesignal_4mu->SetBranchAddress("f_mass4l",&masssignal_4mu);
  //treesignal_4mu->SetBranchAddress("f_jet1_pt",&jet1_pt);
  //treesignal_4mu->SetBranchAddress("f_Ngood",&Ngoodsignal_4mu);
  //treesignal_4mu->SetBranchAddress("f_njets_pass",&jetsignal_4mu);
  //treesignal_4mu->SetBranchAddress("f_totbjet",&bjetsignal_4mu);
  //treesignal_4mu->SetBranchAddress("f_pt3",&pt3signal_4mu);
  
  treesignal_4mu->SetBranchAddress("f_weight",&weightsignal_4mu);
  //TBranch *BDT_shape_MLP = treesignal_4mu->Branch("BDT_shape_MLP", &BDT_shape_fill, "BDT_shape_MLP/F");
  //TBranch *BDT_shape_4_conf = treesignal_4mu->Branch("BDT_shape_4_conf", &BDT_shape_fill_4_conf, "BDT_shape_4_conf/F");
  TBranch *BDT_shape_10_conf = treesignal_4mu->Branch("BDT_shape_10_conf", &BDT_shape_fill_10_conf, "BDT_shape_10_conf/F");
  //TBranch *BDT_shape_11_conf = treesignal_4mu->Branch("BDT_shape_11_conf", &BDT_shape_fill_11_conf, "BDT_shape_11_conf/F");

  //for(auto flag : data){    
    //signal_4mu_4mu  
    for(int i = 0; i < treesignal_4mu->GetEntries(); i++){
      
      treesignal_4mu->GetEntry(i);
      /*
      std::vector<double> inputsignal_4mu = {f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_lept4_ptsignal_4mu, f_METsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu};//jetsignal_4mu}; 
      std::vector<string> parameterssignal_4mu =    {"f_lept1_pt", "f_lept2_pt", "f_lept3_pt", "f_lept4_pt", "f_MET_norm", "f_bdiscjet1", "f_bdiscjet2"};//std::vector<double> inputsignal_4mu_4_conf = {f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_ptjet1signal_4mu, f_ptjet2signal_4mu, f_deltarsignal_4mu, f_METsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu}
      ReadBDT_Grad_200_2 read(parameterssignal_4mu);
      double ciao =      read.GetMvaValue(inputsignal_4mu);
      cout << ciao << "\n";*/
/*
      //configuration 4
      std::vector<double> inputsignal_4mu_4_conf = {f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_lept4_ptsignal_4mu, f_ptjet1signal_4mu, f_ptjet2signal_4mu, f_deltarsignal_4mu, f_METsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu};
      std::vector<string> parameterssignal_4mu_4_conf = {"f_lept1_pt", "f_lept2_pt", "f_lept3_pt", "f_lept4_pt", "f_ptjet1", "f_ptjet2", "f_deltar_norm", "f_MET_norm", "f_bdiscjet1", "f_bdiscjet2"};
      
      ReadBDT_Grad_400_2_4 read(parameterssignal_4mu_4_conf);      
     
      double discriminantsignal_4mu_4_conf = read.GetMvaValue(inputsignal_4mu_4_conf);
      //cout << discriminantsignal_4mu_4_conf << "\n";
      
      //if(f_massjetjetsignal_4mu >= 0){
      BDT_shape_fill_4_conf = discriminantsignal_4mu_4_conf;
      BDT_shape_4_conf->Fill();
      //}*/

    //configuration 10
    std::vector<double> inputsignal_4mu_10_conf = {f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_lept4_ptsignal_4mu, f_massjetjetsignal_4mu, f_ptjet1signal_4mu, f_ptjet2signal_4mu, f_deltarsignal_4mu, f_METsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu};
    std::vector<string> parameterssignal_4mu_10_conf = {"f_lept1_pt", "f_lept2_pt", "f_lept3_pt", "f_lept4_pt", "f_massjetjet", "f_ptjet1", "f_ptjet2", "f_deltar_norm", "f_MET_norm", "f_bdiscjet1", "f_bdiscjet2"}; 

    ReadBDT_Grad_400_2 read(parameterssignal_4mu_10_conf);      
     
    double discriminantsignal_4mu_10_conf = read.GetMvaValue(inputsignal_4mu_10_conf);
    //cout << discriminantsignal_4mu_10_conf << "\n";

    BDT_shape_fill_10_conf = discriminantsignal_4mu_10_conf;
    BDT_shape_10_conf->Fill();
/*
    //configuration 11
    std::vector<double> inputsignal_4mu_11_conf = {f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_lept4_ptsignal_4mu, f_massjetjetsignal_4mu, f_ptjet1signal_4mu, f_ptjet2signal_4mu, f_deltarsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu};
    std::vector<string> parameterssignal_4mu_11_conf = {"f_lept1_pt", "f_lept2_pt", "f_lept3_pt", "f_lept4_pt", "f_massjetjet", "f_ptjet1", "f_ptjet2", "f_deltar_norm", "f_bdiscjet1", "f_bdiscjet2"};

    ReadBDT_Grad_400_2_11 read(parameterssignal_4mu_11_conf);      
     
    double discriminantsignal_4mu_11_conf = read.GetMvaValue(inputsignal_4mu_11_conf);
    //cout << discriminantsignal_4mu_11_conf << "\n";
    BDT_shape_fill_11_conf = discriminantsignal_4mu_11_conf;
    BDT_shape_11_conf->Fill();*/
    }
    //cout << "first print " << "\n";
    //treesignal_4mu->Print();

    treesignal_4mu->Write();

    //treesignal_4mu->Refresh();
    //cout << "second print " << "\n";
    //treesignal_4mu->Print();
    //treesignal_4mu->Refresh();
    flag->Close();
  }
  
}
