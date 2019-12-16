// ************
// run with: root -l readBDT.C
// ************

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
//#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include "TMVAClassification_BDT_Grad_400_2.class.C"


void readBDT(){
using namespace std;
//int main (int argc, char ** argv){
  
  std::ifstream infile_data, infile_mc;
  infile_data.open("samples_2018.txt");

  std::string inputfilename_data;
  std::vector < TFile* > data;

  while(std::getline(infile_data,inputfilename_data)){
    data.push_back(new TFile(static_cast<TString>(inputfilename_data), "update"));
  }
  infile_data.close();

  for(auto flag : data){

  flag->cd();  
  cout << "name " << flag->GetName() << "\n";

  TTree *treesignal_4mu = (TTree*)flag->Get("reducedTree");//use your tree

  float f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_lept4_ptsignal_4mu, f_deltarsignal_4mu, f_METsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu, weightsignal_4mu;
  float BDT_shape_fill;

  //get the branch to feed the BDT
  treesignal_4mu->SetBranchAddress("f_lept1_pt",&f_lept1_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_lept2_pt",&f_lept2_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_lept3_pt",&f_lept3_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_lept4_pt",&f_lept4_ptsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_bdiscjet1",&f_bdiscjet1signal_4mu);
  treesignal_4mu->SetBranchAddress("f_bdiscjet2",&f_bdiscjet2signal_4mu);
  treesignal_4mu->SetBranchAddress("f_deltar_norm",&f_deltarsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_MET_norm",&f_METsignal_4mu);
  treesignal_4mu->SetBranchAddress("f_weight",&weightsignal_4mu);
  TBranch *BDT_shape_4 = treesignal_4mu->Branch("BDT_shape_4", &BDT_shape_fill, "BDT_shape_4/F");

    for(int i = 0; i < treesignal_4mu->GetEntries(); i++){
      
      treesignal_4mu->GetEntry(i);
      //how the TMVAClassification_BDT_Grad_400_2.class.C works 
      //PAY ATTENTION TO THE CORRECT ORDER (is this one DON'T change the order)
      std::vector<double> inputsignal_4mu = {f_lept1_ptsignal_4mu, f_lept2_ptsignal_4mu, f_lept3_ptsignal_4mu, f_lept4_ptsignal_4mu, f_deltarsignal_4mu, f_METsignal_4mu, f_bdiscjet1signal_4mu, f_bdiscjet2signal_4mu};//jetsignal_4mu};
      
      std::vector<string> parameterssignal_4mu =    {"f_lept1_pt", "f_lept2_pt", "f_lept3_pt", "f_lept4_pt", "f_deltar_norm", "f_MET_norm", "f_bdiscjet1", "f_bdiscjet2"};
      
      ReadBDT_Grad_400_2 read(parameterssignal_4mu);      

      //the output of this function is the score of the BDT
      double discriminantsignal_4mu = read.GetMvaValue(inputsignal_4mu);
      
      
      BDT_shape_fill = discriminantsignal_4mu;
      BDT_shape_4->Fill();
    }

    treesignal_4mu->Write();
    flag->Close();
  }
  
}
