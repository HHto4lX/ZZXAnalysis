#include <TH2F.h> 
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

void yield_counter_ZZbb (const std::string &inputList) 
{ 
  std::string outname = "/afs/cern.ch/work/l/lborgono/private/CMSSW_10_2_5_patch1/src/ZZXAnalysis/AnalysisStep/test/Plotter/yieldCounter_" + inputList;
  std::ofstream outFile(outname);

  vector<std::string> inputFile;
  std::ifstream inFile;
  std::string filename;
  float yield = 0;
  std::cout << "Input file: " << inputList.c_str() << "\n";
  inFile.open(inputList.c_str());
  
  while (std::getline(inFile, filename)) inputFile.push_back(filename.c_str());
  std::cout << "Number of input files: " << inputFile.size() << "\n";
  
  for(int i = 0; i<(inputFile.size());i++)
    {
      TFile* file;
      file = new TFile (inputFile[i].c_str(), "READ"); 
   
      TH1F* h_massjj = (TH1F*)file->Get("h1_methodBM_M_4L");
      yield = (Float_t)h_massjj->Integral(0,500);
     
      std::cout << "file " << inputFile[i].c_str() << endl;
      std::cout << "yield " << yield << endl;
      outFile  << yield << "\n";
      
    }
  
}
