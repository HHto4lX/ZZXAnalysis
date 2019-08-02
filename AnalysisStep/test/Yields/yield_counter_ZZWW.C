#include <TH2F.h> 
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

void yield_counter_ZZWW (const std::string &inputList) 
{ 
  std::string outname = "/afs/cern.ch/work/l/lborgono/private/CMSSW_10_2_5_patch1/src/ZZXAnalysis/AnalysisStep/test/Plotter/yieldCounter_" + inputList;
  std::ofstream outFile(outname);

  vector<std::string> inputFile;
  std::ifstream inFile;
  std::string filename;
  float yield1 = 0;
  float yield2 = 0;
  float yield2OS = 0;
  float yield2emu = 0;

  std::cout << "Input file: " << inputList.c_str() << "\n";
  inFile.open(inputList.c_str());
  
  while (std::getline(inFile, filename)) inputFile.push_back(filename.c_str());
  std::cout << "Number of input files: " << inputFile.size() << "\n";
  
  for(int i = 0; i<(inputFile.size());i++)
    {
      TFile* file;
      file = new TFile (inputFile[i].c_str(), "READ"); 
   
      TH1F* y1 = (TH1F*)file->Get("h1_yield_1extraL_4L");
      TH1F* y2 = (TH1F*)file->Get("h1_yield_2extraL_4L");
      TH1F* y2OS = (TH1F*)file->Get("h1_yield_2extraLOS_4L");
      TH1F* y2emu = (TH1F*)file->Get("h1_yield_2extraLemu_4L");

      yield1 = (Float_t)y1->GetBinContent(2);
      yield2 = (Float_t)y2->GetBinContent(2);
      yield2OS = (Float_t)y2OS->GetBinContent(2);
      yield2emu = (Float_t)y2emu->GetBinContent(2);
     
      std::cout << "file " << inputFile[i].c_str() << endl;
      std::cout << "yield1 " << yield1 << endl;
      std::cout << "yield2 " << yield2 << endl;
      std::cout << "yield2OS " << yield2OS << endl;
      std::cout << "yield2emu " << yield2emu << endl;
      outFile << yield1 << "\t" << yield2 << "\t" << yield2OS << "\t" << yield2emu << "\n";
      
    }
  
}
