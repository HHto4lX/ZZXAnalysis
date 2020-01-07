#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TFile.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMath.h"
#include "TSystem.h"
#include "TList.h"
#include <libgen.h>

using namespace std;

class loglike_r_4Lbb{
  
public: 
  loglike_r_4Lbb();
private:
 
};
 
loglike_r_4Lbb::loglike_r_4Lbb(){

 
  gROOT->ForceStyle(); 

//READ FILES
//======================================	
/////////////////////////////Signal_strenght/////////////////////////
  TFile *_file0 = TFile::Open                                                                                                                
    ("/afs/cern.ch/work/l/lborgono/private/HH/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/4LX/Cards_endorsement/Run2_exp_S0.root"); 
  TTree *ttstatonly=(TTree*)_file0->Get("limit");                                     
  TFile *_file1 = TFile::Open                                                                                                                  
    ("/afs/cern.ch/work/l/lborgono/private/HH/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/4LX/Cards_endorsement/Run2_exp_S1.root");                                   
  TTree *tterror=(TTree*)_file1->Get("limit");            
///////////////////////////klambda////////////////////////////////////
 
  // drawing
  gStyle->SetOptFit(101);
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 

  TCanvas *c1 = new TCanvas ("c1", "c1", 600, 600);
  c1->cd();
  c1->SetTicks(1,1);  
  c1->SetLeftMargin(0.14);
  c1->SetRightMargin(0.08);
  //c1->SetFillColor(0);
  //c1->SetBorderMode(0);
  //c1->SetFrameFillStyle(0);
  //c1->SetFrameBorderMode(0);
  //c1->SetTickx(0);
  //c1->SetTicky(0);
  c1->SetLogy(0);
  


//CREATE GRAPHS FROM ttree
//======================================	

  ttstatonly->Draw("2*deltaNLL:r","2*deltaNLL<10. && r != 1.0");
  TGraph *grstatonly = new TGraph(ttstatonly->GetSelectedRows(), ttstatonly->GetV2(), ttstatonly->GetV1());
  
  tterror->Draw("2*deltaNLL:r","2*deltaNLL<10. && r != 1.0");
  TGraph *grerror = new TGraph(tterror->GetSelectedRows(), tterror->GetV2(), tterror->GetV1());
  
  grstatonly->SetName("grstatonly");
  grstatonly->SetLineColor(kBlue);
  grstatonly->SetLineWidth(3);
  grstatonly->SetMarkerColor(kBlue);
  grstatonly->SetMarkerSize(0.);
  grstatonly->SetFillStyle(0.);
  grstatonly->SetMaximum(6);

  grstatonly->GetXaxis()->SetRangeUser(-90, 140);

  grstatonly->GetYaxis()->SetTitleOffset(1.75);
  grstatonly->GetXaxis()->SetTitleOffset(1.35);
  grstatonly->GetXaxis()->SetTitleSize(0.035);
  grstatonly->GetYaxis()->SetTitleSize(0.035);
  grstatonly->GetYaxis()->SetTitle("-2#Delta ln L");
  grstatonly->GetXaxis()->SetTitle("#mu = #sigma_{HH }/#sigma_{HH}^{SM}");  
  // grstatonly->GetXaxis()->SetTitle("k_{#lambda} = #lambda_{obs}/#lambda_{SM}");  
  grstatonly->Draw("ACP"); 
  //grstatonly->Draw("AL"); 
  
  
  grerror->SetName("grerror");
  grerror->SetLineColor(kRed+1);
  grerror->SetLineWidth(3);
  grerror->SetMarkerColor(kRed+1);
  grerror->SetMarkerSize(0.);
  grerror->SetFillStyle(0.);

  grerror->GetXaxis()->SetRangeUser(-90, 140);  

  // grSandB->SetMarkerStyle(22);
  // grSandB->GetYaxis()->SetTitle("2#Delta NLogL");
  // grSandB->GetXaxis()->SetTitle("#Delta #mu/#mu");  
  grerror->Draw("same CP"); 
  // grSandB->Draw("LSAME"); 
  
  
  
  
  //  TLegend* legend = new TLegend( 0.4, 0.5, 0.7, 0.8);
  TLegend* legend = new TLegend( 0.5, 0.7, 0.7, 0.8);
  // legend->SetFillColor(0);
  legend->SetTextSize(0.03);
  legend->SetTextFont(42);
  legend->SetBorderSize(0);
  legend->SetLineColor(0);
  legend->SetShadowColor(10);
  legend->SetHeader("#bf{#it{HH #rightarrow 4lbb}}");
  legend->AddEntry("grstatonly","Stat. only","lp");
  legend->AddEntry("grerror","Stat. + Sys.","lp");
    
  legend->Draw("same");

  //Spaces in the text are important to align! Don't remove them!

  TPaveText *Text_CMS = new TPaveText(0.1, 0.88, 0.93, 0.95, "brNDC");
  Text_CMS->SetTextAlign(31);
  Text_CMS->SetTextSize(0.05); 
  TString text_cms = "CMS                                ";
  Text_CMS->AddText(text_cms);
  Text_CMS->SetFillColor(kWhite);
  Text_CMS->SetFillStyle(0);
  Text_CMS->SetLineStyle(0);
  Text_CMS->SetBorderSize(0);
  Text_CMS->Draw();
  
  
  TPaveText *Text_prel = new TPaveText(0.1, 0.88, 0.93, 0.95, "brNDC");
  Text_prel->SetTextAlign(31);
  Text_prel->SetTextSize(0.035); 
  TString text_prel = "#bf{#it{Preliminary                                               }}";
  Text_prel->AddText(text_prel);
  Text_prel->SetFillColor(kWhite);
  Text_prel->SetFillStyle(0);
  Text_prel->SetLineStyle(0);
  Text_prel->SetBorderSize(0);
  Text_prel->Draw();

  TPaveText *Text_lumi_en = new TPaveText(0.1, 0.88, 1., 0.95, "brNDC");
  Text_lumi_en->SetTextAlign(31);
  Text_lumi_en->SetTextSize(0.035); 
  // TString text_lumi_en = "#bf{59.7 fb^{-1} (13 TeV)}   ";
  TString text_lumi_en = "#bf{137.2 fb^{-1} (13 TeV)}   ";
  Text_lumi_en->AddText(text_lumi_en);
  Text_lumi_en->SetFillColor(kWhite);
  Text_lumi_en->SetFillStyle(0);
  Text_lumi_en->SetLineStyle(0);
  Text_lumi_en->SetBorderSize(0);
  Text_lumi_en->Draw();
  

  
  TLine *line = new TLine(-90,1.,140.5,1.);
  //TLine *line = new TLine(0.68,1.,1.3,1.);
  line->SetLineColor(kBlack);
  //  line->SetLineStyle(9);
  line->SetLineStyle(2);
  line->SetLineWidth(3);
  line->Draw("same");
 
  TLine *line1 = new TLine(-90,3.8,140,3.8);
  //TLine *line1 = new TLine(0.68,3.8,1.3,3.8);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(2);
  line1->SetLineWidth(3);
  line1->Draw("same");


  TLatex sigma1;
  sigma1.SetTextAlign(12);
  sigma1.SetTextSize(0.025);
  //sigma1.DrawLatex(1.25,1.5,"1#sigma");
  sigma1.DrawLatex(115,0.85,"#bf{68% CL}");
  //  sigma1.DrawLatex(210,0.85,"#bf{68% CL}");



  TLatex sigma2;
  sigma2.SetTextAlign(12);
  sigma2.SetTextSize(0.025);
  //sigma1.DrawLatex(1.25,1.5,"1#sigma");
  sigma2.DrawLatex(115,3.65,"#bf{95% CL}");
  //  sigma2.DrawLatex(210,3.65,"#bf{95% CL}");



  c1->RedrawAxis();
  c1->GetFrame()->SetBorderSize( 12 );
  c1->Modified();
  c1->Update();


  c1->Print("plot_limit_4Lbb_Run2.pdf","pdf");
  c1->Print("plot_limit_4Lbb_Run2.png","png");

}                
