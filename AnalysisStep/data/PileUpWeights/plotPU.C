
void plotPU(){

  TFile* fIn2018 = TFile::Open("pu_weights_2018.root");
  
  TH1F* hMC          = (TH1F*)fIn2018->Get("MC_out_of_the_box");
  TH1F* hMC_reweight = (TH1F*)fIn2018->Get("MC_reweighted");
  TH1F* hData2018    = (TH1F*)fIn2018->Get("Data");

  gStyle->SetOptStat(false);  

  hMC_reweight->SetMaximum(0.065);
  hMC_reweight->SetLineColor(kBlue);
  hMC_reweight->SetFillColor(kBlue-9);
  hMC_reweight->SetFillStyle(3345);

  hMC->SetLineColor(kBlack);
  
  hData2018->SetMarkerStyle(21);
  hData2018->SetMarkerSize(0.8);
  hData2018->SetMarkerColor(kRed);


  hMC_reweight->Draw("histo");
  hMC->Draw("samehisto");
  hData2018->Draw("sameP");
  TLegend* leg = new TLegend(0.56,0.71,0.88,0.88);
  leg->AddEntry(hMC,          "MC out-of-the-box",   "f");
  leg->AddEntry(hMC_reweight, "MC reweighed",        "f");
  leg->AddEntry(hData2018,    "Data 2018, 59.74/fb", "p");
  leg->SetLineColor(kWhite);
  leg->Draw();

}
