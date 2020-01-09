void TMVA_CNN_Classification_2017() {

  //Signal directory
  //TString S_DIR = "/lustre/cms/store/user/taliercio/Twiki/histos4mu_25ns_ANok/output_HHTo2B4L_madgraph_pythia8.root";
  TString S_DIR = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_sig_good/output_merged_node_SM.root";
  //TString S_DIR_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_HHTo2B4L_madgraph_pythia8.root";
  //TString S_DIR_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_HHTo2B4L_madgraph_pythia8.root";
  TString S_DIR_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_sig_good/output_merged_node_SM.root";
  TString S_DIR_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_sig_good/output_merged_node_SM.root";

  //Bkg directory 4mu

  TString B_ttH = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_ttZ = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root"; // manca histos4mu_25 che il file era corrotto
  TString B_ggH = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_VBF = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_WminusH = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root";
  TString B_WplusH = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root";
  TString B_ZH = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8.root";
  TString B_ZZ = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_ZZTo4L_13TeV_powheg_pythia8.root";
  TString B_ttW = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root";
  TString B_ZZ_4mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root";
  TString B_ZZ_2mu2tau = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4mu_25ns_bkg2017_new/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root";
  //TString B_ggH = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/ggH.root";

  //Bkg directory 4e
  
  TString B_ttH_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_ttZ_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root";
  TString B_ggH_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_VBF_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_WminusH_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root";
  TString B_WplusH_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root";
  TString B_ZH_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8.root";
  TString B_ZZ_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_ZZTo4L_13TeV_powheg_pythia8.root";
  TString B_ttW_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root";
  TString B_ZZ_4e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root";
  TString B_ZZ_2e2tau = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root";
  TString B_WZZ_e = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos4e_25ns_bkg2017_new/output_WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.root";
 
  //Bkg directory 2e2mu

  TString B_ttH_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_ttZ_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root";
  TString B_ggH_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_VBF_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8.root";
  TString B_WminusH_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root";
  TString B_WplusH_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8.root";
  TString B_ZH_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8.root";
  TString B_ZZ2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_ZZTo4L_13TeV_powheg_pythia8.root";
  TString B_ttW_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root";
  TString B_ZZ_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root";
  TString B_ZZ_2e2tau_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root";
  TString B_ZZ_2mu2tau_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root";
  TString B_WZZ_2e2mu = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.root";
  TString B_tt = "/lustre/cms/store/user/atalierc/HH_signal_2017_new/histos2e2mu_25ns_bkg2017_new/output_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root";

TMVA::Tools::Instance();

//TMVA::PyMethodBase::PyInitialize();
gSystem->Setenv("KERAS_BACKEND","tensorflow");

 bool is_4mu = false;
 bool is_4e = false;
 bool is_2e2mu = true;

 auto outputFile = TFile::Open("MLP_ClassificationOutput_2017_2e2mu_Nomassjet.root", "RECREATE"); //out della bdt
 
 TMVA::Factory factory("TMVAClassification", outputFile,
		       "!V:ROC:!Silent:Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
 

 TMVA::DataLoader * loader = new TMVA::DataLoader("dataset_new_2e2mu_2017_new_new_Nomassjet"); //crea directory con weights and plots , deve contenere le variabili necessarie per il training

 loader->AddVariable("f_lept1_pt",'F');
 loader->AddVariable("f_lept2_pt",'F');
 loader->AddVariable("f_lept3_pt",'F');
 loader->AddVariable("f_lept4_pt",'F');
 /*loader->AddVariable("f_lept1_eta",'F');
   loader->AddVariable("f_lept2_eta",'F');
   loader->AddVariable("f_lept3_eta",'F');
   loader->AddVariable("f_lept4_eta",'F');
   loader->AddVariable("f_lept1_phi",'F');
   loader->AddVariable("f_lept2_phi",'F');
   loader->AddVariable("f_lept3_phi",'F');
   loader->AddVariable("f_lept4_phi",'F');*/
 //loader->AddVariable("f_massjetjet",'F');
 //loader->AddVariable("f_pt3",'F');
 //loader->AddVariable("f_ptjet1",'F');
 //loader->AddVariable("f_ptjet2",'F');
 //loader->AddVariable("f_phijet1",'F');
 //loader->AddVariable("f_phijet2",'F');
 
 
 //loader->AddVariable("f_deltaphi_norm",'F');//tolgo
 loader->AddVariable("f_deltar_norm",'F');//tolgo
 loader->AddVariable("f_MET_norm",'F');
 loader->AddVariable("f_bdiscjet1",'F');//tolgo
 loader->AddVariable("f_bdiscjet2",'F');//tolgo
 //loader->AddVariable("f_etajet1",'F');//attenzione che questo è il phi
 //loader->AddVariable("f_etajet2",'F');
 //loader->AddVariable("f_Nbjets",'I');//f_Nbjets //f_totbjet
 //loader->AddVariable("f_njets_pass",'F');
 //}
  
 TString inputFileName_s = S_DIR; //"Higgs_data.root";
 TString inputFileName_s_e = S_DIR_e; 
 TString inputFileName_s_2e2mu = S_DIR_2e2mu; 



 TString inputFileName_bttH = B_ttH;
 TString inputFileName_bttZ = B_ttZ;
 TString inputFileName_bggH = B_ggH;
 TString inputFileName_bVBF = B_VBF;
 TString inputFileName_bWminusH = B_WminusH;
 TString inputFileName_bWplusH = B_WplusH;
 TString inputFileName_bZH = B_ZH;
 TString inputFileName_bZZ = B_ZZ;
 TString inputFileName_bttW = B_ttW;
 TString inputFileName_bZZ4mu = B_ZZ_4mu;
 TString inputFileName_bZZ2mu2tau = B_ZZ_2mu2tau;
 
 TString inputFileName_bttH_e = B_ttH_e;
 TString inputFileName_bttZ_e = B_ttZ_e;
 TString inputFileName_bggH_e = B_ggH_e;
 TString inputFileName_bVBF_e = B_VBF_e;
 TString inputFileName_bWminusH_e = B_WminusH_e;
 TString inputFileName_bWplusH_e = B_WplusH_e;
 TString inputFileName_bZH_e = B_ZH_e;
 TString inputFileName_bZZ_e  = B_ZZ_e;
 TString inputFileName_bttW_e = B_ttW_e;
 TString inputFileName_bZZ4e = B_ZZ_4e;
 TString inputFileName_bZZ2e2tau = B_ZZ_2e2tau;
 TString inputFileName_bWZZ_e = B_WZZ_e;

 TString inputFileName_bttH_2e2mu = B_ttH_2e2mu;
 TString inputFileName_bttZ_2e2mu = B_ttZ_2e2mu;
 TString inputFileName_bggH_2e2mu = B_ggH_2e2mu;
 TString inputFileName_bVBF_2e2mu = B_VBF_2e2mu;
 TString inputFileName_bWminusH_2e2mu = B_WminusH_2e2mu;
 TString inputFileName_bWplusH_2e2mu = B_WplusH_2e2mu;
 TString inputFileName_bZH_2e2mu = B_ZH_2e2mu;
 TString inputFileName_bZZ_2e2mu  = B_ZZ2e2mu;
 TString inputFileName_bttW_2e2mu = B_ttW_2e2mu;
 TString inputFileName_bZZ2e2mu = B_ZZ_2e2mu;
 TString inputFileName_bZZ2e2tau_2e2mu = B_ZZ_2e2tau_2e2mu;
 TString inputFileName_bZZ2mu2tau_2e2mu = B_ZZ_2mu2tau_2e2mu;
 TString inputFileName_bWZZ_2e2mu = B_WZZ_2e2mu;
 TString inputFileName_btt = B_tt;


 auto inputFile_s = TFile::Open( inputFileName_s );
 auto inputFile_s_e = TFile::Open( inputFileName_s_e );
 auto inputFile_s_2e2mu = TFile::Open( inputFileName_s_2e2mu );


 auto inputFile_bttH = TFile::Open( inputFileName_bttH );
 auto inputFile_bttZ = TFile::Open( inputFileName_bttZ );
 auto inputFile_bggH = TFile::Open( inputFileName_bggH );
 auto inputFile_bVBF = TFile::Open( inputFileName_bVBF );
 auto inputFile_bWminusH = TFile::Open( inputFileName_bWminusH );
 auto inputFile_bWplusH = TFile::Open( inputFileName_bWplusH );
 auto inputFile_bZH = TFile::Open( inputFileName_bZH );
 auto inputFile_bZZ = TFile::Open( inputFileName_bZZ );
 auto inputFile_bttW = TFile::Open( inputFileName_bttW );
 auto inputFile_bZZ4mu = TFile::Open( inputFileName_bZZ4mu );
 auto inputFile_bZZ2mu2tau = TFile::Open( inputFileName_bZZ2mu2tau );

 
 auto inputFile_bttH_e = TFile::Open( inputFileName_bttH_e );
 auto inputFile_bttZ_e = TFile::Open( inputFileName_bttZ_e );
 auto inputFile_bggH_e = TFile::Open( inputFileName_bggH_e );
 auto inputFile_bVBF_e = TFile::Open( inputFileName_bVBF_e );
 auto inputFile_bWminusH_e = TFile::Open( inputFileName_bWminusH_e );
 auto inputFile_bWplusH_e = TFile::Open( inputFileName_bWplusH_e );
 auto inputFile_bZH_e = TFile::Open( inputFileName_bZH_e );
 auto inputFile_bZZ_e = TFile::Open( inputFileName_bZZ_e );
 auto inputFile_bttW_e = TFile::Open( inputFileName_bttW_e );
 auto inputFile_bZZ4e = TFile::Open( inputFileName_bZZ4e );
 auto inputFile_bZZ2e2tau = TFile::Open( inputFileName_bZZ2e2tau );
 auto inputFile_bWZZ_e = TFile::Open( inputFileName_bWZZ_e );
 
 auto inputFile_bttH_2e2mu = TFile::Open( inputFileName_bttH_2e2mu );
 auto inputFile_bttZ_2e2mu = TFile::Open( inputFileName_bttZ_2e2mu );
 auto inputFile_bggH_2e2mu = TFile::Open( inputFileName_bggH_2e2mu );
 auto inputFile_bVBF_2e2mu = TFile::Open( inputFileName_bVBF_2e2mu );
 auto inputFile_bWminusH_2e2mu = TFile::Open( inputFileName_bWminusH_2e2mu );
 auto inputFile_bWplusH_2e2mu = TFile::Open( inputFileName_bWplusH_2e2mu );
 auto inputFile_bZH_2e2mu = TFile::Open( inputFileName_bZH_2e2mu );
 auto inputFile_bZZ_2e2mu = TFile::Open( inputFileName_bZZ_2e2mu );
 auto inputFile_bttW_2e2mu = TFile::Open( inputFileName_bttW_2e2mu );
 auto inputFile_bZZ2e2mu = TFile::Open( inputFileName_bZZ2e2mu );
 auto inputFile_bZZ2e2tau_2e2mu = TFile::Open( inputFileName_bZZ2e2tau_2e2mu );
 auto inputFile_bZZ2mu2tau_2e2mu = TFile::Open( inputFileName_bZZ2mu2tau_2e2mu );
 auto inputFile_bWZZ_2e2mu = TFile::Open( inputFileName_bWZZ_2e2mu );
 auto inputFile_btt = TFile::Open( inputFileName_btt );
 
// --- Register the training and test trees
 TString tree_name = "HZZ4LeptonsAnalysisReduced;1";
 TTree *signalTree     = (TTree*)inputFile_s->Get(tree_name);
 TTree *signalTree_e     = (TTree*)inputFile_s_e->Get(tree_name);
 TTree *signalTree_2e2mu     = (TTree*)inputFile_s_2e2mu->Get(tree_name);


 TTree *backgroundTree_ttH = (TTree*)inputFile_bttH->Get(tree_name);
 TTree *backgroundTree_ttZ = (TTree*)inputFile_bttZ->Get(tree_name);
 TTree *backgroundTree_ggH = (TTree*)inputFile_bggH->Get(tree_name);
 TTree *backgroundTree_VBF = (TTree*)inputFile_bVBF->Get(tree_name);
 TTree *backgroundTree_WminusH = (TTree*)inputFile_bWminusH->Get(tree_name);
 TTree *backgroundTree_WplusH = (TTree*)inputFile_bWplusH->Get(tree_name);
 TTree *backgroundTree_ZH = (TTree*)inputFile_bZH->Get(tree_name);
 TTree *backgroundTree_ZZ = (TTree*)inputFile_bZZ->Get(tree_name);
 TTree *backgroundTree_ttW = (TTree*)inputFile_bttW->Get(tree_name);
 TTree *backgroundTree_ZZ4mu = (TTree*)inputFile_bZZ4mu->Get(tree_name);
 TTree *backgroundTree_ZZ2mu2tau = (TTree*)inputFile_bZZ2mu2tau->Get(tree_name);

 
 TTree *backgroundTree_ttH_e = (TTree*)inputFile_bttH_e->Get(tree_name);
 TTree *backgroundTree_ttZ_e = (TTree*)inputFile_bttZ_e->Get(tree_name);
 TTree *backgroundTree_ggH_e = (TTree*)inputFile_bggH_e->Get(tree_name);
 TTree *backgroundTree_VBF_e = (TTree*)inputFile_bVBF_e->Get(tree_name);
 TTree *backgroundTree_WminusH_e = (TTree*)inputFile_bWminusH_e->Get(tree_name);
 TTree *backgroundTree_WplusH_e = (TTree*)inputFile_bWplusH_e->Get(tree_name);
 TTree *backgroundTree_ZH_e = (TTree*)inputFile_bZH_e->Get(tree_name);
 TTree *backgroundTree_ZZ_e = (TTree*)inputFile_bZZ_e->Get(tree_name);
 TTree *backgroundTree_ttW_e = (TTree*)inputFile_bttW_e->Get(tree_name);
 TTree *backgroundTree_ZZ4e = (TTree*)inputFile_bZZ4e->Get(tree_name);
 TTree *backgroundTree_ZZ2e2tau = (TTree*)inputFile_bZZ2e2tau->Get(tree_name);
 TTree *backgroundTree_WZZ_e = (TTree*)inputFile_bWZZ_e->Get(tree_name);
 
 TTree *backgroundTree_ttH_2e2mu = (TTree*)inputFile_bttH_2e2mu->Get(tree_name);
 TTree *backgroundTree_ttZ_2e2mu = (TTree*)inputFile_bttZ_2e2mu->Get(tree_name);
 TTree *backgroundTree_ggH_2e2mu = (TTree*)inputFile_bggH_2e2mu->Get(tree_name);
 TTree *backgroundTree_VBF_2e2mu = (TTree*)inputFile_bVBF_2e2mu->Get(tree_name);
 TTree *backgroundTree_WminusH_2e2mu = (TTree*)inputFile_bWminusH_2e2mu->Get(tree_name);
 TTree *backgroundTree_WplusH_2e2mu = (TTree*)inputFile_bWplusH_2e2mu->Get(tree_name);
 TTree *backgroundTree_ZH_2e2mu = (TTree*)inputFile_bZH_2e2mu->Get(tree_name);
 TTree *backgroundTree_ZZ_2e2mu = (TTree*)inputFile_bZZ_2e2mu->Get(tree_name);
 TTree *backgroundTree_ttW_2e2mu = (TTree*)inputFile_bttW_2e2mu->Get(tree_name);
 TTree *backgroundTree_ZZ2e2mu = (TTree*)inputFile_bZZ2e2mu->Get(tree_name);
 TTree *backgroundTree_ZZ2e2tau_2e2mu = (TTree*)inputFile_bZZ2e2tau_2e2mu->Get(tree_name);
 TTree *backgroundTree_ZZ2mu2tau_2e2mu = (TTree*)inputFile_bZZ2mu2tau_2e2mu->Get(tree_name);
 TTree *backgroundTree_WZZ_2e2mu = (TTree*)inputFile_bWZZ_2e2mu->Get(tree_name);
 TTree *backgroundTree_btt = (TTree*)inputFile_btt->Get(tree_name);
 
// global event weights per tree (see below for setting event-wise weights)
 Double_t signalWeight     = 1.0;
 Double_t signalWeight_e     = 1.0;
 Double_t signalWeight_2e2mu     = 1.0;

 
 Double_t backgroundWeightttH = 1.0;
 Double_t backgroundWeightttZ = 1.0;
 Double_t backgroundWeightggH = 1.0; 
 Double_t backgroundWeightVBF = 1.0; 
 Double_t backgroundWeightWminusH = 1.0; 
 Double_t backgroundWeightWplusH = 1.0; 
 Double_t backgroundWeightZH = 1.0; 
 Double_t backgroundWeightZZ = 1.0; 
 Double_t backgroundWeightttW = 1.0; 
 Double_t backgroundWeightZZ4mu = 1.0; 
 Double_t backgroundWeightZZ2mu2tau = 1.0; 

 
 Double_t backgroundWeightttH_e = 1.0;
 Double_t backgroundWeightttZ_e = 1.0;
 Double_t backgroundWeightggH_e = 1.0; 
 Double_t backgroundWeightVBF_e = 1.0; 
 Double_t backgroundWeightWminusH_e = 1.0; 
 Double_t backgroundWeightWplusH_e = 1.0; 
 Double_t backgroundWeightZH_e = 1.0; 
 Double_t backgroundWeightZZ_e = 1.0; 
 Double_t backgroundWeightttW_e = 1.0; 
 Double_t backgroundWeightZZ4e = 1.0; 
 Double_t backgroundWeightZZ2e2tau = 1.0;
 Double_t backgroundWeightWZZ_e = 1.0; 
 
 Double_t backgroundWeightttH_2e2mu = 1.0;
 Double_t backgroundWeightttZ_2e2mu = 1.0;
 Double_t backgroundWeightggH_2e2mu = 1.0; 
 Double_t backgroundWeightVBF_2e2mu = 1.0; 
 Double_t backgroundWeightWminusH_2e2mu = 1.0; 
 Double_t backgroundWeightWplusH_2e2mu = 1.0; 
 Double_t backgroundWeightZH_2e2mu = 1.0; 
 Double_t backgroundWeightZZ_2e2mu = 1.0; 
 Double_t backgroundWeightttW_2e2mu = 1.0; 
 Double_t backgroundWeightZZ2e2mu = 1.0; 
 Double_t backgroundWeightZZ2e2tau_2e2mu = 1.0; 
 Double_t backgroundWeightZZ2mu2tau_2e2mu = 1.0; 
 Double_t backgroundWeightWZZ_2e2mu = 1.0;
 Double_t backgroundWeightbtt = 1.0; 
 
// You can add an arbitrary number of signal or background trees

 // controlla che siano valri fisici (nelle ntuple di Angela)
 TCut mycuts = "f_deltaphi_norm > 0 && f_deltar_norm > 0 && f_MET_norm > 0 && f_bdiscjet1 >= 0 && f_bdiscjet2 >=0"; 
 TCut mycutb = "f_deltaphi_norm > 0 && f_deltar_norm > 0 && f_MET_norm > 0 && f_bdiscjet1 >= 0 && f_bdiscjet2 >=0";

 // si prepara al training (carica ntuple di segnale e fondo)
 if(is_4mu == true){ 
 loader->AddSignalTree    ( signalTree,     signalWeight     );

 loader->AddBackgroundTree( backgroundTree_ggH, backgroundWeightggH );
 loader->AddBackgroundTree( backgroundTree_ttH, backgroundWeightttH );
 loader->AddBackgroundTree( backgroundTree_ZZ, backgroundWeightZZ );
 loader->AddBackgroundTree( backgroundTree_ttZ, backgroundWeightttZ );
 loader->AddBackgroundTree( backgroundTree_VBF, backgroundWeightVBF );
 loader->AddBackgroundTree( backgroundTree_WminusH, backgroundWeightWminusH );
 loader->AddBackgroundTree( backgroundTree_WplusH, backgroundWeightWplusH );
 loader->AddBackgroundTree( backgroundTree_ZH, backgroundWeightZH );
 loader->AddBackgroundTree( backgroundTree_ttW, backgroundWeightttW );
 loader->AddBackgroundTree( backgroundTree_ZZ2mu2tau, backgroundWeightZZ2mu2tau );

 loader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                     "nTrain_Signal=6675:nTrain_Background=6675:nTest_Signal=2225:nTest_Background=2225:SplitMode=Random:NormMode=NumEvents:!V" );

 }
 
 if(is_4e == true){

 loader->AddSignalTree    ( signalTree_e,     signalWeight_e     );
 
 loader->AddBackgroundTree( backgroundTree_ggH_e, backgroundWeightggH_e );
 loader->AddBackgroundTree( backgroundTree_ttH_e, backgroundWeightttH_e );
 loader->AddBackgroundTree( backgroundTree_ZZ_e, backgroundWeightZZ_e );
 loader->AddBackgroundTree( backgroundTree_ttZ_e, backgroundWeightttZ_e );
 loader->AddBackgroundTree( backgroundTree_VBF_e, backgroundWeightVBF_e );
 loader->AddBackgroundTree( backgroundTree_WminusH_e, backgroundWeightWminusH_e );
 loader->AddBackgroundTree( backgroundTree_WplusH_e, backgroundWeightWplusH_e );
 loader->AddBackgroundTree( backgroundTree_ZH_e, backgroundWeightZH_e );
 loader->AddBackgroundTree( backgroundTree_ttW_e, backgroundWeightttW_e );
 loader->AddBackgroundTree( backgroundTree_ZZ2e2tau, backgroundWeightZZ2e2tau );
 loader->AddBackgroundTree( backgroundTree_WZZ_e, backgroundWeightWZZ_e );

 loader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                     "nTrain_Signal=4425:nTrain_Background=4425:nTest_Signal=1475:nTest_Background=1475:SplitMode=Random:NormMode=NumEvents:!V" );
 }

 
 if(is_2e2mu == true){
 loader->AddSignalTree    ( signalTree_2e2mu,     signalWeight_2e2mu     );

 loader->AddBackgroundTree( backgroundTree_ggH_2e2mu, backgroundWeightggH_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_ttH_2e2mu, backgroundWeightttH_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_ZZ_2e2mu, backgroundWeightZZ_2e2mu );
 loader->AddBackgroundTree( backgroundTree_ttZ_2e2mu, backgroundWeightttZ_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_VBF_2e2mu, backgroundWeightVBF_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_WminusH_2e2mu, backgroundWeightWminusH_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_WplusH_2e2mu, backgroundWeightWplusH_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_ZH_2e2mu, backgroundWeightZH_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_ttW_2e2mu, backgroundWeightttW_2e2mu );//
 loader->AddBackgroundTree( backgroundTree_ZZ2e2tau_2e2mu, backgroundWeightZZ2e2tau_2e2mu );
 loader->AddBackgroundTree( backgroundTree_ZZ2mu2tau_2e2mu, backgroundWeightZZ2mu2tau_2e2mu );
 loader->AddBackgroundTree( backgroundTree_WZZ_2e2mu, backgroundWeightWZZ_2e2mu );
 loader->AddBackgroundTree( backgroundTree_btt, backgroundWeightbtt );

 loader->PrepareTrainingAndTestTree( mycuts, mycutb,
				     "nTrain_Signal=10875:nTrain_Background=10875:nTest_Signal=3625nTest_Background=3625:SplitMode=Random:NormMode=NumEvents:!V" );
 }
 
 //TCut mycuts = ""; 
 //TCut mycutb = ""; 

 //loader->PrepareTrainingAndTestTree( mycuts, mycutb,
 //                                  "nTrain_Signal=10000:nTrain_Background=10000:nTest_Signal=2000:nTest_Background=2000:SplitMode=Random:NormMode=NumEvents:!V" );//+6000 per ggH 



 bool useBDT = true;  
 bool useMLP = true;
 bool useDNN = false;
 bool useCNN = false;
 bool useKeras = false;

 if(useBDT){
//Boosted Decision Trees
/*
factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_800_2", "!H:!V:NTrees=800:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_400_2", "!H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_200_2","!H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDep\
th=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_100_2","!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDep\
th=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_800_3", "!H:!V:NTrees=800:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_400_3", "!H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_200_3","!H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDep\
th=3" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_100_3","!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDep\
th=3" );





 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_800_2", "!H:!V:NTrees=800:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_400_2", "!H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_200_2", "!H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_100_2", "!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_800_3", "!H:!V:NTrees=800:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_400_3", "!H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_200_3", "!H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

 factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Ada_100_3", "!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=AdaBoost:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
*/

//   factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_800_2", "!H:!V:NTrees=800:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" ); 

   factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_400_2", "!H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" ); //metodo usato per il training 

   //   factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_200_2","!H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=5:MaxDepth=2" );

   //   factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_100_2","!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=5:MaxDepth=2" );

   //  factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_800_3", "!H:!V:NTrees=800:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=5:MaxDepth=3" );

   //  factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_400_3", "!H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=5:MaxDepth=3" );

   //  factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_200_3","!H:!V:NTrees=200:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=5:MaxDepth=3" );

   // factory.BookMethod(loader,TMVA::Types::kBDT, "BDT_Grad_100_3","!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=5:MaxDepth=3" );

 }


 // non usare
if (useDNN) {

     TString layoutString ("Layout=TANH|7,TANH|7,TANH|7,LINEAR");
//
      // Training strategies
      // one can catenate several training strategies
      TString training1("LearningRate=1e-3,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=32,TestRepetitions=5,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.2+0.2+0., Multithreading=False");//TestRepetitions=5

      TString trainingStrategyString ("TrainingStrategy=");
      trainingStrategyString += training1; // + "|" + training2 + "|" + training3;

      // General Options.
      TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                          "WeightInitialization=XAVIERUNIFORM");
      dnnOptions.Append (":"); dnnOptions.Append (layoutString);
      dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

      dnnOptions += ":Architecture=CPU";
      factory.BookMethod(loader, TMVA::Types::kDNN, "DNN_CPU", dnnOptions);

}

if (useCNN) {
    TString inputLayoutString("InputLayout=1|6|5");

// Batch Layout
    TString batchLayoutString("BatchLayout=256|1|64");


TString layoutString("Layout=CONV|10|3|3|1|1|1|1|RELU,CONV|10|3|3|1|1|1|1|RELU,MAXPOOL|2|2|1|1,RESHAPE|FLAT,DENSE|64|TANH,DENSE|1|LINEAR");



   // Training strategies.
   TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                     "ConvergenceSteps=20,BatchSize=256,TestRepetitions=5,"
                     "WeightDecay=1e-4,Regularization=None,"
                     "DropConfig=0.0+0.5+0.5+0.5, Multithreading=False");//TestRepetitions=5

   TString trainingStrategyString ("TrainingStrategy=");
   trainingStrategyString += training0; // + "|" + training1 + "|" + training2;   }

// General Options.
   TString cnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                       "WeightInitialization=XAVIERUNIFORM");

   cnnOptions.Append(":"); cnnOptions.Append(inputLayoutString);
   cnnOptions.Append(":"); cnnOptions.Append(batchLayoutString);
   cnnOptions.Append(":"); cnnOptions.Append(layoutString);
   cnnOptions.Append(":"); cnnOptions.Append(trainingStrategyString);
   cnnOptions.Append(":Architecture=CPU");

   //// New DL (CNN)


  factory.BookMethod(loader, TMVA::Types::kLD, "DL_CNN_CPU", cnnOptions);


}

if (useKeras) {
   factory.BookMethod(loader, TMVA::Types::kPyKeras,
                       "PyKeras","H:!V:VarTransform=None:FilenameModel=model_cnn.h5:"
                       "FilenameTrainedModel=trained_model_cnn.h5:NumEpochs=20:BatchSize=256");
}

 if (useMLP){//multilayer perceptron // rete neurale
   //factory.BookMethod(loader,TMVA::Types::kMLP,"MLP","H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=1000;HiddenLayers=N-23:TestRate=5:!UseRegulator");
   factory.BookMethod(loader,TMVA::Types::kMLP,"MLP_500_N-7","H:!V:NeuronType=sigmoid:VarTransform=Norm:NCycles=100;HiddenLayers=N-7:TestRate=5");//N-4;NCycles=500 N+1,N,N,N
   factory.BookMethod(loader,TMVA::Types::kMLP,"MLP_500_N","H:!V:NeuronType=sigmoid:VarTransform=Norm:NCycles=100;HiddenLayers=N:TestRate=5");
 }

 // la cosa piu impo del codice: traina, testa e cortuisce lo score (file che mi ha passato angela)
factory.TrainAllMethods();

factory.TestAllMethods();
factory.EvaluateAllMethods();

//%jsroot on

auto c1 = factory.GetROCCurve(loader);
// if (!gROOT->IsBatch()) TMVA::TMVAGui("MLP_ClassificationOutput_134.root");
 if (!gROOT->IsBatch()) TMVA::TMVAGui("MLP_ClassificationOutput_2017_2e2mu_Nomassjet.root");
c1->Draw();


// close outputfile to save output file
outputFile->Close();


   }
