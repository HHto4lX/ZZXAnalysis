#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#include "evalEventSF.C"

///////////////////////BEFORE LOOPING ON THE EVENTS, PUT THE FOLLOWING///////////////////////////////////


///////////// SET UP B-TAG CALIBRATION ///////////////
    
    // set up calibration + reader
    cout << "Loading the .csv file..." << endl;
    
    std::string inputCSVfile = "CSVv2_Moriond17_B_H.csv";  
    std::string measType = "iterativefit";
    std::string sysType = "central";
    std::string sysTypeJESUp = "up_jes";
    std::string sysTypeJESDown = "down_jes";
    std::string sysTypeHFUp = "up_hf";
    std::string sysTypeHFDown = "down_hf";
    std::string sysTypeLFUp = "up_lf";
    std::string sysTypeLFDown = "down_lf";
    std::string sysTypehfstats1Up = "up_hfstats1";
    std::string sysTypehfstats1Down = "down_hfstats1";
    std::string sysTypehfstats2Up = "up_hfstats2";
    std::string sysTypehfstats2Down = "down_hfstats2";
    std::string sysTypelfstats1Up = "up_lfstats1";
    std::string sysTypelfstats1Down = "down_lfstats1";
    std::string sysTypelfstats2Up = "up_lfstats2";
    std::string sysTypelfstats2Down = "down_lfstats2";
    std::string sysTypecfErr1Up = "up_cferr1";
    std::string sysTypecfErr1Down = "down_cferr1";
    std::string sysTypecfErr2Up = "up_cferr2";
    std::string sysTypecfErr2Down = "down_cferr2";
    
    BTagCalibration calib("csvv2", inputCSVfile);
    //nominal
    BTagCalibrationReader CSVreader(BTagEntry::OP_RESHAPING, sysType);       
    CSVreader.load(calib, BTagEntry::FLAV_B, measType);
    CSVreader.load(calib, BTagEntry::FLAV_C, measType);
    CSVreader.load(calib, BTagEntry::FLAV_UDSG, measType);
    //up_jes shift
    BTagCalibrationReader CSVreaderJESUp(BTagEntry::OP_RESHAPING, sysTypeJESUp);       
    CSVreaderJESUp.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderJESUp.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderJESUp.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_jes shift
    BTagCalibrationReader CSVreaderJESDown(BTagEntry::OP_RESHAPING, sysTypeJESDown);       
    CSVreaderJESDown.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderJESDown.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderJESDown.load(calib, BTagEntry::FLAV_UDSG, measType);
     //up_hf shift
    BTagCalibrationReader CSVreaderHFUp(BTagEntry::OP_RESHAPING, sysTypeHFUp);       
    CSVreaderHFUp.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderHFUp.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderHFUp.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_hf shift
    BTagCalibrationReader CSVreaderHFDown(BTagEntry::OP_RESHAPING, sysTypeHFDown);       
    CSVreaderHFDown.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderHFDown.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderHFDown.load(calib, BTagEntry::FLAV_UDSG, measType);
    //up_lf shift
    BTagCalibrationReader CSVreaderLFUp(BTagEntry::OP_RESHAPING, sysTypeLFUp);       
    CSVreaderLFUp.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderLFUp.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderLFUp.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_lf shift
    BTagCalibrationReader CSVreaderLFDown(BTagEntry::OP_RESHAPING, sysTypeLFDown);       
    CSVreaderLFDown.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderLFDown.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderLFDown.load(calib, BTagEntry::FLAV_UDSG, measType);
    //up_hfstats1 shift
    BTagCalibrationReader CSVreaderhfstats1Up(BTagEntry::OP_RESHAPING, sysTypehfstats1Up);       
    CSVreaderhfstats1Up.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderhfstats1Up.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderhfstats1Up.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_hfstats1 shift
    BTagCalibrationReader CSVreaderhfstats1Down(BTagEntry::OP_RESHAPING, sysTypehfstats1Down);       
    CSVreaderhfstats1Down.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderhfstats1Down.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderhfstats1Down.load(calib, BTagEntry::FLAV_UDSG, measType);
     //up_lfstats2 shift
    BTagCalibrationReader CSVreaderhfstats2Up(BTagEntry::OP_RESHAPING, sysTypehfstats2Up);       
    CSVreaderhfstats2Up.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderhfstats2Up.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderhfstats2Up.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_lfstats2 shift
    BTagCalibrationReader CSVreaderhfstats2Down(BTagEntry::OP_RESHAPING, sysTypehfstats2Down);       
    CSVreaderhfstats2Down.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderhfstats2Down.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderhfstats2Down.load(calib, BTagEntry::FLAV_UDSG, measType);
    //up_lfstats1 shift
    BTagCalibrationReader CSVreaderlfstats1Up(BTagEntry::OP_RESHAPING, sysTypelfstats1Up);       
    CSVreaderlfstats1Up.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderlfstats1Up.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderlfstats1Up.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_lfstats1 shift
    BTagCalibrationReader CSVreaderlfstats1Down(BTagEntry::OP_RESHAPING, sysTypelfstats1Down);       
    CSVreaderlfstats1Down.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderlfstats1Down.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderlfstats1Down.load(calib, BTagEntry::FLAV_UDSG, measType);
     //up_lfstats2 shift
    BTagCalibrationReader CSVreaderlfstats2Up(BTagEntry::OP_RESHAPING, sysTypelfstats2Up);       
    CSVreaderlfstats2Up.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderlfstats2Up.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderlfstats2Up.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_lfstats2 shift
    BTagCalibrationReader CSVreaderlfstats2Down(BTagEntry::OP_RESHAPING, sysTypelfstats2Down);       
    CSVreaderlfstats2Down.load(calib, BTagEntry::FLAV_B, measType);
    CSVreaderlfstats2Down.load(calib, BTagEntry::FLAV_C, measType);
    CSVreaderlfstats2Down.load(calib, BTagEntry::FLAV_UDSG, measType);
    //up_cferr1 shift
    BTagCalibrationReader CSVreadercfErr1Up(BTagEntry::OP_RESHAPING, sysTypecfErr1Up);       
    CSVreadercfErr1Up.load(calib, BTagEntry::FLAV_B, measType);
    CSVreadercfErr1Up.load(calib, BTagEntry::FLAV_C, measType);
    CSVreadercfErr1Up.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_cferr1 shift
    BTagCalibrationReader CSVreadercfErr1Down(BTagEntry::OP_RESHAPING, sysTypecfErr1Down);       
    CSVreadercfErr1Down.load(calib, BTagEntry::FLAV_B, measType);
    CSVreadercfErr1Down.load(calib, BTagEntry::FLAV_C, measType);
    CSVreadercfErr1Down.load(calib, BTagEntry::FLAV_UDSG, measType);
     //up_cferr2 shift
    BTagCalibrationReader CSVreadercfErr2Up(BTagEntry::OP_RESHAPING, sysTypecfErr2Up);       
    CSVreadercfErr2Up.load(calib, BTagEntry::FLAV_B, measType);
    CSVreadercfErr2Up.load(calib, BTagEntry::FLAV_C, measType);
    CSVreadercfErr2Up.load(calib, BTagEntry::FLAV_UDSG, measType);
    //down_cferr2 shift
    BTagCalibrationReader CSVreadercfErr2Down(BTagEntry::OP_RESHAPING, sysTypecfErr2Down);       
    CSVreadercfErr2Down.load(calib, BTagEntry::FLAV_B, measType);
    CSVreadercfErr2Down.load(calib, BTagEntry::FLAV_C, measType);
    CSVreadercfErr2Down.load(calib, BTagEntry::FLAV_UDSG, measType);
    
    
    
    cout << "Input CSV weight file = " << inputCSVfile << "; measurementType = " << measType << ";" << endl;
    
    
    
    
    //////////////////////////////WHEN IT'S TIME TO FILL THE HISTOGRAMS...////////////////////////////////////////
    
    double * scaleFactors;
    scaleFactors = evalEventSF( myak4jet_nJets, mynJets, myak4jetFlavor, myak4jetEta, myak4jetPt, myak4jetBtag, myjetEtaSub0, myjetPtSub0, myjetBtagSub0, myjetFlavorSub0, myjetEtaSub1, myjetPtSub1, myjetBtagSub1, myjetFlavorSub1, CSVreader, CSVreaderJESUp, CSVreaderJESDown, CSVreaderHFUp, CSVreaderHFDown, CSVreaderLFUp, CSVreaderLFDown, CSVreaderhfstats1Up, CSVreaderhfstats1Down, CSVreaderhfstats2Up, CSVreaderhfstats2Down, CSVreaderlfstats1Up, CSVreaderlfstats1Down, CSVreaderlfstats2Up, CSVreaderlfstats2Down, CSVreadercfErr1Up, CSVreadercfErr1Down, CSVreadercfErr2Up, CSVreadercfErr2Down ); 
    
    
    //fill everything you want like... 
    
    //higgsMass->Fill(mH, scaleFactors[0]*yourEventWeights); 
    //higgsMass_JESUp->Fill(mH, scaleFactors[1]*yourEventWeights);
    //higgsMass_JESDown->Fill(mH, scaleFactors[2]*yourEventWeights); 
    //etc etc...
    
    //and finally...
    
    
    delete scaleFactors;
    
    
    
