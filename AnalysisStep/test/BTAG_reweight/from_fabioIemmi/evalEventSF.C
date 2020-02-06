double * evalEventSF( int nAK4Jets, int nAK8Jets, vector<float> * AK4JetFlavor, vector<float> * AK4JetEta, vector<float> * AK4JetPt, vector<float> * AK4JetBtag, vector<float> * Sub0Eta, vector<float> * Sub0Pt, vector<float> * Sub0Btag,  vector<float> * Sub0Flavor, vector<float> * Sub1Eta, vector<float> * Sub1Pt, vector<float> * Sub1Btag, vector<float> * Sub1Flavor, BTagCalibrationReader CSV_nominal, BTagCalibrationReader CSV_JESUp, BTagCalibrationReader CSV_JESDown, BTagCalibrationReader CSV_HFUp, BTagCalibrationReader CSV_HFDown, BTagCalibrationReader CSV_LFUp, BTagCalibrationReader CSV_LFDown, BTagCalibrationReader CSV_hfstats1Up, BTagCalibrationReader CSV_hfstats1Down, BTagCalibrationReader CSV_hfstats2Up, BTagCalibrationReader CSV_hfstats2Down, BTagCalibrationReader CSV_lfstats1Up, BTagCalibrationReader CSV_lfstats1Down, BTagCalibrationReader CSV_lfstats2Up, BTagCalibrationReader CSV_lfstats2Down, BTagCalibrationReader CSV_cfErr1Up, BTagCalibrationReader CSV_cfErr1Down, BTagCalibrationReader CSV_cfErr2Up, BTagCalibrationReader CSV_cfErr2Down ) {
    
        double * output;
        output = new double[22]; //returned by evalEventSF
        
        double csv_weight_BF_ak4 = 1.0;
        double csv_weight_CF_ak4 = 1.0;
        double csv_weight_LF_ak4 = 1.0;
        double csv_weight_BF_ak8 = 1.0;
        double csv_weight_CF_ak8 = 1.0;
        double csv_weight_LF_ak8 = 1.0;
        double csv_weight_total_ak4 = 1.0;
        double csv_weight_total_ak8 = 1.0;
        double csv_weight_total = 1.0;
        
        //JES up and down
        double csv_weight_BF_ak4JESUp = 1.0;
        double csv_weight_CF_ak4JESUp = 1.0;
        double csv_weight_LF_ak4JESUp = 1.0;
        double csv_weight_BF_ak8JESUp = 1.0;
        double csv_weight_CF_ak8JESUp = 1.0;
        double csv_weight_LF_ak8JESUp = 1.0;
        double csv_weight_total_ak4JESUp = 1.0;
        double csv_weight_total_ak8JESUp = 1.0;
        double csv_weight_total_JESUp = 1.0;
        double csv_weight_BF_ak4JESDown = 1.0;
        double csv_weight_CF_ak4JESDown = 1.0;
        double csv_weight_LF_ak4JESDown = 1.0;
        double csv_weight_BF_ak8JESDown = 1.0;
        double csv_weight_CF_ak8JESDown = 1.0;
        double csv_weight_LF_ak8JESDown = 1.0;
        double csv_weight_total_ak4JESDown = 1.0;
        double csv_weight_total_ak8JESDown = 1.0;
        double csv_weight_total_JESDown = 1.0;
        
        //HF purity up and down
        double csv_weight_BF_ak4HFUp = 1.0;
        double csv_weight_CF_ak4HFUp = 1.0;
        double csv_weight_LF_ak4HFUp = 1.0;
        double csv_weight_BF_ak8HFUp = 1.0;
        double csv_weight_CF_ak8HFUp = 1.0;
        double csv_weight_LF_ak8HFUp = 1.0;
        double csv_weight_total_ak4HFUp = 1.0;
        double csv_weight_total_ak8HFUp = 1.0;
        double csv_weight_total_HFUp = 1.0;
        double csv_weight_BF_ak4HFDown = 1.0;
        double csv_weight_CF_ak4HFDown = 1.0;
        double csv_weight_LF_ak4HFDown = 1.0;
        double csv_weight_BF_ak8HFDown = 1.0;
        double csv_weight_CF_ak8HFDown = 1.0;
        double csv_weight_LF_ak8HFDown = 1.0;
        double csv_weight_total_ak4HFDown = 1.0;
        double csv_weight_total_ak8HFDown = 1.0;
        double csv_weight_total_HFDown = 1.0;
        
        //LF purity up and down
        double csv_weight_BF_ak4LFUp = 1.0;
        double csv_weight_CF_ak4LFUp = 1.0;
        double csv_weight_LF_ak4LFUp = 1.0;
        double csv_weight_BF_ak8LFUp = 1.0;
        double csv_weight_CF_ak8LFUp = 1.0;
        double csv_weight_LF_ak8LFUp = 1.0;
        double csv_weight_total_ak4LFUp = 1.0;
        double csv_weight_total_ak8LFUp = 1.0;
        double csv_weight_total_LFUp = 1.0;
        double csv_weight_BF_ak4LFDown = 1.0;
        double csv_weight_CF_ak4LFDown = 1.0;
        double csv_weight_LF_ak4LFDown = 1.0;
        double csv_weight_BF_ak8LFDown = 1.0;
        double csv_weight_CF_ak8LFDown = 1.0;
        double csv_weight_LF_ak8LFDown = 1.0;
        double csv_weight_total_ak4LFDown = 1.0;
        double csv_weight_total_ak8LFDown = 1.0;
        double csv_weight_total_LFDown = 1.0;
        
        //HF linear and quadratic statistics up and down
        double csv_weight_BF_ak4hfstats1Up = 1.0;
        double csv_weight_CF_ak4hfstats1Up = 1.0;
        double csv_weight_LF_ak4hfstats1Up = 1.0;
        double csv_weight_BF_ak8hfstats1Up = 1.0;
        double csv_weight_CF_ak8hfstats1Up = 1.0;
        double csv_weight_LF_ak8hfstats1Up = 1.0;
        double csv_weight_total_ak4hfstats1Up = 1.0;
        double csv_weight_total_ak8hfstats1Up = 1.0;
        double csv_weight_total_hfstats1Up = 1.0;
        double csv_weight_BF_ak4hfstats1Down = 1.0;
        double csv_weight_CF_ak4hfstats1Down = 1.0;
        double csv_weight_LF_ak4hfstats1Down = 1.0;
        double csv_weight_BF_ak8hfstats1Down = 1.0;
        double csv_weight_CF_ak8hfstats1Down = 1.0;
        double csv_weight_LF_ak8hfstats1Down = 1.0;
        double csv_weight_total_ak4hfstats1Down = 1.0;
        double csv_weight_total_ak8hfstats1Down = 1.0;
        double csv_weight_total_hfstats1Down = 1.0;
        double csv_weight_BF_ak4hfstats2Up = 1.0;
        double csv_weight_CF_ak4hfstats2Up = 1.0;
        double csv_weight_LF_ak4hfstats2Up = 1.0;
        double csv_weight_BF_ak8hfstats2Up = 1.0;
        double csv_weight_CF_ak8hfstats2Up = 1.0;
        double csv_weight_LF_ak8hfstats2Up = 1.0;
        double csv_weight_total_ak4hfstats2Up = 1.0;
        double csv_weight_total_ak8hfstats2Up = 1.0;
        double csv_weight_total_hfstats2Up = 1.0;
        double csv_weight_BF_ak4hfstats2Down = 1.0;
        double csv_weight_CF_ak4hfstats2Down = 1.0;
        double csv_weight_LF_ak4hfstats2Down = 1.0;
        double csv_weight_BF_ak8hfstats2Down = 1.0;
        double csv_weight_CF_ak8hfstats2Down = 1.0;
        double csv_weight_LF_ak8hfstats2Down = 1.0;
        double csv_weight_total_ak4hfstats2Down = 1.0;
        double csv_weight_total_ak8hfstats2Down = 1.0;
        double csv_weight_total_hfstats2Down = 1.0;
        
        //LF linear and quadratic statistics up and down
        double csv_weight_BF_ak4lfstats1Up = 1.0;
        double csv_weight_CF_ak4lfstats1Up = 1.0;
        double csv_weight_LF_ak4lfstats1Up = 1.0;
        double csv_weight_BF_ak8lfstats1Up = 1.0;
        double csv_weight_CF_ak8lfstats1Up = 1.0;
        double csv_weight_LF_ak8lfstats1Up = 1.0;
        double csv_weight_total_ak4lfstats1Up = 1.0;
        double csv_weight_total_ak8lfstats1Up = 1.0;
        double csv_weight_total_lfstats1Up = 1.0;
        double csv_weight_BF_ak4lfstats1Down = 1.0;
        double csv_weight_CF_ak4lfstats1Down = 1.0;
        double csv_weight_LF_ak4lfstats1Down = 1.0;
        double csv_weight_BF_ak8lfstats1Down = 1.0;
        double csv_weight_CF_ak8lfstats1Down = 1.0;
        double csv_weight_LF_ak8lfstats1Down = 1.0;
        double csv_weight_total_ak4lfstats1Down = 1.0;
        double csv_weight_total_ak8lfstats1Down = 1.0;
        double csv_weight_total_lfstats1Down = 1.0;
        double csv_weight_BF_ak4lfstats2Up = 1.0;
        double csv_weight_CF_ak4lfstats2Up = 1.0;
        double csv_weight_LF_ak4lfstats2Up = 1.0;
        double csv_weight_BF_ak8lfstats2Up = 1.0;
        double csv_weight_CF_ak8lfstats2Up = 1.0;
        double csv_weight_LF_ak8lfstats2Up = 1.0;
        double csv_weight_total_ak4lfstats2Up = 1.0;
        double csv_weight_total_ak8lfstats2Up = 1.0;
        double csv_weight_total_lfstats2Up = 1.0;
        double csv_weight_BF_ak4lfstats2Down = 1.0;
        double csv_weight_CF_ak4lfstats2Down = 1.0;
        double csv_weight_LF_ak4lfstats2Down = 1.0;
        double csv_weight_BF_ak8lfstats2Down = 1.0;
        double csv_weight_CF_ak8lfstats2Down = 1.0;
        double csv_weight_LF_ak8lfstats2Down = 1.0;
        double csv_weight_total_ak4lfstats2Down = 1.0;
        double csv_weight_total_ak8lfstats2Down = 1.0;
        double csv_weight_total_lfstats2Down = 1.0;
        
        //CF linear and quadratic statistics up and down
        double csv_weight_BF_ak4cfErr1Up = 1.0;
        double csv_weight_CF_ak4cfErr1Up = 1.0;
        double csv_weight_LF_ak4cfErr1Up = 1.0;
        double csv_weight_BF_ak8cfErr1Up = 1.0;
        double csv_weight_CF_ak8cfErr1Up = 1.0;
        double csv_weight_LF_ak8cfErr1Up = 1.0;
        double csv_weight_total_ak4cfErr1Up = 1.0;
        double csv_weight_total_ak8cfErr1Up = 1.0;
        double csv_weight_total_cfErr1Up = 1.0;
        double csv_weight_BF_ak4cfErr1Down = 1.0;
        double csv_weight_CF_ak4cfErr1Down = 1.0;
        double csv_weight_LF_ak4cfErr1Down = 1.0;
        double csv_weight_BF_ak8cfErr1Down = 1.0;
        double csv_weight_CF_ak8cfErr1Down = 1.0;
        double csv_weight_LF_ak8cfErr1Down = 1.0;
        double csv_weight_total_ak4cfErr1Down = 1.0;
        double csv_weight_total_ak8cfErr1Down = 1.0;
        double csv_weight_total_cfErr1Down = 1.0;
        double csv_weight_BF_ak4cfErr2Up = 1.0;
        double csv_weight_CF_ak4cfErr2Up = 1.0;
        double csv_weight_LF_ak4cfErr2Up = 1.0;
        double csv_weight_BF_ak8cfErr2Up = 1.0;
        double csv_weight_CF_ak8cfErr2Up = 1.0;
        double csv_weight_LF_ak8cfErr2Up = 1.0;
        double csv_weight_total_ak4cfErr2Up = 1.0;
        double csv_weight_total_ak8cfErr2Up = 1.0;
        double csv_weight_total_cfErr2Up = 1.0;
        double csv_weight_BF_ak4cfErr2Down = 1.0;
        double csv_weight_CF_ak4cfErr2Down = 1.0;
        double csv_weight_LF_ak4cfErr2Down = 1.0;
        double csv_weight_BF_ak8cfErr2Down = 1.0;
        double csv_weight_CF_ak8cfErr2Down = 1.0;
        double csv_weight_LF_ak8cfErr2Down = 1.0;
        double csv_weight_total_ak4cfErr2Down = 1.0;
        double csv_weight_total_ak8cfErr2Down = 1.0;
        double csv_weight_total_cfErr2Down = 1.0;
        
        double nL = 0;
        double nC = 0;
        double nB = 0;
         
            //loop over AK4 jets and compute scale factors depending on jet flavor, absolute value of jet eta, jet pt, jet CSV score
         
        if ( nAK4Jets > 0 ) { 
        
            for (int i = 0; i < nAK4Jets; i++) {
                
                if ( fabs(AK4JetFlavor->at(i)) == 5 ) { /*b flavor*/
                    
                    nB++;
                
                    double mycsv_weight_BF_ak4 = CSV_nominal.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4JESUp = CSV_JESUp.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4JESDown = CSV_JESDown.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4HFUp = CSV_HFUp.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4HFDown = CSV_HFDown.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4LFUp = CSV_LFUp.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4LFDown = CSV_LFDown.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4hfstats1Up = CSV_hfstats1Up.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4hfstats1Down = CSV_hfstats1Down.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4hfstats2Up = CSV_hfstats2Up.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4hfstats2Down = CSV_hfstats2Down.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4lfstats1Up = CSV_lfstats1Up.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4lfstats1Down = CSV_lfstats1Down.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4lfstats2Up = CSV_lfstats2Up.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_BF_ak4lfstats2Down = CSV_lfstats2Down.eval(BTagEntry::FLAV_B, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    //no cfErr is computed for b-flavored jets
                    
                    if (mycsv_weight_BF_ak4 != 0) csv_weight_BF_ak4 *= mycsv_weight_BF_ak4;
                    if (mycsv_weight_BF_ak4JESUp != 0) csv_weight_BF_ak4JESUp *= mycsv_weight_BF_ak4JESUp;
                    if (mycsv_weight_BF_ak4JESDown != 0) csv_weight_BF_ak4JESDown *= mycsv_weight_BF_ak4JESDown;
                    if (mycsv_weight_BF_ak4HFUp != 0) csv_weight_BF_ak4HFUp *= mycsv_weight_BF_ak4HFUp;
                    if (mycsv_weight_BF_ak4HFDown != 0) csv_weight_BF_ak4HFDown *= mycsv_weight_BF_ak4HFDown;
                    if (mycsv_weight_BF_ak4LFUp != 0) csv_weight_BF_ak4LFUp *= mycsv_weight_BF_ak4LFUp;
                    if (mycsv_weight_BF_ak4LFDown != 0) csv_weight_BF_ak4LFDown *= mycsv_weight_BF_ak4LFDown;
                    if (mycsv_weight_BF_ak4hfstats1Up != 0) csv_weight_BF_ak4hfstats1Up *= mycsv_weight_BF_ak4hfstats1Up;
                    if (mycsv_weight_BF_ak4hfstats1Down != 0) csv_weight_BF_ak4hfstats1Down *= mycsv_weight_BF_ak4hfstats1Down;
                    if (mycsv_weight_BF_ak4hfstats2Up != 0) csv_weight_BF_ak4hfstats2Up *= mycsv_weight_BF_ak4hfstats2Up;
                    if (mycsv_weight_BF_ak4hfstats2Down != 0) csv_weight_BF_ak4hfstats2Down *= mycsv_weight_BF_ak4hfstats2Down;
                    if (mycsv_weight_BF_ak4lfstats1Up != 0) csv_weight_BF_ak4lfstats1Up *= mycsv_weight_BF_ak4lfstats1Up;
                    if (mycsv_weight_BF_ak4lfstats1Down != 0) csv_weight_BF_ak4lfstats1Down *= mycsv_weight_BF_ak4lfstats1Down;
                    if (mycsv_weight_BF_ak4lfstats2Up != 0) csv_weight_BF_ak4lfstats2Up *= mycsv_weight_BF_ak4lfstats2Up;
                    if (mycsv_weight_BF_ak4lfstats2Down != 0) csv_weight_BF_ak4lfstats2Down *= mycsv_weight_BF_ak4lfstats2Down;
                }
                
                else if ( fabs(AK4JetFlavor->at(i)) == 4 ) { /*c flavor*/
                    
                    nC++;
                
                    double mycsv_weight_CF_ak4 = CSV_nominal.eval(BTagEntry::FLAV_C, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_CF_ak4cfErr1Up = CSV_cfErr1Up.eval(BTagEntry::FLAV_C, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_CF_ak4cfErr1Down = CSV_cfErr1Down.eval(BTagEntry::FLAV_C, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_CF_ak4cfErr2Up = CSV_cfErr2Up.eval(BTagEntry::FLAV_C, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_CF_ak4cfErr2Down = CSV_cfErr2Down.eval(BTagEntry::FLAV_C, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    //for c-flavored jets only the cfErr uncertainty is applied
                    
                    if (mycsv_weight_CF_ak4 != 0) csv_weight_CF_ak4 *= mycsv_weight_CF_ak4;
                    if (mycsv_weight_CF_ak4cfErr1Up != 0) csv_weight_CF_ak4cfErr1Up *= mycsv_weight_CF_ak4cfErr1Up;
                    if (mycsv_weight_CF_ak4cfErr1Down != 0) csv_weight_CF_ak4cfErr1Down *= mycsv_weight_CF_ak4cfErr1Down;
                    if (mycsv_weight_CF_ak4cfErr2Up != 0) csv_weight_CF_ak4cfErr2Up *= mycsv_weight_CF_ak4cfErr2Up;
                    if (mycsv_weight_CF_ak4cfErr2Down != 0) csv_weight_CF_ak4cfErr2Down *= mycsv_weight_CF_ak4cfErr2Down;
                
                }
                
                else  { /*light flavor*/
                    
                    nL++;
                
                    double mycsv_weight_LF_ak4 = CSV_nominal.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4JESUp = CSV_JESUp.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4JESDown = CSV_JESDown.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4HFUp = CSV_HFUp.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4HFDown = CSV_HFDown.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4LFUp = CSV_LFUp.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4LFDown = CSV_LFDown.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4hfstats1Up = CSV_hfstats1Up.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4hfstats1Down = CSV_hfstats1Down.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4hfstats2Up = CSV_hfstats2Up.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4hfstats2Down = CSV_hfstats2Down.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4lfstats1Up = CSV_lfstats1Up.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4lfstats1Down = CSV_lfstats1Down.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4lfstats2Up = CSV_lfstats2Up.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    double mycsv_weight_LF_ak4lfstats2Down = CSV_lfstats2Down.eval(BTagEntry::FLAV_UDSG, fabs(AK4JetEta->at(i)), AK4JetPt->at(i), AK4JetBtag->at(i));
                    //no cfErr is computed for light-flavored jets
                    
                    if (mycsv_weight_LF_ak4 != 0) csv_weight_LF_ak4 *= mycsv_weight_LF_ak4;
                    if (mycsv_weight_LF_ak4JESUp != 0) csv_weight_LF_ak4JESUp *= mycsv_weight_LF_ak4JESUp;
                    if (mycsv_weight_LF_ak4JESDown != 0) csv_weight_LF_ak4JESDown *= mycsv_weight_LF_ak4JESDown;
                    if (mycsv_weight_LF_ak4HFUp != 0) csv_weight_LF_ak4HFUp *= mycsv_weight_LF_ak4HFUp;
                    if (mycsv_weight_LF_ak4HFDown != 0) csv_weight_LF_ak4HFDown *= mycsv_weight_LF_ak4HFDown;
                    if (mycsv_weight_LF_ak4LFUp != 0) csv_weight_LF_ak4LFUp *= mycsv_weight_LF_ak4LFUp;
                    if (mycsv_weight_LF_ak4LFDown != 0) csv_weight_LF_ak4LFDown *= mycsv_weight_LF_ak4LFDown;
                    if (mycsv_weight_LF_ak4hfstats1Up != 0) csv_weight_LF_ak4hfstats1Up *= mycsv_weight_LF_ak4hfstats1Up;
                    if (mycsv_weight_LF_ak4hfstats1Down != 0) csv_weight_LF_ak4hfstats1Down *= mycsv_weight_LF_ak4hfstats1Down;
                    if (mycsv_weight_LF_ak4hfstats2Up != 0) csv_weight_LF_ak4hfstats2Up *= mycsv_weight_LF_ak4hfstats2Up;
                    if (mycsv_weight_LF_ak4hfstats2Down != 0) csv_weight_LF_ak4hfstats2Down *= mycsv_weight_LF_ak4hfstats2Down;
                    if (mycsv_weight_LF_ak4lfstats1Up != 0) csv_weight_LF_ak4lfstats1Up *= mycsv_weight_LF_ak4lfstats1Up;
                    if (mycsv_weight_LF_ak4lfstats1Down != 0) csv_weight_LF_ak4lfstats1Down *= mycsv_weight_LF_ak4lfstats1Down;
                    if (mycsv_weight_LF_ak4lfstats2Up != 0) csv_weight_LF_ak4lfstats2Up *= mycsv_weight_LF_ak4lfstats2Up;
                    if (mycsv_weight_LF_ak4lfstats2Down != 0) csv_weight_LF_ak4lfstats2Down *= mycsv_weight_LF_ak4lfstats2Down;
                
                }
            } // end loop over jets
        }
            
            //compute total ak4 weight
            //NOTE: all the weights from all the systematics are included. Some of them may be 1.0 if no uncertainty is present (e.g. csv_weight_BF_ak4cfErr2Down = 1.0)
            csv_weight_total_ak4 = csv_weight_BF_ak4 * csv_weight_CF_ak4 * csv_weight_LF_ak4;
            csv_weight_total_ak4JESUp = csv_weight_BF_ak4JESUp * csv_weight_CF_ak4JESUp * csv_weight_LF_ak4JESUp;
            csv_weight_total_ak4JESDown = csv_weight_BF_ak4JESDown * csv_weight_CF_ak4JESDown * csv_weight_LF_ak4JESDown;
            csv_weight_total_ak4HFUp = csv_weight_BF_ak4HFUp * csv_weight_CF_ak4HFUp * csv_weight_LF_ak4HFUp;
            csv_weight_total_ak4HFDown = csv_weight_BF_ak4HFDown * csv_weight_CF_ak4HFDown * csv_weight_LF_ak4HFDown;
            csv_weight_total_ak4LFUp = csv_weight_BF_ak4LFUp * csv_weight_CF_ak4LFUp * csv_weight_LF_ak4LFUp;
            csv_weight_total_ak4LFDown = csv_weight_BF_ak4LFDown * csv_weight_CF_ak4LFDown * csv_weight_LF_ak4LFDown;
            csv_weight_total_ak4hfstats1Up = csv_weight_BF_ak4hfstats1Up * csv_weight_CF_ak4hfstats1Up * csv_weight_LF_ak4hfstats1Up;
            csv_weight_total_ak4hfstats1Down = csv_weight_BF_ak4hfstats1Down * csv_weight_CF_ak4hfstats1Down * csv_weight_LF_ak4hfstats1Down;
            csv_weight_total_ak4hfstats2Up = csv_weight_BF_ak4hfstats2Up * csv_weight_CF_ak4hfstats2Up * csv_weight_LF_ak4hfstats2Up;
            csv_weight_total_ak4hfstats2Down = csv_weight_BF_ak4hfstats2Down * csv_weight_CF_ak4hfstats2Down * csv_weight_LF_ak4hfstats2Down;
            csv_weight_total_ak4lfstats1Up = csv_weight_BF_ak4lfstats1Up * csv_weight_CF_ak4lfstats1Up * csv_weight_LF_ak4lfstats1Up;
            csv_weight_total_ak4lfstats1Down = csv_weight_BF_ak4lfstats1Down * csv_weight_CF_ak4lfstats1Down * csv_weight_LF_ak4lfstats1Down;
            csv_weight_total_ak4lfstats2Up = csv_weight_BF_ak4lfstats2Up * csv_weight_CF_ak4lfstats2Up * csv_weight_LF_ak4lfstats2Up;
            csv_weight_total_ak4lfstats2Down = csv_weight_BF_ak4lfstats2Down * csv_weight_CF_ak4lfstats2Down * csv_weight_LF_ak4lfstats2Down;
            csv_weight_total_ak4cfErr1Up = csv_weight_BF_ak4cfErr1Up * csv_weight_CF_ak4cfErr1Up * csv_weight_LF_ak4cfErr1Up;
            csv_weight_total_ak4cfErr1Down = csv_weight_BF_ak4cfErr1Down * csv_weight_CF_ak4cfErr1Down * csv_weight_LF_ak4cfErr1Down;
            csv_weight_total_ak4cfErr2Up = csv_weight_BF_ak4cfErr2Up * csv_weight_CF_ak4cfErr2Up * csv_weight_LF_ak4cfErr2Up;
            csv_weight_total_ak4cfErr2Down = csv_weight_BF_ak4cfErr2Down * csv_weight_CF_ak4cfErr2Down * csv_weight_LF_ak4cfErr2Down;
            
            
            
            //loop over AK8 jets and compute scale factors depending on subjet flavor, absolute value of subjet eta, subjet pt, subjet CSV score
         
            //first subjet
            for (int i = 0; i < nAK8Jets; i++) {
                
                if ( fabs(Sub0Flavor->at(i)) == 5 ) { /*b flavor*/
                    
                    nB++;
                
                    double mycsv_weight_BF_ak8 = CSV_nominal.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8JESUp = CSV_JESUp.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8JESDown = CSV_JESDown.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8HFUp = CSV_HFUp.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8HFDown = CSV_HFDown.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8LFUp = CSV_LFUp.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8LFDown = CSV_LFDown.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats1Up = CSV_hfstats1Up.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats1Down = CSV_hfstats1Down.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats2Up = CSV_hfstats2Up.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats2Down = CSV_hfstats2Down.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats1Up = CSV_lfstats1Up.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats1Down = CSV_lfstats1Down.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats2Up = CSV_lfstats2Up.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats2Down = CSV_lfstats2Down.eval(BTagEntry::FLAV_B, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    //no cfErr is computed for b-flavored subjets
                    
                    if (mycsv_weight_BF_ak8 != 0) csv_weight_BF_ak8 *= mycsv_weight_BF_ak8;
                    if (mycsv_weight_BF_ak8JESUp != 0) csv_weight_BF_ak8JESUp *= mycsv_weight_BF_ak8JESUp;
                    if (mycsv_weight_BF_ak8JESDown != 0) csv_weight_BF_ak8JESDown *= mycsv_weight_BF_ak8JESDown;
                    if (mycsv_weight_BF_ak8HFUp != 0) csv_weight_BF_ak8HFUp *= mycsv_weight_BF_ak8HFUp;
                    if (mycsv_weight_BF_ak8HFDown != 0) csv_weight_BF_ak8HFDown *= mycsv_weight_BF_ak8HFDown;
                    if (mycsv_weight_BF_ak8LFUp != 0) csv_weight_BF_ak8LFUp *= mycsv_weight_BF_ak8LFUp;
                    if (mycsv_weight_BF_ak8LFDown != 0) csv_weight_BF_ak8LFDown *= mycsv_weight_BF_ak8LFDown;
                    if (mycsv_weight_BF_ak8hfstats1Up != 0) csv_weight_BF_ak8hfstats1Up *= mycsv_weight_BF_ak8hfstats1Up;
                    if (mycsv_weight_BF_ak8hfstats1Down != 0) csv_weight_BF_ak8hfstats1Down *= mycsv_weight_BF_ak8hfstats1Down;
                    if (mycsv_weight_BF_ak8hfstats2Up != 0) csv_weight_BF_ak8hfstats2Up *= mycsv_weight_BF_ak8hfstats2Up;
                    if (mycsv_weight_BF_ak8hfstats2Down != 0) csv_weight_BF_ak8hfstats2Down *= mycsv_weight_BF_ak8hfstats2Down;
                    if (mycsv_weight_BF_ak8lfstats1Up != 0) csv_weight_BF_ak8lfstats1Up *= mycsv_weight_BF_ak8lfstats1Up;
                    if (mycsv_weight_BF_ak8lfstats1Down != 0) csv_weight_BF_ak8lfstats1Down *= mycsv_weight_BF_ak8lfstats1Down;
                    if (mycsv_weight_BF_ak8lfstats2Up != 0) csv_weight_BF_ak8lfstats2Up *= mycsv_weight_BF_ak8lfstats2Up;
                    if (mycsv_weight_BF_ak8lfstats2Down != 0) csv_weight_BF_ak8lfstats2Down *= mycsv_weight_BF_ak8lfstats2Down;
                
                }
                
                else if ( fabs(Sub0Flavor->at(i)) == 4 ) { /*c flavor*/
                    
                    nC++;
                
                    double mycsv_weight_CF_ak8 = CSV_nominal.eval(BTagEntry::FLAV_C, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr1Up = CSV_cfErr1Up.eval(BTagEntry::FLAV_C, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr1Down = CSV_cfErr1Down.eval(BTagEntry::FLAV_C, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr2Up = CSV_cfErr2Up.eval(BTagEntry::FLAV_C, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr2Down = CSV_cfErr2Down.eval(BTagEntry::FLAV_C, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    //for c-flavored subjets only the cfErr uncertainty is applied
                    
                    if (mycsv_weight_CF_ak8 != 0) csv_weight_CF_ak8 *= mycsv_weight_CF_ak8;
                    if (mycsv_weight_CF_ak8cfErr1Up != 0) csv_weight_CF_ak8cfErr1Up *= mycsv_weight_CF_ak8cfErr1Up;
                    if (mycsv_weight_CF_ak8cfErr1Down != 0) csv_weight_CF_ak8cfErr1Down *= mycsv_weight_CF_ak8cfErr1Down;
                    if (mycsv_weight_CF_ak8cfErr2Up != 0) csv_weight_CF_ak8cfErr2Up *= mycsv_weight_CF_ak8cfErr2Up;
                    if (mycsv_weight_CF_ak8cfErr2Down != 0) csv_weight_CF_ak8cfErr2Down *= mycsv_weight_CF_ak8cfErr2Down;
                    
                }
                
                else  { /*light flavor*/
                    
                    nL++;
                
                    double mycsv_weight_LF_ak8 = CSV_nominal.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8JESUp = CSV_JESUp.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8JESDown = CSV_JESDown.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8HFUp = CSV_HFUp.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8HFDown = CSV_HFDown.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8LFUp = CSV_LFUp.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8LFDown = CSV_LFDown.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats1Up = CSV_hfstats1Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats1Down = CSV_hfstats1Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats2Up = CSV_hfstats2Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats2Down = CSV_hfstats2Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats1Up = CSV_lfstats1Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats1Down = CSV_lfstats1Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats2Up = CSV_lfstats2Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats2Down = CSV_lfstats2Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub0Eta->at(i)), Sub0Pt->at(i), Sub0Btag->at(i));
                    
                    if (mycsv_weight_LF_ak8 != 0) csv_weight_LF_ak8 *= mycsv_weight_LF_ak8;
                    if (mycsv_weight_LF_ak8JESUp != 0) csv_weight_LF_ak8JESUp *= mycsv_weight_LF_ak8JESUp;
                    if (mycsv_weight_LF_ak8JESDown != 0) csv_weight_LF_ak8JESDown *= mycsv_weight_LF_ak8JESDown;
                    if (mycsv_weight_LF_ak8HFUp != 0) csv_weight_LF_ak8HFUp *= mycsv_weight_LF_ak8HFUp;
                    if (mycsv_weight_LF_ak8HFDown != 0) csv_weight_LF_ak8HFDown *= mycsv_weight_LF_ak8HFDown;
                    if (mycsv_weight_LF_ak8LFUp != 0) csv_weight_LF_ak8LFUp *= mycsv_weight_LF_ak8LFUp;
                    if (mycsv_weight_LF_ak8LFDown != 0) csv_weight_LF_ak8LFDown *= mycsv_weight_LF_ak8LFDown;
                    if (mycsv_weight_LF_ak8hfstats1Up != 0) csv_weight_LF_ak8hfstats1Up *= mycsv_weight_LF_ak8hfstats1Up;
                    if (mycsv_weight_LF_ak8hfstats1Down != 0) csv_weight_LF_ak8hfstats1Down *= mycsv_weight_LF_ak8hfstats1Down;
                    if (mycsv_weight_LF_ak8hfstats2Up != 0) csv_weight_LF_ak8hfstats2Up *= mycsv_weight_LF_ak8hfstats2Up;
                    if (mycsv_weight_LF_ak8hfstats2Down != 0) csv_weight_LF_ak8hfstats2Down *= mycsv_weight_LF_ak8hfstats2Down;
                    if (mycsv_weight_LF_ak8lfstats1Up != 0) csv_weight_LF_ak8lfstats1Up *= mycsv_weight_LF_ak8lfstats1Up;
                    if (mycsv_weight_LF_ak8lfstats1Down != 0) csv_weight_LF_ak8lfstats1Down *= mycsv_weight_LF_ak8lfstats1Down;
                    if (mycsv_weight_LF_ak8lfstats2Up != 0) csv_weight_LF_ak8lfstats2Up *= mycsv_weight_LF_ak8lfstats2Up;
                    if (mycsv_weight_LF_ak8lfstats2Down != 0) csv_weight_LF_ak8lfstats2Down *= mycsv_weight_LF_ak8lfstats2Down;
                
                }
            }
            
            
            //second subjet
            for (int i = 0; i < nAK8Jets; i++) {
                
                if ( fabs(Sub1Flavor->at(i)) == 5 ) { /*b flavor*/
                    
                    nB++;
                
                    double mycsv_weight_BF_ak8 = CSV_nominal.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8JESUp = CSV_JESUp.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8JESDown = CSV_JESDown.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8HFUp = CSV_HFUp.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8HFDown = CSV_HFDown.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8LFUp = CSV_LFUp.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8LFDown = CSV_LFDown.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats1Up = CSV_hfstats1Up.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats1Down = CSV_hfstats1Down.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats2Up = CSV_hfstats2Up.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8hfstats2Down = CSV_hfstats2Down.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats1Up = CSV_lfstats1Up.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats1Down = CSV_lfstats1Down.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats2Up = CSV_lfstats2Up.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_BF_ak8lfstats2Down = CSV_lfstats2Down.eval(BTagEntry::FLAV_B, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    //no cfErr is computed for b-flavored subjets
                    
                    if (mycsv_weight_BF_ak8 != 0) csv_weight_BF_ak8 *= mycsv_weight_BF_ak8;
                    if (mycsv_weight_BF_ak8JESUp != 0) csv_weight_BF_ak8JESUp *= mycsv_weight_BF_ak8JESUp;
                    if (mycsv_weight_BF_ak8JESDown != 0) csv_weight_BF_ak8JESDown *= mycsv_weight_BF_ak8JESDown;
                    if (mycsv_weight_BF_ak8HFUp != 0) csv_weight_BF_ak8HFUp *= mycsv_weight_BF_ak8HFUp;
                    if (mycsv_weight_BF_ak8HFDown != 0) csv_weight_BF_ak8HFDown *= mycsv_weight_BF_ak8HFDown;
                    if (mycsv_weight_BF_ak8LFUp != 0) csv_weight_BF_ak8LFUp *= mycsv_weight_BF_ak8LFUp;
                    if (mycsv_weight_BF_ak8LFDown != 0) csv_weight_BF_ak8LFDown *= mycsv_weight_BF_ak8LFDown;
                    if (mycsv_weight_BF_ak8hfstats1Up != 0) csv_weight_BF_ak8hfstats1Up *= mycsv_weight_BF_ak8hfstats1Up;
                    if (mycsv_weight_BF_ak8hfstats1Down != 0) csv_weight_BF_ak8hfstats1Down *= mycsv_weight_BF_ak8hfstats1Down;
                    if (mycsv_weight_BF_ak8hfstats2Up != 0) csv_weight_BF_ak8hfstats2Up *= mycsv_weight_BF_ak8hfstats2Up;
                    if (mycsv_weight_BF_ak8hfstats2Down != 0) csv_weight_BF_ak8hfstats2Down *= mycsv_weight_BF_ak8hfstats2Down;
                    if (mycsv_weight_BF_ak8lfstats1Up != 0) csv_weight_BF_ak8lfstats1Up *= mycsv_weight_BF_ak8lfstats1Up;
                    if (mycsv_weight_BF_ak8lfstats1Down != 0) csv_weight_BF_ak8lfstats1Down *= mycsv_weight_BF_ak8lfstats1Down;
                    if (mycsv_weight_BF_ak8lfstats2Up != 0) csv_weight_BF_ak8lfstats2Up *= mycsv_weight_BF_ak8lfstats2Up;
                    if (mycsv_weight_BF_ak8lfstats2Down != 0) csv_weight_BF_ak8lfstats2Down *= mycsv_weight_BF_ak8lfstats2Down;
                
                }
                
                else if ( fabs(Sub1Flavor->at(i)) == 4 ) { /*c flavor*/
                    
                    nC++;
                                                 
                    double mycsv_weight_CF_ak8 = CSV_nominal.eval(BTagEntry::FLAV_C, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr1Up = CSV_cfErr1Up.eval(BTagEntry::FLAV_C, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr1Down = CSV_cfErr1Down.eval(BTagEntry::FLAV_C, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr2Up = CSV_cfErr2Up.eval(BTagEntry::FLAV_C, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_CF_ak8cfErr2Down = CSV_cfErr2Down.eval(BTagEntry::FLAV_C, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    //for c-flavored subjets only the cfErr uncertainty is applied
                    
                    if (mycsv_weight_CF_ak8 != 0) csv_weight_CF_ak8 *= mycsv_weight_CF_ak8;
                    if (mycsv_weight_CF_ak8cfErr1Up != 0) csv_weight_CF_ak8cfErr1Up *= mycsv_weight_CF_ak8cfErr1Up;
                    if (mycsv_weight_CF_ak8cfErr1Down != 0) csv_weight_CF_ak8cfErr1Down *= mycsv_weight_CF_ak8cfErr1Down;
                    if (mycsv_weight_CF_ak8cfErr2Up != 0) csv_weight_CF_ak8cfErr2Up *= mycsv_weight_CF_ak8cfErr2Up;
                    if (mycsv_weight_CF_ak8cfErr2Down != 0) csv_weight_CF_ak8cfErr2Down *= mycsv_weight_CF_ak8cfErr2Down;
                
                }
                
                else  { /*light flavor*/
                    
                    nL++;
                
                    double mycsv_weight_LF_ak8 = CSV_nominal.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8JESUp = CSV_JESUp.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8JESDown = CSV_JESDown.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8HFUp = CSV_HFUp.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8HFDown = CSV_HFDown.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8LFUp = CSV_LFUp.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8LFDown = CSV_LFDown.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats1Up = CSV_hfstats1Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats1Down = CSV_hfstats1Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats2Up = CSV_hfstats2Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8hfstats2Down = CSV_hfstats2Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats1Up = CSV_lfstats1Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats1Down = CSV_lfstats1Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats2Up = CSV_lfstats2Up.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    double mycsv_weight_LF_ak8lfstats2Down = CSV_lfstats2Down.eval(BTagEntry::FLAV_UDSG, fabs(Sub1Eta->at(i)), Sub1Pt->at(i), Sub1Btag->at(i));
                    
                    if (mycsv_weight_LF_ak8 != 0) csv_weight_LF_ak8 *= mycsv_weight_LF_ak8;
                    if (mycsv_weight_LF_ak8JESUp != 0) csv_weight_LF_ak8JESUp *= mycsv_weight_LF_ak8JESUp;
                    if (mycsv_weight_LF_ak8JESDown != 0) csv_weight_LF_ak8JESDown *= mycsv_weight_LF_ak8JESDown;
                    if (mycsv_weight_LF_ak8HFUp != 0) csv_weight_LF_ak8HFUp *= mycsv_weight_LF_ak8HFUp;
                    if (mycsv_weight_LF_ak8HFDown != 0) csv_weight_LF_ak8HFDown *= mycsv_weight_LF_ak8HFDown;
                    if (mycsv_weight_LF_ak8LFUp != 0) csv_weight_LF_ak8LFUp *= mycsv_weight_LF_ak8LFUp;
                    if (mycsv_weight_LF_ak8LFDown != 0) csv_weight_LF_ak8LFDown *= mycsv_weight_LF_ak8LFDown;
                    if (mycsv_weight_LF_ak8hfstats1Up != 0) csv_weight_LF_ak8hfstats1Up *= mycsv_weight_LF_ak8hfstats1Up;
                    if (mycsv_weight_LF_ak8hfstats1Down != 0) csv_weight_LF_ak8hfstats1Down *= mycsv_weight_LF_ak8hfstats1Down;
                    if (mycsv_weight_LF_ak8hfstats2Up != 0) csv_weight_LF_ak8hfstats2Up *= mycsv_weight_LF_ak8hfstats2Up;
                    if (mycsv_weight_LF_ak8hfstats2Down != 0) csv_weight_LF_ak8hfstats2Down *= mycsv_weight_LF_ak8hfstats2Down;
                    if (mycsv_weight_LF_ak8lfstats1Up != 0) csv_weight_LF_ak8lfstats1Up *= mycsv_weight_LF_ak8lfstats1Up;
                    if (mycsv_weight_LF_ak8lfstats1Down != 0) csv_weight_LF_ak8lfstats1Down *= mycsv_weight_LF_ak8lfstats1Down;
                    if (mycsv_weight_LF_ak8lfstats2Up != 0) csv_weight_LF_ak8lfstats2Up *= mycsv_weight_LF_ak8lfstats2Up;
                    if (mycsv_weight_LF_ak8lfstats2Down != 0) csv_weight_LF_ak8lfstats2Down *= mycsv_weight_LF_ak8lfstats2Down;
                
                }
                
                }
            
            
            //compute total ak8 weight
            csv_weight_total_ak8 = csv_weight_BF_ak8 * csv_weight_CF_ak8 * csv_weight_LF_ak8;
            csv_weight_total_ak8JESUp = csv_weight_BF_ak8JESUp * csv_weight_CF_ak8JESUp * csv_weight_LF_ak8JESUp;
            csv_weight_total_ak8JESDown = csv_weight_BF_ak8JESDown * csv_weight_CF_ak8JESDown * csv_weight_LF_ak8JESDown;
            csv_weight_total_ak8HFUp = csv_weight_BF_ak8HFUp * csv_weight_CF_ak8HFUp * csv_weight_LF_ak8HFUp;
            csv_weight_total_ak8HFDown = csv_weight_BF_ak8HFDown * csv_weight_CF_ak8HFDown * csv_weight_LF_ak8HFDown;
            csv_weight_total_ak8LFUp = csv_weight_BF_ak8LFUp * csv_weight_CF_ak8LFUp * csv_weight_LF_ak8LFUp;
            csv_weight_total_ak8LFDown = csv_weight_BF_ak8LFDown * csv_weight_CF_ak8LFDown * csv_weight_LF_ak8LFDown;
            csv_weight_total_ak8hfstats1Up = csv_weight_BF_ak8hfstats1Up * csv_weight_CF_ak8hfstats1Up * csv_weight_LF_ak8hfstats1Up;
            csv_weight_total_ak8hfstats1Down = csv_weight_BF_ak8hfstats1Down * csv_weight_CF_ak8hfstats1Down * csv_weight_LF_ak8hfstats1Down;
            csv_weight_total_ak8hfstats2Up = csv_weight_BF_ak8hfstats2Up * csv_weight_CF_ak8hfstats2Up * csv_weight_LF_ak8hfstats2Up;
            csv_weight_total_ak8hfstats2Down = csv_weight_BF_ak8hfstats2Down * csv_weight_CF_ak8hfstats2Down * csv_weight_LF_ak8hfstats2Down;
            csv_weight_total_ak8lfstats1Up = csv_weight_BF_ak8lfstats1Up * csv_weight_CF_ak8lfstats1Up * csv_weight_LF_ak8lfstats1Up;
            csv_weight_total_ak8lfstats1Down = csv_weight_BF_ak8lfstats1Down * csv_weight_CF_ak8lfstats1Down * csv_weight_LF_ak8lfstats1Down;
            csv_weight_total_ak8lfstats2Up = csv_weight_BF_ak8lfstats2Up * csv_weight_CF_ak8lfstats2Up * csv_weight_LF_ak8lfstats2Up;
            csv_weight_total_ak8lfstats2Down = csv_weight_BF_ak8lfstats2Down * csv_weight_CF_ak8lfstats2Down * csv_weight_LF_ak8lfstats2Down;
            csv_weight_total_ak8cfErr1Up = csv_weight_BF_ak8cfErr1Up * csv_weight_CF_ak8cfErr1Up * csv_weight_LF_ak8cfErr1Up;
            csv_weight_total_ak8cfErr1Down = csv_weight_BF_ak8cfErr1Down * csv_weight_CF_ak8cfErr1Down * csv_weight_LF_ak8cfErr1Down;
            csv_weight_total_ak8cfErr2Up = csv_weight_BF_ak8cfErr2Up * csv_weight_CF_ak8cfErr2Up * csv_weight_LF_ak8cfErr2Up;
            csv_weight_total_ak8cfErr2Down = csv_weight_BF_ak8cfErr2Down * csv_weight_CF_ak8cfErr2Down * csv_weight_LF_ak8cfErr2Down;
            
            //compute total weight
            csv_weight_total = csv_weight_total_ak4 * csv_weight_total_ak8;
            csv_weight_total_JESUp = csv_weight_total_ak4JESUp * csv_weight_total_ak8JESUp;
            csv_weight_total_JESDown = csv_weight_total_ak4JESDown * csv_weight_total_ak8JESDown;
            csv_weight_total_HFUp = csv_weight_total_ak4HFUp * csv_weight_total_ak8HFUp;
            csv_weight_total_HFDown = csv_weight_total_ak4HFDown * csv_weight_total_ak8HFDown;
            csv_weight_total_LFUp = csv_weight_total_ak4LFUp * csv_weight_total_ak8LFUp;
            csv_weight_total_LFDown = csv_weight_total_ak4LFDown * csv_weight_total_ak8LFDown;
            csv_weight_total_hfstats1Up = csv_weight_total_ak4hfstats1Up * csv_weight_total_ak8hfstats1Up;
            csv_weight_total_hfstats1Down = csv_weight_total_ak4hfstats1Down * csv_weight_total_ak8hfstats1Down;
            csv_weight_total_hfstats2Up = csv_weight_total_ak4hfstats2Up * csv_weight_total_ak8hfstats2Up;
            csv_weight_total_hfstats2Down = csv_weight_total_ak4hfstats2Down * csv_weight_total_ak8hfstats2Down;
            csv_weight_total_lfstats1Up = csv_weight_total_ak4lfstats1Up * csv_weight_total_ak8lfstats1Up;
            csv_weight_total_lfstats1Down = csv_weight_total_ak4lfstats1Down * csv_weight_total_ak8lfstats1Down;
            csv_weight_total_lfstats2Up = csv_weight_total_ak4lfstats2Up * csv_weight_total_ak8lfstats2Up;
            csv_weight_total_lfstats2Down = csv_weight_total_ak4lfstats2Down * csv_weight_total_ak8lfstats2Down;
            csv_weight_total_cfErr1Up = csv_weight_total_ak4cfErr1Up * csv_weight_total_ak8cfErr1Up;
            csv_weight_total_cfErr1Down = csv_weight_total_ak4cfErr1Down * csv_weight_total_ak8cfErr1Down;
            csv_weight_total_cfErr2Up = csv_weight_total_ak4cfErr2Up * csv_weight_total_ak8cfErr2Up;
            csv_weight_total_cfErr2Down = csv_weight_total_ak4cfErr2Down * csv_weight_total_ak8cfErr2Down;

            output[0] = csv_weight_total;
            output[1] = csv_weight_total_JESUp;
            output[2] = csv_weight_total_JESDown;
            output[3] = csv_weight_total_HFUp;
            output[4] = csv_weight_total_HFDown;
            output[5] = csv_weight_total_LFUp;
            output[6] = csv_weight_total_LFDown;
            output[7] = csv_weight_total_hfstats1Up;
            output[8] = csv_weight_total_hfstats1Down;
            output[9] = csv_weight_total_hfstats2Up;
            output[10] = csv_weight_total_hfstats2Down;
            output[11] = csv_weight_total_lfstats1Up;
            output[12] = csv_weight_total_lfstats1Down;
            output[13] = csv_weight_total_lfstats2Up;
            output[14] = csv_weight_total_lfstats2Down;
            output[15] = csv_weight_total_cfErr1Up;
            output[16] = csv_weight_total_cfErr1Down;
            output[17] = csv_weight_total_cfErr2Up;
            output[18] = csv_weight_total_cfErr2Down;
            output[19] = nL;
            output[20] = nC;
            output[21] = nB;
            
            return output;
            delete output;
            delete AK4JetFlavor;
            delete AK4JetEta;
            delete AK4JetEta;
            delete AK4JetBtag;
            delete Sub0Eta;
            delete Sub0Pt;
            delete Sub0Btag;
            delete Sub0Flavor;
            delete Sub1Eta;
            delete Sub1Pt;
            delete Sub1Btag;
            delete Sub1Flavor;
            
} 
