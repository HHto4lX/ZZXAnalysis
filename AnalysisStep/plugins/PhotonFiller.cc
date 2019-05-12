/** \class PhotonFiller
 *
 *  Preselect photons from all packedPFCandidates
 *
 *  $Date: 2019/03/04 $
 *  \author A. Cappati (Torino)
 *  \author B. Kiani   (Torino)
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
//#include <DataFormats/PatCandidates/interface/PFParticle.h>
//#include <ZZXAnalysis/AnalysisStep/interface/PhotonFwd.h>
//#include <EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h>
//#include <DataFormats/GeometryVector/interface/VectorUtil.h> 
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include <ZZXAnalysis/AnalysisStep/interface/CutSet.h>

//#include <ZZXAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <Math/VectorUtil.h>
#include <TMath.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

class PhotonFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit PhotonFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~PhotonFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<edm::View<pat::Photon> > photonToken;
  edm::EDGetTokenT<double> rhoToken;
  int selectionMode;
  int sampleType;
  int setup;
  bool debug;
};


PhotonFiller::PhotonFiller(const edm::ParameterSet& iConfig) :
  photonToken(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("pikaphotonSrc"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  debug(iConfig.getUntrackedParameter<bool>("debug",false))
{
  
  rhoToken = consumes<double> (edm::InputTag("fixedGridRhoFastjetAll")); 

  produces<pat::PhotonCollection>(); 
}


void
PhotonFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //--- get the photon cand
  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByToken(photonToken, photonHandle);

 
  // Output collections
  auto result = std::make_unique<pat::PhotonCollection>();




  //--- Get the fixedGridRhoFastjetAll value 
  //    (for isolation ID criteria)
  edm::Handle<double> rhoHandle; 
  iEvent.getByToken(rhoToken, rhoHandle);
  double rho = *rhoHandle;


  //--- Define Effective Area vectors
  //    (for isolation ID criteria) - 2017 v2 ID
  //    https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Offline_selection_criteria_AN1 
  //
  //    values for different eta regions: 
  //    |eta|<1, 1<|eta|<1.479, 1.479<|eta|<2, 2<|eta|<2.2, 2.2<|eta|<2.3, 2.3<|eta|<2.4, |eta|>2.4
  double EA_chargedHadrons[7] = {0.0112, 0.0108, 0.0106, 0.01002, 0.0098, 0.0089, 0.0087};
  double EA_neutralHadrons[7] = {0.0668, 0.1054, 0.0786, 0.0233,  0.0078, 0.0028, 0.0137};
  double EA_photons[7]        = {0.1113, 0.0953, 0.0619, 0.0837,  0.1070, 0.1212, 0.1466};





  //----------------------
  // Loop on photons
  //----------------------
  for (auto pht_i = photonHandle->begin(); pht_i != photonHandle->end(); ++pht_i) {
    

    pat::Photon g(*pht_i);
    float g_pt     = g.pt();
    float g_etaSC  = g.superCluster()->eta();
    float g_HoverE = g.hadronicOverEm();
    float g_sigmaIEtaIEta    = g.full5x5_sigmaIetaIeta();
    float g_chargedHadronIso = g.chargedHadronIso();
    float g_neutralHadronIso = g.neutralHadronIso();
    float g_photonIso        = g.photonIso();


    // compute iso rho corrected
    float g_chargedHadronIso_corr = 0.;
    float g_neutralHadronIso_corr = 0.;
    float g_photonIso_corr        = 0.;

    if( fabs(g_etaSC) <= 1. )
    {
      g_chargedHadronIso_corr = max(g_chargedHadronIso - rho*EA_chargedHadrons[0], 0.);
      g_neutralHadronIso_corr = max(g_neutralHadronIso - rho*EA_neutralHadrons[0], 0.);
      g_photonIso_corr        = max(g_photonIso        - rho*EA_photons[0],        0.);
    }
    else if( fabs(g_etaSC) > 1. && fabs(g_etaSC) <= 1.479 )
    {
      g_chargedHadronIso_corr = max(g_chargedHadronIso - rho*EA_chargedHadrons[1], 0.);
      g_neutralHadronIso_corr = max(g_neutralHadronIso - rho*EA_neutralHadrons[1], 0.);
      g_photonIso_corr        = max(g_photonIso        - rho*EA_photons[1],        0.);
    }
    else if( fabs(g_etaSC) > 1.479 && fabs(g_etaSC) <= 2. )
    {
      g_chargedHadronIso_corr = max(g_chargedHadronIso - rho*EA_chargedHadrons[2], 0.);
      g_neutralHadronIso_corr = max(g_neutralHadronIso - rho*EA_neutralHadrons[2], 0.);
      g_photonIso_corr        = max(g_photonIso        - rho*EA_photons[2],        0.);
    }
    else if( fabs(g_etaSC) > 2. && fabs(g_etaSC) <= 2.2 )
    {
      g_chargedHadronIso_corr = max(g_chargedHadronIso - rho*EA_chargedHadrons[3], 0.);
      g_neutralHadronIso_corr = max(g_neutralHadronIso - rho*EA_neutralHadrons[3], 0.);
      g_photonIso_corr        = max(g_photonIso        - rho*EA_photons[3],        0.);
    }
    else if( fabs(g_etaSC) > 2.2 && fabs(g_etaSC) <= 2.3 )
    {
      g_chargedHadronIso_corr = max(g_chargedHadronIso - rho*EA_chargedHadrons[4], 0.);
      g_neutralHadronIso_corr = max(g_neutralHadronIso - rho*EA_neutralHadrons[4], 0.);
      g_photonIso_corr        = max(g_photonIso        - rho*EA_photons[4],        0.);
    }
    else if( fabs(g_etaSC) > 2.3 && fabs(g_etaSC) <= 2.4 )
    {
      g_chargedHadronIso_corr = max(g_chargedHadronIso - rho*EA_chargedHadrons[5], 0.);
      g_neutralHadronIso_corr = max(g_neutralHadronIso - rho*EA_neutralHadrons[5], 0.);
      g_photonIso_corr        = max(g_photonIso        - rho*EA_photons[5],        0.);
    }
    else if( fabs(g_etaSC) > 2.4 )
    {
      g_chargedHadronIso_corr = max(g_chargedHadronIso - rho*EA_chargedHadrons[6], 0.);
      g_neutralHadronIso_corr = max(g_neutralHadronIso - rho*EA_neutralHadrons[6], 0.);
      g_photonIso_corr        = max(g_photonIso        - rho*EA_photons[6],        0.);
    }
    else continue;


    
    
    // -------------
    // Photon LOOSE selection
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Offline_selection_criteria_AN1
    bool photonID_loose = false;
    
    // barrel photons
    if( fabs(g_etaSC) <= 1.479 ){

      if( g_HoverE                < 0.04596 &&
          g_sigmaIEtaIEta         < 0.0106  &&
          g_chargedHadronIso_corr < 1.694   &&
          g_neutralHadronIso_corr < 24.032 + 0.01512*g_pt + 2.259e-05*g_pt*g_pt &&
          g_photonIso_corr        < 2.876  + 0.004017*g_pt 
        )
      {
	photonID_loose = true;
      }

    }
    // endcap photons
    else if( fabs(g_etaSC) > 1.479  &&  fabs(g_etaSC) <= 2.8 ){

      if( g_HoverE                < 0.0590 &&
          g_sigmaIEtaIEta         < 0.0272 &&
          g_chargedHadronIso_corr < 2.089  &&
          g_neutralHadronIso_corr < 19.722 + 0.0117*g_pt + 2.3e-05*g_pt*g_pt &&
          g_photonIso_corr        < 4.162  + 0.0037*g_pt 
        )
      {
	photonID_loose = true;
      }

    }
    // drop photons with eta>2.8
    else { continue; }

    // -------------




    // add variable to the photon collection
    g.addUserFloat("photonID_loose",photonID_loose);
    g.addUserFloat("photon_rho",    rho );
    g.addUserFloat("photon_chargedHadronIso_corr", g_chargedHadronIso_corr );
    g.addUserFloat("photon_neutralHadronIso_corr", g_neutralHadronIso_corr );
    g.addUserFloat("photon_photonIso_corr",        g_photonIso_corr );


    
    result->push_back(g); 
    //result->back().setStatus(0);

  } // end of loop over photon collection



  //--- Reorder photons by pT
  std::sort(result->begin(),result->end(), [](const Photon& g1, const Photon& g2){
    return g1.pt()>g2.pt();
  });

  
  //Put the result in the event
  iEvent.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(PhotonFiller);

