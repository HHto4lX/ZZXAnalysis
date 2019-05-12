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
  //edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfCandToken;
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
  
  //pfCandToken = consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 

  produces<pat::PhotonCollection>(); 
}


void
PhotonFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //--- get the photon cand
  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByToken(photonToken, photonHandle);

  //--- Get the PF cands
  // use pf candidates and then select photon among these
  //edm::Handle<edm::View<pat::PackedCandidate> > pfCands; 
  //iEvent.getByToken(pfCandToken, pfCands);

  // Output collections
  auto result = std::make_unique<pat::PhotonCollection>();



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
   
    
    // -------------
    // Photon LOOSE selection
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2 
    bool photonID_loose = false;
    
    // barrel photons
    if( fabs(g_etaSC) <= 1.479 ){

      if( g_HoverE           < 0.04596 &&
          g_sigmaIEtaIEta    < 0.0106  &&
          g_chargedHadronIso < 1.694   &&
          g_neutralHadronIso < 24.032 + 0.01512*g_pt + 2.259e-05*g_pt*g_pt &&
          g_photonIso        < 2.876  + 0.004017*g_pt 
        )
      {
	photonID_loose = true;
      }

    }
    // endcap photons
    else if( fabs(g_etaSC) > 1.479  &&  fabs(g_etaSC) <= 2.8 ){

      if( g_HoverE           < 0.0590 &&
          g_sigmaIEtaIEta    < 0.0272 &&
          g_chargedHadronIso < 2.089  &&
          g_neutralHadronIso < 19.722 + 0.0117*g_pt + 2.3e-05*g_pt*g_pt &&
          g_photonIso        < 4.162  + 0.0037*g_pt 
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

