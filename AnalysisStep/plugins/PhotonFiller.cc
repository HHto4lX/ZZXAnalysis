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
    

    // Get the candidate as edm::Ptr
    //edm::Ptr<pat::PackedCandidate> g = pfCands->ptrAt(i);
    //edm::Ptr<pat::Photon> g = photonHandle->ptrAt(i);

    pat::Photon g(*pht_i);
    double g_pt  = g.pt();
    double g_eta = g.eta();
    double g_phi = g.phi();
    
    // We only want photons: select photons among all the pf candidates
    //if (g->pdgId()!=22) continue;

    // Photon preselection 
    if (!(g_pt>2. && fabs(g_eta)<2.4)) continue;

    
    
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

