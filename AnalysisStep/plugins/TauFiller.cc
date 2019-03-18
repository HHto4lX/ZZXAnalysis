/** \class TauFiller
 *
 *  Preselect taus from all packedPFCandidates
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
#include <DataFormats/PatCandidates/interface/Tau.h>
//#include <DataFormats/PatCandidates/interface/PFParticle.h>
//#include <ZZXAnalysis/AnalysisStep/interface/TauFwd.h>
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

class TauFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit TauFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~TauFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<edm::View<pat::Tau> > tauToken;
  //edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfCandToken;
  int selectionMode;
  int sampleType;
  int setup;
  bool debug;
};


TauFiller::TauFiller(const edm::ParameterSet& iConfig) :
  tauToken(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("pikatauSrc"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  debug(iConfig.getUntrackedParameter<bool>("debug",false))
{
  
  //pfCandToken = consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 

  produces<pat::TauCollection>(); 
}


void
TauFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //--- get the tau cand
//  edm::Handle<edm::View<pat::Tau> > tauHandle;
  edm::Handle<pat::TauRefVector> tauHandle;
  iEvent.getByToken(tauToken, tauHandle);

  //--- Get the PF cands
  // use pf candidates and then select tau among these
  //edm::Handle<edm::View<pat::PackedCandidate> > pfCands; 
  //iEvent.getByToken(pfCandToken, pfCands);

  // Output collections
  auto result = std::make_unique<pat::TauCollection>();



  //----------------------
  // Loop on taus
  //----------------------
    for (unsigned int itau = 0; itau < tauHandle->size(); ++itau){
    //---Clone the pat::Tau
    pat::Tau g(*((*tauHandle)[itau].get()));
    
    double g_pt  = g.pt();
    double g_eta = g.eta();
    //double g_phi = g.phi();
    
    // We only want taus: select taus among all the pf candidates
    //if (g->pdgId()!=22) continue;

    // Tau preselection 
    if (!(g_pt>2. && fabs(g_eta)<2.4)) continue;

    
    
    result->push_back(g); 
    //result->back().setStatus(0);

  } // end of loop over tau collection



  //--- Reorder taus by pT
  //std::sort(result->begin(),result->end(), [](const Tau& g1, const Tau& g2){
  //  return g1.pt()>g2.pt();
  //});

  
  //Put the result in the event
  iEvent.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauFiller);

