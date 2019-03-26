/** \class TauFiller_basic
 *
 *
 *  $Date: 2019/03/04 $
 *  \author A. Cappati (Torino)
 *  \author B. Kiani   (Torino)
 *
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <DataFormats/PatCandidates/interface/Tau.h>
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

class TauFiller_basic : public edm::EDProducer {
 public:
  /// Constructor
  explicit TauFiller_basic(const edm::ParameterSet&);
    
  /// Destructor
  ~TauFiller_basic(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<edm::View<pat::Tau> > tauToken;
  int selectionMode;
  int sampleType;
  int setup;
  bool debug;
};


TauFiller_basic::TauFiller_basic(const edm::ParameterSet& iConfig) :
  tauToken(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("pikatauSrc"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  debug(iConfig.getUntrackedParameter<bool>("debug",false))
{
  
  produces<pat::TauCollection>(); 
}


void
TauFiller_basic::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //--- get the tau cand
  edm::Handle<edm::View<pat::Tau> > tauHandle;
  //edm::Handle<pat::TauRefVector> tauHandle;
  iEvent.getByToken(tauToken, tauHandle);

  
  // Output collections
  auto result = std::make_unique<pat::TauCollection>();



  //----------------------
  // Loop on taus
  //----------------------
  for (auto itau = tauHandle->begin(); itau != tauHandle->end(); ++itau){


    pat::Tau t(*itau);
    double t_pt  = t.pt();
    double t_eta = t.eta();
    //double t_phi = t.phi();
    

    // Tau preselection 
    if (!(t_pt>0.5 && fabs(t_eta)<2.4)) continue;

    
    
    result->push_back(t); 
    //result->back().setStatus(0);

  } // end of loop over tau collection



  //--- Reorder taus by pT
  //    object type is BaseTau from DataFormats/PatCandidates/interface/Tau.h class definition:  class Tau : public Lepton<reco::BaseTau> 
  std::sort(result->begin(),result->end(), [](const BaseTau& t1, const BaseTau& t2){    
   return t1.pt()>t2.pt();
  });

  
  //Put the result in the event
  iEvent.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauFiller_basic);

