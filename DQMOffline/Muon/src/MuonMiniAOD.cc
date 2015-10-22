#include "DQMOffline/Muon/interface/MuonMiniAOD.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>
#include "TMath.h"
using namespace std;
using namespace edm;



MuonMiniAOD::MuonMiniAOD(const edm::ParameterSet& pSet) {
  parameters = pSet;

  // the services:
  
  theMuonCollectionLabel_ = consumes<edm::View<pat::Muon> >  (parameters.getParameter<edm::InputTag>("MuonCollection"));
  theVertexLabel_          = consumes<reco::VertexCollection>(parameters.getParameter<edm::InputTag>("VertexLabel"));
  theBeamSpotLabel_        = mayConsume<reco::BeamSpot>      (parameters.getParameter<edm::InputTag>("BeamSpotLabel"));

}


MuonMiniAOD::~MuonMiniAOD() { }
void MuonMiniAOD::bookHistograms(DQMStore::IBooker & ibooker,
				      edm::Run const & /*iRun*/,
				      edm::EventSetup const & /* iSetup */){
    
  ibooker.cd();
  ibooker.setCurrentFolder("Muons_miniAOD/MuonMiniAOD");
  
  tightMuons = ibooker.book1D("tightMuons","Tight Muons",2,1,3);
  tightMuons -> setBinLabel(1,"OK");
  tightMuons -> setBinLabel(2,"No OK");

  mediumMuons = ibooker.book1D("mediumMuons","Medium Muons",2,1,3);
  mediumMuons -> setBinLabel(1,"OK");
  mediumMuons -> setBinLabel(2,"No OK");

  looseMuons = ibooker.book1D("looseMuons","Loose Muons",2,1,3);
  looseMuons -> setBinLabel(1,"OK");
  looseMuons -> setBinLabel(2,"No OK");

  softMuons = ibooker.book1D("softMuons","Soft Muons",2,1,3);
  softMuons -> setBinLabel(1,"OK");
  softMuons -> setBinLabel(2,"No OK");

  highPtMuons = ibooker.book1D("highPtMuons","High Pt Muons",2,1,3);
  highPtMuons -> setBinLabel(1,"OK");
  highPtMuons -> setBinLabel(2,"No OK");

}


void MuonMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  LogTrace(metname)<<"[MuonMiniAOD] Analyze the mu";
  
 
  // Take the muon container
  edm::Handle<edm::View<pat::Muon> > muons; 
  iEvent.getByToken(theMuonCollectionLabel_,muons);
 
  //Vertex information
  edm::Handle<reco::VertexCollection> vertex;
  iEvent.getByToken(theVertexLabel_, vertex);  


  if(!muons.isValid()) return;

  reco::Vertex::Point posVtx;
  reco::Vertex::Error errVtx;
  unsigned int theIndexOfThePrimaryVertex = 999.;
  if (!vertex.isValid()) {
    LogTrace(metname) << "[EfficiencyAnalyzer] Could not find vertex collection" << std::endl;
    for (unsigned int ind=0; ind<vertex->size(); ++ind) {
      if ( (*vertex)[ind].isValid() && !((*vertex)[ind].isFake()) ) {
	theIndexOfThePrimaryVertex = ind;
	break;
      }
    }
  }

  if (theIndexOfThePrimaryVertex<100) {
    posVtx = ((*vertex)[theIndexOfThePrimaryVertex]).position();
    errVtx = ((*vertex)[theIndexOfThePrimaryVertex]).error();
  }   
  else {
    LogInfo("RecoMuonValidator") << "reco::PrimaryVertex not found, use BeamSpot position instead\n";
    
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken(theBeamSpotLabel_,recoBeamSpotHandle);
    reco::BeamSpot bs = *recoBeamSpotHandle;
    
    posVtx = bs.position();
    errVtx(0,0) = bs.BeamWidthX();
    errVtx(1,1) = bs.BeamWidthY();
    errVtx(2,2) = bs.sigmaZ();
  }
    


  const reco::Vertex thePrimaryVertex(posVtx,errVtx);

  for (edm::View<pat::Muon>::const_iterator muon1 = muons->begin(); muon1 != muons->end(); ++muon1){
    if (muon::isTightMuon(*muon1,thePrimaryVertex) == muon1->isTightMuon(thePrimaryVertex))
      tightMuons->Fill(1);
    else
      tightMuons->Fill(2);
  
    if (muon::isMediumMuon(*muon1) == muon1->isMediumMuon())
      mediumMuons->Fill(1);
    else
      mediumMuons->Fill(2);

    if (muon::isLooseMuon(*muon1) == muon1->isLooseMuon())
      looseMuons->Fill(1);
    else
      looseMuons->Fill(2);
    
    if (muon::isSoftMuon(*muon1,thePrimaryVertex) == muon1->isSoftMuon(thePrimaryVertex))
      softMuons->Fill(1);
    else
      softMuons->Fill(2);
    
    if (muon::isHighPtMuon(*muon1,thePrimaryVertex) == muon1->isHighPtMuon(thePrimaryVertex))
      highPtMuons->Fill(1);
    else
      highPtMuons->Fill(2);
  }
}
