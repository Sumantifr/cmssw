// -*- C++ -*-
//
// Package:    PhysicsTools/PatAlgos
// Class:      LeptonTagInfoCollectionProducer
//
/**\class LeptonTagInfoCollectionProducer LeptonTagInfoCollectionProducer.cc PhysicsTools/PatAlgos/plugins/PNETLeptonProducer.cc


*/
//
// Original Author:  Sergio Sanchez Cruz
//         Created:  Mon, 15 May 2023 08:32:03 GMT
//
//

#include"PhysicsTools/PatAlgos/interface/LeptonTagInfoCollectionProducer.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/NanoAOD/interface/MatchingUtils.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"


namespace pat {
  template <typename T> LeptonTagInfoCollectionProducer<T>::LeptonTagInfoCollectionProducer(const edm::ParameterSet &iConfig) : 
    src_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("src"))),
    sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
    pv_token_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    lepton_varsPSet_(iConfig.getParameter<edm::ParameterSet>("leptonVars")),
    pf_varsPSet_(iConfig.getParameter<edm::ParameterSet>("pfVars")),
    sv_varsPSet_(iConfig.getParameter<edm::ParameterSet>("svVars"))
  {
    std::cout << "Check1" << std::endl;
    produces<LeptonTagInfoCollection<T>>();
    parse_vars_into(lepton_varsPSet_, lepton_vars_);
    parse_vars_into(pf_varsPSet_    , pf_vars_);
    parse_vars_into(sv_varsPSet_    , sv_vars_);
    
  }
  
  template <typename T>
  void LeptonTagInfoCollectionProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
  {
    std::cout << "Check2" << std::endl;
    auto src = iEvent.getHandle(src_token_);
    iEvent.getByToken(sv_token_, svs_);
    iEvent.getByToken(pv_token_, pvs_);

    auto output_info = std::make_unique<LeptonTagInfoCollection<T>>();
    
    
    for (size_t ilep=0; ilep < src->size(); ilep++){
      const auto &lep = (*src)[ilep];
      edm::RefToBase<T> lep_ref(src, ilep);
      btagbtvdeep::DeepBoostedJetFeatures features;
      std::cout << "Filling leptons" << std::endl;
      fill_lepton_features( lep, features );
      std::cout << "Filling pf" << std::endl;
      fill_pf_features( lep, features );
      std::cout << "Filling svs" << std::endl;
      fill_sv_features( lep, features );


      output_info->emplace_back(features, lep_ref);

    }
    iEvent.put(std::move(output_info));
    
  }

  template <typename T>
  template <typename T2>
  void LeptonTagInfoCollectionProducer<T>::parse_vars_into(const edm::ParameterSet &varsPSet, std::vector<std::unique_ptr<varWithName<T2>>>& vars)
  {
    
    for (const std::string &vname : varsPSet.getParameterNamesForType<std::string>()) {
      const std::string &func = varsPSet.getParameter<std::string>(vname);
      std::cout << "Parsing var " << vname << func << std::endl;
      vars.push_back(std::make_unique<varWithName<T2>>(vname, StringObjectFunction<T2, true>(func)));
    }
  }
    
  template <typename T>
  void LeptonTagInfoCollectionProducer<T>::fill_lepton_features(const T& lep, btagbtvdeep::DeepBoostedJetFeatures& features){
    for (auto& var : lepton_vars_){
      std::cout << "Fillin var " << var->first << " " << var->second(lep) << std::endl;
      features.add(var->first);
      features.reserve(var->first,1);
      features.fill(var->first, var->second(lep));
    }
  }

  template <typename T>
  void LeptonTagInfoCollectionProducer<T>::fill_pf_features(const T& lep, btagbtvdeep::DeepBoostedJetFeatures& features){
    auto jet = dynamic_cast<const pat::Jet*>(&(*lep.userCand("jetForLepJetVar"))); 

    std::vector<edm::Ptr<reco::Candidate> > daughters;
    if (jet)
      daughters=jet->daughterPtrVector();


    for (auto& var : pf_vars_){
      std::cout << "Fillin var " << var->first << " for " << daughters.size() << std::endl;
      features.add(var->first);
      features.reserve(var->first,daughters.size());
      for(const auto _d : daughters) {
	const auto d = dynamic_cast<const pat::PackedCandidate*>(_d.get());
	features.fill(var->first, var->second(*d));
      }
    }

    // afaik these need to be hardcoded because I cannot put userFloats to pat::packedCandidates
    features.add("PF_phi_rel");
    features.reserve("PF_phi_rel",daughters.size());
    features.add("PF_eta_rel");
    features.reserve("PF_eta_rel",daughters.size());

    for(const auto _d : daughters) {
      const auto d = dynamic_cast<const pat::PackedCandidate*>(_d.get());
      features.fill("PF_phi_rel", deltaPhi(lep.phi(), d->phi()));
      features.fill("PF_eta_rel", lep.eta()-d->eta());
    }

  }

  template <typename T>
  void LeptonTagInfoCollectionProducer<T>::fill_sv_features(const T& lep, btagbtvdeep::DeepBoostedJetFeatures& features){

    auto jet = dynamic_cast<const pat::Jet*>(&(*lep.userCand("jetForLepJetVar"))); 
    std::vector<size_t> jetSVs;

    if (jet){
      for (size_t isv =0; isv < svs_->size(); ++isv){
	if (matchByCommonSourceCandidatePtr(*jet, svs_->at(isv))){
	  jetSVs.push_back(isv);
	  break;
	}
      }
    }

    for (auto& var : sv_vars_){
      std::cout << "Filling " << var->first << "for " << jetSVs.size() << std::endl;
      
      features.add(var->first);
      features.reserve(var->first, jetSVs.size());
      for (auto& isv : jetSVs)
	features.fill(var->first, var->second(svs_->at(isv)));
    }


    // afaik these need to be hardcoded 
    const auto& PV0 = pvs_->front();
    VertexDistance3D vdist;
    VertexDistanceXY vdistXY;

    features.add("MuonSV_dlenSig");
    features.reserve("MuonSV_dlenSig", jetSVs.size());
    features.add("MuonSV_dxy");
    features.reserve("MuonSV_dxy", jetSVs.size());


    for (auto& isv : jetSVs){
      auto sv = svs_->at(isv);
      Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
      features.fill("MuonSV_dlenSig",dl.significance());
      Measurement1D d2d = vdistXY.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
      features.fill("MuonSV_dxy",d2d.value());
    }

    
  }

  typedef LeptonTagInfoCollectionProducer<Muon> MuonInfoCollectionProducer;
  typedef LeptonTagInfoCollectionProducer<Electron> ElectronInfoCollectionProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonInfoCollectionProducer);
DEFINE_FWK_MODULE(ElectronInfoCollectionProducer);

}
