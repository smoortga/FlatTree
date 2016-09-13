#ifndef MCTRUTH_H
#define MCTRUTH_H

#include <string>
#include <iostream>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FlatTree/FlatTreeProducer/interface/FlatTree.hh"

#define DEFVAL -666

class MCTruth
{
 public:
   
   MCTruth() {};
   
   void Init(FlatTree &tree);

   void fillGenParticles(const edm::Event& iEvent,
			 const edm::EventSetup& iSetup,
			 FlatTree& tree,
			 const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillGenPV(const edm::Event& iEvent,
		  const edm::EventSetup& iSetup,
		  FlatTree& tree,
		  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   bool doMatch(const edm::Event& iEvent,
		const edm::EventSetup& iSetup,
		const edm::Handle<std::vector<reco::GenParticle> >& GenParticles,
		reco::GenParticle *genp,
		float &drMin,
		float pt, float eta, float phi, int pdgId);
   
   void p4toTLV(reco::Particle::LorentzVector vp4,
		TLorentzVector& tlv);
   
   const reco::GenParticle* getMother(const reco::GenParticle&);
};

#endif
