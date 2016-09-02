#include "FlatTree/FlatTreeProducer/interface/MCTruth.hh"

void MCTruth::fillGenParticles(const edm::Event& iEvent,
			       const edm::EventSetup& iSetup,
			       FlatTree& tree,
			       const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;

   int gen_n = 0;
   
   std::vector<float> gen_pt;
   std::vector<float> gen_eta;
   std::vector<float> gen_phi;
   std::vector<float> gen_m;
   std::vector<float> gen_E;
   std::vector<int> gen_id;
   std::vector<int> gen_status;
   std::vector<int> gen_charge;
   std::vector<int> gen_index;
   std::vector<int> gen_mother_index;
   std::vector<int> gen_daughter_n;
   std::vector<std::vector<int> > gen_daughter_index;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	float ptGen = mcp->pt();
	float etaGen = mcp->eta();
	float phiGen = mcp->phi();
	float mGen = mcp->mass();
	float EGen = mcp->energy();
	int idGen = mcp->pdgId();
	int statusGen = mcp->status();
	int chargeGen = mcp->charge();
	int indexGen = gen_n;

	const reco::GenParticle* mom = getMother(*mcp);

	reco::GenParticleCollection genParticlesCollection_m = *GenParticles;
	reco::GenParticleCollection::const_iterator genParticleSrc_m;
	
	int mother_index = 0;
	for(genParticleSrc_m = genParticlesCollection_m.begin();
	    genParticleSrc_m != genParticlesCollection_m.end();
	    genParticleSrc_m++)
	  {
	     reco::GenParticle *mcp_m = &(const_cast<reco::GenParticle&>(*genParticleSrc_m));
	     if( fabs(mcp_m->pt()-mom->pt()) < 10E-6 && fabs(mcp_m->eta()-mom->eta()) < 10E-6 )
	       {
		  break;
	       }		       
	     mother_index++;
	  }		  
	
	int daughter_n = 0;
	std::vector<int> daughter_index;
	
	const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	  {
	     if( idr->isAvailable() ) 
	       {		       
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();

		  reco::GenParticleCollection genParticlesCollection_s = *GenParticles;
		  reco::GenParticleCollection::const_iterator genParticleSrc_s;
		  
		  int index = 0;
		  for(genParticleSrc_s = genParticlesCollection_s.begin();
		      genParticleSrc_s != genParticlesCollection_s.end();
		      genParticleSrc_s++)
		    {
		       reco::GenParticle *mcp_s = &(const_cast<reco::GenParticle&>(*genParticleSrc_s));
		       if( fabs(mcp_s->pt()-(*d).pt()) < 10E-6 && fabs(mcp_s->eta()-(*d).eta()) < 10E-6 )
			 {
			    break;
			 }		       
		       index++;
		    }		  
		  
		  daughter_index.push_back(index);
		  daughter_n++;
	       }
	  }	

	gen_pt.push_back(ptGen);
	gen_eta.push_back(etaGen);
	gen_phi.push_back(phiGen);
	gen_m.push_back(mGen);
	gen_E.push_back(EGen);
	gen_id.push_back(idGen);
	gen_status.push_back(statusGen);
	gen_charge.push_back(chargeGen);
	gen_index.push_back(indexGen);
	gen_mother_index.push_back(mother_index);
	gen_daughter_n.push_back(daughter_n);
	gen_daughter_index.push_back(daughter_index);
	
	gen_n++;
     }
   
   tree.gen_n = gen_n;
   tree.gen_pt = gen_pt;
   tree.gen_eta = gen_eta;
   tree.gen_phi = gen_phi;
   tree.gen_m = gen_m;
   tree.gen_E = gen_E;
   tree.gen_status = gen_status;
   tree.gen_id = gen_id;
   tree.gen_charge = gen_charge;
   tree.gen_index = gen_index;
   tree.gen_mother_index = gen_mother_index;
   tree.gen_daughter_n = gen_daughter_n;
   tree.gen_daughter_index = gen_daughter_index;
}

void MCTruth::fillGenPV(const edm::Event& iEvent,
			const edm::EventSetup& iSetup,
			FlatTree& tree,
			const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;

   float gen_PVz = -666;
   for( size_t i=0;i<GenParticles->size();++i )
     {	
	const reco::GenParticle & genIt = (*GenParticles)[i];
	
	int status = genIt.status();
	if( (status>=21 && status<=29) || status==3 ) 
	  {
	     gen_PVz = genIt.vz();
	  }	
     }
   tree.gen_PVz = gen_PVz;
}

bool MCTruth::doMatch(const edm::Event& iEvent,
		      const edm::EventSetup& iSetup,
		      const edm::Handle<std::vector<reco::GenParticle> >& GenParticles,
		      reco::GenParticle *genp,
		      float &drMin,
		      float pt, float eta, float phi, int pdgId)
{
   bool foundMatch = 0;
   
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;

   float drmin = 0.2;
   float ptRatMin = 0.5;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	float ptGen = mcp->pt();
	float etaGen = mcp->eta();
	float phiGen = mcp->phi();
	int idGen = mcp->pdgId();
	int statusGen = mcp->status();

	if( statusGen != 1 && statusGen != 3 ) continue;
	if( abs(pdgId) != abs(idGen) ) continue;
	
	float dr = GetDeltaR(eta,phi,etaGen,phiGen);
	float ptRat = (pt > 0.) ? fabs(pt-ptGen)/pt : 10E+10;
	
	if( dr < drmin && ptRat < ptRatMin )
	  {
	     drmin = dr;
	     foundMatch = 1;
	     genp = mcp;
	  }	
     }
   
   drMin = drmin;
   
   return foundMatch;
}

void MCTruth::Init(FlatTree &tree)
{
   tree.gen_pt.clear();
   tree.gen_eta.clear();
   tree.gen_phi.clear();
   tree.gen_m.clear();
   tree.gen_status.clear();
   tree.gen_id.clear();
   tree.gen_charge.clear();
   tree.gen_index.clear();
   tree.gen_mother_index.clear();
   tree.gen_daughter_n.clear();
   tree.gen_daughter_index.clear();
}

reco::GenParticle* MCTruth::getUnique(const reco::GenParticle* p,bool verbose)
{
   reco::GenParticle *pcur = const_cast<reco::GenParticle*>(p);
   
   if( verbose )
     {	
	std::cout << "---------b--------" << std::endl;
	std::cout << "INITIAL = " << pcur->pdgId() << " " << pcur->status() << std::endl;
     }
   
   while( 1 )
     {
	bool foundDupl = false;

	const reco::GenParticleRefVector& daughterRefs = pcur->daughterRefVector();
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr)
	  {	     
	     if( idr->isAvailable() )
	       {		  
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();

		  if( d )
		    {			       
////		       if( fabs(d->pdgId()) != 15 && d->status() == 2 ) continue;
//		       std::cout << d->pdgId() << " " << d->status() << std::endl;
		       if( verbose )
			 {		  
			    std::cout << "current: " << d->pdgId() << " " << d->status() << std::endl;
			    std::cout << "pcur: " << pcur->pdgId() << " " << pcur->status() << std::endl;
			 }
		       
		       if( d->pdgId() == pcur->pdgId() )
			 {
			    pcur = const_cast<reco::GenParticle*>(d);
			    foundDupl = true;
		       
			    if( verbose )
			      {		       
				 std::cout << "Found duplicate, switch to it" << std::endl;
				 std::cout << "Number of children = " << pcur->numberOfDaughters() << std::endl;
			      }		  
			 }
		    }
		  else break; // the world is fcked up in this case
	       }
	     else break;
	  }
	
	if( !foundDupl ) break;
     }   
   
   if( verbose )
     {   
	std::cout << "FINAL: id=" << pcur->pdgId() << " status=" << pcur->status() << 
	  " daughters=" << pcur->numberOfDaughters() << std::endl;
	
	const reco::GenParticleRefVector& daughterRefs = pcur->daughterRefVector();
	int ip = 0;
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr)
	  {	     
	     if( idr->isAvailable() )
	       {		  
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();
		  std::cout << "daughter #" << ip << ": id=" << d->pdgId() << std::endl;
		  ip++;
	       }
	  }	
	std::cout << "---------e--------" << std::endl;
     }
      
   return pcur;
}

void MCTruth::p4toTLV(reco::Particle::LorentzVector vp4,TLorentzVector& tlv)
{
   return tlv.SetPxPyPzE(vp4.px(),vp4.py(),vp4.pz(),vp4.energy());
}

const reco::GenParticle* MCTruth::getMother(const reco::GenParticle &part)
{
   const reco::GenParticle *mom = &part;
   while( mom->numberOfMothers() > 0 )
     {	     
	for( unsigned int j=0;j<mom->numberOfMothers();++j )
	  {		  
	     mom = dynamic_cast<const reco::GenParticle*>(mom->mother(j));
	     if( mom->pdgId() != part.pdgId() )
	       {		       
		  return mom;
	       }
	  }	     
     }
   
   return mom;
}
