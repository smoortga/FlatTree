#include <memory>
#include <iostream>
#include <sstream>

#include "TRegexp.h"
#include "TString.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "FlatTree/FlatTreeProducer/interface/tinyxml2.h"

#include "FlatTree/FlatTreeProducer/interface/FlatTree.hh"
#include "FlatTree/FlatTreeProducer/interface/MCTruth.hh"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TFile.h"
#include "TThreadSlots.h"
#include "TROOT.h"
#include "Compression.h"

using namespace tinyxml2;

class FlatTreeProducer : public edm::EDAnalyzer
{
 public:

   explicit FlatTreeProducer(const edm::ParameterSet&);
   ~FlatTreeProducer();

   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   bool foundTrigger(const std::string& name) const;

 private:

   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

   virtual void beginRun(const edm::Run&, const edm::EventSetup&);
   virtual void endRun(const edm::Run&, const edm::EventSetup&);

   FlatTree* ftree;
   const edm::Service<TFileService> fs;

   TH1D* hcount;
   TH1D* hweight;

   XMLDocument xmlconf;

   std::string dataFormat_;
   bool isData_;
   bool applyMETFilters_;
   bool fillMCScaleWeight_;
   bool fillPUInfo_;
   int nPdf_;

   HLTConfigProvider hltConfig_;
   HLTPrescaleProvider hltPrescale_;

   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
   edm::EDGetTokenT<edm::TriggerResults> triggerBitsPAT_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
   edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

   edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
   edm::EDGetTokenT<pat::ElectronCollection> electronPATToken_;
   edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronToken_;
   edm::EDGetTokenT<pat::MuonCollection> muonToken_;
   edm::EDGetTokenT<pat::TauCollection> tauToken_;
   edm::EDGetTokenT<pat::JetCollection> jetToken_;
   edm::EDGetTokenT<edm::View<pat::Jet> > viewJetToken_;
   edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
   edm::EDGetTokenT<std::vector<pat::MET> > metTokenAOD_;
   edm::EDGetTokenT<pat::METCollection> metTokenPuppi_;
   edm::EDGetTokenT<pat::METCollection> metTokenNoHF_;
   edm::EDGetTokenT<pat::METCollection> metTokenMINIAOD_;
   edm::EDGetTokenT<double> rhoToken_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
   edm::EDGetTokenT<LHEEventProduct> LHEEventProductToken_;
   edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_;
   edm::EDGetTokenT<reco::ConversionCollection> hConversionsToken_;

   edm::EDGetTokenT<bool> badMuonFilterToken_;
   edm::EDGetTokenT<bool> badChargedCandidateFilterToken_;

   edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoCBIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseCBIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumCBIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleTightCBIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPCBIdMapToken_;
   
   edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumMVAIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleTightMVAIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
   edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

   edm::EDGetTokenT<double> metSigToken_;
   edm::EDGetTokenT<math::Error<2>::type> metCovToken_;

   std::vector<std::string> filterTriggerNames_;

   JetCorrectionUncertainty *jecUnc;
};

FlatTreeProducer::FlatTreeProducer(const edm::ParameterSet& iConfig):
hltPrescale_(iConfig,consumesCollector(),*this)
{
   // ########################
   // #  Create output tree  #
   // ########################
   //
   TFile& f = fs->file();
   f.SetCompressionAlgorithm(ROOT::kZLIB);
   f.SetCompressionLevel(9);
   ftree = new FlatTree(fs->make<TTree>("tree","tree"));

   // #############################################################
   // #  Read parameters from python file and get consume tokens  #
   // #############################################################
   //
   dataFormat_           = iConfig.getParameter<std::string>("dataFormat");
   nPdf_                 = iConfig.getParameter<int>("nPDF");
   fillMCScaleWeight_    = iConfig.getParameter<bool>("fillMCScaleWeight");
   fillPUInfo_           = iConfig.getParameter<bool>("fillPUInfo");
   isData_               = iConfig.getParameter<bool>("isData");
   applyMETFilters_      = iConfig.getParameter<bool>("applyMETFilters");
   triggerBits_          = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT")));
   triggerBitsPAT_       = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("PAT")));
   triggerObjects_       = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"));
   triggerPrescales_     = consumes<pat::PackedTriggerPrescales>(edm::InputTag(std::string("patTrigger")));
   vertexToken_          = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexInput"));
   electronPATToken_     = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronPATInput"));
   electronToken_        = consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronInput"));
   muonToken_            = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonInput"));
   tauToken_             = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauInput"));
   jetToken_             = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetInput"));
   viewJetToken_         = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetInput"));
   genJetToken_          = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetInput"));
   metTokenAOD_          = consumes<std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("metInput"));
   metTokenMINIAOD_      = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metInput"));
   metTokenNoHF_         = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metNoHFInput"));
   rhoToken_             = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInput"));
   genParticlesToken_    = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInput"));
   genEventInfoToken_    = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoInput"));
   LHEEventProductToken_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductInput"));
   puInfoToken_          = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoInput"));
   bsToken_              = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bsInput"));
   hConversionsToken_    = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("hConversionsInput"));

   badMuonFilterToken_   = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"));
   badChargedCandidateFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"));

   eleVetoCBIdMapToken_    = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoCBIdMap"));
   eleLooseCBIdMapToken_   = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseCBIdMap"));
   eleMediumCBIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumCBIdMap"));
   eleTightCBIdMapToken_   = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightCBIdMap"));
   eleHEEPCBIdMapToken_    = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPCBIdMap"));

   eleMediumMVAIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumMVAIdMap"));
   eleTightMVAIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightMVAIdMap"));
   mvaValuesMapToken_      = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"));
   mvaCategoriesMapToken_  = consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"));

   filterTriggerNames_     = iConfig.getUntrackedParameter<std::vector<std::string> >("filterTriggerNames");

   metSigToken_            = consumes<double>(iConfig.getParameter<edm::InputTag>("metSigInput"));
   metCovToken_            = consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("metCovInput"));

   // #########################
   // #  Read XML config file #
   // #########################
   //
   std::string confFile = iConfig.getParameter<std::string>("confFile");
   int buffersize = iConfig.getParameter<int>("bufferSize");
   if (buffersize <= 0) buffersize = 32000;
   
   xmlconf.LoadFile("conf.xml");
   XMLElement* tElement = xmlconf.FirstChildElement("variables")->FirstChildElement("var");

   for( XMLElement* child=tElement;child!=0;child=child->NextSiblingElement() )
     {
        std::string vname = child->ToElement()->Attribute("name");
        std::string vsave = child->ToElement()->Attribute("save");
        bool bsave = atoi(vsave.c_str());
        if( child->ToElement()->Attribute("mc") )
	  {
	     std::string vmc = child->ToElement()->Attribute("mc");
	     bool mc =  atoi(vmc.c_str());
	     if( isData_ && mc ) bsave = 0; // force the exclusion of mc-related variables when running on data
	  }

        ftree->conf.insert(std::make_pair(vname,bsave));
     }
   
   ftree->CreateBranches(buffersize);
   
   // ###############################
   // #  Add count & weight histos  #
   // ###############################
   //
   hcount = fs->make<TH1D>("hcount","hcount",1,0.,1.);
   hweight = fs->make<TH1D>("hweight","hweight",1,0.,1.);
}

FlatTreeProducer::~FlatTreeProducer()
{
}

void FlatTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   hcount->SetBinContent(1,hcount->GetBinContent(1)+1);

   ftree->Init();

   // Initial-state info
   edm::Handle<GenEventInfoProduct> genEventInfo;
   if( !isData_ ) iEvent.getByToken(genEventInfoToken_,genEventInfo);

   // LHE
   edm::Handle<LHEEventProduct> lheEventProduct;
   if( !isData_ && fillMCScaleWeight_ ) iEvent.getByToken(LHEEventProductToken_,lheEventProduct);

   // Gen particles
   edm::Handle<reco::GenParticleCollection> genParticlesHandle;
   if( !isData_ ) iEvent.getByToken(genParticlesToken_,genParticlesHandle);

   // Beamspot
   edm::Handle<reco::BeamSpot> bsHandle;
   iEvent.getByToken(bsToken_, bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();

   // Primary vertex
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken_,vertices);

   // Triggers
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_,triggerBits);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   edm::Handle<bool> badMuonFilter;
   edm::Handle<bool> badChargedCandidateFilter;

   iEvent.getByToken(badMuonFilterToken_,badMuonFilter);
   iEvent.getByToken(badChargedCandidateFilterToken_,badChargedCandidateFilter);

   edm::Handle<edm::TriggerResults> triggerBitsPAT;
   iEvent.getByToken(triggerBitsPAT_,triggerBitsPAT);
   edm::TriggerNames namesPAT;
   if( triggerBitsPAT.isValid() )
     {
        namesPAT = iEvent.triggerNames(*triggerBitsPAT);
     }
   
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);

   // Pile-up
   edm::Handle<std::vector< PileupSummaryInfo> > pileupInfo;
   if( !isData_ && fillPUInfo_ ) iEvent.getByToken(puInfoToken_,pileupInfo);

   // Rho info
   edm::Handle<double> rhoPtr;
   iEvent.getByToken(rhoToken_,rhoPtr);

   // MET significance
   edm::Handle<double> metSigPtr;
   iEvent.getByToken(metSigToken_,metSigPtr);
   std::vector<edm::Handle<double> > doubleVec;
   iEvent.getManyByType(doubleVec);
   for(unsigned int i=0;i<doubleVec.size();i++)
     {
        if(doubleVec[i].provenance()->moduleLabel() == "METSignificance")
	  metSigPtr = doubleVec[i];
     }
   if( !metSigPtr.isValid() or metSigPtr.failedToGet() )
     std::cerr << " Fail to access METSignificance branch " << std::endl;

   // MET covariance
   edm::Handle<math::Error<2>::type> metCovPtr;
   iEvent.getByToken(metCovToken_,metCovPtr);
   ftree->met_cov00 = (*metCovPtr)(0,0);
   ftree->met_cov10 = (*metCovPtr)(1,0);
   ftree->met_cov01 = (*metCovPtr)(0,1);
   ftree->met_cov11 = (*metCovPtr)(1,1);

   // Jets
   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_,jets);

   edm::Handle<edm::View<pat::Jet> > view_jets;
   iEvent.getByToken(viewJetToken_,view_jets);

   // GenJets
   edm::Handle<reco::GenJetCollection> genJets;
   if( !isData_ ) iEvent.getByToken(genJetToken_,genJets);
   edm::Handle<reco::JetFlavourMatchingCollection> genJetFlavourMatching;
   if( !isData_ )
     {
        //	iEvent.getByLabel("genJetFlavour",genJetFlavourMatching);
     }

   // Muons
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_,muons);

   // Electrons
   edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(electronToken_,electrons);

   edm::Handle<pat::ElectronCollection> electronsPAT;
   iEvent.getByToken(electronPATToken_,electronsPAT);

   edm::Handle<edm::ValueMap<bool> > veto_cbid_decisions;
   edm::Handle<edm::ValueMap<bool> > loose_cbid_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_cbid_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_cbid_decisions;
   edm::Handle<edm::ValueMap<bool> > heep_cbid_decisions;

   edm::Handle<edm::ValueMap<bool> > medium_mvaid_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_mvaid_decisions;

   iEvent.getByToken(eleVetoCBIdMapToken_,veto_cbid_decisions);
   iEvent.getByToken(eleLooseCBIdMapToken_,loose_cbid_decisions);
   iEvent.getByToken(eleMediumCBIdMapToken_,medium_cbid_decisions);
   iEvent.getByToken(eleTightCBIdMapToken_,tight_cbid_decisions);
   iEvent.getByToken(eleHEEPCBIdMapToken_,heep_cbid_decisions);

   iEvent.getByToken(eleMediumMVAIdMapToken_,medium_mvaid_decisions);
   iEvent.getByToken(eleTightMVAIdMapToken_,tight_mvaid_decisions);

   edm::Handle<edm::ValueMap<float> > mvaValues;
   edm::Handle<edm::ValueMap<int> > mvaCategories;
   iEvent.getByToken(mvaValuesMapToken_,mvaValues);
   iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

   // Taus
   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_,taus);

   // Conversions info
   edm::Handle<reco::ConversionCollection> hConversions;
   if( dataFormat_ != "AOD" ) iEvent.getByToken(hConversionsToken_,hConversions);

   ftree->ev_run = iEvent.id().run();
   ftree->ev_id = iEvent.id().event();
   ftree->ev_lumi = iEvent.id().luminosityBlock();

   float mc_weight = 1.;

   if( genEventInfo.isValid() )
     {
        float wGen = genEventInfo->weight();
        mc_weight = (wGen > 0 ) ? 1. : -1.;

        ftree->mc_id = genEventInfo->signalProcessID();
        ftree->mc_f1 = genEventInfo->pdf()->id.first;
        ftree->mc_f2 = genEventInfo->pdf()->id.second;
        ftree->mc_x1 = genEventInfo->pdf()->x.first;
        ftree->mc_x2 = genEventInfo->pdf()->x.second;
        ftree->mc_scale = genEventInfo->pdf()->scalePDF;
        if( genEventInfo->binningValues().size() > 0 ) ftree->mc_ptHat = genEventInfo->binningValues()[0];
     }
   
   if(! lheEventProduct.failedToGet())
     {
        if( !isData_ && fillMCScaleWeight_ )
	  {
	     ftree->weight_originalXWGTUP = lheEventProduct->originalXWGTUP();
	     
	     if( lheEventProduct->weights().size() > 0 )
	       {
		  ftree->weight_scale_muF0p5 = (genEventInfo->weight())*(lheEventProduct->weights()[2].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 0.5 | muR = 1
		  ftree->weight_scale_muF2   = (genEventInfo->weight())*(lheEventProduct->weights()[1].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 2   | muR = 1
		  ftree->weight_scale_muR0p5 = (genEventInfo->weight())*(lheEventProduct->weights()[6].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 1   | muR = 0.5
		  ftree->weight_scale_muR2   = (genEventInfo->weight())*(lheEventProduct->weights()[3].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 1   | muR = 2
	       }
	     
	     int nPdfAll = lheEventProduct->weights().size();
	     if( nPdf_ < nPdfAll && nPdf_ >= 0 ) nPdfAll = nPdf_;
	     for( int w=0;w<nPdfAll;w++ )
	       {
		  const LHEEventProduct::WGT& wgt = lheEventProduct->weights().at(w);
		  ftree->mc_pdfweights.push_back(wgt.wgt);
		  ftree->mc_pdfweightIds.push_back(wgt.id);
	       }
	  }
     }
   
   ftree->mc_weight = mc_weight;

   hweight->SetBinContent(1,hweight->GetBinContent(1)+ftree->mc_weight);

   if( !isData_ && fillPUInfo_)
     {
        ftree->mc_pu_Npvi = pileupInfo->size();
        for(std::vector<PileupSummaryInfo>::const_iterator pvi=pileupInfo->begin();
	    pvi!=pileupInfo->end();pvi++)
	  {
	     signed int n_bc = pvi->getBunchCrossing();
	     ftree->mc_pu_BunchCrossing.push_back(n_bc);
	     if( n_bc == 0 )
	       {
		  ftree->mc_pu_intime_NumInt = pvi->getPU_NumInteractions();
		  ftree->mc_pu_trueNumInt = pvi->getTrueNumInteractions();
	       }
	     else if( n_bc == -1 ) ftree->mc_pu_before_npu = pvi->getPU_NumInteractions();
	     else if( n_bc == +1 ) ftree->mc_pu_after_npu  = pvi->getPU_NumInteractions();
	     
	     std::vector<float> mc_pu_zpositions;
	     std::vector<float> mc_pu_sumpT_lowpT;
	     std::vector<float> mc_pu_sumpT_highpT;
	     std::vector<int> mc_pu_ntrks_lowpT;
	     std::vector<int> mc_pu_ntrks_highpT;

	     ftree->mc_pu_Nzpositions.push_back(pvi->getPU_zpositions().size());
	     for( unsigned int ipu=0;ipu<pvi->getPU_zpositions().size();ipu++ )
	       {
		  mc_pu_zpositions.push_back((pvi->getPU_zpositions())[ipu]);
		  mc_pu_sumpT_lowpT.push_back((pvi->getPU_sumpT_lowpT())[ipu]);
		  mc_pu_sumpT_highpT.push_back((pvi->getPU_sumpT_highpT())[ipu]);
		  mc_pu_ntrks_lowpT.push_back((pvi->getPU_ntrks_lowpT())[ipu]);
		  mc_pu_ntrks_highpT.push_back((pvi->getPU_ntrks_highpT())[ipu]);
	       }

	     ftree->mc_pu_zpositions.push_back(mc_pu_zpositions);
	     ftree->mc_pu_sumpT_lowpT.push_back(mc_pu_sumpT_lowpT);
	     ftree->mc_pu_sumpT_highpT.push_back(mc_pu_sumpT_highpT);
	     ftree->mc_pu_ntrks_lowpT.push_back(mc_pu_ntrks_lowpT);
	     ftree->mc_pu_ntrks_highpT.push_back(mc_pu_ntrks_highpT);
	  }
    }
   
   bool do_gen_all = ftree->doWrite("gen_all");
   
   MCTruth *mc_truth = new MCTruth();
   
   if( do_gen_all &&
       !isData_ )
     {
        mc_truth->Init(*ftree);
        mc_truth->fillGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
     }
   if( !isData_ )
     {
        mc_truth->Init(*ftree);
        mc_truth->fillGenPV(iEvent,iSetup,*ftree,genParticlesHandle);
     }
   
   bool passMETFilters = 1;
   bool pass_HBHENoiseFilter = 1;
   bool pass_HBHENoiseIsoFilter = 1;
   bool pass_EcalDeadCellTriggerPrimitiveFilter = 1;
   bool pass_goodVertices = 1;
   bool pass_eeBadScFilter = 1;
   bool pass_globalTightHalo2016Filter = 1;

   bool pass_badMuonFilter = *badMuonFilter;
   bool pass_badChargedCandidateFilter = *badChargedCandidateFilter;
   
   if( triggerBitsPAT.isValid() )
     {
        for (unsigned int i = 0, n = triggerBitsPAT->size(); i < n; ++i)
	  {
	     std::string triggerName = namesPAT.triggerName(i);
	     
	     bool isFired = (triggerBitsPAT->accept(i) ? true : false);
	     
	     if( strcmp(triggerName.c_str(),"Flag_HBHENoiseFilter") == 0 )
	       {
		  if( !isFired ) pass_HBHENoiseFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_HBHENoiseIsoFilter") == 0 )
	       {
		  if( !isFired ) pass_HBHENoiseIsoFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_EcalDeadCellTriggerPrimitiveFilter") == 0 )
	       {
		  if( !isFired ) pass_EcalDeadCellTriggerPrimitiveFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_goodVertices") == 0 )
	       {
		  if( !isFired ) pass_goodVertices = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_eeBadScFilter") == 0 )
	       {
		  if( !isFired ) pass_eeBadScFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_globalTightHalo2016Filter") == 0 )
	       {
		  if( !isFired ) pass_globalTightHalo2016Filter = 0;
	       }
	  }
     }
   
   passMETFilters = (pass_HBHENoiseFilter &&
		     pass_HBHENoiseIsoFilter &&
		     pass_EcalDeadCellTriggerPrimitiveFilter &&
		     pass_goodVertices &&
		     pass_eeBadScFilter &&
		     pass_globalTightHalo2016Filter &&
		     pass_badMuonFilter &&
		     pass_badChargedCandidateFilter);

   //std::cout << "\n === TRIGGER PATHS === " << std::endl;
   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
     {
        //std::cout << "[" << i << "] " << (triggerBits->accept(i) ? "1" : "0") << "  " << names.triggerName(i)  << std::endl;

        std::string triggerName = names.triggerName(i);

        if( !foundTrigger(triggerName) ) continue;

        ftree->trigger.push_back(i);
        ftree->trigger_name.push_back(triggerName);
        ftree->trigger_pass.push_back(triggerBits->accept(i) ? true : false);
        ftree->trigger_prescale.push_back(triggerPrescales->getPrescaleForIndex(i));

        float HLTprescale = 1.;
        float L1prescale = 1.;
	
        if( isData_ )
	  {
	     std::pair<std::vector<std::pair<std::string,int> >,int> detailedPrescaleInfo =
	       hltPrescale_.prescaleValuesInDetail(iEvent,iSetup,triggerName);

	     HLTprescale = triggerPrescales.isValid() ? detailedPrescaleInfo.second : -1;

	     std::vector <int> l1prescalevals;
	     for( size_t varind=0;varind<detailedPrescaleInfo.first.size();varind++ )
	       {
		  l1prescalevals.push_back(detailedPrescaleInfo.first.at(varind).second);
	       }
	     
	     // find and save minimum l1 prescale of any ORed L1 that seeds the HLT
	     std::vector<int>::iterator result = std::min_element(std::begin(l1prescalevals),
								  std::end(l1prescalevals));
	     size_t minind = std::distance(std::begin(l1prescalevals), result);
	     // sometimes there are no L1s associated with a HLT.
	     // In that case, this branch stores -1 for the l1prescale
	     L1prescale = minind < l1prescalevals.size() ? l1prescalevals.at(minind) : -1;
	     
	     //	     std::cout << HLTprescale << " " << L1prescale << std::endl;
	  }
	
        ftree->trigger_HLTprescale.push_back(HLTprescale);
        ftree->trigger_L1prescale.push_back(L1prescale);
     }

   //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
   for (pat::TriggerObjectStandAlone obj : *triggerObjects)
     {
        // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
	
        // Trigger object basic informations (pt, eta, phi)
        //std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;

        ftree->triggerobject_pt.push_back(obj.pt());
        ftree->triggerobject_eta.push_back(obj.eta());
        ftree->triggerobject_phi.push_back(obj.phi());

        // Trigger object collection
        //std::cout << "\t   Collection: " << obj.collection() << std::endl;
        ftree->triggerobject_collection.push_back(obj.collection());

        // Trigger object type IDs
        ftree->triggerobject_filterIds_n.push_back(obj.filterIds().size());
        for (unsigned h = 0; h < obj.filterIds().size(); ++h)
	  {
	     ftree->triggerobject_isTriggerL1Mu.push_back(obj.filterIds()[h] == -81 ? true : false);
	     ftree->triggerobject_isTriggerL1NoIsoEG.push_back(obj.filterIds()[h] == -82 ? true : false);
	     ftree->triggerobject_isTriggerL1IsoEG.push_back(obj.filterIds()[h] == -83 ? true : false);
	     ftree->triggerobject_isTriggerL1CenJet.push_back(obj.filterIds()[h] == -84 ? true : false);
	     ftree->triggerobject_isTriggerL1ForJet.push_back(obj.filterIds()[h] == -85 ? true : false);
	     ftree->triggerobject_isTriggerL1TauJet.push_back(obj.filterIds()[h] == -86 ? true : false);
	     ftree->triggerobject_isTriggerL1ETM.push_back(obj.filterIds()[h] == -87 ? true : false);
	     ftree->triggerobject_isTriggerL1ETT.push_back(obj.filterIds()[h] == -88 ? true : false);
	     ftree->triggerobject_isTriggerL1HTT.push_back(obj.filterIds()[h] == -89 ? true : false);
	     ftree->triggerobject_isTriggerL1HTM.push_back(obj.filterIds()[h] == -90 ? true : false);
	     ftree->triggerobject_isTriggerL1JetCounts.push_back(obj.filterIds()[h] == -91 ? true : false);
	     ftree->triggerobject_isTriggerL1HfBitCounts.push_back(obj.filterIds()[h] == -92 ? true : false);
	     ftree->triggerobject_isTriggerL1HfRingEtSums.push_back(obj.filterIds()[h] == -93 ? true : false);
	     ftree->triggerobject_isTriggerL1TechTrig.push_back(obj.filterIds()[h] == -94 ? true : false);
	     ftree->triggerobject_isTriggerL1Castor.push_back(obj.filterIds()[h] == -95 ? true : false);
	     ftree->triggerobject_isTriggerL1BPTX.push_back(obj.filterIds()[h] == -96 ? true : false);
	     ftree->triggerobject_isTriggerL1GtExternal.push_back(obj.filterIds()[h] == -97 ? true : false);

	     ftree->triggerobject_isHLT_TriggerPhoton.push_back(obj.filterIds()[h] == 81 ? true : false);
	     ftree->triggerobject_isHLT_TriggerElectron.push_back(obj.filterIds()[h] == 82 ? true : false);
	     ftree->triggerobject_isHLT_TriggerMuon.push_back(obj.filterIds()[h] == 83 ? true : false);
	     ftree->triggerobject_isHLT_TriggerTau.push_back(obj.filterIds()[h] == 84 ? true : false);
	     ftree->triggerobject_isHLT_TriggerJet.push_back(obj.filterIds()[h] == 85 ? true : false);
	     ftree->triggerobject_isHLT_TriggerBJet.push_back(obj.filterIds()[h] == 86 ? true : false);
	     ftree->triggerobject_isHLT_TriggerMET.push_back(obj.filterIds()[h] == 87 ? true : false);
	     ftree->triggerobject_isHLT_TriggerTET.push_back(obj.filterIds()[h] == 88 ? true : false);
	     ftree->triggerobject_isHLT_TriggerTHT.push_back(obj.filterIds()[h] == 89 ? true : false);
	     ftree->triggerobject_isHLT_TriggerMHT.push_back(obj.filterIds()[h] == 90 ? true : false);
	     ftree->triggerobject_isHLT_TriggerTrack.push_back(obj.filterIds()[h] == 91 ? true : false);
	     ftree->triggerobject_isHLT_TriggerCluster.push_back(obj.filterIds()[h] == 92 ? true : false);
	     ftree->triggerobject_isHLT_TriggerMETSig.push_back(obj.filterIds()[h] == 93 ? true : false);
	     ftree->triggerobject_isHLT_TriggerELongit.push_back(obj.filterIds()[h] == 94 ? true : false);
	     ftree->triggerobject_isHLT_TriggerMHTSig.push_back(obj.filterIds()[h] == 95 ? true : false);
	     ftree->triggerobject_isHLT_TriggerHLongit.push_back(obj.filterIds()[h] == 96 ? true : false);

	     ftree->triggerobject_filterIds.push_back(obj.filterIds()[h]);
	  }

        // Trigger object filter
        ftree->triggerobject_filterLabels_n.push_back(obj.filterLabels().size());
        for (unsigned h = 0; h < obj.filterLabels().size(); ++h)
	  {
	     //std::cout << "FilterLabel: " << obj.filterLabels()[h] << std::endl;
	     ftree->triggerobject_filterLabels.push_back(obj.filterLabels()[h]);
	  }

        //std::cout << std::endl;
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);

        // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
        // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
        // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
        //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
        ftree->triggerobject_pathNamesAll_n.push_back(pathNamesAll.size());
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
	  {
	     bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
	     bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
	     bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
	     bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
	     //std::cout << "   " << pathNamesAll[h];
	     //if (isBoth) std::cout << "(L,3)" << std::endl;
	     //if (isL3 && !isBoth) std::cout << "(*,3)" << std::endl;
	     //if (isLF && !isBoth) std::cout << "(L,*)" << std::endl;
	     //if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)" << std::endl;

	     ftree->triggerobject_pathNamesAll.push_back(pathNamesAll[h]);
	     ftree->triggerobject_pathNamesAll_isBoth.push_back(isBoth);
	     ftree->triggerobject_pathNamesAll_isL3.push_back(isL3);
	     ftree->triggerobject_pathNamesAll_isLF.push_back(isLF);
	     ftree->triggerobject_pathNamesAll_isNone.push_back(isNone);
	  }
	
        ftree->triggerobject_n = ftree->triggerobject_pt.size();
     }
   
   // =========== END OF TRIGGER ==========
   //
   reco::Vertex *primVtx = NULL;

   ftree->nvertex = int(vertices->size());

   if( ! vertices->empty() )
     {
        const reco::Vertex &PV = vertices->front();
        primVtx = (reco::Vertex*)&PV;
     }
   
   if( primVtx )
     {
        ftree->pv_x = primVtx->position().x();
        ftree->pv_y = primVtx->position().y();
        ftree->pv_z = primVtx->position().z();
	
        ftree->pv_xError = primVtx->xError();
        ftree->pv_yError = primVtx->yError();
        ftree->pv_zError = primVtx->zError();

        ftree->pv_ndof = primVtx->chi2();
        ftree->pv_ndof = primVtx->ndof();
        ftree->pv_rho = primVtx->position().Rho();
        ftree->pv_isFake = primVtx->isFake();
     }
   
   ftree->ev_rho = *rhoPtr;
   
   ftree->met_sig = *metSigPtr;

   // MET
   edm::Handle<std::vector<pat::MET> > metAOD;
   edm::Handle<pat::METCollection> metMINIAOD;
   if( dataFormat_ != "AOD" )
     {
        iEvent.getByToken(metTokenMINIAOD_,metMINIAOD);
        const pat::MET &metv = metMINIAOD->front();
	
        ftree->met_px    = metv.px();
        ftree->met_py    = metv.py();
        ftree->met_pt    = metv.pt();
        ftree->met_phi   = metv.phi();
        ftree->met_sumet = metv.sumEt();

        if( !isData_ )
	  {
	     ftree->metGen_px    = metv.genMET()->px();
	     ftree->metGen_py    = metv.genMET()->py();
	     ftree->metGen_pt    = metv.genMET()->pt();
	     ftree->metGen_phi   = metv.genMET()->phi();
	     ftree->metGen_sumet = metv.genMET()->sumEt();

	     ftree->metGen_NeutralEMEt  = metv.genMET()->NeutralEMEt();
	     ftree->metGen_ChargedEMEt  = metv.genMET()->ChargedEMEt();
	     ftree->metGen_NeutralHadEt = metv.genMET()->NeutralHadEt();
	     ftree->metGen_ChargedHadEt = metv.genMET()->ChargedHadEt();
	     ftree->metGen_MuonEt       = metv.genMET()->MuonEt();
	     ftree->metGen_InvisibleEt  = metv.genMET()->InvisibleEt();
	  }
	
        ftree->met_uncorrectedPt    = metv.uncorPt();
        ftree->met_uncorrectedPhi   = metv.uncorPhi();
        ftree->met_uncorrectedSumEt = metv.uncorSumEt();

        ftree->met_caloMETPt    = metv.caloMETPt();
        ftree->met_caloMETPhi   = metv.caloMETPhi();
        ftree->met_caloMETSumEt = metv.caloMETSumEt();

        pat::MET::METCorrectionLevel level = pat::MET::METCorrectionLevel::Type1;

        ftree->met_shiftedPx_JetEnUp = metv.shiftedPx(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedPx_JetEnDown = metv.shiftedPx(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedPx_JetResUp = metv.shiftedPx(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedPx_JetResDown = metv.shiftedPx(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedPx_MuonEnUp = metv.shiftedPx(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedPx_MuonEnDown = metv.shiftedPx(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedPx_ElectronEnUp = metv.shiftedPx(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedPx_ElectronEnDown = metv.shiftedPx(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedPx_TauEnUp = metv.shiftedPx(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedPx_TauEnDown = metv.shiftedPx(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedPx_UnclusteredEnUp = metv.shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedPx_UnclusteredEnDown = metv.shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedPx_NoShift = metv.shiftedPx(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedPx_PhotonEnUp = metv.shiftedPx(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedPx_PhotonEnDown = metv.shiftedPx(pat::MET::METUncertainty::PhotonEnDown,level);

        ftree->met_shiftedPy_JetEnUp = metv.shiftedPy(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedPy_JetEnDown = metv.shiftedPy(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedPy_JetResUp = metv.shiftedPy(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedPy_JetResDown = metv.shiftedPy(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedPy_MuonEnUp = metv.shiftedPy(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedPy_MuonEnDown = metv.shiftedPy(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedPy_ElectronEnUp = metv.shiftedPy(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedPy_ElectronEnDown = metv.shiftedPy(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedPy_TauEnUp = metv.shiftedPy(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedPy_TauEnDown = metv.shiftedPy(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedPy_UnclusteredEnUp = metv.shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedPy_UnclusteredEnDown = metv.shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedPy_NoShift = metv.shiftedPy(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedPy_PhotonEnUp = metv.shiftedPy(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedPy_PhotonEnDown = metv.shiftedPy(pat::MET::METUncertainty::PhotonEnDown,level);

        ftree->met_shiftedPhi_JetEnUp = metv.shiftedPhi(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedPhi_JetEnDown = metv.shiftedPhi(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedPhi_JetResUp = metv.shiftedPhi(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedPhi_JetResDown = metv.shiftedPhi(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedPhi_MuonEnUp = metv.shiftedPhi(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedPhi_MuonEnDown = metv.shiftedPhi(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedPhi_ElectronEnUp = metv.shiftedPhi(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedPhi_ElectronEnDown = metv.shiftedPhi(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedPhi_TauEnUp = metv.shiftedPhi(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedPhi_TauEnDown = metv.shiftedPhi(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedPhi_UnclusteredEnUp = metv.shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedPhi_UnclusteredEnDown = metv.shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedPhi_NoShift = metv.shiftedPhi(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedPhi_PhotonEnUp = metv.shiftedPhi(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedPhi_PhotonEnDown = metv.shiftedPhi(pat::MET::METUncertainty::PhotonEnDown,level);

        ftree->met_shiftedSumEt_JetEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedSumEt_JetEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedSumEt_JetResUp = metv.shiftedSumEt(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedSumEt_JetResDown = metv.shiftedSumEt(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedSumEt_MuonEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedSumEt_MuonEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedSumEt_ElectronEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedSumEt_ElectronEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedSumEt_TauEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedSumEt_TauEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedSumEt_UnclusteredEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedSumEt_UnclusteredEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedSumEt_NoShift = metv.shiftedSumEt(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedSumEt_PhotonEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedSumEt_PhotonEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::PhotonEnDown,level);
     }
   else
     {
        iEvent.getByToken(metTokenAOD_,metAOD);
        const pat::MET &metv = metAOD->front();
	
        ftree->met_pt = metv.pt();
        ftree->met_phi = metv.phi();
        ftree->met_sumet = metv.sumEt();
     }
   
   // MET no HF
   edm::Handle<pat::METCollection> metNoHF;
   try
     {
        iEvent.getByToken(metTokenNoHF_,metNoHF);
     }
   catch (...)
     {;
     }
   if( metNoHF.isValid() )
     {
        const pat::MET &met = metNoHF->front();
	
        ftree->metNoHF_pt    = met.pt();
        ftree->metNoHF_phi   = met.phi();
        ftree->metNoHF_sumet = met.sumEt();
     }

   int nElec = electrons->size();
   int nElecPass = 0;

   for(int ie=0;ie<nElec;ie++)
     {
        const pat::Electron& elec = electronsPAT->at(ie);

        // Skimming electrons with pT < 5 GeV.
        if (elec.pt() < 5) continue;

        ftree->el_pt.push_back(elec.pt());
        ftree->el_eta.push_back(elec.eta());
        ftree->el_phi.push_back(elec.phi());
        ftree->el_m.push_back(elec.mass());
        ftree->el_E.push_back(elec.energy());
        ftree->el_id.push_back(elec.pdgId());
        ftree->el_charge.push_back(elec.charge());

        ftree->el_isGsfCtfScPixChargeConsistent.push_back(elec.isGsfCtfScPixChargeConsistent());
        ftree->el_isGsfScPixChargeConsistent.push_back(elec.isGsfScPixChargeConsistent());
        ftree->el_hadronicOverEm.push_back(elec.hadronicOverEm());

        // IP
        const reco::GsfTrackRef gsfTrack = elec.gsfTrack();
        bool hasGsfTrack = ( gsfTrack.isNonnull() );
        ftree->el_hasGsfTrack.push_back(hasGsfTrack);
        ftree->el_gsfTrack_d0.push_back((hasGsfTrack) ? gsfTrack->d0() : -666);
        ftree->el_gsfTrack_z0.push_back((hasGsfTrack) ? gsfTrack->dz() : -666);
        ftree->el_gsfTrack_d0Error.push_back((hasGsfTrack) ? gsfTrack->d0Error() : -666);
        ftree->el_gsfTrack_z0Error.push_back((hasGsfTrack) ? gsfTrack->dzError() : -666);
        ftree->el_gsfTrack_PV_dxy.push_back((hasGsfTrack) ? gsfTrack->dxy(primVtx->position()) : -666);
        ftree->el_gsfTrack_PV_dz.push_back((hasGsfTrack) ? gsfTrack->dz(primVtx->position()) : -666);
        ftree->el_gsfTrack_RP_dxy.push_back((hasGsfTrack) ? gsfTrack->dxy(gsfTrack->referencePoint()) : -666);
        ftree->el_gsfTrack_RP_dz.push_back((hasGsfTrack) ? gsfTrack->dz(gsfTrack->referencePoint()) : -666);
        ftree->el_gsfTrack_BS_dxy.push_back((hasGsfTrack) ? gsfTrack->dxy(beamspot.position()) : -666);
        ftree->el_gsfTrack_BS_dz.push_back((hasGsfTrack) ? gsfTrack->dz(beamspot.position()) : -666);
        ftree->el_gsfTrack_dxyError.push_back((hasGsfTrack) ? gsfTrack->dxyError() : -666);
        ftree->el_gsfTrack_dzError.push_back((hasGsfTrack) ? gsfTrack->dzError() : -666);
        ftree->el_gsfTrack_normalizedChi2.push_back((hasGsfTrack) ? gsfTrack->normalizedChi2() : -666);

        ftree->el_ip3d.push_back(elec.dB(pat::Electron::PV3D));
        ftree->el_ip3dErr.push_back(elec.edB(pat::Electron::PV3D));
        ftree->el_ip2d.push_back(elec.dB(pat::Electron::PV2D));
        ftree->el_ip2dErr.push_back(elec.edB(pat::Electron::PV2D));
        ftree->el_ip3dBS.push_back(elec.dB(pat::Electron::BS3D));
        ftree->el_ip3dBSErr.push_back(elec.edB(pat::Electron::BS3D));
        ftree->el_ip2dBS.push_back(elec.dB(pat::Electron::BS2D));
        ftree->el_ip2dBSErr.push_back(elec.edB(pat::Electron::BS2D));

        // Energy cluster
        ftree->el_superCluster_eta.push_back(elec.superCluster()->eta());
        ftree->el_superCluster_phi.push_back(elec.superCluster()->phi());
        ftree->el_superCluster_energy.push_back(elec.superCluster()->energy());
        ftree->el_superCluster_rawEnergy.push_back(elec.superCluster()->rawEnergy());
        ftree->el_superCluster_preshowerEnergy.push_back(elec.superCluster()->preshowerEnergy());
        ftree->el_superCluster_etaWidth.push_back(elec.superCluster()->etaWidth());
        ftree->el_superCluster_phiWidth.push_back(elec.superCluster()->phiWidth());
        ftree->el_superCluster_preshowerEnergyPlane1.push_back(elec.superCluster()->preshowerEnergyPlane1());
        ftree->el_superCluster_preshowerEnergyPlane2.push_back(elec.superCluster()->preshowerEnergyPlane2());
        ftree->el_superCluster_positionR.push_back(elec.superCluster()->position().R());

        ftree->el_basicClustersSize.push_back(elec.basicClustersSize());
        ftree->el_e1x5.push_back(elec.e1x5());
        ftree->el_e5x5.push_back(elec.e5x5());
        ftree->el_e2x5Max.push_back(elec.e2x5Max());
        ftree->el_sigmaEtaEta.push_back(elec.sigmaEtaEta());
        ftree->el_sigmaIetaIeta.push_back(elec.sigmaIetaIeta());
        ftree->el_sigmaIphiIphi.push_back(elec.sigmaIphiIphi());
        ftree->el_sigmaIetaIphi.push_back(elec.sigmaIetaIphi());
        ftree->el_full5x5_sigmaIphiIphi.push_back(elec.full5x5_sigmaIphiIphi());
        ftree->el_full5x5_sigmaEtaEta.push_back(elec.full5x5_sigmaEtaEta());
        ftree->el_full5x5_sigmaIetaIeta.push_back(elec.full5x5_sigmaIetaIeta());
        ftree->el_full5x5_sigmaIetaIphi.push_back(elec.full5x5_sigmaIetaIphi());
        ftree->el_full5x5_r9.push_back(elec.full5x5_r9());
        ftree->el_full5x5_e1x5.push_back(elec.full5x5_e1x5());
        ftree->el_full5x5_e5x5.push_back(elec.full5x5_e5x5());
        ftree->el_full5x5_e2x5Max.push_back(elec.full5x5_e2x5Max());

        double OneMinusE1x5E5x5 = (elec.e5x5() != 0.) ? 1.-(elec.e1x5()/elec.e5x5()) : -1.;
        double full5x5_OneMinusE1x5E5x5 = (elec.full5x5_e5x5() != 0.) ? 1.-(elec.full5x5_e1x5()/elec.full5x5_e5x5()) : -1.;

        ftree->el_full5x5_OneMinusE1x5E5x5.push_back(full5x5_OneMinusE1x5E5x5);
        ftree->el_OneMinusE1x5E5x5.push_back(OneMinusE1x5E5x5);

        double IoEmIoP = (1.0/elec.ecalEnergy())-(1.0/elec.p());
        ftree->el_IoEmIoP.push_back(IoEmIoP);
        ftree->el_eleEoPout.push_back(elec.eEleClusterOverPout());
        double PreShowerOverRaw = elec.superCluster()->preshowerEnergy()/elec.superCluster()->rawEnergy();
        ftree->el_PreShowerOverRaw.push_back(PreShowerOverRaw);
        double ooEmooP = fabs(1.0/elec.ecalEnergy()-elec.eSuperClusterOverP()/elec.ecalEnergy());
        ftree->el_ooEmooP.push_back(ooEmooP);

        // Track hits
        const reco::HitPattern& pattern = gsfTrack->hitPattern();
        ftree->el_numberOfLostHits.push_back(pattern.numberOfLostHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS));
        ftree->el_expectedMissingInnerHits.push_back(pattern.numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
        ftree->el_numberOfHits.push_back(pattern.numberOfHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS));

        ftree->el_expectedMissingOuterHits.push_back(pattern.numberOfHits(reco::HitPattern::MISSING_OUTER_HITS));
        ftree->el_numberOfValidPixelHits.push_back(pattern.numberOfValidPixelHits());
        ftree->el_numberOfLostPixelHits.push_back(pattern.numberOfLostPixelHits(reco::HitPattern::TRACK_HITS));
        ftree->el_trackerLayersWithMeasurement.push_back(pattern.trackerLayersWithMeasurement());
        ftree->el_pixelLayersWithMeasurement.push_back(pattern.pixelLayersWithMeasurement());
        ftree->el_numberOfValidStripLayersWithMonoAndStereo.push_back(pattern.numberOfValidStripLayersWithMonoAndStereo());
        ftree->el_trackerLayersWithoutMeasurement.push_back(pattern.trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));

        ftree->el_numberOfValidHits.push_back((hasGsfTrack) ? gsfTrack->numberOfValidHits() : -666);
        ftree->el_numberOfLostHitsDefault.push_back((hasGsfTrack) ? gsfTrack->numberOfLostHits() : -666);

        ftree->el_fbrem.push_back(elec.fbrem());

        ftree->el_deltaEtaSuperClusterTrackAtVtx.push_back(elec.deltaEtaSuperClusterTrackAtVtx());
        ftree->el_deltaPhiSuperClusterTrackAtVtx.push_back(elec.deltaPhiSuperClusterTrackAtVtx());
        ftree->el_deltaEtaSeedClusterTrackAtCalo.push_back(elec.deltaEtaSeedClusterTrackAtCalo());
        ftree->el_deltaPhiSeedClusterTrackAtCalo.push_back(elec.deltaPhiSeedClusterTrackAtCalo());
        ftree->el_eSuperClusterOverP.push_back(elec.eSuperClusterOverP());

        const auto el = electrons->ptrAt(ie);
        ftree->el_mvaNonTrigV0.push_back((*mvaValues)[el]);
        ftree->el_mvaNonTrigCat.push_back((*mvaCategories)[el]);

        ftree->el_mediumMVAId.push_back((*medium_mvaid_decisions)[el]);
        ftree->el_tightMVAId.push_back((*tight_mvaid_decisions)[el]);

        ftree->el_vetoCBId.push_back((*veto_cbid_decisions)[el]);
        ftree->el_looseCBId.push_back((*loose_cbid_decisions)[el]);
        ftree->el_mediumCBId.push_back((*medium_cbid_decisions)[el]);
        ftree->el_tightCBId.push_back((*tight_cbid_decisions)[el]);
        ftree->el_heepCBId.push_back((*heep_cbid_decisions)[el]);		

        ftree->el_ecalEnergy.push_back(elec.ecalEnergy());
        ftree->el_correctedEcalEnergy.push_back(elec.correctedEcalEnergy());
        ftree->el_correctedEcalEnergyError.push_back(elec.correctedEcalEnergyError());
        ftree->el_trackMomentumError.push_back(elec.trackMomentumError());

        ftree->el_hcalOverEcal.push_back(elec.hcalOverEcal());
        ftree->el_hcalOverEcalBc.push_back(elec.hcalOverEcalBc());
        ftree->el_hcalDepth1OverEcal.push_back(elec.hcalDepth1OverEcal());
        ftree->el_hcalDepth2OverEcal.push_back(elec.hcalDepth2OverEcal());
        ftree->el_eSeedClusterOverPout.push_back(elec.eSeedClusterOverPout());
        ftree->el_eSeedClusterOverP.push_back(elec.eSeedClusterOverP());
        ftree->el_eEleClusterOverPout.push_back(elec.eEleClusterOverPout());
        ftree->el_deltaEtaEleClusterTrackAtCalo.push_back(elec.deltaEtaEleClusterTrackAtCalo());
        ftree->el_deltaPhiEleClusterTrackAtCalo.push_back(elec.deltaPhiEleClusterTrackAtCalo());

        ftree->el_neutralHadronIso.push_back(elec.neutralHadronIso());
        ftree->el_chargedHadronIso.push_back(elec.chargedHadronIso());
        ftree->el_puChargedHadronIso.push_back(elec.puChargedHadronIso());
        ftree->el_ecalIso.push_back(elec.ecalIso());
        ftree->el_hcalIso.push_back(elec.hcalIso());
        ftree->el_particleIso.push_back(elec.particleIso());
        ftree->el_photonIso.push_back(elec.photonIso());
        ftree->el_trackIso.push_back(elec.trackIso());

        ftree->el_pfIso_sumChargedHadronPt.push_back(elec.pfIsolationVariables().sumChargedHadronPt);
        ftree->el_pfIso_sumNeutralHadronEt.push_back(elec.pfIsolationVariables().sumNeutralHadronEt);
        ftree->el_pfIso_sumPhotonEt.push_back(elec.pfIsolationVariables().sumPhotonEt);
        ftree->el_pfIso_sumPUPt.push_back(elec.pfIsolationVariables().sumPUPt);

        ftree->el_ecalPFClusterIso.push_back(elec.ecalPFClusterIso());
        ftree->el_hcalPFClusterIso.push_back(elec.hcalPFClusterIso());

        ftree->el_dr03EcalRecHitSumEt.push_back(elec.dr03EcalRecHitSumEt());
        ftree->el_dr03HcalTowerSumEt.push_back(elec.dr03HcalTowerSumEt());
        ftree->el_dr03HcalDepth1TowerSumEt.push_back(elec.dr03HcalDepth1TowerSumEt());
        ftree->el_dr03HcalDepth2TowerSumEt.push_back(elec.dr03HcalDepth2TowerSumEt());
        ftree->el_dr03TkSumPt.push_back(elec.dr03TkSumPt());

        ftree->el_dr04EcalRecHitSumEt.push_back(elec.dr04EcalRecHitSumEt());
        ftree->el_dr04HcalTowerSumEt.push_back(elec.dr04HcalTowerSumEt());
        ftree->el_dr04HcalDepth1TowerSumEt.push_back(elec.dr04HcalDepth1TowerSumEt());
        ftree->el_dr04HcalDepth2TowerSumEt.push_back(elec.dr04HcalDepth2TowerSumEt());
        ftree->el_dr04TkSumPt.push_back(elec.dr04TkSumPt());

        ftree->el_vx.push_back(elec.vx());
        ftree->el_vy.push_back(elec.vy());
        ftree->el_vz.push_back(elec.vz());

        ftree->el_hasGsfTrack.push_back(hasGsfTrack);

        ftree->el_passConversionVeto.push_back(elec.passConversionVeto());

        if( !isData_ )
	  {
	     // Internal matching
	     reco::GenParticle *genp = new reco::GenParticle();

	     float drmin;
	     bool hasMCMatch = mc_truth->doMatch(iEvent,iSetup,genParticlesHandle,genp,drmin,
						 elec.pt(),elec.eta(),elec.phi(),elec.pdgId());
	     ftree->el_hasMCMatch.push_back(hasMCMatch);
	     if( hasMCMatch )
	       {
		  ftree->el_gen_pt.push_back(genp->pt());
		  ftree->el_gen_eta.push_back(genp->eta());
		  ftree->el_gen_phi.push_back(genp->phi());
		  ftree->el_gen_m.push_back(genp->mass());
		  ftree->el_gen_status.push_back(genp->status());
		  ftree->el_gen_id.push_back(genp->pdgId());
		  ftree->el_gen_charge.push_back(genp->charge());
		  ftree->el_gen_dr.push_back(drmin);
	       }
	     else
	       {
		  ftree->el_gen_pt.push_back(-666);
		  ftree->el_gen_eta.push_back(-666);
		  ftree->el_gen_phi.push_back(-666);
		  ftree->el_gen_m.push_back(-666);
		  ftree->el_gen_status.push_back(-666);
		  ftree->el_gen_id.push_back(-666);
		  ftree->el_gen_charge.push_back(-666);
		  ftree->el_gen_dr.push_back(-666);
	       }
	     delete genp;

	     // PAT matching
	     const reco::GenParticle *genpPAT = elec.genParticle();
	     bool hasMCMatchPAT = (genpPAT != 0);
	     ftree->el_hasMCMatchPAT.push_back(hasMCMatchPAT);
	     if( hasMCMatchPAT )
	       {
		  ftree->el_genPAT_pt.push_back(genpPAT->pt());
		  ftree->el_genPAT_eta.push_back(genpPAT->eta());
		  ftree->el_genPAT_phi.push_back(genpPAT->phi());
		  ftree->el_genPAT_m.push_back(genpPAT->mass());
		  ftree->el_genPAT_status.push_back(genpPAT->status());
		  ftree->el_genPAT_id.push_back(genpPAT->pdgId());
		  ftree->el_genPAT_charge.push_back(genpPAT->charge());
	       }
	     else
	       {
		  ftree->el_genPAT_pt.push_back(-666);
		  ftree->el_genPAT_eta.push_back(-666);
		  ftree->el_genPAT_phi.push_back(-666);
		  ftree->el_genPAT_m.push_back(-666);
		  ftree->el_genPAT_status.push_back(-666);
		  ftree->el_genPAT_id.push_back(-666);
		  ftree->el_genPAT_charge.push_back(-666);
	       }
	  }
	
	bool allowCkfMatch = true;
        float lxyMin = 2.0;
        float probMin = 1e-6;
        uint nHitsBeforeVtxMax = 0;

        if( dataFormat_ != "AOD" )
	  {
	     bool matchConv = 0;
	     if( &beamspot ) matchConv = ConversionTools::hasMatchedConversion(elec,hConversions,beamspot.position(),allowCkfMatch,lxyMin,probMin,nHitsBeforeVtxMax);
	     ftree->el_hasMatchedConversion.push_back(matchConv);
	  }
	
	if( ftree->el_looseCBId.back() ) nElecPass++;
     }   
   ftree->el_n = ftree->el_pt.size();

   int nMuon = muons->size();
   int nMuonPass = 0;
   
   for(int im=0;im<nMuon;im++)
     {
        const pat::Muon& muon = muons->at(im);

        // Skimming muons with pT < 5 GeV.
        if (muon.pt() < 5) continue;
	
        ftree->mu_pt.push_back(muon.pt());
        ftree->mu_eta.push_back(muon.eta());
        ftree->mu_phi.push_back(muon.phi());
        ftree->mu_m.push_back(muon.mass());
        ftree->mu_E.push_back(muon.energy());
        ftree->mu_id.push_back(muon.pdgId());
        ftree->mu_charge.push_back(muon.charge());

        // IP
        ftree->mu_ip3d.push_back(muon.dB(pat::Muon::PV3D));
        ftree->mu_ip3dErr.push_back(muon.edB(pat::Muon::PV3D));
        ftree->mu_ip2d.push_back(muon.dB(pat::Muon::PV2D));
        ftree->mu_ip2dErr.push_back(muon.edB(pat::Muon::PV2D));
        ftree->mu_ip3dBS.push_back(muon.dB(pat::Muon::BS3D));
        ftree->mu_ip3dBSErr.push_back(muon.edB(pat::Muon::BS3D));
        ftree->mu_ip2dBS.push_back(muon.dB(pat::Muon::BS2D));
        ftree->mu_ip2dBSErr.push_back(muon.edB(pat::Muon::BS2D));

        const reco::MuonQuality combQuality = muon.combinedQuality();
        ftree->mu_combinedQuality_chi2LocalPosition.push_back(combQuality.chi2LocalPosition);
        ftree->mu_combinedQuality_trkKink.push_back(combQuality.trkKink);

        ftree->mu_numberOfMatches.push_back(muon.isMatchesValid() ? muon.numberOfMatches() : -666);
        ftree->mu_numberOfMatchedStations.push_back(muon.numberOfMatchedStations());

        // GlobalTrack
        const reco::TrackRef globalTrack = muon.globalTrack();
        bool hasGlobalTrack = globalTrack.isNonnull();

        ftree->mu_hasGlobalTrack.push_back(hasGlobalTrack);
        ftree->mu_globalTrack_d0.push_back((hasGlobalTrack) ? globalTrack->d0() : -666);
        ftree->mu_globalTrack_z0.push_back((hasGlobalTrack) ? globalTrack->dz() : -666);
        ftree->mu_globalTrack_d0Error.push_back((hasGlobalTrack) ? globalTrack->d0Error() : -666);
        ftree->mu_globalTrack_z0Error.push_back((hasGlobalTrack) ? globalTrack->dzError() : -666);
        ftree->mu_globalTrack_PV_dxy.push_back((hasGlobalTrack) ? globalTrack->dxy(primVtx->position()) : -666);
        ftree->mu_globalTrack_PV_dz.push_back((hasGlobalTrack) ? globalTrack->dz(primVtx->position()) : -666);
        ftree->mu_globalTrack_RP_dxy.push_back((hasGlobalTrack) ? globalTrack->dxy(globalTrack->referencePoint()) : -666);
        ftree->mu_globalTrack_RP_dz.push_back((hasGlobalTrack) ? globalTrack->dz(globalTrack->referencePoint()) : -666);
        ftree->mu_globalTrack_BS_dxy.push_back((hasGlobalTrack) ? globalTrack->dxy(beamspot.position()) : -666);
        ftree->mu_globalTrack_BS_dz.push_back((hasGlobalTrack) ? globalTrack->dz(beamspot.position()) : -666);
        ftree->mu_globalTrack_dxyError.push_back((hasGlobalTrack) ? globalTrack->dxyError() : -666);
        ftree->mu_globalTrack_dzError.push_back((hasGlobalTrack) ? globalTrack->dzError() : -666);
        ftree->mu_globalTrack_normalizedChi2.push_back((hasGlobalTrack) ? globalTrack->normalizedChi2() : -666);
        ftree->mu_globalTrack_numberOfValidHits.push_back((hasGlobalTrack) ? globalTrack->numberOfValidHits() : -666);
        ftree->mu_globalTrack_numberOfValidMuonHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfValidMuonHits() : -666);
        ftree->mu_globalTrack_numberOfLostHits.push_back((hasGlobalTrack) ? globalTrack->numberOfLostHits() : -666);
        ftree->mu_globalTrack_pt.push_back((hasGlobalTrack) ? globalTrack->pt() : -666);
        ftree->mu_globalTrack_eta.push_back((hasGlobalTrack) ? globalTrack->eta() : -666);
        ftree->mu_globalTrack_phi.push_back((hasGlobalTrack) ? globalTrack->phi() : -666);
        ftree->mu_globalTrack_ptError.push_back((hasGlobalTrack) ? globalTrack->ptError() : -666);
        ftree->mu_globalTrack_etaError.push_back((hasGlobalTrack) ? globalTrack->etaError() : -666);
        ftree->mu_globalTrack_phiError.push_back((hasGlobalTrack) ? globalTrack->phiError() : -666);
        ftree->mu_globalTrack_vx.push_back((hasGlobalTrack) ? globalTrack->vx() : -666);
        ftree->mu_globalTrack_vy.push_back((hasGlobalTrack) ? globalTrack->vy() : -666);
        ftree->mu_globalTrack_vz.push_back((hasGlobalTrack) ? globalTrack->vz() : -666);
        ftree->mu_globalTrack_qoverp.push_back((hasGlobalTrack) ? globalTrack->qoverp() : -666);
        ftree->mu_globalTrack_qoverpError.push_back((hasGlobalTrack) ? globalTrack->qoverpError() : -666);
        ftree->mu_globalTrack_charge.push_back((hasGlobalTrack) ? globalTrack->charge() : -666);
        ftree->mu_globalTrack_trackerLayersWithMeasurement.push_back((hasGlobalTrack) ? globalTrack->hitPattern().trackerLayersWithMeasurement() : -666);
        ftree->mu_globalTrack_pixelLayersWithMeasurement.push_back((hasGlobalTrack) ? globalTrack->hitPattern().pixelLayersWithMeasurement() : -666);
        ftree->mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo() : -666);
        ftree->mu_globalTrack_trackerLayersWithoutMeasurement.push_back((hasGlobalTrack) ? globalTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_globalTrack_numberOfValidPixelHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfValidPixelHits() : -666);
        ftree->mu_globalTrack_numberOfLostPixelHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_globalTrack_numberOfInnerHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) : -666);
        ftree->mu_globalTrack_numberOfOuterHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_OUTER_HITS) : -666);
        ftree->mu_globalTrack_validFraction.push_back((hasGlobalTrack) ? globalTrack->validFraction() : -666);

        // BestTrack
        ftree->mu_bestTrackType.push_back(muon.muonBestTrackType());
        const reco::TrackRef bestTrack = muon.muonBestTrack();
        bool hasBestTrack = bestTrack.isNonnull();

        ftree->mu_hasBestTrack.push_back(hasBestTrack);
        ftree->mu_bestTrack_d0.push_back((hasBestTrack) ? bestTrack->d0() : -666);
        ftree->mu_bestTrack_z0.push_back((hasBestTrack) ? bestTrack->dz() : -666);
        ftree->mu_bestTrack_d0Error.push_back((hasBestTrack) ? bestTrack->d0Error() : -666);
        ftree->mu_bestTrack_z0Error.push_back((hasBestTrack) ? bestTrack->dzError() : -666);
        ftree->mu_bestTrack_PV_dxy.push_back((hasBestTrack) ? bestTrack->dxy(primVtx->position()) : -666);
        ftree->mu_bestTrack_PV_dz.push_back((hasBestTrack) ? bestTrack->dz(primVtx->position()) : -666);
        ftree->mu_bestTrack_RP_dxy.push_back((hasBestTrack) ? bestTrack->dxy(bestTrack->referencePoint()) : -666);
        ftree->mu_bestTrack_RP_dz.push_back((hasBestTrack) ? bestTrack->dz(bestTrack->referencePoint()) : -666);
        ftree->mu_bestTrack_BS_dxy.push_back((hasBestTrack) ? bestTrack->dxy(beamspot.position()) : -666);
        ftree->mu_bestTrack_BS_dz.push_back((hasBestTrack) ? bestTrack->dz(beamspot.position()) : -666);
        ftree->mu_bestTrack_dxyError.push_back((hasBestTrack) ? bestTrack->dxyError() : -666);
        ftree->mu_bestTrack_dzError.push_back((hasBestTrack) ? bestTrack->dzError() : -666);
        ftree->mu_bestTrack_normalizedChi2.push_back((hasBestTrack) ? bestTrack->normalizedChi2() : -666);
        ftree->mu_bestTrack_numberOfValidHits.push_back((hasBestTrack) ? bestTrack->numberOfValidHits() : -666);
        ftree->mu_bestTrack_numberOfLostHits.push_back((hasBestTrack) ? bestTrack->numberOfLostHits() : -666);
        ftree->mu_bestTrack_pt.push_back((hasBestTrack) ? bestTrack->pt() : -666);
        ftree->mu_bestTrack_eta.push_back((hasBestTrack) ? bestTrack->eta() : -666);
        ftree->mu_bestTrack_phi.push_back((hasBestTrack) ? bestTrack->phi() : -666);
        ftree->mu_bestTrack_ptError.push_back((hasBestTrack) ? bestTrack->ptError() : -666);
        ftree->mu_bestTrack_etaError.push_back((hasBestTrack) ? bestTrack->etaError() : -666);
        ftree->mu_bestTrack_phiError.push_back((hasBestTrack) ? bestTrack->phiError() : -666);
        ftree->mu_bestTrack_vx.push_back((hasBestTrack) ? bestTrack->vx() : -666);
        ftree->mu_bestTrack_vy.push_back((hasBestTrack) ? bestTrack->vy() : -666);
        ftree->mu_bestTrack_vz.push_back((hasBestTrack) ? bestTrack->vz() : -666);
        ftree->mu_bestTrack_qoverp.push_back((hasBestTrack) ? bestTrack->qoverp() : -666);
        ftree->mu_bestTrack_qoverpError.push_back((hasBestTrack) ? bestTrack->qoverpError() : -666);
        ftree->mu_bestTrack_charge.push_back((hasBestTrack) ? bestTrack->charge() : -666);
        ftree->mu_bestTrack_trackerLayersWithMeasurement.push_back((hasBestTrack) ? bestTrack->hitPattern().trackerLayersWithMeasurement() : -666);
        ftree->mu_bestTrack_pixelLayersWithMeasurement.push_back((hasBestTrack) ? bestTrack->hitPattern().pixelLayersWithMeasurement() : -666);
        ftree->mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo() : -666);
        ftree->mu_bestTrack_trackerLayersWithoutMeasurement.push_back((hasBestTrack) ? bestTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_bestTrack_numberOfValidPixelHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfValidPixelHits() : -666);
        ftree->mu_bestTrack_numberOfLostPixelHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_bestTrack_numberOfInnerHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) : -666);
        ftree->mu_bestTrack_numberOfOuterHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_OUTER_HITS) : -666);
        ftree->mu_bestTrack_validFraction.push_back((hasBestTrack) ? bestTrack->validFraction() : -666);

        // InnerTrack
        const reco::TrackRef innerTrack = muon.innerTrack();
        bool hasInnerTrack = innerTrack.isNonnull();

        ftree->mu_hasInnerTrack.push_back(hasInnerTrack);
        ftree->mu_innerTrack_d0.push_back((hasInnerTrack) ? innerTrack->d0() : -666);
        ftree->mu_innerTrack_z0.push_back((hasInnerTrack) ? innerTrack->dz() : -666);
        ftree->mu_innerTrack_d0Error.push_back((hasInnerTrack) ? innerTrack->d0Error() : -666);
        ftree->mu_innerTrack_z0Error.push_back((hasInnerTrack) ? innerTrack->dzError() : -666);
        ftree->mu_innerTrack_PV_dxy.push_back((hasInnerTrack) ? innerTrack->dxy(primVtx->position()) : -666);
        ftree->mu_innerTrack_PV_dz.push_back((hasInnerTrack) ? innerTrack->dz(primVtx->position()) : -666);
        ftree->mu_innerTrack_RP_dxy.push_back((hasInnerTrack) ? innerTrack->dxy(bestTrack->referencePoint()) : -666);
        ftree->mu_innerTrack_RP_dz.push_back((hasInnerTrack) ? innerTrack->dz(bestTrack->referencePoint()) : -666);
        ftree->mu_innerTrack_BS_dxy.push_back((hasInnerTrack) ? innerTrack->dxy(beamspot.position()) : -666);
        ftree->mu_innerTrack_BS_dz.push_back((hasInnerTrack) ? innerTrack->dz(beamspot.position()) : -666);
        ftree->mu_innerTrack_dxyError.push_back((hasInnerTrack) ? innerTrack->dxyError() : -666);
        ftree->mu_innerTrack_dzError.push_back((hasInnerTrack) ? innerTrack->dzError() : -666);
        ftree->mu_innerTrack_normalizedChi2.push_back((hasInnerTrack) ? innerTrack->normalizedChi2() : -666);
        ftree->mu_innerTrack_numberOfValidHits.push_back((hasInnerTrack) ? innerTrack->numberOfValidHits() : -666);
        ftree->mu_innerTrack_numberOfLostHits.push_back((hasInnerTrack) ? innerTrack->numberOfLostHits() : -666);
        ftree->mu_innerTrack_pt.push_back((hasInnerTrack) ? innerTrack->pt() : -666);
        ftree->mu_innerTrack_eta.push_back((hasInnerTrack) ? innerTrack->eta() : -666);
        ftree->mu_innerTrack_phi.push_back((hasInnerTrack) ? innerTrack->phi() : -666);
        ftree->mu_innerTrack_ptError.push_back((hasInnerTrack) ? innerTrack->ptError() : -666);
        ftree->mu_innerTrack_etaError.push_back((hasInnerTrack) ? innerTrack->etaError() : -666);
        ftree->mu_innerTrack_phiError.push_back((hasInnerTrack) ? innerTrack->phiError() : -666);
        ftree->mu_innerTrack_vx.push_back((hasInnerTrack) ? innerTrack->vx() : -666);
        ftree->mu_innerTrack_vy.push_back((hasInnerTrack) ? innerTrack->vy() : -666);
        ftree->mu_innerTrack_vz.push_back((hasInnerTrack) ? innerTrack->vz() : -666);
        ftree->mu_innerTrack_qoverp.push_back((hasInnerTrack) ? innerTrack->qoverp() : -666);
        ftree->mu_innerTrack_qoverpError.push_back((hasInnerTrack) ? innerTrack->qoverpError() : -666);
        ftree->mu_innerTrack_charge.push_back((hasInnerTrack) ? innerTrack->charge() : -666);
        ftree->mu_innerTrack_trackerLayersWithMeasurement.push_back((hasInnerTrack) ? innerTrack->hitPattern().trackerLayersWithMeasurement() : -666);
        ftree->mu_innerTrack_pixelLayersWithMeasurement.push_back((hasInnerTrack) ? innerTrack->hitPattern().pixelLayersWithMeasurement() : -666);
        ftree->mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo() : -666);
        ftree->mu_innerTrack_trackerLayersWithoutMeasurement.push_back((hasInnerTrack) ? innerTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_innerTrack_numberOfValidPixelHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfValidPixelHits() : -666);
        ftree->mu_innerTrack_numberOfLostPixelHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_innerTrack_numberOfInnerHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) : -666);
        ftree->mu_innerTrack_numberOfOuterHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_OUTER_HITS) : -666);
        ftree->mu_innerTrack_validFraction.push_back((hasInnerTrack) ? innerTrack->validFraction() : -666);

        // PF Isolation
        reco::MuonPFIsolation pfR03 = muon.pfIsolationR03();
        ftree->mu_pfIso03_sumChargedHadronPt.push_back(pfR03.sumChargedHadronPt);
        ftree->mu_pfIso03_sumChargedParticlePt.push_back(pfR03.sumChargedParticlePt);
        ftree->mu_pfIso03_sumNeutralHadronEt.push_back(pfR03.sumNeutralHadronEt);
        ftree->mu_pfIso03_sumNeutralHadronEtHighThreshold.push_back(pfR03.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfIso03_sumPhotonEt.push_back(pfR03.sumPhotonEt);
        ftree->mu_pfIso03_sumPhotonEtHighThreshold.push_back(pfR03.sumPhotonEtHighThreshold);
        ftree->mu_pfIso03_sumPUPt.push_back(pfR03.sumPUPt);

        reco::MuonPFIsolation pfR04 = muon.pfIsolationR04();
        ftree->mu_pfIso04_sumChargedHadronPt.push_back(pfR04.sumChargedHadronPt);
        ftree->mu_pfIso04_sumChargedParticlePt.push_back(pfR04.sumChargedParticlePt);
        ftree->mu_pfIso04_sumNeutralHadronEt.push_back(pfR04.sumNeutralHadronEt);
        ftree->mu_pfIso04_sumNeutralHadronEtHighThreshold.push_back(pfR04.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfIso04_sumPhotonEt.push_back(pfR04.sumPhotonEt);
        ftree->mu_pfIso04_sumPhotonEtHighThreshold.push_back(pfR04.sumPhotonEtHighThreshold);
        ftree->mu_pfIso04_sumPUPt.push_back(pfR04.sumPUPt);

        reco::MuonPFIsolation pfMeanR03 = muon.pfMeanDRIsoProfileR03();
        ftree->mu_pfMeanIso03_sumChargedHadronPt.push_back(pfMeanR03.sumChargedHadronPt);
        ftree->mu_pfMeanIso03_sumChargedParticlePt.push_back(pfMeanR03.sumChargedParticlePt);
        ftree->mu_pfMeanIso03_sumNeutralHadronEt.push_back(pfMeanR03.sumNeutralHadronEt);
        ftree->mu_pfMeanIso03_sumNeutralHadronEtHighThreshold.push_back(pfMeanR03.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfMeanIso03_sumPhotonEt.push_back(pfMeanR03.sumPhotonEt);
        ftree->mu_pfMeanIso03_sumPhotonEtHighThreshold.push_back(pfMeanR03.sumPhotonEtHighThreshold);
        ftree->mu_pfMeanIso03_sumPUPt.push_back(pfMeanR03.sumPUPt);

        reco::MuonPFIsolation pfSumR03 = muon.pfSumDRIsoProfileR03();
        ftree->mu_pfSumIso03_sumChargedHadronPt.push_back(pfSumR03.sumChargedHadronPt);
        ftree->mu_pfSumIso03_sumChargedParticlePt.push_back(pfSumR03.sumChargedParticlePt);
        ftree->mu_pfSumIso03_sumNeutralHadronEt.push_back(pfSumR03.sumNeutralHadronEt);
        ftree->mu_pfSumIso03_sumNeutralHadronEtHighThreshold.push_back(pfSumR03.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfSumIso03_sumPhotonEt.push_back(pfSumR03.sumPhotonEt);
        ftree->mu_pfSumIso03_sumPhotonEtHighThreshold.push_back(pfSumR03.sumPhotonEtHighThreshold);
        ftree->mu_pfSumIso03_sumPUPt.push_back(pfSumR03.sumPUPt);

        reco::MuonPFIsolation pfMeanR04 = muon.pfMeanDRIsoProfileR04();
        ftree->mu_pfMeanIso04_sumChargedHadronPt.push_back(pfMeanR04.sumChargedHadronPt);
        ftree->mu_pfMeanIso04_sumChargedParticlePt.push_back(pfMeanR04.sumChargedParticlePt);
        ftree->mu_pfMeanIso04_sumNeutralHadronEt.push_back(pfMeanR04.sumNeutralHadronEt);
        ftree->mu_pfMeanIso04_sumNeutralHadronEtHighThreshold.push_back(pfMeanR04.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfMeanIso04_sumPhotonEt.push_back(pfMeanR04.sumPhotonEt);
        ftree->mu_pfMeanIso04_sumPhotonEtHighThreshold.push_back(pfMeanR04.sumPhotonEtHighThreshold);
        ftree->mu_pfMeanIso04_sumPUPt.push_back(pfMeanR04.sumPUPt);

        reco::MuonPFIsolation pfSumR04 = muon.pfSumDRIsoProfileR04();
        ftree->mu_pfSumIso04_sumChargedHadronPt.push_back(pfSumR04.sumChargedHadronPt);
        ftree->mu_pfSumIso04_sumChargedParticlePt.push_back(pfSumR04.sumChargedParticlePt);
        ftree->mu_pfSumIso04_sumNeutralHadronEt.push_back(pfSumR04.sumNeutralHadronEt);
        ftree->mu_pfSumIso04_sumNeutralHadronEtHighThreshold.push_back(pfSumR04.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfSumIso04_sumPhotonEt.push_back(pfSumR04.sumPhotonEt);
        ftree->mu_pfSumIso04_sumPhotonEtHighThreshold.push_back(pfSumR04.sumPhotonEtHighThreshold);
        ftree->mu_pfSumIso04_sumPUPt.push_back(pfSumR04.sumPUPt);

        ftree->mu_neutralHadronIso.push_back(muon.neutralHadronIso());
        ftree->mu_chargedHadronIso.push_back(muon.chargedHadronIso());
        ftree->mu_puChargedHadronIso.push_back(muon.puChargedHadronIso());
        ftree->mu_ecalIso.push_back(muon.ecalIso());
        ftree->mu_hcalIso.push_back(muon.hcalIso());
        ftree->mu_photonIso.push_back(muon.photonIso());
        ftree->mu_trackIso.push_back(muon.trackIso());

        // ID
        ftree->mu_isGlobalMuon.push_back(muon.isGlobalMuon());
        ftree->mu_isTrackerMuon.push_back(muon.isTrackerMuon());
        ftree->mu_isStandAloneMuon.push_back(muon.isStandAloneMuon());
        ftree->mu_isCaloMuon.push_back(muon.isCaloMuon());
        ftree->mu_isPFMuon.push_back(muon.isPFMuon());
        ftree->mu_isRPCMuon.push_back(muon.isRPCMuon());

        ftree->mu_isLooseMuon.push_back(muon.isLooseMuon());
        ftree->mu_isMediumMuon.push_back(muon.isMediumMuon());

        bool isTightMuon = 0;
        if( primVtx ) isTightMuon = muon.isTightMuon(*primVtx);
        ftree->mu_isTightMuon.push_back(isTightMuon);
        bool isSoftMuon = 0;
        if( primVtx ) isSoftMuon = muon.isSoftMuon(*primVtx);
        ftree->mu_isSoftMuon.push_back(isSoftMuon);
        bool isHighPtMuon = 0;
        if( primVtx ) isHighPtMuon = muon.isHighPtMuon(*primVtx);
        ftree->mu_isHighPtMuon.push_back(isHighPtMuon);

        ftree->mu_type.push_back(muon.type());

        ftree->mu_caloCompatibility.push_back(muon.caloCompatibility());
        ftree->mu_segmentCompatibility.push_back(muon.segmentCompatibility());

        // vertex
        ftree->mu_vx.push_back(muon.vx());
        ftree->mu_vy.push_back(muon.vy());
        ftree->mu_vz.push_back(muon.vz());

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis#Muon_identification
        ftree->mu_isGoodMuon_AllGlobalMuons.push_back(muon::isGoodMuon(muon,muon::AllGlobalMuons));
        ftree->mu_isGoodMuon_AllStandAloneMuons.push_back(muon::isGoodMuon(muon,muon::AllStandAloneMuons));
        ftree->mu_isGoodMuon_AllTrackerMuons.push_back(muon::isGoodMuon(muon,muon::AllTrackerMuons));
        ftree->mu_isGoodMuon_TrackerMuonArbitrated.push_back(muon::isGoodMuon(muon,muon::TrackerMuonArbitrated));
        ftree->mu_isGoodMuon_AllArbitrated.push_back(muon::isGoodMuon(muon,muon::AllArbitrated));
        ftree->mu_isGoodMuon_GlobalMuonPromptTight.push_back(muon::isGoodMuon(muon,muon::GlobalMuonPromptTight));
        ftree->mu_isGoodMuon_TMLastStationLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationLoose));
        ftree->mu_isGoodMuon_TMLastStationTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationTight));
        ftree->mu_isGoodMuon_TM2DCompatibilityLoose.push_back(muon::isGoodMuon(muon,muon::TM2DCompatibilityLoose));
        ftree->mu_isGoodMuon_TM2DCompatibilityTight.push_back(muon::isGoodMuon(muon,muon::TM2DCompatibilityTight));
        ftree->mu_isGoodMuon_TMOneStationLoose.push_back(muon::isGoodMuon(muon,muon::TMOneStationLoose));
        ftree->mu_isGoodMuon_TMOneStationTight.push_back(muon::isGoodMuon(muon,muon::TMOneStationTight));
        ftree->mu_isGoodMuon_TMLastStationOptimizedLowPtLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedLowPtLoose));
        ftree->mu_isGoodMuon_TMLastStationOptimizedLowPtTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedLowPtTight));
        ftree->mu_isGoodMuon_GMTkChiCompatibility.push_back(muon::isGoodMuon(muon,muon::GMTkChiCompatibility));
        ftree->mu_isGoodMuon_GMStaChiCompatibility.push_back(muon::isGoodMuon(muon,muon::GMStaChiCompatibility));
        ftree->mu_isGoodMuon_GMTkKinkTight.push_back(muon::isGoodMuon(muon,muon::GMTkKinkTight));
        ftree->mu_isGoodMuon_TMLastStationAngLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationAngLoose));
        ftree->mu_isGoodMuon_TMLastStationAngTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationAngTight));
        ftree->mu_isGoodMuon_TMOneStationAngLoose.push_back(muon::isGoodMuon(muon,muon::TMOneStationAngLoose));
        ftree->mu_isGoodMuon_TMOneStationAngTight.push_back(muon::isGoodMuon(muon,muon::TMOneStationAngTight));
        ftree->mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedBarrelLowPtLoose));
        ftree->mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedBarrelLowPtTight));

        bool energyIsValid = muon.isEnergyValid();
        ftree->mu_calEnergy_em.push_back(energyIsValid ? muon.calEnergy().em : -666.);
        ftree->mu_calEnergy_had.push_back(energyIsValid ? muon.calEnergy().had : -666.);
        ftree->mu_calEnergy_ho.push_back(energyIsValid ? muon.calEnergy().ho : -666.);
        ftree->mu_calEnergy_emS9.push_back(energyIsValid ? muon.calEnergy().emS9 : -666.);
        ftree->mu_calEnergy_hadS9.push_back(energyIsValid ? muon.calEnergy().hadS9 : -666.);
        ftree->mu_calEnergy_hoS9.push_back(energyIsValid ? muon.calEnergy().hoS9 : -666.);
        ftree->mu_calEnergy_emS25.push_back(energyIsValid ? muon.calEnergy().emS25 : -666.);
        ftree->mu_calEnergy_emMax.push_back(energyIsValid ? muon.calEnergy().emMax : -666.);
        ftree->mu_calEnergy_hadMax.push_back(energyIsValid ? muon.calEnergy().hadMax : -666.);
        ftree->mu_calEnergy_ecal_time.push_back(energyIsValid ? muon.calEnergy().ecal_time : -666.);
        ftree->mu_calEnergy_hcal_time.push_back(energyIsValid ? muon.calEnergy().hcal_time : -666.);
        ftree->mu_calEnergy_ecal_rawId.push_back(energyIsValid ? muon.calEnergy().ecal_id.rawId() : -666.);
        ftree->mu_calEnergy_hcal_rawId.push_back(energyIsValid ? muon.calEnergy().hcal_id.rawId() : -666.);

        bool isoIsValid = muon.isIsolationValid();

        ftree->mu_isolationR03_trackerVetoPt.push_back(isoIsValid ? muon.isolationR03().trackerVetoPt : -666.);
        ftree->mu_isolationR03_emVetoEt.push_back(isoIsValid ? muon.isolationR03().emVetoEt : -666.);
        ftree->mu_isolationR03_hadVetoEt.push_back(isoIsValid ? muon.isolationR03().hadVetoEt : -666.);
        ftree->mu_isolationR03_hoVetoEt.push_back(isoIsValid ? muon.isolationR03().hoVetoEt : -666.);
        ftree->mu_isolationR03_sumPt.push_back(isoIsValid ? muon.isolationR03().sumPt : -666.);
        ftree->mu_isolationR03_emEt.push_back(isoIsValid ? muon.isolationR03().emEt : -666.);
        ftree->mu_isolationR03_hadEt.push_back(isoIsValid ? muon.isolationR03().hadEt : -666.);
        ftree->mu_isolationR03_hoEt.push_back(isoIsValid ? muon.isolationR03().hoEt : -666.);
        ftree->mu_isolationR03_nTracks.push_back(isoIsValid ? muon.isolationR03().nTracks : -666.);
        ftree->mu_isolationR03_nJets.push_back(isoIsValid ? muon.isolationR03().nJets : -666.);

        ftree->mu_isolationR05_trackerVetoPt.push_back(isoIsValid ? muon.isolationR05().trackerVetoPt : -666.);
        ftree->mu_isolationR05_emVetoEt.push_back(isoIsValid ? muon.isolationR05().emVetoEt : -666.);
        ftree->mu_isolationR05_hadVetoEt.push_back(isoIsValid ? muon.isolationR05().hadVetoEt : -666.);
        ftree->mu_isolationR05_hoVetoEt.push_back(isoIsValid ? muon.isolationR05().hoVetoEt : -666.);
        ftree->mu_isolationR05_sumPt.push_back(isoIsValid ? muon.isolationR05().sumPt : -666.);
        ftree->mu_isolationR05_emEt.push_back(isoIsValid ? muon.isolationR05().emEt : -666.);
        ftree->mu_isolationR05_hadEt.push_back(isoIsValid ? muon.isolationR05().hadEt : -666.);
        ftree->mu_isolationR05_hoEt.push_back(isoIsValid ? muon.isolationR05().hoEt : -666.);
        ftree->mu_isolationR05_nTracks.push_back(isoIsValid ? muon.isolationR05().nTracks : -666.);
        ftree->mu_isolationR05_nJets.push_back(isoIsValid ? muon.isolationR05().nJets : -666.);

        if( !isData_ )
	  {
	     // Internal matching
	     reco::GenParticle *genp = new reco::GenParticle();

	     float drmin;
	     bool hasMCMatch = mc_truth->doMatch(iEvent,iSetup,genParticlesHandle,genp,drmin,
						 muon.pt(),muon.eta(),muon.phi(),muon.pdgId());
	     ftree->mu_hasMCMatch.push_back(hasMCMatch);
	     if( hasMCMatch )
	       {
		  ftree->mu_gen_pt.push_back(genp->pt());
		  ftree->mu_gen_eta.push_back(genp->eta());
		  ftree->mu_gen_phi.push_back(genp->phi());
		  ftree->mu_gen_m.push_back(genp->mass());
		  ftree->mu_gen_status.push_back(genp->status());
		  ftree->mu_gen_id.push_back(genp->pdgId());
		  ftree->mu_gen_charge.push_back(genp->charge());
		  ftree->mu_gen_dr.push_back(drmin);
	       }
	     else
	       {
		  ftree->mu_gen_pt.push_back(-666);
		  ftree->mu_gen_eta.push_back(-666);
		  ftree->mu_gen_phi.push_back(-666);
		  ftree->mu_gen_m.push_back(-666);
		  ftree->mu_gen_status.push_back(-666);
		  ftree->mu_gen_id.push_back(-666);
		  ftree->mu_gen_charge.push_back(-666);
		  ftree->mu_gen_dr.push_back(-666);
	       }
	     delete genp;

	     // PAT matching
	     const reco::GenParticle *genpPAT = muon.genParticle();
	     bool hasMCMatchPAT = (genpPAT != 0);
	     ftree->mu_hasMCMatchPAT.push_back(hasMCMatchPAT);
	     if( hasMCMatchPAT )
	       {
		  ftree->mu_genPAT_pt.push_back(genpPAT->pt());
		  ftree->mu_genPAT_eta.push_back(genpPAT->eta());
		  ftree->mu_genPAT_phi.push_back(genpPAT->phi());
		  ftree->mu_genPAT_m.push_back(genpPAT->mass());
		  ftree->mu_genPAT_status.push_back(genpPAT->status());
		  ftree->mu_genPAT_id.push_back(genpPAT->pdgId());
		  ftree->mu_genPAT_charge.push_back(genpPAT->charge());
	       }
	     else
	       {
		  ftree->mu_genPAT_pt.push_back(-666);
		  ftree->mu_genPAT_eta.push_back(-666);
		  ftree->mu_genPAT_phi.push_back(-666);
		  ftree->mu_genPAT_m.push_back(-666);
		  ftree->mu_genPAT_status.push_back(-666);
		  ftree->mu_genPAT_id.push_back(-666);
		  ftree->mu_genPAT_charge.push_back(-666);
	       }
	  }
	
	if( ftree->mu_isTightMuon.back() ) nMuonPass++;
     }
   ftree->mu_n = ftree->mu_pt.size();
   
   int nTau = taus->size();
   for(int it=0;it<nTau;it++)
     {
        const pat::Tau& tau = taus->at(it);

        // Skimming taus with pT < 5 GeV. (should do nothing for miniAOD where pT > 18 GeV is applied)
        if (tau.pt() < 5) continue;

        ftree->tau_pt.push_back(tau.pt());
        ftree->tau_eta.push_back(tau.eta());
        ftree->tau_phi.push_back(tau.phi());
        ftree->tau_m.push_back(tau.mass());
        ftree->tau_E.push_back(tau.energy());
        ftree->tau_id.push_back(tau.pdgId());
        ftree->tau_charge.push_back(tau.charge());

        // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Taus
        //
        float tau_leadingTrackPt = -666;
        int tau_leadingTrackCharge = -666;
        float tau_leadingTrackDz = -666;
        float tau_leadingTrackDxy = -666;

        ftree->tau_hasLeadChargedHadrCand.push_back(tau.leadChargedHadrCand().isNonnull());

        if( tau.leadChargedHadrCand().isNonnull() )
	  {
	     pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
	     tau_leadingTrackPt = packedLeadTauCand->pt();
	     tau_leadingTrackCharge = packedLeadTauCand->charge();
	     tau_leadingTrackDz = packedLeadTauCand->dz();
	     tau_leadingTrackDxy = packedLeadTauCand->dxy();
	  }
        ftree->tau_leadingTrackPt.push_back(tau_leadingTrackPt);
        ftree->tau_leadingTrackCharge.push_back(tau_leadingTrackCharge);
        ftree->tau_leadingTrackDz.push_back(tau_leadingTrackDz);
        ftree->tau_leadingTrackDxy.push_back(tau_leadingTrackDxy);

        ftree->tau_decayMode.push_back(tau.decayMode());
        ftree->tau_decayModeFindingOldDMs.push_back(tau.tauID("decayModeFinding")); // For Eric
        ftree->tau_decayModeFindingNewDMs.push_back(tau.tauID("decayModeFindingNewDMs"));

        ftree->tau_puCorrPtSum.push_back(tau.tauID("puCorrPtSum"));
        ftree->tau_neutralIsoPtSum.push_back(tau.tauID("neutralIsoPtSum"));
        ftree->tau_chargedIsoPtSum.push_back(tau.tauID("chargedIsoPtSum"));
        ftree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));

        ftree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
        ftree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")); // For Eric
        ftree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));

        ftree->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"));
        ftree->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"));
        ftree->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"));
        ftree->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"));

        ftree->tau_againstMuonLoose3.push_back(tau.tauID("againstMuonLoose3"));
        ftree->tau_againstMuonTight3.push_back(tau.tauID("againstMuonTight3"));

        ftree->tau_pfEssential_jet_pt.push_back(tau.pfEssential().p4Jet_.pt());
        ftree->tau_pfEssential_jet_eta.push_back(tau.pfEssential().p4Jet_.eta());
        ftree->tau_pfEssential_jet_phi.push_back(tau.pfEssential().p4Jet_.phi());
        ftree->tau_pfEssential_jet_m.push_back(tau.pfEssential().p4Jet_.mass());

        ftree->tau_pfEssential_jetCorr_pt.push_back(tau.pfEssential().p4CorrJet_.pt());
        ftree->tau_pfEssential_jetCorr_eta.push_back(tau.pfEssential().p4CorrJet_.eta());
        ftree->tau_pfEssential_jetCorr_phi.push_back(tau.pfEssential().p4CorrJet_.phi());
        ftree->tau_pfEssential_jetCorr_m.push_back(tau.pfEssential().p4CorrJet_.mass());

        float tau_pfEssential_sv_x = -666;
        float tau_pfEssential_sv_y = -666;
        float tau_pfEssential_sv_z = -666;

        ftree->tau_pfEssential_hasSV.push_back(tau.pfEssential().sv_.isNonnull());

        if( tau.pfEssential().sv_.isNonnull() )
	  {
	     tau_pfEssential_sv_x = tau.pfEssential().sv_->x();
	     tau_pfEssential_sv_y = tau.pfEssential().sv_->y();
	     tau_pfEssential_sv_z = tau.pfEssential().sv_->z();
	  }

        ftree->tau_pfEssential_sv_x.push_back(tau_pfEssential_sv_x);
        ftree->tau_pfEssential_sv_y.push_back(tau_pfEssential_sv_y);
        ftree->tau_pfEssential_sv_z.push_back(tau_pfEssential_sv_z);

        ftree->tau_pfEssential_flightLengthSig.push_back(tau.pfEssential().flightLengthSig_);
        ftree->tau_pfEssential_dxy.push_back(tau.pfEssential().dxy_);
        ftree->tau_pfEssential_dxy_error.push_back(tau.pfEssential().dxy_error_);
        ftree->tau_pfEssential_dxy_Sig.push_back(tau.pfEssential().dxy_Sig_);
    }
   ftree->tau_n = ftree->tau_pt.size();

   // ##########################
   // #       _      _         #
   // #      | | ___| |_ ___   #
   // #   _  | |/ _ \ __/ __|  #
   // #  | |_| |  __/ |_\__ \  #
   // #   \___/ \___|\__|___/  #
   // #                        #
   // ##########################
   //
   // Jets
   //
   int nJet = jets->size();
   ftree->jet_n = nJet;
   for(int ij=0;ij<nJet;ij++)
     {
        const pat::Jet& jet = jets->at(ij);
	
        //std::cout<< "jet pt"  << jet.pt() << " jet eta: " << jet.eta() << " jet phi " << jet.phi() << "jet E " << jet.energy() << " jet mass " << jet.mass() << std::endl;

        bool debug = false;

        if(debug)
	  {
	     std::cout << "jet["       << ij                      << "]: "
	       << " pt= "            << jet.pt()
		 << " eta= "           << jet.eta()
		   << " phi= "           << jet.phi()               << std::endl;
	  }
	
        ftree->jet_pt.push_back(jet.pt());
        ftree->jet_eta.push_back(jet.eta());
        ftree->jet_phi.push_back(jet.phi());
        ftree->jet_m.push_back(jet.mass());
        ftree->jet_E.push_back(jet.energy());
        ftree->jet_JBP.push_back(jet.bDiscriminator("pfJetBProbabilityBJetTags"));
        ftree->jet_JP.push_back(jet.bDiscriminator("pfJetProbabilityBJetTags"));
        ftree->jet_TCHP.push_back(jet.bDiscriminator("pfTrackCountingHighPurBJetTags"));
        ftree->jet_TCHE.push_back(jet.bDiscriminator("pfTrackCountingHighEffBJetTags"));
        ftree->jet_SSVHE.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
        ftree->jet_SSVHP.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
        ftree->jet_CMVA.push_back(jet.bDiscriminator("pfCombinedMVABJetTags"));

        ftree->jet_chargedMultiplicity.push_back(jet.chargedMultiplicity());
        ftree->jet_neutralMultiplicity.push_back(jet.neutralMultiplicity());
        ftree->jet_chargedHadronMultiplicity.push_back(jet.chargedHadronMultiplicity());

        ftree->jet_jecFactorUncorrected.push_back(jet.jecFactor("Uncorrected"));
        ftree->jet_jecFactorL1FastJet.push_back(jet.jecFactor("L1FastJet"));
        ftree->jet_jecFactorL2Relative.push_back(jet.jecFactor("L2Relative"));
        ftree->jet_jecFactorL3Absolute.push_back(jet.jecFactor("L3Absolute"));

        if(debug)
	  {
	     std::cout << "jet["       << ij                      << "]: "
	       << " Uncorrected= "         << jet.jecFactor("Uncorrected")
		 << " L1FastJet  = "         << jet.jecFactor("L1FastJet")
		   << " L2Relative = "         << jet.jecFactor("L2Relative")
		     << " L3Absolute = "         << jet.jecFactor("L3Absolute")
		       << std::endl;
	  }

        ftree->jet_jetArea.push_back(jet.jetArea());

        jecUnc->setJetEta(fabs(jet.eta()));
        jecUnc->setJetPt(jet.pt());

        ftree->jet_Unc.push_back(jecUnc->getUncertainty(true));

        ftree->jet_ntrk.push_back(jet.associatedTracks().size());
        //	std::cout << jet.hasTagInfo("pfInclusiveSecondaryVertexFinderTagInfos") << std::endl;

        float CSVIVF = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        ftree->jet_CSVv2.push_back(CSVIVF);

        ftree->jet_cMVAv2.push_back(jet.bDiscriminator("pfCombinedMVAV2BJetTags"));
        ftree->jet_CharmCvsL.push_back(jet.bDiscriminator("pfCombinedCvsLJetTags"));
        ftree->jet_CharmCvsB.push_back(jet.bDiscriminator("pfCombinedCvsBJetTags"));

        ftree->jet_partonFlavour.push_back(jet.partonFlavour());
        ftree->jet_hadronFlavour.push_back(jet.hadronFlavour());
	
        ftree->jet_neutralHadronEnergy.push_back(jet.neutralHadronEnergy());
        ftree->jet_neutralEmEnergy.push_back(jet.neutralEmEnergy());
        ftree->jet_chargedHadronEnergy.push_back(jet.chargedHadronEnergy());
        ftree->jet_chargedEmEnergy.push_back(jet.chargedEmEnergy());
        ftree->jet_electronEnergy.push_back(jet.electronEnergy());
        ftree->jet_muonEnergy.push_back(jet.muonEnergy());
        ftree->jet_photonEnergy.push_back(jet.photonEnergy());

        if( jet.hasUserFloat("pileupJetId:fullDiscriminant") )
	  ftree->jet_pileupJetId.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
        else
	  ftree->jet_pileupJetId.push_back(-666.);

        // Jet ID
        float NHF = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF = jet.chargedHadronEnergyFraction();
        float MUF = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
        float CHM = jet.chargedMultiplicity();
        float eta = jet.eta();
        bool looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta)>2.4);
        bool tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || fabs(eta)>2.4);

        ftree->jet_neutralHadronEnergyFraction.push_back(jet.neutralHadronEnergyFraction());
        ftree->jet_neutralEmEnergyFraction.push_back(jet.neutralEmEnergyFraction());
        ftree->jet_chargedHadronEnergyFraction.push_back(jet.chargedHadronEnergyFraction());
        ftree->jet_muonEnergyFraction.push_back(jet.muonEnergyFraction());
        ftree->jet_chargedEmEnergyFraction.push_back(jet.chargedEmEnergyFraction());

        ftree->jet_looseJetID.push_back(looseJetID);
        ftree->jet_tightJetID.push_back(tightJetID);

        const reco::GenJet* genJet = jet.genJet();
        bool hasGenInfo = (genJet);
        ftree->jet_hasGenJet.push_back(hasGenInfo);

        float gen_jet_pt = -666;
        float gen_jet_eta = -666;
        float gen_jet_phi = -666;
        float gen_jet_m = -666;
        float gen_jet_E = -666;
        int gen_jet_status = -666;
        int gen_jet_id = -666;

        if( hasGenInfo )
	  {
	     gen_jet_pt = genJet->pt();
	     gen_jet_eta = genJet->eta();
	     gen_jet_phi = genJet->phi();
	     gen_jet_m = genJet->mass();
	     gen_jet_E = genJet->energy();
	     gen_jet_status = genJet->status();
	     gen_jet_id = genJet->pdgId();
	  }
	
        ftree->jet_genJet_pt.push_back(gen_jet_pt);
        ftree->jet_genJet_eta.push_back(gen_jet_eta);
        ftree->jet_genJet_phi.push_back(gen_jet_phi);
        ftree->jet_genJet_m.push_back(gen_jet_m);
        ftree->jet_genJet_E.push_back(gen_jet_E);

        ftree->jet_genJet_status.push_back(gen_jet_status);
        ftree->jet_genJet_id.push_back(gen_jet_id);

        const reco::GenParticle* genParton = (!isData_) ? jet.genParton() : 0;
        bool hasGenPartonInfo = (genParton);
        ftree->jet_hasGenParton.push_back(hasGenPartonInfo);

        float gen_parton_pt = -666;
        float gen_parton_eta = -666;
        float gen_parton_phi = -666;
        float gen_parton_m = -666;
        float gen_parton_E = -666;
        int gen_parton_status = -666;
        int gen_parton_id = -666;

        if( hasGenPartonInfo )
	  {
	     gen_parton_pt = genParton->pt();
	     gen_parton_eta = genParton->eta();
	     gen_parton_phi = genParton->phi();
	     gen_parton_m = genParton->mass();
	     gen_parton_E = genParton->energy();
	     gen_parton_status = genParton->status();
	     gen_parton_id = genParton->pdgId();
	  }

        ftree->jet_genParton_pt.push_back(gen_parton_pt);
        ftree->jet_genParton_eta.push_back(gen_parton_eta);
        ftree->jet_genParton_phi.push_back(gen_parton_phi);
        ftree->jet_genParton_m.push_back(gen_parton_m);
        ftree->jet_genParton_E.push_back(gen_parton_E);

        ftree->jet_genParton_status.push_back(gen_parton_status);
        ftree->jet_genParton_id.push_back(gen_parton_id);
     }

   // GenJets
   //
   if( !isData_ )
     {
        int nGenJet = genJets->size();
        ftree->genJet_n = nGenJet;
        for(int ij=0;ij<nGenJet;ij++)
	  {
	     const reco::GenJet& genJet = genJets->at(ij);

	     ftree->genJet_pt.push_back(genJet.pt());
	     ftree->genJet_eta.push_back(genJet.eta());
	     ftree->genJet_phi.push_back(genJet.phi());
	     ftree->genJet_m.push_back(genJet.mass());
	     ftree->genJet_E.push_back(genJet.energy());

	     ftree->genJet_emEnergy.push_back(genJet.emEnergy());
	     ftree->genJet_hadEnergy.push_back(genJet.hadEnergy());
	     ftree->genJet_invisibleEnergy.push_back(genJet.invisibleEnergy());
	     ftree->genJet_auxiliaryEnergy.push_back(genJet.auxiliaryEnergy());

	     //	     RefToBase<reco::Jet> jetRef(RefToBaseProd<reco::Jet>(genJets),ij);
	     //	     int genJet_flavour = (*genJetFlavourMatching)[jetRef].getFlavour();
	     //	     ftree->genJet_flavour.push_back(genJet_flavour);
	     ftree->genJet_flavour.push_back(-666); // FIXME
	  }
     }
   
   if( (applyMETFilters_ && passMETFilters) || !applyMETFilters_ ) // MET filters
     {
	if( nElecPass+nMuonPass > 0 ) // skim
	  {	     
	     ftree->tree->Fill();
	  }
     }
   
   delete mc_truth;
}

// ------------ method called when starting to processes a run  ------------
void FlatTreeProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
   bool changed = true;
   if(!hltConfig_.init(iRun, iSetup, "HLT", changed))
     std::cout << "Warning, didn't find HLTConfigProvider with label "
     << "HLT" << " in run " << iRun.run() << std::endl;

   if(!hltPrescale_.init(iRun, iSetup, "HLT", changed))
     std::cout << "Warning, didn't find HLTPrescaleProvider with label "
     << "HLT" << " in run " << iRun.run() << std::endl;

    /*   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;

         iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
    //   iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];

    jecUnc = new JetCorrectionUncertainty(JetCorPar);*/

   const char* cmssw_base = std::getenv("CMSSW_BASE");
   std::string JECUncertaintyPath = std::string(cmssw_base)+"/src/FlatTree/FlatTreeProducer/data/jecFiles/Fall15_25nsV2_MC/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
   jecUnc = new JetCorrectionUncertainty(JECUncertaintyPath.c_str());
}

// ------------ method called when ending the processing of a run  ------------
void FlatTreeProducer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
   delete jecUnc;
}

void FlatTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

bool FlatTreeProducer::foundTrigger(const std::string& name) const
{
   for( unsigned int i=0;i<filterTriggerNames_.size();++i )
     {
        TString pattern(filterTriggerNames_[i]);
        pattern.ToLower();
        TRegexp reg(Form("%s",pattern.Data()),true);
        TString sname(name);
        sname.ToLower();
        if( sname.Index(reg) >= 0 ) return true;
     }

   return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlatTreeProducer);

