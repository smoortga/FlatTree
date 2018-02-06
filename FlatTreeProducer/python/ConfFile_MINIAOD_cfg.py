import FWCore.ParameterSet.Config as cms

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run on real data')
options.register('applyMETFilters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')
options.register('applyJEC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply JEC corrections')
options.register('fillMCScaleWeight',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Fill PDF weights')
options.register('fillPUInfo',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Fill PU info')
options.register('nPDF', -1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "nPDF")
options.register('confFile', 'conf.xml', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Flattree variables configuration")
options.register('bufferSize', 32000, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Buffer size for branches of the flat tree")
options.parseArguments()

##########################
#  Global configuration  #
##########################

process = cms.Process("FlatTree")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.suppressWarning = cms.untracked.vstring(["JetPtMismatchAtLowPt","NullTransverseMomentum"])

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

if options.isData:
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v5'
else:    
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    
corName="Summer16_23Sep2016V3_MC"
corTag="JetCorrectorParametersCollection_"+corName
if options.isData:
    corName="Summer16_23Sep2016AllV3_DATA"
    corTag="JetCorrectorParametersCollection_"+corName
dBFile=corName+".db"

process.load("CondCore.CondDB.CondDB_cfi")
process.jec = cms.ESSource("PoolDBESSource",
                           DBParameters = cms.PSet(
                           messageLevel = cms.untracked.int32(0)
                           ),
                           timetype = cms.string('runnumber'),
                           toGet = cms.VPSet(
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK4PF"),
                                    label  = cms.untracked.string('AK4PF')
                                    ),
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK4PFchs"),
                                    label  = cms.untracked.string('AK4PFchs')
                                    ),
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK8PF"),
                                    label  = cms.untracked.string('AK8PF')
                                    ),
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK8PFchs"),
                                    label  = cms.untracked.string('AK8PFchs')
                                    ),
                           ),
                           connect = cms.string("sqlite_file:"+dBFile)
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')                           

process.load('Configuration.StandardSequences.Services_cff')
if not options.isData:
    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                           calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                           engineName = cms.untracked.string('TRandom3'),
                                           ),
                                           calibratedPatPhotons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                           engineName = cms.untracked.string('TRandom3'),
                                           ),
    )

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

corList = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
if options.isData:
    corList = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

# Re-apply JEC to AK4
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', corList, 'None')
)

jetsNameAK4="selectedUpdatedPatJetsUpdatedJEC"
#jetsNameAK4="slimmedJets"

########################
#  Additional modules  #
########################

if not options.isData:
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

    from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
    process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
            particles = "prunedGenParticles"
    )

    from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
    process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
            jets = "slimmedGenJets"
    )

    from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
    process.matchGenBHadron = matchGenBHadron.clone(
            genParticles = "prunedGenParticles",
            jetFlavourInfos = "genJetFlavourInfos"
    )

    from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
    process.matchGenCHadron = matchGenCHadron.clone(
            genParticles = "prunedGenParticles",
            jetFlavourInfos = "genJetFlavourInfos"
    )

# egamma
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2

my_id_modules = [
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff'
]

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#####################
# MET Significance  #
#####################
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")
from RecoMET.METProducers.testInputFiles_cff import recoMETtestInputFiles

###########
#  Input  #
###########

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"), # WARNING / FIXME for test only !
    fileNames = cms.untracked.vstring(
#    '/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/00000/68C5F6A4-1999-E611-B338-02163E013B09.root'
    '/store/mc/RunIISummer16MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/0015BB42-9BAA-E611-8C7F-0CC47A7E0196.root'
        )
)

############
#  Output  #
############

process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root"))

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True)	 # needed for ak10 computation (JMEAnalysis/JetToolbox)
)
        
#############################
#  Flat Tree configuration  #
#############################

process.FlatTree = cms.EDAnalyzer('FlatTreeProducer',

                  dataFormat        = cms.string("MINIAOD"),

                  bufferSize        = cms.int32(options.bufferSize),
                  confFile          = cms.string(options.confFile),

                  isData            = cms.bool(options.isData),
                  applyMETFilters   = cms.bool(options.applyMETFilters),
                  fillMCScaleWeight = cms.bool(options.fillMCScaleWeight),
                  fillPUInfo	    = cms.bool(options.fillPUInfo),
                  nPDF              = cms.int32(options.nPDF),
                  
                  vertexInput              = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  electronInput            = cms.InputTag("slimmedElectrons"),
                  #electronPATInput         = cms.InputTag("slimmedElectrons"),
                  electronPATInput         = cms.InputTag("calibratedPatElectrons"),

#                  eleVetoCBIdMap           = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
#                  eleLooseCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
#                  eleMediumCBIdMap         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
#                  eleTightCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                  eleVetoCBIdMap           = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                  eleLooseCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                  eleMediumCBIdMap         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                  eleTightCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
                  eleHEEPCBIdMap           = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),

                  eleMediumMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                  eleTightMVAIdMap         = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                  mvaValuesMap             = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                  mvaCategoriesMap         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),

                  BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
                  BadChargedCandidateFilter  = cms.InputTag("BadChargedCandidateFilter",""),
                  
                  filterTriggerNames       = cms.untracked.vstring(
#                  "*"
                  # "HLT_Ele35_WPLoose_Gsf_v*",
#                   "HLT_Ele27_WPTight_Gsf_v*",
#                   "HLT_Ele32_eta2p1_WPTight_Gsf_v*",
#                   "HLT_IsoMu22_v*",
#                   "HLT_IsoTkMu22_v*",
#                   "HLT_IsoMu24_v*",
#                   "HLT_IsoTkMu24_v*"
                    #MuonEG
                    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                    #DoubleEG
                    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                    #DoubleMuon
                    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"
                  ),
                  
                  muonInput                = cms.InputTag("slimmedMuons"),
                  tauInput                 = cms.InputTag("slimmedTaus"),
                  jetInput                 = cms.InputTag(jetsNameAK4),
                  genJetInput              = cms.InputTag("slimmedGenJets"),
                  jetFlavorMatchTokenInput = cms.InputTag("jetFlavourMatch"),
                  metInput                 = cms.InputTag("slimmedMETs"),
                  metPuppiInput            = cms.InputTag("slimmedMETsPuppi"),
                  metNoHFInput             = cms.InputTag("slimmedMETsNoHF"),
                  metSigInput              = cms.InputTag("METSignificance"),
                  metCovInput              = cms.InputTag("METSignificance","METCovariance"),
                  rhoInput                 = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                  genParticlesInput        = cms.InputTag("prunedGenParticles"),
                  genEventInfoInput        = cms.InputTag("generator"),
                  LHEEventProductInput     = cms.InputTag("externalLHEProducer"),
                  bsInput                  = cms.InputTag("offlineBeamSpot"),
                  pfcandsInput             = cms.InputTag("packedPFCandidates"),
                  hConversionsInput        = cms.InputTag("reducedEgamma","reducedConversions"),
                  puInfoInput		   = cms.InputTag("slimmedAddPileupInfo"),
#                  puInfoInput		   = cms.InputTag("addPileupInfo"),
                  objects                  = cms.InputTag("selectedPatTrigger"),
                  
                  genTTXJets                    = cms.InputTag("slimmedGenJets"),
                  genTTXBHadJetIndex            = cms.InputTag("matchGenBHadron","genBHadJetIndex"),
                  genTTXBHadFlavour             = cms.InputTag("matchGenBHadron","genBHadFlavour"),
                  genTTXBHadFromTopWeakDecay    = cms.InputTag("matchGenBHadron","genBHadFromTopWeakDecay"),
                  genTTXBHadPlusMothers         = cms.InputTag("matchGenBHadron","genBHadPlusMothers"),
                  genTTXBHadPlusMothersIndices  = cms.InputTag("matchGenBHadron","genBHadPlusMothersIndices"),
                  genTTXBHadIndex               = cms.InputTag("matchGenBHadron","genBHadIndex"),
                  genTTXBHadLeptonHadronIndex   = cms.InputTag("matchGenBHadron","genBHadLeptonHadronIndex"),
                  genTTXBHadLeptonViaTau        = cms.InputTag("matchGenBHadron","genBHadLeptonViaTau"),
                  genTTXCHadJetIndex            = cms.InputTag("matchGenCHadron","genCHadJetIndex"),
                  genTTXCHadFlavour             = cms.InputTag("matchGenCHadron","genCHadFlavour"),
                  genTTXCHadFromTopWeakDecay    = cms.InputTag("matchGenCHadron","genCHadFromTopWeakDecay"),
                  genTTXCHadBHadronId           = cms.InputTag("matchGenCHadron","genCHadBHadronId")
)

##########
#  Path  #
##########

process.p = cms.Path(
                    process.calibratedPatElectrons+
                    process.electronMVAValueMapProducer+
                    process.egmGsfElectronIDSequence+
                    process.METSignificance+
                    process.BadChargedCandidateFilter+
                    process.BadPFMuonFilter+
                    process.FlatTree
                    )
