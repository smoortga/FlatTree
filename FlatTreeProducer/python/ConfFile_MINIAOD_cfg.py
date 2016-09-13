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

# remove verbose from patTrigger due to missing L1 prescales for some trigger paths
#process.MessageLogger.suppressWarning.append('patTrigger')
#process.MessageLogger.cerr.FwkJob.limit=1
#process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

if options.isData:
    process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8'
else:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'

corName="Spring16_25nsV6_DATA"
corTag="JetCorrectorParametersCollection_"+corName
#if options.isData:
#    corName="Fall15_25nsV2_DATA"
#    corTag="JetCorrectorParametersCollection_"+corName
dBFile=corName+".db"

if options.isData:
    process.load("CondCore.CondDB.CondDB_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
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

# egamma
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
my_id_modules = [
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V2_cff',
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff'
]

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer - added unwillingly by Xavier
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
                '/store/mc/RunIISpring16MiniAODv1/ttHToNonbb_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/50000/0ADF7BAE-0914-E611-B788-0025905A6068.root'
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

                  eleVetoCBIdMap           = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                  eleLooseCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                  eleMediumCBIdMap         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                  eleTightCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                  eleHEEPCBIdMap           = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),

                  eleMediumMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                  eleTightMVAIdMap         = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                  mvaValuesMap             = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                  mvaCategoriesMap         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),

                  BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
                  BadChargedCandidateFilter  = cms.InputTag("BadChargedCandidateFilter",""),
                  
                  filterTriggerNames       = cms.untracked.vstring(
#                  "*"
                  "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*",
                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
                  "HLT_DoubleMu33NoFiltersNoVtx_v*",
                  "HLT_DoubleMu38NoFiltersNoVtx_v*",
                  "HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
                  "HLT_Ele22_eta2p1_WPTight_Gsf_v*",
                  "HLT_Ele22_eta2p1_WP75_Gsf_v*",
                  "HLT_Ele23_WPLoose_Gsf_v*",
                  "HLT_Ele23_WP75_Gsf_v*",
                  "HLT_Ele27_eta2p1_WP75_Gsf_v*",
                  "HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v*",
                  "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
                  "HLT_Ele27_eta2p1_WPTight_Gsf_v*",
                  "HLT_Ele32_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v*",
                  "HLT_Ele32_eta2p1_WPLoose_Gsf_v*",
                  "HLT_Ele32_eta2p1_WPTight_Gsf_v*",
                  "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
                  "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*",
                  "HLT_IsoMu17_eta2p1_v*",
                  "HLT_DoubleIsoMu17_eta2p1_v*",
                  "HLT_IsoMu18_v*",
                  "HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v*",
                  "HLT_IsoMu20_v*",
                  "HLT_IsoMu20_eta2p1_v*",
                  "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v*",
                  "HLT_IsoMu24_eta2p1_v*",
                  "HLT_IsoMu27_v*",
                  "HLT_IsoTkMu20_v*",
                  "HLT_IsoTkMu20_eta2p1_v*",
                  "HLT_IsoTkMu24_eta2p1_v*",
                  "HLT_IsoTkMu27_v*",
                  "HLT_Mu17_Mu8_v*",
                  "HLT_Mu17_Mu8_DZ_v*",
                  "HLT_Mu17_Mu8_SameSign_DZ_v*",
                  "HLT_Mu20_Mu10_v*",
                  "HLT_Mu20_Mu10_DZ_v*",
                  "HLT_Mu20_Mu10_SameSign_DZ_v*",
                  "HLT_Mu17_TkMu8_DZ_v*",
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
                  "HLT_Mu27_TkMu8_v*",
                  "HLT_Mu30_TkMu11_v*",
                  "HLT_Mu40_TkMu11_v*",
                  "HLT_Mu20_v*",
                  "HLT_TkMu20_v*",
                  "HLT_Mu24_eta2p1_v*",
                  "HLT_TkMu24_eta2p1_v*",
                  "HLT_Mu27_v*",
                  "HLT_TkMu27_v*",
                  "HLT_Mu50_v*",
                  "HLT_Mu55_v*",
                  "HLT_Mu45_eta2p1_v*",
                  "HLT_Mu50_eta2p1_v*",
                  "HLT_Mu8_TrkIsoVVL_v*",
                  "HLT_Mu17_TrkIsoVVL_v*",
                  "HLT_Mu24_TrkIsoVVL_v*",
                  "HLT_Mu34_TrkIsoVVL_v*",
                  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                  "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*",
                  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*",
                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                  "HLT_Mu8_v*",
                  "HLT_Mu17_v*",
                  "HLT_Mu24_v*",
                  "HLT_Mu34_v*",
                  "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele18_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v*",
                  "HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v*"
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
                  objects                  = cms.InputTag("selectedPatTrigger")
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
