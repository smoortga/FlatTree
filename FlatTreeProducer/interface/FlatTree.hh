#ifndef FLATTREE_H
#define FLATTREE_H

#include <TTree.h>
#include <TLorentzVector.h>
#include <string>
#include <iostream>
#include <vector>

#define DEFVAL -666

#include "FlatTree/FlatTreeProducer/interface/Helper.hh"

class FlatTree
{
 public:

   FlatTree(TTree* _tree);

   TTree* tree;

   std::map<std::string, bool> conf;

   void Init();
   void CreateBranches(int buffersize);
   bool doWrite(const std::string& name);

   int ev_run;
   int ev_id;
   int ev_lumi;
   float ev_rho;

   float met_px;
   float met_py;
   float met_pt;
   float met_phi;
   float met_sumet;

   float metGen_px;
   float metGen_py;
   float metGen_pt;
   float metGen_phi;
   float metGen_sumet;

   float metGen_NeutralEMEt;
   float metGen_ChargedEMEt;
   float metGen_NeutralHadEt;
   float metGen_ChargedHadEt;
   float metGen_MuonEt;
   float metGen_InvisibleEt;

   float met_uncorrectedPt;
   float met_uncorrectedPhi;
   float met_uncorrectedSumEt;

   float met_caloMETPt;
   float met_caloMETPhi;
   float met_caloMETSumEt;

   float met_shiftedPx_JetEnUp;
   float met_shiftedPx_JetEnDown;
   float met_shiftedPx_JetResUp;
   float met_shiftedPx_JetResDown;
   float met_shiftedPx_MuonEnUp;
   float met_shiftedPx_MuonEnDown;
   float met_shiftedPx_ElectronEnUp;
   float met_shiftedPx_ElectronEnDown;
   float met_shiftedPx_TauEnUp;
   float met_shiftedPx_TauEnDown;
   float met_shiftedPx_UnclusteredEnUp;
   float met_shiftedPx_UnclusteredEnDown;
   float met_shiftedPx_NoShift;
   float met_shiftedPx_PhotonEnUp;
   float met_shiftedPx_PhotonEnDown;

   float met_shiftedPy_JetEnUp;
   float met_shiftedPy_JetEnDown;
   float met_shiftedPy_JetResUp;
   float met_shiftedPy_JetResDown;
   float met_shiftedPy_MuonEnUp;
   float met_shiftedPy_MuonEnDown;
   float met_shiftedPy_ElectronEnUp;
   float met_shiftedPy_ElectronEnDown;
   float met_shiftedPy_TauEnUp;
   float met_shiftedPy_TauEnDown;
   float met_shiftedPy_UnclusteredEnUp;
   float met_shiftedPy_UnclusteredEnDown;
   float met_shiftedPy_NoShift;
   float met_shiftedPy_PhotonEnUp;
   float met_shiftedPy_PhotonEnDown;

   float met_shiftedPhi_JetEnUp;
   float met_shiftedPhi_JetEnDown;
   float met_shiftedPhi_JetResUp;
   float met_shiftedPhi_JetResDown;
   float met_shiftedPhi_MuonEnUp;
   float met_shiftedPhi_MuonEnDown;
   float met_shiftedPhi_ElectronEnUp;
   float met_shiftedPhi_ElectronEnDown;
   float met_shiftedPhi_TauEnUp;
   float met_shiftedPhi_TauEnDown;
   float met_shiftedPhi_UnclusteredEnUp;
   float met_shiftedPhi_UnclusteredEnDown;
   float met_shiftedPhi_NoShift;
   float met_shiftedPhi_PhotonEnUp;
   float met_shiftedPhi_PhotonEnDown;

   float met_shiftedSumEt_JetEnUp;
   float met_shiftedSumEt_JetEnDown;
   float met_shiftedSumEt_JetResUp;
   float met_shiftedSumEt_JetResDown;
   float met_shiftedSumEt_MuonEnUp;
   float met_shiftedSumEt_MuonEnDown;
   float met_shiftedSumEt_ElectronEnUp;
   float met_shiftedSumEt_ElectronEnDown;
   float met_shiftedSumEt_TauEnUp;
   float met_shiftedSumEt_TauEnDown;
   float met_shiftedSumEt_UnclusteredEnUp;
   float met_shiftedSumEt_UnclusteredEnDown;
   float met_shiftedSumEt_NoShift;
   float met_shiftedSumEt_PhotonEnUp;
   float met_shiftedSumEt_PhotonEnDown;

   float metNoHF_pt;
   float metNoHF_phi;
   float metNoHF_sumet;

   int   pv_n;

   float pv_x;
   float pv_y;
   float pv_z;

   float pv_xError;
   float pv_yError;
   float pv_zError;

   float pv_chi2;
   float pv_ndof;
   float pv_rho;
   int   pv_isFake;

   double met_sig;

   double met_cov00;
   double met_cov10;
   double met_cov01;
   double met_cov11;

   float mc_weight;
   int   mc_id;
   int   mc_f1;
   int   mc_f2;
   float mc_x1;
   float mc_x2;
   float mc_scale;
   float mc_ptHat;

   float weight_originalXWGTUP;
   float weight_scale_muF0p5;
   float weight_scale_muF2;
   float weight_scale_muR0p5;
   float weight_scale_muR2;

   std::vector<float> mc_pdfweights;
   std::vector<std::string> mc_pdfweightIds;

   int mc_pu_intime_NumInt;
   int mc_pu_trueNumInt;
   int mc_pu_before_npu;
   int mc_pu_after_npu;

   int mc_pu_Npvi;
   std::vector<int> mc_pu_Nzpositions;
   std::vector<int> mc_pu_BunchCrossing;
   std::vector<std::vector<float> > mc_pu_zpositions;
   std::vector<std::vector<float> > mc_pu_sumpT_lowpT;
   std::vector<std::vector<float> > mc_pu_sumpT_highpT;
   std::vector<std::vector<int> > mc_pu_ntrks_lowpT;
   std::vector<std::vector<int> > mc_pu_ntrks_highpT;

   // Trigger
   //
   int                       trigger_n;
   std::vector<int>          trigger;
   std::vector<std::string>  trigger_name;
   std::vector<bool>         trigger_pass;
   std::vector<int>          trigger_prescale;
   std::vector<int>          trigger_HLTprescale;
   std::vector<int>          trigger_L1prescale;

   // trigger object general informations
   int                       triggerobject_n;
   std::vector<float>        triggerobject_pt;
   std::vector<float>        triggerobject_eta;
   std::vector<float>        triggerobject_phi;

   std::vector<std::string>  triggerobject_collection;

   // filter Ids...
   std::vector<int>     triggerobject_filterIds_n;
   std::vector<int>     triggerobject_filterIds;

   std::vector<bool>    triggerobject_isTriggerL1Mu;            //-81
   std::vector<bool>    triggerobject_isTriggerL1NoIsoEG;       //-82
   std::vector<bool>    triggerobject_isTriggerL1IsoEG;
   std::vector<bool>    triggerobject_isTriggerL1CenJet;
   std::vector<bool>    triggerobject_isTriggerL1ForJet;
   std::vector<bool>    triggerobject_isTriggerL1TauJet;
   std::vector<bool>    triggerobject_isTriggerL1ETM;
   std::vector<bool>    triggerobject_isTriggerL1ETT;
   std::vector<bool>    triggerobject_isTriggerL1HTT;
   std::vector<bool>    triggerobject_isTriggerL1HTM;
   std::vector<bool>    triggerobject_isTriggerL1JetCounts;
   std::vector<bool>    triggerobject_isTriggerL1HfBitCounts;
   std::vector<bool>    triggerobject_isTriggerL1HfRingEtSums;
   std::vector<bool>    triggerobject_isTriggerL1TechTrig;
   std::vector<bool>    triggerobject_isTriggerL1Castor;
   std::vector<bool>    triggerobject_isTriggerL1BPTX;
   std::vector<bool>    triggerobject_isTriggerL1GtExternal;

   std::vector<bool>    triggerobject_isHLT_TriggerPhoton;      //+81
   std::vector<bool>    triggerobject_isHLT_TriggerElectron;    //+82
   std::vector<bool>    triggerobject_isHLT_TriggerMuon;
   std::vector<bool>    triggerobject_isHLT_TriggerTau;
   std::vector<bool>    triggerobject_isHLT_TriggerJet;
   std::vector<bool>    triggerobject_isHLT_TriggerBJet;
   std::vector<bool>    triggerobject_isHLT_TriggerMET;
   std::vector<bool>    triggerobject_isHLT_TriggerTET;
   std::vector<bool>    triggerobject_isHLT_TriggerTHT;
   std::vector<bool>    triggerobject_isHLT_TriggerMHT;
   std::vector<bool>    triggerobject_isHLT_TriggerTrack;
   std::vector<bool>    triggerobject_isHLT_TriggerCluster;
   std::vector<bool>    triggerobject_isHLT_TriggerMETSig;
   std::vector<bool>    triggerobject_isHLT_TriggerELongit;
   std::vector<bool>    triggerobject_isHLT_TriggerMHTSig;
   std::vector<bool>    triggerobject_isHLT_TriggerHLongit;

   // filters label...
   std::vector<int>           triggerobject_filterLabels_n;
   std::vector<std::string>   triggerobject_filterLabels;

   // paths names and status
   std::vector<int>           triggerobject_pathNamesAll_n;
   std::vector<std::string>   triggerobject_pathNamesAll;
   std::vector<bool>          triggerobject_pathNamesAll_isBoth;
   std::vector<bool>          triggerobject_pathNamesAll_isL3;
   std::vector<bool>          triggerobject_pathNamesAll_isLF;
   std::vector<bool>          triggerobject_pathNamesAll_isNone;

   int nvertex;

   // Electrons
   //
   int el_n;
   std::vector<float> el_pt;
   std::vector<float> el_eta;
   std::vector<float> el_phi;
   std::vector<float> el_m;
   std::vector<float> el_E;
   std::vector<int> el_id;
   std::vector<int> el_charge;

   std::vector<int> el_passConversionVeto;
   std::vector<int> el_isGsfCtfScPixChargeConsistent;
   std::vector<int> el_isGsfScPixChargeConsistent;

   std::vector<float> el_ecalEnergy;
   std::vector<float> el_correctedEcalEnergy;
   std::vector<float> el_correctedEcalEnergyError;
   std::vector<float> el_trackMomentumError;

   std::vector<float> el_ip3d;
   std::vector<float> el_ip3dErr;
   std::vector<float> el_ip2d;
   std::vector<float> el_ip2dErr;
   std::vector<float> el_ip3dBS;
   std::vector<float> el_ip3dBSErr;
   std::vector<float> el_ip2dBS;
   std::vector<float> el_ip2dBSErr;

   std::vector<float> el_neutralHadronIso;
   std::vector<float> el_chargedHadronIso;
   std::vector<float> el_puChargedHadronIso;
   std::vector<float> el_ecalIso;
   std::vector<float> el_hcalIso;
   std::vector<float> el_particleIso;
   std::vector<float> el_photonIso;
   std::vector<float> el_trackIso;

   std::vector<float> el_ecalPFClusterIso;
   std::vector<float> el_hcalPFClusterIso;

   std::vector<float> el_pfIso_sumChargedHadronPt;
   std::vector<float> el_pfIso_sumNeutralHadronEt;
   std::vector<float> el_pfIso_sumPhotonEt;
   std::vector<float> el_pfIso_sumPUPt;

   std::vector<float> el_dr03EcalRecHitSumEt;
   std::vector<float> el_dr03HcalTowerSumEt;
   std::vector<float> el_dr03HcalDepth1TowerSumEt;
   std::vector<float> el_dr03HcalDepth2TowerSumEt;
   std::vector<float> el_dr03TkSumPt;

   std::vector<float> el_dr04EcalRecHitSumEt;
   std::vector<float> el_dr04HcalTowerSumEt;
   std::vector<float> el_dr04HcalDepth1TowerSumEt;
   std::vector<float> el_dr04HcalDepth2TowerSumEt;
   std::vector<float> el_dr04TkSumPt;

   std::vector<float> el_hcalOverEcal;
   std::vector<float> el_hcalOverEcalBc;
   std::vector<float> el_hcalDepth1OverEcal;
   std::vector<float> el_hcalDepth2OverEcal;
   std::vector<float> el_eSeedClusterOverPout;
   std::vector<float> el_eSeedClusterOverP;
   std::vector<float> el_eEleClusterOverPout;
   std::vector<float> el_deltaEtaEleClusterTrackAtCalo;
   std::vector<float> el_deltaPhiEleClusterTrackAtCalo;

   std::vector<float> el_vx;
   std::vector<float> el_vy;
   std::vector<float> el_vz;

   std::vector<bool> el_hasGsfTrack;
   std::vector<float> el_gsfTrack_d0;
   std::vector<float> el_gsfTrack_z0;
   std::vector<float> el_gsfTrack_d0Error;
   std::vector<float> el_gsfTrack_z0Error;
   std::vector<float> el_gsfTrack_PV_dxy;
   std::vector<float> el_gsfTrack_PV_dz;
   std::vector<float> el_gsfTrack_RP_dxy;
   std::vector<float> el_gsfTrack_RP_dz;
   std::vector<float> el_gsfTrack_BS_dxy;
   std::vector<float> el_gsfTrack_BS_dz;
   std::vector<float> el_gsfTrack_dxyError;
   std::vector<float> el_gsfTrack_dzError;
   std::vector<float> el_gsfTrack_normalizedChi2;

   std::vector<int> el_numberOfHits;
   std::vector<int> el_numberOfValidHits;

   std::vector<int> el_expectedMissingOuterHits;
   std::vector<int> el_numberOfValidPixelHits;
   std::vector<int> el_numberOfLostPixelHits;
   std::vector<int> el_trackerLayersWithMeasurement;
   std::vector<int> el_pixelLayersWithMeasurement;
   std::vector<int> el_numberOfValidStripLayersWithMonoAndStereo;
   std::vector<int> el_trackerLayersWithoutMeasurement;

   std::vector<float> el_superCluster_eta;
   std::vector<float> el_superCluster_phi;
   std::vector<float> el_superCluster_energy;
   std::vector<float> el_superCluster_rawEnergy;
   std::vector<float> el_superCluster_preshowerEnergy;
   std::vector<float> el_superCluster_etaWidth;
   std::vector<float> el_superCluster_phiWidth;
   std::vector<float> el_superCluster_preshowerEnergyPlane1;
   std::vector<float> el_superCluster_preshowerEnergyPlane2;
   std::vector<float> el_superCluster_positionR;

   std::vector<int> el_basicClustersSize;
   std::vector<float> el_e1x5;
   std::vector<float> el_e5x5;
   std::vector<float> el_e2x5Max;
   std::vector<float> el_sigmaEtaEta;
   std::vector<float> el_sigmaIetaIeta;
   std::vector<float> el_sigmaIphiIphi;
   std::vector<float> el_sigmaIetaIphi;
   std::vector<float> el_full5x5_sigmaIphiIphi;
   std::vector<float> el_full5x5_sigmaEtaEta;
   std::vector<float> el_full5x5_sigmaIetaIeta;
   std::vector<float> el_full5x5_sigmaIetaIphi;
   std::vector<float> el_full5x5_r9;
   std::vector<float> el_full5x5_e1x5;
   std::vector<float> el_full5x5_e5x5;
   std::vector<float> el_full5x5_e2x5Max;

   std::vector<float> el_hadronicOverEm;
   std::vector<int> el_numberOfLostHits;
   std::vector<int> el_numberOfLostHitsDefault;

   std::vector<float> el_fbrem;
   std::vector<float> el_kf_normalizedChi2;
   std::vector<float> el_gsf_normalizedChi2;
   std::vector<float> el_deltaEtaSuperClusterTrackAtVtx;
   std::vector<float> el_deltaPhiSuperClusterTrackAtVtx;
   std::vector<float> el_deltaEtaSeedClusterTrackAtCalo;
   std::vector<float> el_deltaPhiSeedClusterTrackAtCalo;
   std::vector<float> el_superClusterEtaWidth;
   std::vector<float> el_superClusterPhiWidth;
   std::vector<float> el_full5x5_OneMinusE1x5E5x5;
   std::vector<float> el_OneMinusE1x5E5x5;
   std::vector<float> el_eSuperClusterOverP;
   std::vector<float> el_IoEmIoP;
   std::vector<float> el_ooEmooP;
   std::vector<float> el_eleEoPout;
   std::vector<float> el_PreShowerOverRaw;

   std::vector<float> el_mvaNonTrigV0;
   std::vector<float> el_mvaNonTrigCat;

   std::vector<bool> el_vetoCBId;
   std::vector<bool> el_looseCBId;
   std::vector<bool> el_mediumCBId;
   std::vector<bool> el_tightCBId;
   std::vector<bool> el_heepCBId;

   std::vector<bool> el_mediumMVAId;
   std::vector<bool> el_tightMVAId;

   std::vector<int> el_hasMCMatch;
   std::vector<float> el_gen_pt;
   std::vector<float> el_gen_eta;
   std::vector<float> el_gen_phi;
   std::vector<float> el_gen_m;
   std::vector<float> el_gen_E;
   std::vector<int> el_gen_status;
   std::vector<int> el_gen_id;
   std::vector<int> el_gen_charge;
   std::vector<float> el_gen_dr;

   std::vector<int> el_hasMCMatchPAT;
   std::vector<float> el_genPAT_pt;
   std::vector<float> el_genPAT_eta;
   std::vector<float> el_genPAT_phi;
   std::vector<float> el_genPAT_m;
   std::vector<float> el_genPAT_E;
   std::vector<int> el_genPAT_status;
   std::vector<int> el_genPAT_id;
   std::vector<int> el_genPAT_charge;

   std::vector<bool> el_hasMatchedConversion;
   std::vector<int> el_expectedMissingInnerHits;

   // Muons
   //
   int mu_n;
   std::vector<float> mu_pt;
   std::vector<float> mu_eta;
   std::vector<float> mu_phi;
   std::vector<float> mu_m;
   std::vector<float> mu_E;
   std::vector<int> mu_id;
   std::vector<int> mu_charge;

   std::vector<float> mu_ip3d;
   std::vector<float> mu_ip3dErr;
   std::vector<float> mu_ip2d;
   std::vector<float> mu_ip2dErr;
   std::vector<float> mu_ip3dBS;
   std::vector<float> mu_ip3dBSErr;
   std::vector<float> mu_ip2dBS;
   std::vector<float> mu_ip2dBSErr;

   std::vector<float> mu_neutralHadronIso;
   std::vector<float> mu_chargedHadronIso;
   std::vector<float> mu_puChargedHadronIso;
   std::vector<float> mu_ecalIso;
   std::vector<float> mu_hcalIso;
   std::vector<float> mu_photonIso;
   std::vector<float> mu_trackIso;

   std::vector<float> mu_pfIso03_sumChargedHadronPt;
   std::vector<float> mu_pfIso03_sumChargedParticlePt;
   std::vector<float> mu_pfIso03_sumNeutralHadronEt;
   std::vector<float> mu_pfIso03_sumNeutralHadronEtHighThreshold;
   std::vector<float> mu_pfIso03_sumPhotonEt;
   std::vector<float> mu_pfIso03_sumPhotonEtHighThreshold;
   std::vector<float> mu_pfIso03_sumPUPt;

   std::vector<float> mu_pfIso04_sumChargedHadronPt;
   std::vector<float> mu_pfIso04_sumChargedParticlePt;
   std::vector<float> mu_pfIso04_sumNeutralHadronEt;
   std::vector<float> mu_pfIso04_sumNeutralHadronEtHighThreshold;
   std::vector<float> mu_pfIso04_sumPhotonEt;
   std::vector<float> mu_pfIso04_sumPhotonEtHighThreshold;
   std::vector<float> mu_pfIso04_sumPUPt;

   std::vector<float> mu_pfMeanIso03_sumChargedHadronPt;
   std::vector<float> mu_pfMeanIso03_sumChargedParticlePt;
   std::vector<float> mu_pfMeanIso03_sumNeutralHadronEt;
   std::vector<float> mu_pfMeanIso03_sumNeutralHadronEtHighThreshold;
   std::vector<float> mu_pfMeanIso03_sumPhotonEt;
   std::vector<float> mu_pfMeanIso03_sumPhotonEtHighThreshold;
   std::vector<float> mu_pfMeanIso03_sumPUPt;

   std::vector<float> mu_pfSumIso03_sumChargedHadronPt;
   std::vector<float> mu_pfSumIso03_sumChargedParticlePt;
   std::vector<float> mu_pfSumIso03_sumNeutralHadronEt;
   std::vector<float> mu_pfSumIso03_sumNeutralHadronEtHighThreshold;
   std::vector<float> mu_pfSumIso03_sumPhotonEt;
   std::vector<float> mu_pfSumIso03_sumPhotonEtHighThreshold;
   std::vector<float> mu_pfSumIso03_sumPUPt;

   std::vector<float> mu_pfMeanIso04_sumChargedHadronPt;
   std::vector<float> mu_pfMeanIso04_sumChargedParticlePt;
   std::vector<float> mu_pfMeanIso04_sumNeutralHadronEt;
   std::vector<float> mu_pfMeanIso04_sumNeutralHadronEtHighThreshold;
   std::vector<float> mu_pfMeanIso04_sumPhotonEt;
   std::vector<float> mu_pfMeanIso04_sumPhotonEtHighThreshold;
   std::vector<float> mu_pfMeanIso04_sumPUPt;

   std::vector<float> mu_pfSumIso04_sumChargedHadronPt;
   std::vector<float> mu_pfSumIso04_sumChargedParticlePt;
   std::vector<float> mu_pfSumIso04_sumNeutralHadronEt;
   std::vector<float> mu_pfSumIso04_sumNeutralHadronEtHighThreshold;
   std::vector<float> mu_pfSumIso04_sumPhotonEt;
   std::vector<float> mu_pfSumIso04_sumPhotonEtHighThreshold;
   std::vector<float> mu_pfSumIso04_sumPUPt;

   std::vector<int> mu_isGlobalMuon;
   std::vector<int> mu_isTrackerMuon;
   std::vector<int> mu_isStandAloneMuon;
   std::vector<int> mu_isCaloMuon;
   std::vector<int> mu_isPFMuon;
   std::vector<int> mu_isRPCMuon;

   std::vector<float> mu_vx;
   std::vector<float> mu_vy;
   std::vector<float> mu_vz;

   std::vector<int> mu_numberOfMatches;
   std::vector<int> mu_numberOfMatchedStations;

   std::vector<float> mu_combinedQuality_chi2LocalPosition;
   std::vector<float> mu_combinedQuality_trkKink;

   std::vector<float> mu_segmentCompatibility;
   std::vector<float> mu_caloCompatibility;

   std::vector<bool> mu_isLooseMuon;
   std::vector<bool> mu_isMediumMuon;
   std::vector<bool> mu_isTightMuon;
   std::vector<bool> mu_isSoftMuon;
   std::vector<bool> mu_isHighPtMuon;

   std::vector<bool> mu_isGoodMuon_AllGlobalMuons;
   std::vector<bool> mu_isGoodMuon_AllStandAloneMuons;
   std::vector<bool> mu_isGoodMuon_AllTrackerMuons;
   std::vector<bool> mu_isGoodMuon_TrackerMuonArbitrated;
   std::vector<bool> mu_isGoodMuon_AllArbitrated;
   std::vector<bool> mu_isGoodMuon_GlobalMuonPromptTight;
   std::vector<bool> mu_isGoodMuon_TMLastStationLoose;
   std::vector<bool> mu_isGoodMuon_TMLastStationTight;
   std::vector<bool> mu_isGoodMuon_TM2DCompatibilityLoose;
   std::vector<bool> mu_isGoodMuon_TM2DCompatibilityTight;
   std::vector<bool> mu_isGoodMuon_TMOneStationLoose;
   std::vector<bool> mu_isGoodMuon_TMOneStationTight;
   std::vector<bool> mu_isGoodMuon_TMLastStationOptimizedLowPtLoose;
   std::vector<bool> mu_isGoodMuon_TMLastStationOptimizedLowPtTight;
   std::vector<bool> mu_isGoodMuon_GMTkChiCompatibility;
   std::vector<bool> mu_isGoodMuon_GMStaChiCompatibility;
   std::vector<bool> mu_isGoodMuon_GMTkKinkTight;
   std::vector<bool> mu_isGoodMuon_TMLastStationAngLoose;
   std::vector<bool> mu_isGoodMuon_TMLastStationAngTight;
   std::vector<bool> mu_isGoodMuon_TMOneStationAngLoose;
   std::vector<bool> mu_isGoodMuon_TMOneStationAngTight;
   std::vector<bool> mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose;
   std::vector<bool> mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;

   std::vector<float> mu_calEnergy_em;
   std::vector<float> mu_calEnergy_had;
   std::vector<float> mu_calEnergy_ho;
   std::vector<float> mu_calEnergy_emS9;
   std::vector<float> mu_calEnergy_hadS9;
   std::vector<float> mu_calEnergy_hoS9;
   std::vector<float> mu_calEnergy_emS25;
   std::vector<float> mu_calEnergy_emMax;
   std::vector<float> mu_calEnergy_hadMax;
   std::vector<float> mu_calEnergy_ecal_time;
   std::vector<float> mu_calEnergy_hcal_time;
   std::vector<int> mu_calEnergy_ecal_rawId;
   std::vector<int> mu_calEnergy_hcal_rawId;

   std::vector<float> mu_isolationR03_trackerVetoPt;
   std::vector<float> mu_isolationR03_emVetoEt;
   std::vector<float> mu_isolationR03_hadVetoEt;
   std::vector<float> mu_isolationR03_hoVetoEt;
   std::vector<float> mu_isolationR03_sumPt;
   std::vector<float> mu_isolationR03_emEt;
   std::vector<float> mu_isolationR03_hadEt;
   std::vector<float> mu_isolationR03_hoEt;
   std::vector<int> mu_isolationR03_nTracks;
   std::vector<int> mu_isolationR03_nJets;

   std::vector<float> mu_isolationR05_trackerVetoPt;
   std::vector<float> mu_isolationR05_emVetoEt;
   std::vector<float> mu_isolationR05_hadVetoEt;
   std::vector<float> mu_isolationR05_hoVetoEt;
   std::vector<float> mu_isolationR05_sumPt;
   std::vector<float> mu_isolationR05_emEt;
   std::vector<float> mu_isolationR05_hadEt;
   std::vector<float> mu_isolationR05_hoEt;
   std::vector<int> mu_isolationR05_nTracks;
   std::vector<int> mu_isolationR05_nJets;

   std::vector<int> mu_hasGlobalTrack;
   std::vector<float> mu_globalTrack_d0;
   std::vector<float> mu_globalTrack_z0;
   std::vector<float> mu_globalTrack_d0Error;
   std::vector<float> mu_globalTrack_z0Error;
   std::vector<float> mu_globalTrack_PV_dxy;
   std::vector<float> mu_globalTrack_PV_dz;
   std::vector<float> mu_globalTrack_RP_dxy;
   std::vector<float> mu_globalTrack_RP_dz;
   std::vector<float> mu_globalTrack_BS_dxy;
   std::vector<float> mu_globalTrack_BS_dz;
   std::vector<float> mu_globalTrack_dxyError;
   std::vector<float> mu_globalTrack_dzError;
   std::vector<float> mu_globalTrack_normalizedChi2;
   std::vector<int> mu_globalTrack_numberOfValidHits;
   std::vector<int> mu_globalTrack_numberOfValidMuonHits;
   std::vector<int> mu_globalTrack_numberOfLostHits;
   std::vector<float> mu_globalTrack_pt;
   std::vector<float> mu_globalTrack_eta;
   std::vector<float> mu_globalTrack_phi;
   std::vector<float> mu_globalTrack_ptError;
   std::vector<float> mu_globalTrack_etaError;
   std::vector<float> mu_globalTrack_phiError;
   std::vector<float> mu_globalTrack_vx;
   std::vector<float> mu_globalTrack_vy;
   std::vector<float> mu_globalTrack_vz;
   std::vector<float> mu_globalTrack_qoverp;
   std::vector<float> mu_globalTrack_qoverpError;
   std::vector<int> mu_globalTrack_charge;
   std::vector<int> mu_globalTrack_trackerLayersWithMeasurement;
   std::vector<int> mu_globalTrack_pixelLayersWithMeasurement;
   std::vector<int> mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo;
   std::vector<int> mu_globalTrack_trackerLayersWithoutMeasurement;
   std::vector<int> mu_globalTrack_numberOfValidPixelHits;
   std::vector<int> mu_globalTrack_numberOfLostPixelHits;
   std::vector<int> mu_globalTrack_numberOfInnerHits;
   std::vector<int> mu_globalTrack_numberOfOuterHits;
   std::vector<float> mu_globalTrack_validFraction;

   std::vector<int> mu_bestTrackType;
   std::vector<int> mu_hasBestTrack;
   std::vector<float> mu_bestTrack_d0;
   std::vector<float> mu_bestTrack_z0;
   std::vector<float> mu_bestTrack_d0Error;
   std::vector<float> mu_bestTrack_z0Error;
   std::vector<float> mu_bestTrack_PV_dxy;
   std::vector<float> mu_bestTrack_PV_dz;
   std::vector<float> mu_bestTrack_RP_dxy;
   std::vector<float> mu_bestTrack_RP_dz;
   std::vector<float> mu_bestTrack_BS_dxy;
   std::vector<float> mu_bestTrack_BS_dz;
   std::vector<float> mu_bestTrack_dxyError;
   std::vector<float> mu_bestTrack_dzError;
   std::vector<float> mu_bestTrack_normalizedChi2;
   std::vector<int> mu_bestTrack_numberOfValidHits;
   std::vector<int> mu_bestTrack_numberOfLostHits;
   std::vector<float> mu_bestTrack_pt;
   std::vector<float> mu_bestTrack_eta;
   std::vector<float> mu_bestTrack_phi;
   std::vector<float> mu_bestTrack_ptError;
   std::vector<float> mu_bestTrack_etaError;
   std::vector<float> mu_bestTrack_phiError;
   std::vector<float> mu_bestTrack_vx;
   std::vector<float> mu_bestTrack_vy;
   std::vector<float> mu_bestTrack_vz;
   std::vector<float> mu_bestTrack_qoverp;
   std::vector<float> mu_bestTrack_qoverpError;
   std::vector<int> mu_bestTrack_charge;
   std::vector<int> mu_bestTrack_trackerLayersWithMeasurement;
   std::vector<int> mu_bestTrack_pixelLayersWithMeasurement;
   std::vector<int> mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo;
   std::vector<int> mu_bestTrack_trackerLayersWithoutMeasurement;
   std::vector<int> mu_bestTrack_numberOfValidPixelHits;
   std::vector<int> mu_bestTrack_numberOfLostPixelHits;
   std::vector<int> mu_bestTrack_numberOfInnerHits;
   std::vector<int> mu_bestTrack_numberOfOuterHits;
   std::vector<float> mu_bestTrack_validFraction;

   std::vector<int> mu_hasInnerTrack;
   std::vector<float> mu_innerTrack_d0;
   std::vector<float> mu_innerTrack_z0;
   std::vector<float> mu_innerTrack_d0Error;
   std::vector<float> mu_innerTrack_z0Error;
   std::vector<float> mu_innerTrack_PV_dxy;
   std::vector<float> mu_innerTrack_PV_dz;
   std::vector<float> mu_innerTrack_RP_dxy;
   std::vector<float> mu_innerTrack_RP_dz;
   std::vector<float> mu_innerTrack_BS_dxy;
   std::vector<float> mu_innerTrack_BS_dz;
   std::vector<float> mu_innerTrack_dxyError;
   std::vector<float> mu_innerTrack_dzError;
   std::vector<float> mu_innerTrack_normalizedChi2;
   std::vector<int> mu_innerTrack_numberOfValidHits;
   std::vector<int> mu_innerTrack_numberOfLostHits;
   std::vector<float> mu_innerTrack_pt;
   std::vector<float> mu_innerTrack_eta;
   std::vector<float> mu_innerTrack_phi;
   std::vector<float> mu_innerTrack_ptError;
   std::vector<float> mu_innerTrack_etaError;
   std::vector<float> mu_innerTrack_phiError;
   std::vector<float> mu_innerTrack_vx;
   std::vector<float> mu_innerTrack_vy;
   std::vector<float> mu_innerTrack_vz;
   std::vector<float> mu_innerTrack_qoverp;
   std::vector<float> mu_innerTrack_qoverpError;
   std::vector<int> mu_innerTrack_charge;
   std::vector<int> mu_innerTrack_trackerLayersWithMeasurement;
   std::vector<int> mu_innerTrack_pixelLayersWithMeasurement;
   std::vector<int> mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo;
   std::vector<int> mu_innerTrack_trackerLayersWithoutMeasurement;
   std::vector<int> mu_innerTrack_numberOfValidPixelHits;
   std::vector<int> mu_innerTrack_numberOfLostPixelHits;
   std::vector<int> mu_innerTrack_numberOfInnerHits;
   std::vector<int> mu_innerTrack_numberOfOuterHits;
   std::vector<float> mu_innerTrack_validFraction;

   std::vector<int> mu_type;

   std::vector<int> mu_hasMCMatch;
   std::vector<float> mu_gen_pt;
   std::vector<float> mu_gen_eta;
   std::vector<float> mu_gen_phi;
   std::vector<float> mu_gen_m;
   std::vector<int> mu_gen_status;
   std::vector<int> mu_gen_id;
   std::vector<int> mu_gen_charge;
   std::vector<float> mu_gen_dr;

   std::vector<int> mu_hasMCMatchPAT;
   std::vector<float> mu_genPAT_pt;
   std::vector<float> mu_genPAT_eta;
   std::vector<float> mu_genPAT_phi;
   std::vector<float> mu_genPAT_m;
   std::vector<int> mu_genPAT_status;
   std::vector<int> mu_genPAT_id;
   std::vector<int> mu_genPAT_charge;

   // Taus
   //
   int tau_n;
   std::vector<float> tau_pt;
   std::vector<float> tau_eta;
   std::vector<float> tau_phi;
   std::vector<float> tau_m;
   std::vector<float> tau_E;
   std::vector<int> tau_id;
   std::vector<int> tau_charge;

   std::vector<bool> tau_hasLeadChargedHadrCand;
   std::vector<float> tau_leadingTrackPt;
   std::vector<int> tau_leadingTrackCharge;
   std::vector<float> tau_leadingTrackDz;
   std::vector<float> tau_leadingTrackDxy;

   std::vector<int> tau_decayMode;
   std::vector<float> tau_decayModeFindingOldDMs;
   std::vector<float> tau_decayModeFindingNewDMs;

   std::vector<float> tau_puCorrPtSum;
   std::vector<float> tau_neutralIsoPtSum;
   std::vector<float> tau_chargedIsoPtSum;
   std::vector<float> tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;

   std::vector<float> tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<float> tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<float> tau_byTightCombinedIsolationDeltaBetaCorr3Hits;

   std::vector<float> tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   std::vector<float> tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   std::vector<float> tau_byTightIsolationMVArun2v1DBdR03oldDMwLT;
   std::vector<float> tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;

   std::vector<float> tau_againstMuonLoose3;
   std::vector<float> tau_againstMuonTight3;

   std::vector<float> tau_pfEssential_jet_pt;
   std::vector<float> tau_pfEssential_jet_eta;
   std::vector<float> tau_pfEssential_jet_phi;
   std::vector<float> tau_pfEssential_jet_m;

   std::vector<float> tau_pfEssential_jetCorr_pt;
   std::vector<float> tau_pfEssential_jetCorr_eta;
   std::vector<float> tau_pfEssential_jetCorr_phi;
   std::vector<float> tau_pfEssential_jetCorr_m;

   std::vector<bool> tau_pfEssential_hasSV;
   std::vector<float> tau_pfEssential_sv_x;
   std::vector<float> tau_pfEssential_sv_y;
   std::vector<float> tau_pfEssential_sv_z;

   std::vector<float> tau_pfEssential_flightLengthSig;
   std::vector<float> tau_pfEssential_dxy;
   std::vector<float> tau_pfEssential_dxy_error;
   std::vector<float> tau_pfEssential_dxy_Sig;

   // Jets
   //
   int jet_n;
   std::vector<float> jet_pt;
   std::vector<float> jet_eta;
   std::vector<float> jet_phi;
   std::vector<float> jet_m;
   std::vector<float> jet_E;

   std::vector<int> jet_ntrk;

   std::vector<float> jet_JBP;
   std::vector<float> jet_JP;
   std::vector<float> jet_TCHP;
   std::vector<float> jet_TCHE;
   std::vector<float> jet_SSVHE;
   std::vector<float> jet_SSVHP;
   std::vector<float> jet_CMVA;
   std::vector<float> jet_CSVv2;
   std::vector<float> jet_cMVAv2;
   std::vector<float> jet_CharmCvsL;
   std::vector<float> jet_CharmCvsB;
   std::vector<int> jet_partonFlavour;
   std::vector<int> jet_hadronFlavour;

   std::vector<float> jet_neutralHadronEnergy;
   std::vector<float> jet_neutralEmEnergy;
   std::vector<float> jet_chargedHadronEnergy;
   std::vector<float> jet_chargedEmEnergy;
   std::vector<float> jet_electronEnergy;
   std::vector<float> jet_muonEnergy;
   std::vector<float> jet_photonEnergy;

   std::vector<float> jet_charge;
   std::vector<float> jet_chargeVec;
   std::vector<int> jet_chargedMultiplicity;
   std::vector<int> jet_neutralMultiplicity;
   std::vector<int> jet_chargedHadronMultiplicity;

   std::vector<float> jet_jetArea;

   std::vector<float> jet_jecFactorUncorrected;
   std::vector<float> jet_jecFactorL1FastJet;
   std::vector<float> jet_jecFactorL2Relative;
   std::vector<float> jet_jecFactorL3Absolute;

   std::vector<float> jet_neutralHadronEnergyFraction;
   std::vector<float> jet_neutralEmEnergyFraction;
   std::vector<float> jet_chargedHadronEnergyFraction;
   std::vector<float> jet_muonEnergyFraction;
   std::vector<float> jet_chargedEmEnergyFraction;

   std::vector<float> jet_Unc;

   std::vector<float> jet_pileupJetId;

   std::vector<bool> jet_looseJetID;
   std::vector<bool> jet_tightJetID;

   std::vector<bool> jet_hasGenJet;
   std::vector<float> jet_genJet_pt;
   std::vector<float> jet_genJet_eta;
   std::vector<float> jet_genJet_phi;
   std::vector<float> jet_genJet_m;
   std::vector<float> jet_genJet_E;
   std::vector<int> jet_genJet_status;
   std::vector<int> jet_genJet_id;

   std::vector<bool> jet_hasGenParton;
   std::vector<float> jet_genParton_pt;
   std::vector<float> jet_genParton_eta;
   std::vector<float> jet_genParton_phi;
   std::vector<float> jet_genParton_m;
   std::vector<float> jet_genParton_E;
   std::vector<int> jet_genParton_status;
   std::vector<int> jet_genParton_id;

   // GenJets
   //
   int genJet_n;
   std::vector<float> genJet_pt;
   std::vector<float> genJet_eta;
   std::vector<float> genJet_phi;
   std::vector<float> genJet_m;
   std::vector<float> genJet_E;
   std::vector<float> genJet_emEnergy;
   std::vector<float> genJet_hadEnergy;
   std::vector<float> genJet_invisibleEnergy;
   std::vector<float> genJet_auxiliaryEnergy;
   std::vector<int>   genJet_flavour;

   // gen
   float gen_PVz;
   int gen_n;
   std::vector<float> gen_pt;
   std::vector<float> gen_eta;
   std::vector<float> gen_phi;
   std::vector<float> gen_m;
   std::vector<float> gen_E;
   std::vector<int> gen_id;
   std::vector<int> gen_charge;
   std::vector<int> gen_status;
   std::vector<int> gen_index;
   std::vector<int> gen_mother_index;
   std::vector<int> gen_daughter_n;
   std::vector<std::vector<int> > gen_daughter_index;
   
   int genTTX_id;
};

#endif
