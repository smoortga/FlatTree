#include "FlatTree/FlatTreeProducer/interface/FlatTree.hh"

FlatTree::FlatTree(TTree* _tree)
{
   tree = _tree;
}

void FlatTree::Init()
{
   ev_run = DEFVAL;
   ev_id = DEFVAL;
   ev_lumi = DEFVAL;
   ev_rho = DEFVAL;

   met_px = DEFVAL;
   met_py = DEFVAL;
   met_pt = DEFVAL;
   met_phi = DEFVAL;
   met_sumet = DEFVAL;
   met_sig = DEFVAL;

   met_cov00 = DEFVAL;
   met_cov10 = DEFVAL;
   met_cov01 = DEFVAL;
   met_cov11 = DEFVAL;

   metGen_px = DEFVAL;
   metGen_py = DEFVAL;
   metGen_pt = DEFVAL;
   metGen_phi = DEFVAL;
   metGen_sumet = DEFVAL;

   metGen_NeutralEMEt = DEFVAL;
   metGen_ChargedEMEt = DEFVAL;
   metGen_NeutralHadEt = DEFVAL;
   metGen_ChargedHadEt = DEFVAL;
   metGen_MuonEt = DEFVAL;
   metGen_InvisibleEt = DEFVAL;

   met_uncorrectedPt = DEFVAL;
   met_uncorrectedPhi = DEFVAL;
   met_uncorrectedSumEt = DEFVAL;

   met_caloMETPt = DEFVAL;
   met_caloMETPhi = DEFVAL;
   met_caloMETSumEt = DEFVAL;

   met_shiftedPx_JetEnUp = DEFVAL;
   met_shiftedPx_JetEnDown = DEFVAL;
   met_shiftedPx_JetResUp = DEFVAL;
   met_shiftedPx_JetResDown = DEFVAL;
   met_shiftedPx_MuonEnUp = DEFVAL;
   met_shiftedPx_MuonEnDown = DEFVAL;
   met_shiftedPx_ElectronEnUp = DEFVAL;
   met_shiftedPx_ElectronEnDown = DEFVAL;
   met_shiftedPx_TauEnUp = DEFVAL;
   met_shiftedPx_TauEnDown = DEFVAL;
   met_shiftedPx_UnclusteredEnUp = DEFVAL;
   met_shiftedPx_UnclusteredEnDown = DEFVAL;
   met_shiftedPx_NoShift = DEFVAL;
   met_shiftedPx_PhotonEnUp = DEFVAL;
   met_shiftedPx_PhotonEnDown = DEFVAL;

   met_shiftedPy_JetEnUp = DEFVAL;
   met_shiftedPy_JetEnDown = DEFVAL;
   met_shiftedPy_JetResUp = DEFVAL;
   met_shiftedPy_JetResDown = DEFVAL;
   met_shiftedPy_MuonEnUp = DEFVAL;
   met_shiftedPy_MuonEnDown = DEFVAL;
   met_shiftedPy_ElectronEnUp = DEFVAL;
   met_shiftedPy_ElectronEnDown = DEFVAL;
   met_shiftedPy_TauEnUp = DEFVAL;
   met_shiftedPy_TauEnDown = DEFVAL;
   met_shiftedPy_UnclusteredEnUp = DEFVAL;
   met_shiftedPy_UnclusteredEnDown = DEFVAL;
   met_shiftedPy_NoShift = DEFVAL;
   met_shiftedPy_PhotonEnUp = DEFVAL;
   met_shiftedPy_PhotonEnDown = DEFVAL;

   met_shiftedPhi_JetEnUp = DEFVAL;
   met_shiftedPhi_JetEnDown = DEFVAL;
   met_shiftedPhi_JetResUp = DEFVAL;
   met_shiftedPhi_JetResDown = DEFVAL;
   met_shiftedPhi_MuonEnUp = DEFVAL;
   met_shiftedPhi_MuonEnDown = DEFVAL;
   met_shiftedPhi_ElectronEnUp = DEFVAL;
   met_shiftedPhi_ElectronEnDown = DEFVAL;
   met_shiftedPhi_TauEnUp = DEFVAL;
   met_shiftedPhi_TauEnDown = DEFVAL;
   met_shiftedPhi_UnclusteredEnUp = DEFVAL;
   met_shiftedPhi_UnclusteredEnDown = DEFVAL;
   met_shiftedPhi_NoShift = DEFVAL;
   met_shiftedPhi_PhotonEnUp = DEFVAL;
   met_shiftedPhi_PhotonEnDown = DEFVAL;

   met_shiftedSumEt_JetEnUp = DEFVAL;
   met_shiftedSumEt_JetEnDown = DEFVAL;
   met_shiftedSumEt_JetResUp = DEFVAL;
   met_shiftedSumEt_JetResDown = DEFVAL;
   met_shiftedSumEt_MuonEnUp = DEFVAL;
   met_shiftedSumEt_MuonEnDown = DEFVAL;
   met_shiftedSumEt_ElectronEnUp = DEFVAL;
   met_shiftedSumEt_ElectronEnDown = DEFVAL;
   met_shiftedSumEt_TauEnUp = DEFVAL;
   met_shiftedSumEt_TauEnDown = DEFVAL;
   met_shiftedSumEt_UnclusteredEnUp = DEFVAL;
   met_shiftedSumEt_UnclusteredEnDown = DEFVAL;
   met_shiftedSumEt_NoShift = DEFVAL;
   met_shiftedSumEt_PhotonEnUp = DEFVAL;
   met_shiftedSumEt_PhotonEnDown = DEFVAL;

   metNoHF_pt     = DEFVAL;
   metNoHF_phi    = DEFVAL;
   metNoHF_sumet  = DEFVAL;

   nvertex = DEFVAL;

   pv_n = DEFVAL;
   pv_x = DEFVAL;
   pv_y = DEFVAL;
   pv_z = DEFVAL;

   pv_xError = DEFVAL;
   pv_yError = DEFVAL;
   pv_zError = DEFVAL;

   pv_chi2 = DEFVAL;
   pv_ndof = DEFVAL;
   pv_rho = DEFVAL;
   pv_isFake = DEFVAL;

   mc_weight = DEFVAL;
   mc_id = DEFVAL;
   mc_f1 = DEFVAL;
   mc_f2 = DEFVAL;
   mc_x1 = DEFVAL;
   mc_x2 = DEFVAL;
   mc_scale = DEFVAL;
   mc_ptHat = DEFVAL;

   weight_originalXWGTUP = DEFVAL;
   weight_scale_muF0p5 = DEFVAL;
   weight_scale_muF2   = DEFVAL;
   weight_scale_muR0p5 = DEFVAL;
   weight_scale_muR2   = DEFVAL;
   mc_pdfweights.clear();
   mc_pdfweightIds.clear();

   mc_pu_intime_NumInt = DEFVAL;
   mc_pu_trueNumInt = DEFVAL;
   mc_pu_before_npu = DEFVAL;
   mc_pu_after_npu = DEFVAL;

   mc_pu_Npvi = DEFVAL;
   mc_pu_Nzpositions.clear();
   mc_pu_BunchCrossing.clear();
   for(unsigned int i=0;i<mc_pu_zpositions.size();i++) mc_pu_zpositions[i].clear();
   mc_pu_zpositions.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_lowpT.size();i++) mc_pu_sumpT_lowpT[i].clear();
   mc_pu_sumpT_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_highpT.size();i++) mc_pu_sumpT_highpT[i].clear();
   mc_pu_sumpT_highpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_lowpT.size();i++) mc_pu_ntrks_lowpT[i].clear();
   mc_pu_ntrks_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_highpT.size();i++) mc_pu_ntrks_highpT[i].clear();
   mc_pu_ntrks_highpT.clear();

   trigger_n = 0;
   trigger.clear();
   trigger_name.clear();
   trigger_pass.clear();
   trigger_prescale.clear();
   trigger_HLTprescale.clear();
   trigger_L1prescale.clear();

   triggerobject_n = 0;
   triggerobject_pt.clear();
   triggerobject_eta.clear();
   triggerobject_phi.clear();

   triggerobject_collection.clear();

   triggerobject_filterIds_n.clear();
   triggerobject_filterIds.clear();

   triggerobject_isTriggerL1Mu.clear();
   triggerobject_isTriggerL1NoIsoEG.clear();
   triggerobject_isTriggerL1IsoEG.clear();
   triggerobject_isTriggerL1CenJet.clear();
   triggerobject_isTriggerL1ForJet.clear();
   triggerobject_isTriggerL1TauJet.clear();
   triggerobject_isTriggerL1ETM.clear();
   triggerobject_isTriggerL1ETT.clear();
   triggerobject_isTriggerL1HTT.clear();
   triggerobject_isTriggerL1HTM.clear();
   triggerobject_isTriggerL1JetCounts.clear();
   triggerobject_isTriggerL1HfBitCounts.clear();
   triggerobject_isTriggerL1HfRingEtSums.clear();
   triggerobject_isTriggerL1TechTrig.clear();
   triggerobject_isTriggerL1Castor.clear();
   triggerobject_isTriggerL1BPTX.clear();
   triggerobject_isTriggerL1GtExternal.clear();

   triggerobject_isHLT_TriggerPhoton.clear();
   triggerobject_isHLT_TriggerElectron.clear();
   triggerobject_isHLT_TriggerMuon.clear();
   triggerobject_isHLT_TriggerTau.clear();
   triggerobject_isHLT_TriggerJet.clear();
   triggerobject_isHLT_TriggerBJet.clear();
   triggerobject_isHLT_TriggerMET.clear();
   triggerobject_isHLT_TriggerTET.clear();
   triggerobject_isHLT_TriggerTHT.clear();
   triggerobject_isHLT_TriggerMHT.clear();
   triggerobject_isHLT_TriggerTrack.clear();
   triggerobject_isHLT_TriggerCluster.clear();
   triggerobject_isHLT_TriggerMETSig.clear();
   triggerobject_isHLT_TriggerELongit.clear();
   triggerobject_isHLT_TriggerMHTSig.clear();
   triggerobject_isHLT_TriggerHLongit.clear();

   triggerobject_filterLabels_n.clear();
   triggerobject_filterLabels.clear();

   triggerobject_pathNamesAll_n.clear();
   triggerobject_pathNamesAll.clear();
   triggerobject_pathNamesAll_isBoth.clear();
   triggerobject_pathNamesAll_isL3.clear();
   triggerobject_pathNamesAll_isLF.clear();
   triggerobject_pathNamesAll_isNone.clear();

   el_n = 0;
   el_pt.clear();
   el_eta.clear();
   el_phi.clear();
   el_m.clear();
   el_E.clear();
   el_id.clear();
   el_charge.clear();

   el_isGsfCtfScPixChargeConsistent.clear();
   el_isGsfScPixChargeConsistent.clear();
   el_passConversionVeto.clear();

   el_ip3d.clear();
   el_ip3dErr.clear();
   el_ip2d.clear();
   el_ip2dErr.clear();
   el_ip3dBS.clear();
   el_ip3dBSErr.clear();
   el_ip2dBS.clear();
   el_ip2dBSErr.clear();

   el_ecalEnergy.clear();
   el_correctedEcalEnergy.clear();
   el_correctedEcalEnergyError.clear();
   el_trackMomentumError.clear();

   el_neutralHadronIso.clear();
   el_chargedHadronIso.clear();
   el_puChargedHadronIso.clear();
   el_ecalIso.clear();
   el_hcalIso.clear();
   el_particleIso.clear();
   el_photonIso.clear();
   el_trackIso.clear();

   el_ecalPFClusterIso.clear();
   el_hcalPFClusterIso.clear();

   el_dr03EcalRecHitSumEt.clear();
   el_dr03HcalTowerSumEt.clear();
   el_dr03HcalDepth1TowerSumEt.clear();
   el_dr03HcalDepth2TowerSumEt.clear();
   el_dr03TkSumPt.clear();

   el_dr04EcalRecHitSumEt.clear();
   el_dr04HcalTowerSumEt.clear();
   el_dr04HcalDepth1TowerSumEt.clear();
   el_dr04HcalDepth2TowerSumEt.clear();
   el_dr04TkSumPt.clear();

   el_hcalOverEcal.clear();
   el_hcalOverEcalBc.clear();
   el_hcalDepth1OverEcal.clear();
   el_hcalDepth2OverEcal.clear();
   el_eSeedClusterOverPout.clear();
   el_eSeedClusterOverP.clear();
   el_eEleClusterOverPout.clear();
   el_deltaEtaEleClusterTrackAtCalo.clear();
   el_deltaPhiEleClusterTrackAtCalo.clear();

   el_pfIso_sumChargedHadronPt.clear();
   el_pfIso_sumNeutralHadronEt.clear();
   el_pfIso_sumPhotonEt.clear();
   el_pfIso_sumPUPt.clear();

   el_vx.clear();
   el_vy.clear();
   el_vz.clear();

   el_hasGsfTrack.clear();
   el_gsfTrack_d0.clear();
   el_gsfTrack_z0.clear();
   el_gsfTrack_d0Error.clear();
   el_gsfTrack_z0Error.clear();
   el_gsfTrack_PV_dxy.clear();
   el_gsfTrack_PV_dz.clear();
   el_gsfTrack_RP_dxy.clear();
   el_gsfTrack_RP_dz.clear();
   el_gsfTrack_BS_dxy.clear();
   el_gsfTrack_BS_dz.clear();
   el_gsfTrack_dxyError.clear();
   el_gsfTrack_dzError.clear();
   el_gsfTrack_normalizedChi2.clear();

   el_superCluster_eta.clear();
   el_superCluster_phi.clear();
   el_superCluster_energy.clear();
   el_superCluster_rawEnergy.clear();
   el_superCluster_preshowerEnergy.clear();
   el_superCluster_etaWidth.clear();
   el_superCluster_phiWidth.clear();
   el_superCluster_preshowerEnergyPlane1.clear();
   el_superCluster_preshowerEnergyPlane2.clear();
   el_superCluster_positionR.clear();

   el_basicClustersSize.clear();
   el_e1x5.clear();
   el_e5x5.clear();
   el_e2x5Max.clear();
   el_sigmaEtaEta.clear();
   el_sigmaIetaIeta.clear();
   el_sigmaIphiIphi.clear();
   el_sigmaIetaIphi.clear();
   el_full5x5_sigmaIphiIphi.clear();
   el_full5x5_sigmaEtaEta.clear();
   el_full5x5_sigmaIetaIeta.clear();
   el_full5x5_sigmaIetaIphi.clear();
   el_full5x5_r9.clear();
   el_full5x5_e1x5.clear();
   el_full5x5_e5x5.clear();
   el_full5x5_e2x5Max.clear();

   el_numberOfHits.clear();
   el_numberOfValidHits.clear();

   el_expectedMissingOuterHits.clear();
   el_numberOfValidPixelHits.clear();
   el_numberOfLostPixelHits.clear();
   el_trackerLayersWithMeasurement.clear();
   el_pixelLayersWithMeasurement.clear();
   el_numberOfValidStripLayersWithMonoAndStereo.clear();
   el_trackerLayersWithoutMeasurement.clear();

   el_hadronicOverEm.clear();
   el_numberOfLostHits.clear();
   el_numberOfLostHitsDefault.clear();

   el_fbrem.clear();
   el_kf_normalizedChi2.clear();
   el_gsf_normalizedChi2.clear();
   el_deltaEtaSuperClusterTrackAtVtx.clear();
   el_deltaPhiSuperClusterTrackAtVtx.clear();
   el_deltaEtaSeedClusterTrackAtCalo.clear();
   el_deltaPhiSeedClusterTrackAtCalo.clear();
   el_full5x5_OneMinusE1x5E5x5.clear();
   el_OneMinusE1x5E5x5.clear();
   el_eSuperClusterOverP.clear();
   el_IoEmIoP.clear();
   el_ooEmooP.clear();
   el_eleEoPout.clear();
   el_PreShowerOverRaw.clear();

   el_mvaNonTrigV0.clear();
   el_mvaNonTrigCat.clear();

   el_vetoCBId.clear();
   el_looseCBId.clear();
   el_mediumCBId.clear();
   el_tightCBId.clear();
   el_heepCBId.clear();

   el_mediumMVAId.clear();
   el_tightMVAId.clear();

   el_hasMCMatch.clear();
   el_gen_pt.clear();
   el_gen_eta.clear();
   el_gen_phi.clear();
   el_gen_m.clear();
   el_gen_status.clear();
   el_gen_id.clear();
   el_gen_charge.clear();
   el_gen_dr.clear();

   el_hasMCMatchPAT.clear();
   el_genPAT_pt.clear();
   el_genPAT_eta.clear();
   el_genPAT_phi.clear();
   el_genPAT_m.clear();
   el_genPAT_status.clear();
   el_genPAT_id.clear();
   el_genPAT_charge.clear();

   el_hasMatchedConversion.clear();
   el_expectedMissingInnerHits.clear();

   mu_n = 0;
   mu_pt.clear();
   mu_eta.clear();
   mu_phi.clear();
   mu_m.clear();
   mu_E.clear();
   mu_id.clear();
   mu_charge.clear();

   mu_ip3d.clear();
   mu_ip3dErr.clear();
   mu_ip2d.clear();
   mu_ip2dErr.clear();
   mu_ip3dBS.clear();
   mu_ip3dBSErr.clear();
   mu_ip2dBS.clear();
   mu_ip2dBSErr.clear();

   mu_neutralHadronIso.clear();
   mu_chargedHadronIso.clear();
   mu_puChargedHadronIso.clear();
   mu_ecalIso.clear();
   mu_hcalIso.clear();
   mu_photonIso.clear();
   mu_trackIso.clear();

   mu_pfIso03_sumChargedHadronPt.clear();
   mu_pfIso03_sumChargedParticlePt.clear();
   mu_pfIso03_sumNeutralHadronEt.clear();
   mu_pfIso03_sumNeutralHadronEtHighThreshold.clear();
   mu_pfIso03_sumPhotonEt.clear();
   mu_pfIso03_sumPhotonEtHighThreshold.clear();
   mu_pfIso03_sumPUPt.clear();

   mu_pfIso04_sumChargedHadronPt.clear();
   mu_pfIso04_sumChargedParticlePt.clear();
   mu_pfIso04_sumNeutralHadronEt.clear();
   mu_pfIso04_sumNeutralHadronEtHighThreshold.clear();
   mu_pfIso04_sumPhotonEt.clear();
   mu_pfIso04_sumPhotonEtHighThreshold.clear();
   mu_pfIso04_sumPUPt.clear();

   mu_pfMeanIso03_sumChargedHadronPt.clear();
   mu_pfMeanIso03_sumChargedParticlePt.clear();
   mu_pfMeanIso03_sumNeutralHadronEt.clear();
   mu_pfMeanIso03_sumNeutralHadronEtHighThreshold.clear();
   mu_pfMeanIso03_sumPhotonEt.clear();
   mu_pfMeanIso03_sumPhotonEtHighThreshold.clear();
   mu_pfMeanIso03_sumPUPt.clear();

   mu_pfSumIso03_sumChargedHadronPt.clear();
   mu_pfSumIso03_sumChargedParticlePt.clear();
   mu_pfSumIso03_sumNeutralHadronEt.clear();
   mu_pfSumIso03_sumNeutralHadronEtHighThreshold.clear();
   mu_pfSumIso03_sumPhotonEt.clear();
   mu_pfSumIso03_sumPhotonEtHighThreshold.clear();
   mu_pfSumIso03_sumPUPt.clear();

   mu_pfMeanIso04_sumChargedHadronPt.clear();
   mu_pfMeanIso04_sumChargedParticlePt.clear();
   mu_pfMeanIso04_sumNeutralHadronEt.clear();
   mu_pfMeanIso04_sumNeutralHadronEtHighThreshold.clear();
   mu_pfMeanIso04_sumPhotonEt.clear();
   mu_pfMeanIso04_sumPhotonEtHighThreshold.clear();
   mu_pfMeanIso04_sumPUPt.clear();

   mu_pfSumIso04_sumChargedHadronPt.clear();
   mu_pfSumIso04_sumChargedParticlePt.clear();
   mu_pfSumIso04_sumNeutralHadronEt.clear();
   mu_pfSumIso04_sumNeutralHadronEtHighThreshold.clear();
   mu_pfSumIso04_sumPhotonEt.clear();
   mu_pfSumIso04_sumPhotonEtHighThreshold.clear();
   mu_pfSumIso04_sumPUPt.clear();

   mu_isGlobalMuon.clear();
   mu_isTrackerMuon.clear();
   mu_isStandAloneMuon.clear();
   mu_isCaloMuon.clear();
   mu_isPFMuon.clear();
   mu_isRPCMuon.clear();

   mu_vx.clear();
   mu_vy.clear();
   mu_vz.clear();

   mu_numberOfMatches.clear();
   mu_numberOfMatchedStations.clear();

   mu_segmentCompatibility.clear();
   mu_caloCompatibility.clear();

   mu_combinedQuality_chi2LocalPosition.clear();
   mu_combinedQuality_trkKink.clear();

   mu_isLooseMuon.clear();
   mu_isMediumMuon.clear();
   mu_isTightMuon.clear();
   mu_isSoftMuon.clear();
   mu_isHighPtMuon.clear();

   mu_isGoodMuon_AllGlobalMuons.clear();
   mu_isGoodMuon_AllStandAloneMuons.clear();
   mu_isGoodMuon_AllTrackerMuons.clear();
   mu_isGoodMuon_TrackerMuonArbitrated.clear();
   mu_isGoodMuon_AllArbitrated.clear();
   mu_isGoodMuon_GlobalMuonPromptTight.clear();
   mu_isGoodMuon_TMLastStationLoose.clear();
   mu_isGoodMuon_TMLastStationTight.clear();
   mu_isGoodMuon_TM2DCompatibilityLoose.clear();
   mu_isGoodMuon_TM2DCompatibilityTight.clear();
   mu_isGoodMuon_TMOneStationLoose.clear();
   mu_isGoodMuon_TMOneStationTight.clear();
   mu_isGoodMuon_TMLastStationOptimizedLowPtLoose.clear();
   mu_isGoodMuon_TMLastStationOptimizedLowPtTight.clear();
   mu_isGoodMuon_GMTkChiCompatibility.clear();
   mu_isGoodMuon_GMStaChiCompatibility.clear();
   mu_isGoodMuon_GMTkKinkTight.clear();
   mu_isGoodMuon_TMLastStationAngLoose.clear();
   mu_isGoodMuon_TMLastStationAngTight.clear();
   mu_isGoodMuon_TMOneStationAngLoose.clear();
   mu_isGoodMuon_TMOneStationAngTight.clear();
   mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose.clear();
   mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.clear();

   mu_calEnergy_em.clear();
   mu_calEnergy_had.clear();
   mu_calEnergy_ho.clear();
   mu_calEnergy_emS9.clear();
   mu_calEnergy_hadS9.clear();
   mu_calEnergy_hoS9.clear();
   mu_calEnergy_emS25.clear();
   mu_calEnergy_emMax.clear();
   mu_calEnergy_hadMax.clear();
   mu_calEnergy_ecal_time.clear();
   mu_calEnergy_hcal_time.clear();
   mu_calEnergy_ecal_rawId.clear();
   mu_calEnergy_hcal_rawId.clear();

   mu_isolationR03_trackerVetoPt.clear();
   mu_isolationR03_emVetoEt.clear();
   mu_isolationR03_hadVetoEt.clear();
   mu_isolationR03_hoVetoEt.clear();
   mu_isolationR03_sumPt.clear();
   mu_isolationR03_emEt.clear();
   mu_isolationR03_hadEt.clear();
   mu_isolationR03_hoEt.clear();
   mu_isolationR03_nTracks.clear();
   mu_isolationR03_nJets.clear();

   mu_isolationR05_trackerVetoPt.clear();
   mu_isolationR05_emVetoEt.clear();
   mu_isolationR05_hadVetoEt.clear();
   mu_isolationR05_hoVetoEt.clear();
   mu_isolationR05_sumPt.clear();
   mu_isolationR05_emEt.clear();
   mu_isolationR05_hadEt.clear();
   mu_isolationR05_hoEt.clear();
   mu_isolationR05_nTracks.clear();
   mu_isolationR05_nJets.clear();

   mu_hasGlobalTrack.clear();
   mu_globalTrack_d0.clear();
   mu_globalTrack_z0.clear();
   mu_globalTrack_d0Error.clear();
   mu_globalTrack_z0Error.clear();
   mu_globalTrack_RP_dxy.clear();
   mu_globalTrack_RP_dz.clear();
   mu_globalTrack_PV_dxy.clear();
   mu_globalTrack_PV_dz.clear();
   mu_globalTrack_BS_dxy.clear();
   mu_globalTrack_BS_dz.clear();
   mu_globalTrack_dxyError.clear();
   mu_globalTrack_dzError.clear();
   mu_globalTrack_normalizedChi2.clear();
   mu_globalTrack_numberOfValidHits.clear();
   mu_globalTrack_numberOfValidMuonHits.clear();
   mu_globalTrack_numberOfLostHits.clear();
   mu_globalTrack_pt.clear();
   mu_globalTrack_eta.clear();
   mu_globalTrack_phi.clear();
   mu_globalTrack_ptError.clear();
   mu_globalTrack_etaError.clear();
   mu_globalTrack_phiError.clear();
   mu_globalTrack_vx.clear();
   mu_globalTrack_vy.clear();
   mu_globalTrack_vz.clear();
   mu_globalTrack_qoverp.clear();
   mu_globalTrack_qoverpError.clear();
   mu_globalTrack_charge.clear();
   mu_globalTrack_trackerLayersWithMeasurement.clear();
   mu_globalTrack_pixelLayersWithMeasurement.clear();
   mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo.clear();
   mu_globalTrack_trackerLayersWithoutMeasurement.clear();
   mu_globalTrack_numberOfValidPixelHits.clear();
   mu_globalTrack_numberOfLostPixelHits.clear();
   mu_globalTrack_numberOfInnerHits.clear();
   mu_globalTrack_numberOfOuterHits.clear();
   mu_globalTrack_validFraction.clear();

   mu_bestTrackType.clear();
   mu_hasBestTrack.clear();
   mu_bestTrack_d0.clear();
   mu_bestTrack_z0.clear();
   mu_bestTrack_d0Error.clear();
   mu_bestTrack_z0Error.clear();
   mu_bestTrack_RP_dxy.clear();
   mu_bestTrack_RP_dz.clear();
   mu_bestTrack_PV_dxy.clear();
   mu_bestTrack_PV_dz.clear();
   mu_bestTrack_BS_dxy.clear();
   mu_bestTrack_BS_dz.clear();
   mu_bestTrack_dxyError.clear();
   mu_bestTrack_dzError.clear();
   mu_bestTrack_normalizedChi2.clear();
   mu_bestTrack_numberOfValidHits.clear();
   mu_bestTrack_numberOfLostHits.clear();
   mu_bestTrack_pt.clear();
   mu_bestTrack_eta.clear();
   mu_bestTrack_phi.clear();
   mu_bestTrack_ptError.clear();
   mu_bestTrack_etaError.clear();
   mu_bestTrack_phiError.clear();
   mu_bestTrack_vx.clear();
   mu_bestTrack_vy.clear();
   mu_bestTrack_vz.clear();
   mu_bestTrack_qoverp.clear();
   mu_bestTrack_qoverpError.clear();
   mu_bestTrack_charge.clear();
   mu_bestTrack_trackerLayersWithMeasurement.clear();
   mu_bestTrack_pixelLayersWithMeasurement.clear();
   mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo.clear();
   mu_bestTrack_trackerLayersWithoutMeasurement.clear();
   mu_bestTrack_numberOfValidPixelHits.clear();
   mu_bestTrack_numberOfLostPixelHits.clear();
   mu_bestTrack_numberOfInnerHits.clear();
   mu_bestTrack_numberOfOuterHits.clear();
   mu_bestTrack_validFraction.clear();

   mu_hasInnerTrack.clear();
   mu_innerTrack_d0.clear();
   mu_innerTrack_z0.clear();
   mu_innerTrack_d0Error.clear();
   mu_innerTrack_z0Error.clear();
   mu_innerTrack_RP_dxy.clear();
   mu_innerTrack_RP_dz.clear();
   mu_innerTrack_PV_dxy.clear();
   mu_innerTrack_PV_dz.clear();
   mu_innerTrack_BS_dxy.clear();
   mu_innerTrack_BS_dz.clear();
   mu_innerTrack_dxyError.clear();
   mu_innerTrack_dzError.clear();
   mu_innerTrack_normalizedChi2.clear();
   mu_innerTrack_numberOfValidHits.clear();
   mu_innerTrack_numberOfLostHits.clear();
   mu_innerTrack_pt.clear();
   mu_innerTrack_eta.clear();
   mu_innerTrack_phi.clear();
   mu_innerTrack_ptError.clear();
   mu_innerTrack_etaError.clear();
   mu_innerTrack_phiError.clear();
   mu_innerTrack_vx.clear();
   mu_innerTrack_vy.clear();
   mu_innerTrack_vz.clear();
   mu_innerTrack_qoverp.clear();
   mu_innerTrack_qoverpError.clear();
   mu_innerTrack_charge.clear();
   mu_innerTrack_trackerLayersWithMeasurement.clear();
   mu_innerTrack_pixelLayersWithMeasurement.clear();
   mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo.clear();
   mu_innerTrack_trackerLayersWithoutMeasurement.clear();
   mu_innerTrack_numberOfValidPixelHits.clear();
   mu_innerTrack_numberOfLostPixelHits.clear();
   mu_innerTrack_numberOfInnerHits.clear();
   mu_innerTrack_numberOfOuterHits.clear();
   mu_innerTrack_validFraction.clear();

   mu_type.clear();

   mu_hasMCMatch.clear();
   mu_gen_pt.clear();
   mu_gen_eta.clear();
   mu_gen_phi.clear();
   mu_gen_m.clear();
   mu_gen_status.clear();
   mu_gen_id.clear();
   mu_gen_charge.clear();
   mu_gen_dr.clear();

   mu_hasMCMatchPAT.clear();
   mu_genPAT_pt.clear();
   mu_genPAT_eta.clear();
   mu_genPAT_phi.clear();
   mu_genPAT_m.clear();
   mu_genPAT_status.clear();
   mu_genPAT_id.clear();
   mu_genPAT_charge.clear();

   tau_n = 0;
   tau_pt.clear();
   tau_eta.clear();
   tau_phi.clear();
   tau_m.clear();
   tau_E.clear();
   tau_id.clear();
   tau_charge.clear();

   tau_hasLeadChargedHadrCand.clear();
   tau_leadingTrackPt.clear();
   tau_leadingTrackCharge.clear();
   tau_leadingTrackDz.clear();
   tau_leadingTrackDxy.clear();

   tau_decayMode.clear();
   tau_decayModeFindingOldDMs.clear();
   tau_decayModeFindingNewDMs.clear();

   tau_puCorrPtSum.clear();
   tau_neutralIsoPtSum.clear();
   tau_chargedIsoPtSum.clear();
   tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();

   tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
   tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
   tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();

   tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
   tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
   tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
   tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.clear();

   tau_againstMuonLoose3.clear();
   tau_againstMuonTight3.clear();

   tau_pfEssential_jet_pt.clear();
   tau_pfEssential_jet_eta.clear();
   tau_pfEssential_jet_phi.clear();
   tau_pfEssential_jet_m.clear();

   tau_pfEssential_jetCorr_pt.clear();
   tau_pfEssential_jetCorr_eta.clear();
   tau_pfEssential_jetCorr_phi.clear();
   tau_pfEssential_jetCorr_m.clear();

   tau_pfEssential_hasSV.clear();
   tau_pfEssential_sv_x.clear();
   tau_pfEssential_sv_y.clear();
   tau_pfEssential_sv_z.clear();

   tau_pfEssential_flightLengthSig.clear();
   tau_pfEssential_dxy.clear();
   tau_pfEssential_dxy_error.clear();
   tau_pfEssential_dxy_Sig.clear();

   jet_n = 0;
   jet_pt.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_m.clear();
   jet_E.clear();

   jet_ntrk.clear();

   jet_JBP.clear();
   jet_JP.clear();
   jet_TCHP.clear();
   jet_TCHE.clear();
   jet_SSVHE.clear();
   jet_SSVHP.clear();
   jet_CMVA.clear();
   jet_CSVv2.clear();
   jet_cMVAv2.clear();
   jet_CharmCvsL.clear();
   jet_CharmCvsB.clear();
   jet_partonFlavour.clear();
   jet_hadronFlavour.clear();

   jet_neutralHadronEnergy.clear();
   jet_neutralEmEnergy.clear();
   jet_chargedHadronEnergy.clear();
   jet_chargedEmEnergy.clear();
   jet_electronEnergy.clear();
   jet_muonEnergy.clear();
   jet_photonEnergy.clear();

   jet_charge.clear();
   jet_chargeVec.clear();
   jet_chargedMultiplicity.clear();
   jet_neutralMultiplicity.clear();
   jet_chargedHadronMultiplicity.clear();

   jet_jetArea.clear();

   jet_jecFactorUncorrected.clear();
   jet_jecFactorL1FastJet.clear();
   jet_jecFactorL2Relative.clear();
   jet_jecFactorL3Absolute.clear();

   jet_neutralHadronEnergyFraction.clear();
   jet_neutralEmEnergyFraction.clear();
   jet_chargedHadronEnergyFraction.clear();
   jet_muonEnergyFraction.clear();
   jet_chargedEmEnergyFraction.clear();

   jet_Unc.clear();

   jet_pileupJetId.clear();

   jet_looseJetID.clear();
   jet_tightJetID.clear();

   jet_hasGenJet.clear();
   jet_genJet_pt.clear();
   jet_genJet_eta.clear();
   jet_genJet_phi.clear();
   jet_genJet_m.clear();
   jet_genJet_E.clear();
   jet_genJet_status.clear();
   jet_genJet_id.clear();

   jet_hasGenParton.clear();
   jet_genParton_pt.clear();
   jet_genParton_eta.clear();
   jet_genParton_phi.clear();
   jet_genParton_m.clear();
   jet_genParton_E.clear();
   jet_genParton_status.clear();
   jet_genParton_id.clear();

   //------------------------
   //  GenJet collection
   //------------------------

   genJet_n = 0;
   genJet_pt.clear();
   genJet_eta.clear();
   genJet_phi.clear();
   genJet_m.clear();
   genJet_E.clear();
   genJet_emEnergy.clear();
   genJet_hadEnergy.clear();
   genJet_invisibleEnergy.clear();
   genJet_auxiliaryEnergy.clear();
   genJet_flavour.clear();
}

void FlatTree::CreateBranches(int buffersize = 32000)
{
   if( doWrite("ev_run") ) tree->Branch("ev_run", &ev_run, "ev_run/I", buffersize);
   if( doWrite("ev_id") ) tree->Branch("ev_id", &ev_id, "ev_id/I", buffersize);
   if( doWrite("ev_lumi") ) tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/I", buffersize);
   if( doWrite("ev_rho") ) tree->Branch("ev_rho", &ev_rho, "ev_rho/F", buffersize);

   if( doWrite("met_px") ) tree->Branch("met_px", &met_px, "met_px/F", buffersize);
   if( doWrite("met_py") ) tree->Branch("met_py", &met_py, "met_py/F", buffersize);
   if( doWrite("met_pt") ) tree->Branch("met_pt", &met_pt, "met_pt/F", buffersize);
   if( doWrite("met_phi") ) tree->Branch("met_phi", &met_phi, "met_phi/F", buffersize);
   if( doWrite("met_sumet") ) tree->Branch("met_sumet", &met_sumet, "met_sumet/F", buffersize);
   if( doWrite("met_sig") ) tree->Branch("met_sig", &met_sig, "met_sig/D", buffersize);

   if( doWrite("met_cov00") ) tree->Branch("met_cov00", &met_cov00, "met_cov00/D", buffersize);
   if( doWrite("met_cov10") ) tree->Branch("met_cov10", &met_cov10, "met_cov10/D", buffersize);
   if( doWrite("met_cov01") ) tree->Branch("met_cov01", &met_cov01, "met_cov01/D", buffersize);
   if( doWrite("met_cov11") ) tree->Branch("met_cov11", &met_cov11, "met_cov11/D", buffersize);

   if( doWrite("metGen_px") ) tree->Branch("metGen_px", &metGen_px, "metGen_px/F", buffersize);
   if( doWrite("metGen_py") ) tree->Branch("metGen_py", &metGen_py, "metGen_py/F", buffersize);
   if( doWrite("metGen_pt") ) tree->Branch("metGen_pt", &metGen_pt, "metGen_pt/F", buffersize);
   if( doWrite("metGen_phi") ) tree->Branch("metGen_phi", &metGen_phi, "metGen_phi/F", buffersize);
   if( doWrite("metGen_sumet") ) tree->Branch("metGen_sumet", &metGen_sumet, "metGen_sumet/F", buffersize);

   if( doWrite("metGen_NeutralEMEt") ) tree->Branch("metGen_NeutralEMEt", &metGen_NeutralEMEt, "metGen_NeutralEMEt/F", buffersize);
   if( doWrite("metGen_ChargedEMEt") ) tree->Branch("metGen_ChargedEMEt", &metGen_ChargedEMEt, "metGen_ChargedEMEt/F", buffersize);
   if( doWrite("metGen_NeutralHadEt") ) tree->Branch("metGen_NeutralHadEt", &metGen_NeutralHadEt, "metGen_NeutralHadEt/F", buffersize);
   if( doWrite("metGen_ChargedHadEt") ) tree->Branch("metGen_ChargedHadEt", &metGen_ChargedHadEt, "metGen_ChargedHadEt/F", buffersize);
   if( doWrite("metGen_MuonEt") ) tree->Branch("metGen_MuonEt", &metGen_MuonEt, "metGen_MuonEt/F", buffersize);
   if( doWrite("metGen_InvisibleEt") ) tree->Branch("metGen_InvisibleEt", &metGen_InvisibleEt, "metGen_InvisibleEt/F", buffersize);

   if( doWrite("met_uncorrectedPt") ) tree->Branch("met_uncorrectedPt", &met_uncorrectedPt, "met_uncorrectedPt/F", buffersize);
   if( doWrite("met_uncorrectedPhi") ) tree->Branch("met_uncorrectedPhi", &met_uncorrectedPhi, "met_uncorrectedPhi/F", buffersize);
   if( doWrite("met_uncorrectedSumEt") ) tree->Branch("met_uncorrectedSumEt", &met_uncorrectedSumEt, "met_uncorrectedSumEt/F", buffersize);

   if( doWrite("met_caloMETPt") ) tree->Branch("met_caloMETPt", &met_caloMETPt, "met_caloMETPt/F", buffersize);
   if( doWrite("met_caloMETPhi") ) tree->Branch("met_caloMETPhi", &met_caloMETPhi, "met_caloMETPhi/F", buffersize);
   if( doWrite("met_caloMETSumEt") ) tree->Branch("met_caloMETSumEt", &met_caloMETSumEt, "met_caloMETSumEt/F", buffersize);

   if( doWrite("met_shiftedPx_JetEnUp") ) tree->Branch("met_shiftedPx_JetEnUp", &met_shiftedPx_JetEnUp, "met_shiftedPx_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_JetEnDown") ) tree->Branch("met_shiftedPx_JetEnDown", &met_shiftedPx_JetEnDown, "met_shiftedPx_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_JetResUp") ) tree->Branch("met_shiftedPx_JetResUp", &met_shiftedPx_JetResUp, "met_shiftedPx_JetResUp/F", buffersize);
   if( doWrite("met_shiftedPx_JetResDown") ) tree->Branch("met_shiftedPx_JetResDown", &met_shiftedPx_JetResDown, "met_shiftedPx_JetResDown/F", buffersize);
   if( doWrite("met_shiftedPx_MuonEnUp") ) tree->Branch("met_shiftedPx_MuonEnUp", &met_shiftedPx_MuonEnUp, "met_shiftedPx_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_MuonEnDown") ) tree->Branch("met_shiftedPx_MuonEnDown", &met_shiftedPx_MuonEnDown, "met_shiftedPx_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_ElectronEnUp") ) tree->Branch("met_shiftedPx_ElectronEnUp", &met_shiftedPx_ElectronEnUp, "met_shiftedPx_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_ElectronEnDown") ) tree->Branch("met_shiftedPx_ElectronEnDown", &met_shiftedPx_ElectronEnDown, "met_shiftedPx_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_TauEnUp") ) tree->Branch("met_shiftedPx_TauEnUp", &met_shiftedPx_TauEnUp, "met_shiftedPx_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_TauEnDown") ) tree->Branch("met_shiftedPx_TauEnDown", &met_shiftedPx_TauEnDown, "met_shiftedPx_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_UnclusteredEnUp") ) tree->Branch("met_shiftedPx_UnclusteredEnUp", &met_shiftedPx_UnclusteredEnUp, "met_shiftedPx_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_UnclusteredEnDown") ) tree->Branch("met_shiftedPx_UnclusteredEnDown", &met_shiftedPx_UnclusteredEnDown, "met_shiftedPx_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_NoShift") ) tree->Branch("met_shiftedPx_NoShift", &met_shiftedPx_NoShift, "met_shiftedPx_NoShift/F", buffersize);
   if( doWrite("met_shiftedPx_PhotonEnUp") ) tree->Branch("met_shiftedPx_PhotonEnUp", &met_shiftedPx_PhotonEnUp, "met_shiftedPx_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_PhotonEnDown") ) tree->Branch("met_shiftedPx_PhotonEnDown", &met_shiftedPx_PhotonEnDown, "met_shiftedPx_PhotonEnDown/F", buffersize);

   if( doWrite("met_shiftedPy_JetEnUp") ) tree->Branch("met_shiftedPy_JetEnUp", &met_shiftedPy_JetEnUp, "met_shiftedPy_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_JetEnDown") ) tree->Branch("met_shiftedPy_JetEnDown", &met_shiftedPy_JetEnDown, "met_shiftedPy_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_JetResUp") ) tree->Branch("met_shiftedPy_JetResUp", &met_shiftedPy_JetResUp, "met_shiftedPy_JetResUp/F", buffersize);
   if( doWrite("met_shiftedPy_JetResDown") ) tree->Branch("met_shiftedPy_JetResDown", &met_shiftedPy_JetResDown, "met_shiftedPy_JetResDown/F", buffersize);
   if( doWrite("met_shiftedPy_MuonEnUp") ) tree->Branch("met_shiftedPy_MuonEnUp", &met_shiftedPy_MuonEnUp, "met_shiftedPy_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_MuonEnDown") ) tree->Branch("met_shiftedPy_MuonEnDown", &met_shiftedPy_MuonEnDown, "met_shiftedPy_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_ElectronEnUp") ) tree->Branch("met_shiftedPy_ElectronEnUp", &met_shiftedPy_ElectronEnUp, "met_shiftedPy_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_ElectronEnDown") ) tree->Branch("met_shiftedPy_ElectronEnDown", &met_shiftedPy_ElectronEnDown, "met_shiftedPy_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_TauEnUp") ) tree->Branch("met_shiftedPy_TauEnUp", &met_shiftedPy_TauEnUp, "met_shiftedPy_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_TauEnDown") ) tree->Branch("met_shiftedPy_TauEnDown", &met_shiftedPy_TauEnDown, "met_shiftedPy_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_UnclusteredEnUp") ) tree->Branch("met_shiftedPy_UnclusteredEnUp", &met_shiftedPy_UnclusteredEnUp, "met_shiftedPy_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_UnclusteredEnDown") ) tree->Branch("met_shiftedPy_UnclusteredEnDown", &met_shiftedPy_UnclusteredEnDown, "met_shiftedPy_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_NoShift") ) tree->Branch("met_shiftedPy_NoShift", &met_shiftedPy_NoShift, "met_shiftedPy_NoShift/F", buffersize);
   if( doWrite("met_shiftedPy_PhotonEnUp") ) tree->Branch("met_shiftedPy_PhotonEnUp", &met_shiftedPy_PhotonEnUp, "met_shiftedPy_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_PhotonEnDown") ) tree->Branch("met_shiftedPy_PhotonEnDown", &met_shiftedPy_PhotonEnDown, "met_shiftedPy_PhotonEnDown/F", buffersize);

   if( doWrite("met_shiftedPhi_JetEnUp") ) tree->Branch("met_shiftedPhi_JetEnUp", &met_shiftedPhi_JetEnUp, "met_shiftedPhi_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_JetEnDown") ) tree->Branch("met_shiftedPhi_JetEnDown", &met_shiftedPhi_JetEnDown, "met_shiftedPhi_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_JetResUp") ) tree->Branch("met_shiftedPhi_JetResUp", &met_shiftedPhi_JetResUp, "met_shiftedPhi_JetResUp/F", buffersize);
   if( doWrite("met_shiftedPhi_JetResDown") ) tree->Branch("met_shiftedPhi_JetResDown", &met_shiftedPhi_JetResDown, "met_shiftedPhi_JetResDown/F", buffersize);
   if( doWrite("met_shiftedPhi_MuonEnUp") ) tree->Branch("met_shiftedPhi_MuonEnUp", &met_shiftedPhi_MuonEnUp, "met_shiftedPhi_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_MuonEnDown") ) tree->Branch("met_shiftedPhi_MuonEnDown", &met_shiftedPhi_MuonEnDown, "met_shiftedPhi_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_ElectronEnUp") ) tree->Branch("met_shiftedPhi_ElectronEnUp", &met_shiftedPhi_ElectronEnUp, "met_shiftedPhi_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_ElectronEnDown") ) tree->Branch("met_shiftedPhi_ElectronEnDown", &met_shiftedPhi_ElectronEnDown, "met_shiftedPhi_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_TauEnUp") ) tree->Branch("met_shiftedPhi_TauEnUp", &met_shiftedPhi_TauEnUp, "met_shiftedPhi_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_TauEnDown") ) tree->Branch("met_shiftedPhi_TauEnDown", &met_shiftedPhi_TauEnDown, "met_shiftedPhi_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_UnclusteredEnUp") ) tree->Branch("met_shiftedPhi_UnclusteredEnUp", &met_shiftedPhi_UnclusteredEnUp, "met_shiftedPhi_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_UnclusteredEnDown") ) tree->Branch("met_shiftedPhi_UnclusteredEnDown", &met_shiftedPhi_UnclusteredEnDown, "met_shiftedPhi_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_NoShift") ) tree->Branch("met_shiftedPhi_NoShift", &met_shiftedPhi_NoShift, "met_shiftedPhi_NoShift/F", buffersize);
   if( doWrite("met_shiftedPhi_PhotonEnUp") ) tree->Branch("met_shiftedPhi_PhotonEnUp", &met_shiftedPhi_PhotonEnUp, "met_shiftedPhi_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_PhotonEnDown") ) tree->Branch("met_shiftedPhi_PhotonEnDown", &met_shiftedPhi_PhotonEnDown, "met_shiftedPhi_PhotonEnDown/F", buffersize);

   if( doWrite("met_shiftedSumEt_JetEnUp") ) tree->Branch("met_shiftedSumEt_JetEnUp", &met_shiftedSumEt_JetEnUp, "met_shiftedSumEt_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_JetEnDown") ) tree->Branch("met_shiftedSumEt_JetEnDown", &met_shiftedSumEt_JetEnDown, "met_shiftedSumEt_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_JetResUp") ) tree->Branch("met_shiftedSumEt_JetResUp", &met_shiftedSumEt_JetResUp, "met_shiftedSumEt_JetResUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_JetResDown") ) tree->Branch("met_shiftedSumEt_JetResDown", &met_shiftedSumEt_JetResDown, "met_shiftedSumEt_JetResDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_MuonEnUp") ) tree->Branch("met_shiftedSumEt_MuonEnUp", &met_shiftedSumEt_MuonEnUp, "met_shiftedSumEt_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_MuonEnDown") ) tree->Branch("met_shiftedSumEt_MuonEnDown", &met_shiftedSumEt_MuonEnDown, "met_shiftedSumEt_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_ElectronEnUp") ) tree->Branch("met_shiftedSumEt_ElectronEnUp", &met_shiftedSumEt_ElectronEnUp, "met_shiftedSumEt_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_ElectronEnDown") ) tree->Branch("met_shiftedSumEt_ElectronEnDown", &met_shiftedSumEt_ElectronEnDown, "met_shiftedSumEt_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_TauEnUp") ) tree->Branch("met_shiftedSumEt_TauEnUp", &met_shiftedSumEt_TauEnUp, "met_shiftedSumEt_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_TauEnDown") ) tree->Branch("met_shiftedSumEt_TauEnDown", &met_shiftedSumEt_TauEnDown, "met_shiftedSumEt_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_UnclusteredEnUp") ) tree->Branch("met_shiftedSumEt_UnclusteredEnUp", &met_shiftedSumEt_UnclusteredEnUp, "met_shiftedSumEt_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_UnclusteredEnDown") ) tree->Branch("met_shiftedSumEt_UnclusteredEnDown", &met_shiftedSumEt_UnclusteredEnDown, "met_shiftedSumEt_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_NoShift") ) tree->Branch("met_shiftedSumEt_NoShift", &met_shiftedSumEt_NoShift, "met_shiftedSumEt_NoShift/F", buffersize);
   if( doWrite("met_shiftedSumEt_PhotonEnUp") ) tree->Branch("met_shiftedSumEt_PhotonEnUp", &met_shiftedSumEt_PhotonEnUp, "met_shiftedSumEt_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_PhotonEnDown") ) tree->Branch("met_shiftedSumEt_PhotonEnDown", &met_shiftedSumEt_PhotonEnDown, "met_shiftedSumEt_PhotonEnDown/F", buffersize);

   if( doWrite("metNoHF_pt") )     tree->Branch("metNoHF_pt",    &metNoHF_pt,    "metNoHF_pt/F",    buffersize);
   if( doWrite("metNoHF_phi") )    tree->Branch("metNoHF_phi",   &metNoHF_phi,   "metNoHF_phi/F",   buffersize);
   if( doWrite("metNoHF_sumet") )  tree->Branch("metNoHF_sumet", &metNoHF_sumet, "metNoHF_sumet/F", buffersize);

   if( doWrite("pv_n") ) tree->Branch("pv_n", &pv_n, "pv_n/I", buffersize);
   if( doWrite("pv_x") ) tree->Branch("pv_x", &pv_x, "pv_x/F", buffersize);
   if( doWrite("pv_y") ) tree->Branch("pv_y", &pv_y, "pv_y/F", buffersize);
   if( doWrite("pv_z") ) tree->Branch("pv_z", &pv_z, "pv_z/F", buffersize);

   if( doWrite("pv_xError") ) tree->Branch("pv_xError", &pv_xError, "pv_xError/F", buffersize);
   if( doWrite("pv_yError") ) tree->Branch("pv_yError", &pv_yError, "pv_yError/F", buffersize);
   if( doWrite("pv_zError") ) tree->Branch("pv_zError", &pv_zError, "pv_zError/F", buffersize);

   if( doWrite("pv_chi2") ) tree->Branch("pv_chi2", &pv_chi2, "pv_chi2/F", buffersize);
   if( doWrite("pv_ndof") ) tree->Branch("pv_ndof", &pv_ndof, "pv_ndof/I", buffersize);
   if( doWrite("pv_rho") ) tree->Branch("pv_rho", &pv_rho, "pv_rho/F", buffersize);
   if( doWrite("pv_isFake") ) tree->Branch("pv_isFake", &pv_isFake, "pv_isFake/I", buffersize);

   if( doWrite("mc_weight") ) tree->Branch("mc_weight", &mc_weight, "mc_weight/F", buffersize);
   if( doWrite("mc_id") ) tree->Branch("mc_id", &mc_id, "mc_id/I", buffersize);
   if( doWrite("mc_f1") ) tree->Branch("mc_f1", &mc_f1, "mc_f1/I", buffersize);
   if( doWrite("mc_f2") ) tree->Branch("mc_f2", &mc_f2, "mc_f2/I", buffersize);
   if( doWrite("mc_x1") ) tree->Branch("mc_x1", &mc_x1, "mc_x1/F", buffersize);
   if( doWrite("mc_x2") ) tree->Branch("mc_x2", &mc_x2, "mc_x2/F", buffersize);
   if( doWrite("mc_scale") ) tree->Branch("mc_scale", &mc_scale, "mc_scale/F", buffersize);
   if( doWrite("mc_ptHat") ) tree->Branch("mc_ptHat", &mc_ptHat, "mc_ptHat/F", buffersize);

   if( doWrite("weight_originalXWGTUP") ) tree->Branch("weight_originalXWGTUP", &weight_originalXWGTUP, "weight_originalXWGTUP/F", buffersize);
   if( doWrite("weight_scale_muF0p5") ) tree->Branch("weight_scale_muF0p5", &weight_scale_muF0p5, "weight_scale_muF0p5/F", buffersize);
   if( doWrite("weight_scale_muF2"  ) ) tree->Branch("weight_scale_muF2",   &weight_scale_muF2,   "weight_scale_muF2/F", buffersize);
   if( doWrite("weight_scale_muR0p5") ) tree->Branch("weight_scale_muR0p5", &weight_scale_muR0p5, "weight_scale_muR0p5/F", buffersize);
   if( doWrite("weight_scale_muR2"  ) ) tree->Branch("weight_scale_muR2",   &weight_scale_muR2,   "weight_scale_muR2/F", buffersize);
   if( doWrite("mc_pdfweights") ) tree->Branch("mc_pdfweights", "std::vector<float>", &mc_pdfweights, buffersize);
   if( doWrite("mc_pdfweightIds") ) tree->Branch("mc_pdfweightIds", "std::vector<std::string>", &mc_pdfweightIds, buffersize);

   if( doWrite("mc_pu_intime_NumInt") ) tree->Branch("mc_pu_intime_NumInt", &mc_pu_intime_NumInt, "mc_pu_intime_NumInt/I", buffersize);
   if( doWrite("mc_pu_trueNumInt") ) tree->Branch("mc_pu_trueNumInt", &mc_pu_trueNumInt, "mc_pu_trueNumInt/I", buffersize);
   if( doWrite("mc_pu_before_npu") ) tree->Branch("mc_pu_before_npu", &mc_pu_before_npu, "mc_pu_before_npu/I", buffersize);
   if( doWrite("mc_pu_after_npu") ) tree->Branch("mc_pu_after_npu", &mc_pu_after_npu, "mc_pu_after_npu/I", buffersize);

   if( doWrite("mc_pu_Npvi") ) tree->Branch("mc_pu_Npvi", &mc_pu_Npvi, "mc_pu_Npvi/I", buffersize);
   if( doWrite("mc_pu_Nzpositions") ) tree->Branch("mc_pu_Nzpositions", "std::vector<int>", &mc_pu_Nzpositions, buffersize);
   if( doWrite("mc_pu_BunchCrossing") ) tree->Branch("mc_pu_BunchCrossing", "std::vector<int>", &mc_pu_BunchCrossing, buffersize );
   if( doWrite("mc_pu_zpositions") ) tree->Branch("mc_pu_zpositions", "std::vector<std::vector<float> >", &mc_pu_zpositions, buffersize);
   if( doWrite("mc_pu_sumpT_lowpT") ) tree->Branch("mc_pu_sumpT_lowpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_lowpT, buffersize);
   if( doWrite("mc_pu_sumpT_highpT") ) tree->Branch("mc_pu_sumpT_highpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_highpT, buffersize);
   if( doWrite("mc_pu_ntrks_lowpT") ) tree->Branch("mc_pu_ntrks_lowpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_lowpT, buffersize);
   if( doWrite("mc_pu_ntrks_highpT") ) tree->Branch("mc_pu_ntrks_highpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_highpT, buffersize);

   if( doWrite("trigger_n") ) tree->Branch("trigger_n", &trigger_n, "trigger_n/I", buffersize);
   if( doWrite("trigger") ) tree->Branch("trigger", "std::vector<int>", &trigger, buffersize);
   if( doWrite("trigger_name") ) tree->Branch("trigger_name", "std::vector<string>", &trigger_name, buffersize);
   if( doWrite("trigger_pass") ) tree->Branch("trigger_pass", "std::vector<bool>", &trigger_pass, buffersize);
   if( doWrite("trigger_prescale") ) tree->Branch("trigger_prescale", "std::vector<int>", &trigger_prescale, buffersize);
   if( doWrite("trigger_HLTprescale") ) tree->Branch("trigger_HLTprescale", "std::vector<int>", &trigger_HLTprescale, buffersize);
   if( doWrite("trigger_L1prescale") ) tree->Branch("trigger_L1prescale", "std::vector<int>", &trigger_L1prescale, buffersize);

   if( doWrite("triggerobject_n") ) tree->Branch("triggerobject_n", &triggerobject_n, "triggerobject_n/I",buffersize);
   if( doWrite("triggerobject_pt") ) tree->Branch("triggerobject_pt", "std::vector<float>", &triggerobject_pt, buffersize);
   if( doWrite("triggerobject_eta") ) tree->Branch("triggerobject_eta", "std::vector<float>", &triggerobject_eta, buffersize);
   if( doWrite("triggerobject_phi") ) tree->Branch("triggerobject_phi", "std::vector<float>", &triggerobject_phi, buffersize);

   if( doWrite("triggerobject_collection") ) tree->Branch("triggerobject_collection", "std::vector<std::string>", &triggerobject_collection, buffersize);

   if( doWrite("triggerobject_filterIds_n") ) tree->Branch("triggerobject_filterIds_n", "std::vector<int>", &triggerobject_filterIds_n, buffersize);
   if( doWrite("triggerobject_filterIds") ) tree->Branch("triggerobject_filterIds", "std::vector<int>", &triggerobject_filterIds, buffersize);

   if( doWrite("triggerobject_isTriggerL1Mu") ) tree->Branch("triggerobject_isTriggerL1Mu", "std::vector<bool>", &triggerobject_isTriggerL1Mu, buffersize);
   if( doWrite("triggerobject_isTriggerL1NoIsoEG") ) tree->Branch("triggerobject_isTriggerL1NoIsoEG", "std::vector<bool>", &triggerobject_isTriggerL1NoIsoEG, buffersize);
   if( doWrite("triggerobject_isTriggerL1IsoEG") ) tree->Branch("triggerobject_isTriggerL1IsoEG", "std::vector<bool>", &triggerobject_isTriggerL1IsoEG, buffersize);
   if( doWrite("triggerobject_isTriggerL1CenJet") ) tree->Branch("triggerobject_isTriggerL1CenJet", "std::vector<bool>", &triggerobject_isTriggerL1CenJet, buffersize);
   if( doWrite("triggerobject_isTriggerL1ForJet") ) tree->Branch("triggerobject_isTriggerL1ForJet", "std::vector<bool>", &triggerobject_isTriggerL1ForJet, buffersize);
   if( doWrite("triggerobject_isTriggerL1TauJet") ) tree->Branch("triggerobject_isTriggerL1TauJet", "std::vector<bool>", &triggerobject_isTriggerL1TauJet, buffersize);
   if( doWrite("triggerobject_isTriggerL1ETM") ) tree->Branch("triggerobject_isTriggerL1ETM", "std::vector<bool>", &triggerobject_isTriggerL1ETM, buffersize);
   if( doWrite("triggerobject_isTriggerL1ETT") ) tree->Branch("triggerobject_isTriggerL1ETT", "std::vector<bool>", &triggerobject_isTriggerL1ETT, buffersize);
   if( doWrite("triggerobject_isTriggerL1HTT") ) tree->Branch("triggerobject_isTriggerL1HTT", "std::vector<bool>", &triggerobject_isTriggerL1HTT, buffersize);
   if( doWrite("triggerobject_isTriggerL1HTM") ) tree->Branch("triggerobject_isTriggerL1HTM", "std::vector<bool>", &triggerobject_isTriggerL1HTM, buffersize);
   if( doWrite("triggerobject_isTriggerL1JetCounts") ) tree->Branch("triggerobject_isTriggerL1JetCounts", "std::vector<bool>", &triggerobject_isTriggerL1JetCounts, buffersize);
   if( doWrite("triggerobject_isTriggerL1HfBitCounts") ) tree->Branch("triggerobject_isTriggerL1HfBitCounts", "std::vector<bool>", &triggerobject_isTriggerL1HfBitCounts, buffersize);
   if( doWrite("triggerobject_isTriggerL1HfRingEtSums") ) tree->Branch("triggerobject_isTriggerL1HfRingEtSums", "std::vector<bool>", &triggerobject_isTriggerL1HfRingEtSums, buffersize);
   if( doWrite("triggerobject_isTriggerL1TechTrig") ) tree->Branch("triggerobject_isTriggerL1TechTrig", "std::vector<bool>", &triggerobject_isTriggerL1TechTrig, buffersize);
   if( doWrite("triggerobject_isTriggerL1Castor") ) tree->Branch("triggerobject_isTriggerL1Castor", "std::vector<bool>", &triggerobject_isTriggerL1Castor, buffersize);
   if( doWrite("triggerobject_isTriggerL1BPTX") ) tree->Branch("triggerobject_isTriggerL1BPTX", "std::vector<bool>", &triggerobject_isTriggerL1BPTX, buffersize);
   if( doWrite("triggerobject_isTriggerL1GtExternal") ) tree->Branch("triggerobject_isTriggerL1GtExternal", "std::vector<bool>", &triggerobject_isTriggerL1GtExternal, buffersize);

   if( doWrite("triggerobject_isHLT_TriggerPhoton") ) tree->Branch("triggerobject_isHLT_TriggerPhoton", "std::vector<bool>", &triggerobject_isHLT_TriggerPhoton, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerElectron") ) tree->Branch("triggerobject_isHLT_TriggerElectron", "std::vector<bool>", &triggerobject_isHLT_TriggerElectron, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMuon") ) tree->Branch("triggerobject_isHLT_TriggerMuon", "std::vector<bool>", &triggerobject_isHLT_TriggerMuon, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTau") ) tree->Branch("triggerobject_isHLT_TriggerTau", "std::vector<bool>", &triggerobject_isHLT_TriggerTau, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerJet") ) tree->Branch("triggerobject_isHLT_TriggerJet", "std::vector<bool>", &triggerobject_isHLT_TriggerJet, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerBJet") ) tree->Branch("triggerobject_isHLT_TriggerBJet", "std::vector<bool>", &triggerobject_isHLT_TriggerBJet, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMET") ) tree->Branch("triggerobject_isHLT_TriggerMET", "std::vector<bool>", &triggerobject_isHLT_TriggerMET, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTET") ) tree->Branch("triggerobject_isHLT_TriggerTET", "std::vector<bool>", &triggerobject_isHLT_TriggerTET, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTHT") ) tree->Branch("triggerobject_isHLT_TriggerTHT", "std::vector<bool>", &triggerobject_isHLT_TriggerTHT, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMHT") ) tree->Branch("triggerobject_isHLT_TriggerMHT", "std::vector<bool>", &triggerobject_isHLT_TriggerMHT, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerTrack") ) tree->Branch("triggerobject_isHLT_TriggerTrack", "std::vector<bool>", &triggerobject_isHLT_TriggerTrack, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerCluster") ) tree->Branch("triggerobject_isHLT_TriggerCluster", "std::vector<bool>", &triggerobject_isHLT_TriggerCluster, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMETSig") ) tree->Branch("triggerobject_isHLT_TriggerMETSig", "std::vector<bool>", &triggerobject_isHLT_TriggerMETSig, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerELongit") ) tree->Branch("triggerobject_isHLT_TriggerELongit", "std::vector<bool>", &triggerobject_isHLT_TriggerELongit, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerMHTSig") ) tree->Branch("triggerobject_isHLT_TriggerMHTSig", "std::vector<bool>", &triggerobject_isHLT_TriggerMHTSig, buffersize);
   if( doWrite("triggerobject_isHLT_TriggerHLongit") ) tree->Branch("triggerobject_isHLT_TriggerHLongit", "std::vector<bool>", &triggerobject_isHLT_TriggerHLongit, buffersize);

   if( doWrite("triggerobject_filterLabels_n") ) tree->Branch("triggerobject_filterLabels_n", "std::vector<int>", &triggerobject_filterLabels_n, buffersize);
   if( doWrite("triggerobject_filterLabels") ) tree->Branch("triggerobject_filterLabels", "std::vector<string>", &triggerobject_filterLabels, buffersize);

   if( doWrite("triggerobject_pathNamesAll_n") ) tree->Branch("triggerobject_pathNamesAll_n", "std::vector<int>", &triggerobject_pathNamesAll_n, buffersize);
   if( doWrite("triggerobject_pathNamesAll") ) tree->Branch("triggerobject_pathNamesAll", "std::vector<string>", &triggerobject_pathNamesAll, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isBoth") ) tree->Branch("triggerobject_pathNamesAll_isBoth", "std::vector<bool>", &triggerobject_pathNamesAll_isBoth, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isL3") ) tree->Branch("triggerobject_pathNamesAll_isL3", "std::vector<bool>", &triggerobject_pathNamesAll_isL3, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isLF") ) tree->Branch("triggerobject_pathNamesAll_isLF", "std::vector<bool>", &triggerobject_pathNamesAll_isLF, buffersize);
   if( doWrite("triggerobject_pathNamesAll_isNone") ) tree->Branch("triggerobject_pathNamesAll_isNone", "std::vector<bool>", &triggerobject_pathNamesAll_isNone, buffersize);

   if( doWrite("nvertex") ) tree->Branch("nvertex", &nvertex, "nvertex/I", buffersize);

   if( doWrite("el_n") ) tree->Branch("el_n", &el_n, "el_n/I", buffersize);
   if( doWrite("el_pt") ) tree->Branch("el_pt", "std::vector<float>", &el_pt, buffersize);
   if( doWrite("el_eta") ) tree->Branch("el_eta", "std::vector<float>", &el_eta, buffersize);
   if( doWrite("el_phi") ) tree->Branch("el_phi", "std::vector<float>", &el_phi, buffersize);
   if( doWrite("el_m") ) tree->Branch("el_m", "std::vector<float>", &el_m, buffersize);
   if( doWrite("el_E") ) tree->Branch("el_E", "std::vector<float>", &el_E, buffersize);
   if( doWrite("el_id") ) tree->Branch("el_id", "std::vector<int>", &el_id, buffersize);
   if( doWrite("el_charge") ) tree->Branch("el_charge", "std::vector<int>", &el_charge, buffersize);

   if( doWrite("el_passConversionVeto") ) tree->Branch("el_passConversionVeto", "std::vector<int>", &el_passConversionVeto, buffersize);
   if( doWrite("el_isGsfCtfScPixChargeConsistent") ) tree->Branch("el_isGsfCtfScPixChargeConsistent", "std::vector<int>", &el_isGsfCtfScPixChargeConsistent, buffersize);
   if( doWrite("el_isGsfScPixChargeConsistent") ) tree->Branch("el_isGsfScPixChargeConsistent", "std::vector<int>", &el_isGsfScPixChargeConsistent, buffersize);

   if( doWrite("el_ecalEnergy") ) tree->Branch("el_ecalEnergy", "std::vector<float>", &el_ecalEnergy, buffersize);
   if( doWrite("el_correctedEcalEnergy") ) tree->Branch("el_correctedEcalEnergy", "std::vector<float>", &el_correctedEcalEnergy, buffersize);
   if( doWrite("el_correctedEcalEnergyError") ) tree->Branch("el_correctedEcalEnergyError", "std::vector<float>", &el_correctedEcalEnergyError, buffersize);
   if( doWrite("el_trackMomentumError") ) tree->Branch("el_trackMomentumError", "std::vector<float>", &el_trackMomentumError, buffersize);

   if( doWrite("el_ip3d") ) tree->Branch("el_ip3d", "std::vector<float>", &el_ip3d, buffersize);
   if( doWrite("el_ip3dErr") ) tree->Branch("el_ip3dErr", "std::vector<float>", &el_ip3dErr, buffersize);
   if( doWrite("el_ip2d") ) tree->Branch("el_ip2d", "std::vector<float>", &el_ip2d, buffersize);
   if( doWrite("el_ip2dErr") ) tree->Branch("el_ip2dErr", "std::vector<float>", &el_ip2dErr, buffersize);
   if( doWrite("el_ip3dBS") ) tree->Branch("el_ip3dBS", "std::vector<float>", &el_ip3dBS, buffersize);
   if( doWrite("el_ip3dBSErr") ) tree->Branch("el_ip3dBSErr", "std::vector<float>", &el_ip3dBSErr, buffersize);
   if( doWrite("el_ip2dBS") ) tree->Branch("el_ip2dBS", "std::vector<float>", &el_ip2dBS, buffersize);
   if( doWrite("el_ip2dBSErr") ) tree->Branch("el_ip2dBSErr", "std::vector<float>", &el_ip2dBSErr, buffersize);

   if( doWrite("el_neutralHadronIso") ) tree->Branch("el_neutralHadronIso", "std::vector<float>", &el_neutralHadronIso, buffersize);
   if( doWrite("el_chargedHadronIso") ) tree->Branch("el_chargedHadronIso", "std::vector<float>", &el_chargedHadronIso, buffersize);
   if( doWrite("el_puChargedHadronIso") ) tree->Branch("el_puChargedHadronIso", "std::vector<float>", &el_puChargedHadronIso, buffersize);
   if( doWrite("el_ecalIso") ) tree->Branch("el_ecalIso", "std::vector<float>", &el_ecalIso, buffersize);
   if( doWrite("el_hcalIso") ) tree->Branch("el_hcalIso", "std::vector<float>", &el_hcalIso, buffersize);
   if( doWrite("el_particleIso") ) tree->Branch("el_particleIso", "std::vector<float>", &el_particleIso, buffersize);
   if( doWrite("el_photonIso") ) tree->Branch("el_photonIso", "std::vector<float>", &el_photonIso, buffersize);
   if( doWrite("el_trackIso") ) tree->Branch("el_trackIso", "std::vector<float>", &el_trackIso, buffersize);

   if( doWrite("el_ecalPFClusterIso") ) tree->Branch("el_ecalPFClusterIso", "std::vector<float>", &el_ecalPFClusterIso, buffersize);
   if( doWrite("el_hcalPFClusterIso") ) tree->Branch("el_hcalPFClusterIso", "std::vector<float>", &el_hcalPFClusterIso, buffersize);

   if( doWrite("el_dr03EcalRecHitSumEt") ) tree->Branch("el_dr03EcalRecHitSumEt", "std::vector<float>", &el_dr03EcalRecHitSumEt, buffersize);
   if( doWrite("el_dr03HcalTowerSumEt") ) tree->Branch("el_dr03HcalTowerSumEt", "std::vector<float>", &el_dr03HcalTowerSumEt, buffersize);
   if( doWrite("el_dr03HcalDepth1TowerSumEt") ) tree->Branch("el_dr03HcalDepth1TowerSumEt", "std::vector<float>", &el_dr03HcalDepth1TowerSumEt, buffersize);
   if( doWrite("el_dr03HcalDepth2TowerSumEt") ) tree->Branch("el_dr03HcalDepth2TowerSumEt", "std::vector<float>", &el_dr03HcalDepth2TowerSumEt, buffersize);
   if( doWrite("el_dr03TkSumPt") ) tree->Branch("el_dr03TkSumPt", "std::vector<float>", &el_dr03TkSumPt, buffersize);

   if( doWrite("el_dr04EcalRecHitSumEt") ) tree->Branch("el_dr04EcalRecHitSumEt", "std::vector<float>", &el_dr04EcalRecHitSumEt, buffersize);
   if( doWrite("el_dr04HcalTowerSumEt") ) tree->Branch("el_dr04HcalTowerSumEt", "std::vector<float>", &el_dr04HcalTowerSumEt, buffersize);
   if( doWrite("el_dr04HcalDepth1TowerSumEt") ) tree->Branch("el_dr04HcalDepth1TowerSumEt", "std::vector<float>", &el_dr04HcalDepth1TowerSumEt, buffersize);
   if( doWrite("el_dr04HcalDepth2TowerSumEt") ) tree->Branch("el_dr04HcalDepth2TowerSumEt", "std::vector<float>", &el_dr04HcalDepth2TowerSumEt, buffersize);
   if( doWrite("el_dr04TkSumPt") ) tree->Branch("el_dr04TkSumPt", "std::vector<float>", &el_dr04TkSumPt, buffersize);

   if( doWrite("el_hcalOverEcal") ) tree->Branch("el_hcalOverEcal", "std::vector<float>", &el_hcalOverEcal, buffersize);
   if( doWrite("el_hcalOverEcalBc") ) tree->Branch("el_hcalOverEcalBc", "std::vector<float>", &el_hcalOverEcalBc, buffersize);
   if( doWrite("el_hcalDepth1OverEcal") ) tree->Branch("el_hcalDepth1OverEcal", "std::vector<float>", &el_hcalDepth1OverEcal, buffersize);
   if( doWrite("el_hcalDepth2OverEcal") ) tree->Branch("el_hcalDepth2OverEcal", "std::vector<float>", &el_hcalDepth2OverEcal, buffersize);
   if( doWrite("el_eSeedClusterOverPout") ) tree->Branch("el_eSeedClusterOverPout", "std::vector<float>", &el_eSeedClusterOverPout, buffersize);
   if( doWrite("el_eSeedClusterOverP") ) tree->Branch("el_eSeedClusterOverP", "std::vector<float>", &el_eSeedClusterOverP, buffersize);
   if( doWrite("el_eEleClusterOverPout") ) tree->Branch("el_eEleClusterOverPout", "std::vector<float>", &el_eEleClusterOverPout, buffersize);
   if( doWrite("el_deltaEtaEleClusterTrackAtCalo") ) tree->Branch("el_deltaEtaEleClusterTrackAtCalo", "std::vector<float>", &el_deltaEtaEleClusterTrackAtCalo, buffersize);
   if( doWrite("el_deltaPhiEleClusterTrackAtCalo") ) tree->Branch("el_deltaPhiEleClusterTrackAtCalo", "std::vector<float>", &el_deltaPhiEleClusterTrackAtCalo, buffersize);

   if( doWrite("el_pfIso_sumChargedHadronPt") ) tree->Branch("el_pfIso_sumChargedHadronPt", "std::vector<float>", &el_pfIso_sumChargedHadronPt, buffersize);
   if( doWrite("el_pfIso_sumNeutralHadronEt") ) tree->Branch("el_pfIso_sumNeutralHadronEt", "std::vector<float>", &el_pfIso_sumNeutralHadronEt, buffersize);
   if( doWrite("el_pfIso_sumPhotonEt") ) tree->Branch("el_pfIso_sumPhotonEt", "std::vector<float>", &el_pfIso_sumPhotonEt, buffersize);
   if( doWrite("el_pfIso_sumPUPt") ) tree->Branch("el_pfIso_sumPUPt", "std::vector<float>", &el_pfIso_sumPUPt, buffersize);

   if( doWrite("el_vx") ) tree->Branch("el_vx", "std::vector<float>", &el_vx, buffersize);
   if( doWrite("el_vy") ) tree->Branch("el_vy", "std::vector<float>", &el_vy, buffersize);
   if( doWrite("el_vz") ) tree->Branch("el_vz", "std::vector<float>", &el_vz, buffersize);

   if( doWrite("el_hasGsfTrack") ) tree->Branch("el_hasGsfTrack", "std::vector<bool>", &el_hasGsfTrack, buffersize);
   if( doWrite("el_gsfTrack_d0") ) tree->Branch("el_gsfTrack_d0", "std::vector<float>", &el_gsfTrack_d0, buffersize);
   if( doWrite("el_gsfTrack_z0") ) tree->Branch("el_gsfTrack_z0", "std::vector<float>", &el_gsfTrack_z0, buffersize);
   if( doWrite("el_gsfTrack_d0Error") ) tree->Branch("el_gsfTrack_d0Error", "std::vector<float>", &el_gsfTrack_d0Error, buffersize);
   if( doWrite("el_gsfTrack_z0Error") ) tree->Branch("el_gsfTrack_z0Error", "std::vector<float>", &el_gsfTrack_z0Error, buffersize);
   if( doWrite("el_gsfTrack_PV_dxy") ) tree->Branch("el_gsfTrack_PV_dxy", "std::vector<float>", &el_gsfTrack_PV_dxy, buffersize);
   if( doWrite("el_gsfTrack_PV_dz") ) tree->Branch("el_gsfTrack_PV_dz", "std::vector<float>", &el_gsfTrack_PV_dz, buffersize);
   if( doWrite("el_gsfTrack_RP_dxy") ) tree->Branch("el_gsfTrack_RP_dxy", "std::vector<float>", &el_gsfTrack_RP_dxy, buffersize);
   if( doWrite("el_gsfTrack_RP_dz") ) tree->Branch("el_gsfTrack_RP_dz", "std::vector<float>", &el_gsfTrack_RP_dz, buffersize);
   if( doWrite("el_gsfTrack_BS_dxy") ) tree->Branch("el_gsfTrack_BS_dxy", "std::vector<float>", &el_gsfTrack_BS_dxy, buffersize);
   if( doWrite("el_gsfTrack_BS_dz") ) tree->Branch("el_gsfTrack_BS_dz", "std::vector<float>", &el_gsfTrack_BS_dz, buffersize);
   if( doWrite("el_gsfTrack_dxyError") ) tree->Branch("el_gsfTrack_dxyError", "std::vector<float>", &el_gsfTrack_dxyError, buffersize);
   if( doWrite("el_gsfTrack_dzError") ) tree->Branch("el_gsfTrack_dzError", "std::vector<float>", &el_gsfTrack_dzError, buffersize);
   if( doWrite("el_gsfTrack_normalizedChi2") ) tree->Branch("el_gsfTrack_normalizedChi2", "std::vector<float>", &el_gsfTrack_normalizedChi2, buffersize);

   if( doWrite("el_superCluster_eta") ) tree->Branch("el_superCluster_eta", "std::vector<float>", &el_superCluster_eta, buffersize);
   if( doWrite("el_superCluster_phi") ) tree->Branch("el_superCluster_phi", "std::vector<float>", &el_superCluster_phi, buffersize);
   if( doWrite("el_superCluster_energy") ) tree->Branch("el_superCluster_energy", "std::vector<float>", &el_superCluster_energy, buffersize);
   if( doWrite("el_superCluster_rawEnergy") ) tree->Branch("el_superCluster_rawEnergy", "std::vector<float>", &el_superCluster_rawEnergy, buffersize);
   if( doWrite("el_superCluster_preshowerEnergy") ) tree->Branch("el_superCluster_preshowerEnergy", "std::vector<float>", &el_superCluster_preshowerEnergy, buffersize);
   if( doWrite("el_superCluster_etaWidth") ) tree->Branch("el_superCluster_etaWidth", "std::vector<float>", &el_superCluster_etaWidth, buffersize);
   if( doWrite("el_superCluster_phiWidth") ) tree->Branch("el_superCluster_phiWidth", "std::vector<float>", &el_superCluster_phiWidth, buffersize);
   if( doWrite("el_superCluster_preshowerEnergyPlane1") ) tree->Branch("el_superCluster_preshowerEnergyPlane1", "std::vector<float>", &el_superCluster_preshowerEnergyPlane1, buffersize);
   if( doWrite("el_superCluster_preshowerEnergyPlane2") ) tree->Branch("el_superCluster_preshowerEnergyPlane2", "std::vector<float>", &el_superCluster_preshowerEnergyPlane2, buffersize);
   if( doWrite("el_superCluster_positionR") ) tree->Branch("el_superCluster_positionR", "std::vector<float>", &el_superCluster_positionR, buffersize);

   if( doWrite("el_basicClustersSize") ) tree->Branch("el_basicClustersSize", "std::vector<int>", &el_basicClustersSize, buffersize);
   if( doWrite("el_e1x5") ) tree->Branch("el_e1x5", "std::vector<float>", &el_e1x5, buffersize);
   if( doWrite("el_e5x5") ) tree->Branch("el_e5x5", "std::vector<float>", &el_e5x5, buffersize);
   if( doWrite("el_e2x5Max") ) tree->Branch("el_e2x5Max", "std::vector<float>", &el_e2x5Max, buffersize);
   if( doWrite("el_sigmaEtaEta") ) tree->Branch("el_sigmaEtaEta", "std::vector<float>", &el_sigmaEtaEta, buffersize);
   if( doWrite("el_sigmaIetaIeta") ) tree->Branch("el_sigmaIetaIeta", "std::vector<float>", &el_sigmaIetaIeta, buffersize);
   if( doWrite("el_sigmaIphiIphi") ) tree->Branch("el_sigmaIphiIphi", "std::vector<float>", &el_sigmaIphiIphi, buffersize);
   if( doWrite("el_sigmaIetaIphi") ) tree->Branch("el_sigmaIetaIphi", "std::vector<float>", &el_sigmaIetaIphi, buffersize);
   if( doWrite("el_full5x5_sigmaIphiIphi") ) tree->Branch("el_full5x5_sigmaIphiIphi", "std::vector<float>", &el_full5x5_sigmaIphiIphi, buffersize);
   if( doWrite("el_full5x5_sigmaEtaEta") ) tree->Branch("el_full5x5_sigmaEtaEta", "std::vector<float>", &el_full5x5_sigmaEtaEta, buffersize);
   if( doWrite("el_full5x5_sigmaIetaIeta") ) tree->Branch("el_full5x5_sigmaIetaIeta", "std::vector<float>", &el_full5x5_sigmaIetaIeta, buffersize);
   if( doWrite("el_full5x5_sigmaIetaIphi") ) tree->Branch("el_full5x5_sigmaIetaIphi", "std::vector<float>", &el_full5x5_sigmaIetaIphi, buffersize);
   if( doWrite("el_full5x5_r9") ) tree->Branch("el_full5x5_r9", "std::vector<float>", &el_full5x5_r9, buffersize);
   if( doWrite("el_full5x5_e1x5") ) tree->Branch("el_full5x5_e1x5", "std::vector<float>", &el_full5x5_e1x5, buffersize);
   if( doWrite("el_full5x5_e5x5") ) tree->Branch("el_full5x5_e5x5", "std::vector<float>", &el_full5x5_e5x5, buffersize);
   if( doWrite("el_full5x5_e2x5Max") ) tree->Branch("el_full5x5_e2x5Max", "std::vector<float>", &el_full5x5_e2x5Max, buffersize);

   if( doWrite("el_expectedMissingOuterHits") ) tree->Branch("el_expectedMissingOuterHits", "std::vector<int>", &el_expectedMissingOuterHits, buffersize);
   if( doWrite("el_numberOfValidPixelHits") ) tree->Branch("el_numberOfValidPixelHits", "std::vector<int>", &el_numberOfValidPixelHits, buffersize);
   if( doWrite("el_numberOfLostPixelHits") ) tree->Branch("el_numberOfLostPixelHits", "std::vector<int>", &el_numberOfLostPixelHits, buffersize);
   if( doWrite("el_trackerLayersWithMeasurement") ) tree->Branch("el_trackerLayersWithMeasurement", "std::vector<int>", &el_trackerLayersWithMeasurement, buffersize);
   if( doWrite("el_pixelLayersWithMeasurement") ) tree->Branch("el_pixelLayersWithMeasurement", "std::vector<int>", &el_pixelLayersWithMeasurement, buffersize);
   if( doWrite("el_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("el_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &el_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("el_trackerLayersWithoutMeasurement") ) tree->Branch("el_trackerLayersWithoutMeasurement", "std::vector<int>", &el_trackerLayersWithoutMeasurement, buffersize);

   if( doWrite("el_numberOfHits") ) tree->Branch("el_numberOfHits", "std::vector<int>", &el_numberOfHits, buffersize);
   if( doWrite("el_numberOfValidHits") ) tree->Branch("el_numberOfValidHits", "std::vector<int>", &el_numberOfValidHits, buffersize);

   if( doWrite("el_hadronicOverEm") ) tree->Branch("el_hadronicOverEm", "std::vector<float>", &el_hadronicOverEm, buffersize);
   if( doWrite("el_numberOfLostHits") ) tree->Branch("el_numberOfLostHits", "std::vector<int>", &el_numberOfLostHits, buffersize);
   if( doWrite("el_numberOfLostHitsDefault") ) tree->Branch("el_numberOfLostHitsDefault", "std::vector<int>", &el_numberOfLostHitsDefault, buffersize);

   if( doWrite("el_fbrem") ) tree->Branch("el_fbrem", "std::vector<float>", &el_fbrem, buffersize);
   if( doWrite("el_kf_normalizedChi2") ) tree->Branch("el_kf_normalizedChi2", "std::vector<float>", &el_kf_normalizedChi2, buffersize);
   if( doWrite("el_gsf_normalizedChi2") ) tree->Branch("el_gsf_normalizedChi2", "std::vector<float>", &el_gsf_normalizedChi2, buffersize);
   if( doWrite("el_deltaEtaSuperClusterTrackAtVtx") ) tree->Branch("el_deltaEtaSuperClusterTrackAtVtx", "std::vector<float>", &el_deltaEtaSuperClusterTrackAtVtx, buffersize);
   if( doWrite("el_deltaPhiSuperClusterTrackAtVtx") ) tree->Branch("el_deltaPhiSuperClusterTrackAtVtx", "std::vector<float>", &el_deltaPhiSuperClusterTrackAtVtx, buffersize);
   if( doWrite("el_deltaEtaSeedClusterTrackAtCalo") ) tree->Branch("el_deltaEtaSeedClusterTrackAtCalo", "std::vector<float>", &el_deltaEtaSeedClusterTrackAtCalo, buffersize);
   if( doWrite("el_deltaPhiSeedClusterTrackAtCalo") ) tree->Branch("el_deltaPhiSeedClusterTrackAtCalo", "std::vector<float>", &el_deltaPhiSeedClusterTrackAtCalo, buffersize);
   if( doWrite("el_full5x5_OneMinusE1x5E5x5") ) tree->Branch("el_full5x5_OneMinusE1x5E5x5", "std::vector<float>", &el_full5x5_OneMinusE1x5E5x5, buffersize);
   if( doWrite("el_OneMinusE1x5E5x5") ) tree->Branch("el_OneMinusE1x5E5x5", "std::vector<float>", &el_OneMinusE1x5E5x5, buffersize);
   if( doWrite("el_eSuperClusterOverP") ) tree->Branch("el_eSuperClusterOverP", "std::vector<float>", &el_eSuperClusterOverP, buffersize);
   if( doWrite("el_IoEmIoP") ) tree->Branch("el_IoEmIoP", "std::vector<float>", &el_IoEmIoP, buffersize);
   if( doWrite("el_ooEmooP") ) tree->Branch("el_ooEmooP", "std::vector<float>", &el_ooEmooP, buffersize);
   if( doWrite("el_eleEoPout") ) tree->Branch("el_eleEoPout", "std::vector<float>", &el_eleEoPout, buffersize);
   if( doWrite("el_PreShowerOverRaw") ) tree->Branch("el_PreShowerOverRaw", "std::vector<float>", &el_PreShowerOverRaw, buffersize);

   if( doWrite("el_mvaNonTrigV0") ) tree->Branch("el_mvaNonTrigV0", "std::vector<float>", &el_mvaNonTrigV0, buffersize);
   if( doWrite("el_mvaNonTrigCat") ) tree->Branch("el_mvaNonTrigCat", "std::vector<float>", &el_mvaNonTrigCat, buffersize);

   if( doWrite("el_vetoCBId") ) tree->Branch("el_vetoCBId", "std::vector<bool>", &el_vetoCBId, buffersize);
   if( doWrite("el_looseCBId") ) tree->Branch("el_looseCBId", "std::vector<bool>", &el_looseCBId, buffersize);
   if( doWrite("el_mediumCBId") ) tree->Branch("el_mediumCBId", "std::vector<bool>", &el_mediumCBId, buffersize);
   if( doWrite("el_tightCBId") ) tree->Branch("el_tightCBId", "std::vector<bool>", &el_tightCBId, buffersize);
   if( doWrite("el_heepCBId") ) tree->Branch("el_heepCBId", "std::vector<bool>", &el_heepCBId, buffersize);
   if( doWrite("el_mediumMVAId") ) tree->Branch("el_mediumMVAId", "std::vector<bool>", &el_mediumMVAId, buffersize);
   if( doWrite("el_tightMVAId") ) tree->Branch("el_tightMVAId", "std::vector<bool>", &el_tightMVAId, buffersize);

   if( doWrite("el_hasMCMatch") ) tree->Branch("el_hasMCMatch", "std::vector<int>", &el_hasMCMatch, buffersize);
   if( doWrite("el_gen_pt") ) tree->Branch("el_gen_pt", "std::vector<float>", &el_gen_pt, buffersize);
   if( doWrite("el_gen_eta") ) tree->Branch("el_gen_eta", "std::vector<float>", &el_gen_eta, buffersize);
   if( doWrite("el_gen_phi") ) tree->Branch("el_gen_phi", "std::vector<float>", &el_gen_phi, buffersize);
   if( doWrite("el_gen_m") ) tree->Branch("el_gen_m", "std::vector<float>", &el_gen_m, buffersize);
   if( doWrite("el_gen_status") ) tree->Branch("el_gen_status", "std::vector<int>", &el_gen_status, buffersize);
   if( doWrite("el_gen_id") ) tree->Branch("el_gen_id", "std::vector<int>", &el_gen_id, buffersize);
   if( doWrite("el_gen_charge") ) tree->Branch("el_gen_charge", "std::vector<int>", &el_gen_charge, buffersize);
   if( doWrite("el_gen_dr") ) tree->Branch("el_gen_dr", "std::vector<float>", &el_gen_dr, buffersize);

   if( doWrite("el_hasMCMatchPAT") ) tree->Branch("el_hasMCMatchPAT", "std::vector<int>", &el_hasMCMatchPAT, buffersize);
   if( doWrite("el_genPAT_pt") ) tree->Branch("el_genPAT_pt", "std::vector<float>", &el_genPAT_pt, buffersize);
   if( doWrite("el_genPAT_eta") ) tree->Branch("el_genPAT_eta", "std::vector<float>", &el_genPAT_eta, buffersize);
   if( doWrite("el_genPAT_phi") ) tree->Branch("el_genPAT_phi", "std::vector<float>", &el_genPAT_phi, buffersize);
   if( doWrite("el_genPAT_m") ) tree->Branch("el_genPAT_m", "std::vector<float>", &el_genPAT_m, buffersize);
   if( doWrite("el_genPAT_status") ) tree->Branch("el_genPAT_status", "std::vector<int>", &el_genPAT_status, buffersize);
   if( doWrite("el_genPAT_id") ) tree->Branch("el_genPAT_id", "std::vector<int>", &el_genPAT_id, buffersize);
   if( doWrite("el_genPAT_charge") ) tree->Branch("el_genPAT_charge", "std::vector<int>", &el_genPAT_charge, buffersize);

   if( doWrite("el_hasMatchedConversion") ) tree->Branch("el_hasMatchedConversion", "std::vector<bool>", &el_hasMatchedConversion, buffersize);
   if( doWrite("el_expectedMissingInnerHits") ) tree->Branch("el_expectedMissingInnerHits", "std::vector<int>", &el_expectedMissingInnerHits, buffersize);

   if( doWrite("mu_n") ) tree->Branch("mu_n", &mu_n, "mu_n/I", buffersize);
   if( doWrite("mu_pt") ) tree->Branch("mu_pt", "std::vector<float>", &mu_pt, buffersize);
   if( doWrite("mu_eta") ) tree->Branch("mu_eta", "std::vector<float>", &mu_eta, buffersize);
   if( doWrite("mu_phi") ) tree->Branch("mu_phi", "std::vector<float>", &mu_phi, buffersize);
   if( doWrite("mu_m") ) tree->Branch("mu_m", "std::vector<float>", &mu_m, buffersize);
   if( doWrite("mu_E") ) tree->Branch("mu_E", "std::vector<float>", &mu_E, buffersize);
   if( doWrite("mu_id") ) tree->Branch("mu_id", "std::vector<int>", &mu_id, buffersize);
   if( doWrite("mu_charge") ) tree->Branch("mu_charge", "std::vector<int>", &mu_charge, buffersize);

   if( doWrite("mu_ip3d") ) tree->Branch("mu_ip3d", "std::vector<float>", &mu_ip3d, buffersize);
   if( doWrite("mu_ip3dErr") ) tree->Branch("mu_ip3dErr", "std::vector<float>", &mu_ip3dErr, buffersize);
   if( doWrite("mu_ip2d") ) tree->Branch("mu_ip2d", "std::vector<float>", &mu_ip2d, buffersize);
   if( doWrite("mu_ip2dErr") ) tree->Branch("mu_ip2dErr", "std::vector<float>", &mu_ip2dErr, buffersize);
   if( doWrite("mu_ip3dBS") ) tree->Branch("mu_ip3dBS", "std::vector<float>", &mu_ip3dBS, buffersize);
   if( doWrite("mu_ip3dBSErr") ) tree->Branch("mu_ip3dBSErr", "std::vector<float>", &mu_ip3dBSErr, buffersize);
   if( doWrite("mu_ip2dBS") ) tree->Branch("mu_ip2dBS", "std::vector<float>", &mu_ip2dBS, buffersize);
   if( doWrite("mu_ip2dBSErr") ) tree->Branch("mu_ip2dBSErr", "std::vector<float>", &mu_ip2dBSErr, buffersize);

   if( doWrite("mu_neutralHadronIso") ) tree->Branch("mu_neutralHadronIso", "std::vector<float>", &mu_neutralHadronIso, buffersize);
   if( doWrite("mu_chargedHadronIso") ) tree->Branch("mu_chargedHadronIso", "std::vector<float>", &mu_chargedHadronIso, buffersize);
   if( doWrite("mu_puChargedHadronIso") ) tree->Branch("mu_puChargedHadronIso", "std::vector<float>", &mu_puChargedHadronIso, buffersize);
   if( doWrite("mu_ecalIso") ) tree->Branch("mu_ecalIso", "std::vector<float>", &mu_ecalIso, buffersize);
   if( doWrite("mu_hcalIso") ) tree->Branch("mu_hcalIso", "std::vector<float>", &mu_hcalIso, buffersize);
   if( doWrite("mu_photonIso") ) tree->Branch("mu_photonIso", "std::vector<float>", &mu_photonIso, buffersize);
   if( doWrite("mu_trackIso") ) tree->Branch("mu_trackIso", "std::vector<float>", &mu_trackIso, buffersize);

   if( doWrite("mu_pfIso03_sumChargedHadronPt") ) tree->Branch("mu_pfIso03_sumChargedHadronPt", "std::vector<float>", &mu_pfIso03_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfIso03_sumChargedParticlePt") ) tree->Branch("mu_pfIso03_sumChargedParticlePt", "std::vector<float>", &mu_pfIso03_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfIso03_sumNeutralHadronEt") ) tree->Branch("mu_pfIso03_sumNeutralHadronEt", "std::vector<float>", &mu_pfIso03_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfIso03_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfIso03_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfIso03_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso03_sumPhotonEt") ) tree->Branch("mu_pfIso03_sumPhotonEt", "std::vector<float>", &mu_pfIso03_sumPhotonEt, buffersize);
   if( doWrite("mu_pfIso03_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfIso03_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfIso03_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso03_sumPUPt") ) tree->Branch("mu_pfIso03_sumPUPt", "std::vector<float>", &mu_pfIso03_sumPUPt, buffersize);

   if( doWrite("mu_pfIso04_sumChargedHadronPt") ) tree->Branch("mu_pfIso04_sumChargedHadronPt", "std::vector<float>", &mu_pfIso04_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfIso04_sumChargedParticlePt") ) tree->Branch("mu_pfIso04_sumChargedParticlePt", "std::vector<float>", &mu_pfIso04_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfIso04_sumNeutralHadronEt") ) tree->Branch("mu_pfIso04_sumNeutralHadronEt", "std::vector<float>", &mu_pfIso04_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfIso04_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfIso04_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfIso04_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso04_sumPhotonEt") ) tree->Branch("mu_pfIso04_sumPhotonEt", "std::vector<float>", &mu_pfIso04_sumPhotonEt, buffersize);
   if( doWrite("mu_pfIso04_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfIso04_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfIso04_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso04_sumPUPt") ) tree->Branch("mu_pfIso04_sumPUPt", "std::vector<float>", &mu_pfIso04_sumPUPt, buffersize);

   if( doWrite("mu_pfMeanIso03_sumChargedHadronPt") ) tree->Branch("mu_pfMeanIso03_sumChargedHadronPt", "std::vector<float>", &mu_pfMeanIso03_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumChargedParticlePt") ) tree->Branch("mu_pfMeanIso03_sumChargedParticlePt", "std::vector<float>", &mu_pfMeanIso03_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumNeutralHadronEt") ) tree->Branch("mu_pfMeanIso03_sumNeutralHadronEt", "std::vector<float>", &mu_pfMeanIso03_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfMeanIso03_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfMeanIso03_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso03_sumPhotonEt") ) tree->Branch("mu_pfMeanIso03_sumPhotonEt", "std::vector<float>", &mu_pfMeanIso03_sumPhotonEt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfMeanIso03_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfMeanIso03_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso03_sumPUPt") ) tree->Branch("mu_pfMeanIso03_sumPUPt", "std::vector<float>", &mu_pfMeanIso03_sumPUPt, buffersize);

   if( doWrite("mu_pfSumIso03_sumChargedHadronPt") ) tree->Branch("mu_pfSumIso03_sumChargedHadronPt", "std::vector<float>", &mu_pfSumIso03_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfSumIso03_sumChargedParticlePt") ) tree->Branch("mu_pfSumIso03_sumChargedParticlePt", "std::vector<float>", &mu_pfSumIso03_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfSumIso03_sumNeutralHadronEt") ) tree->Branch("mu_pfSumIso03_sumNeutralHadronEt", "std::vector<float>", &mu_pfSumIso03_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfSumIso03_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfSumIso03_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfSumIso03_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso03_sumPhotonEt") ) tree->Branch("mu_pfSumIso03_sumPhotonEt", "std::vector<float>", &mu_pfSumIso03_sumPhotonEt, buffersize);
   if( doWrite("mu_pfSumIso03_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfSumIso03_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfSumIso03_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso03_sumPUPt") ) tree->Branch("mu_pfSumIso03_sumPUPt", "std::vector<float>", &mu_pfSumIso03_sumPUPt, buffersize);

   if( doWrite("mu_pfMeanIso04_sumChargedHadronPt") ) tree->Branch("mu_pfMeanIso04_sumChargedHadronPt", "std::vector<float>", &mu_pfMeanIso04_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumChargedParticlePt") ) tree->Branch("mu_pfMeanIso04_sumChargedParticlePt", "std::vector<float>", &mu_pfMeanIso04_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumNeutralHadronEt") ) tree->Branch("mu_pfMeanIso04_sumNeutralHadronEt", "std::vector<float>", &mu_pfMeanIso04_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfMeanIso04_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfMeanIso04_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso04_sumPhotonEt") ) tree->Branch("mu_pfMeanIso04_sumPhotonEt", "std::vector<float>", &mu_pfMeanIso04_sumPhotonEt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfMeanIso04_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfMeanIso04_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso04_sumPUPt") ) tree->Branch("mu_pfMeanIso04_sumPUPt", "std::vector<float>", &mu_pfMeanIso04_sumPUPt, buffersize);

   if( doWrite("mu_pfSumIso04_sumChargedHadronPt") ) tree->Branch("mu_pfSumIso04_sumChargedHadronPt", "std::vector<float>", &mu_pfSumIso04_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfSumIso04_sumChargedParticlePt") ) tree->Branch("mu_pfSumIso04_sumChargedParticlePt", "std::vector<float>", &mu_pfSumIso04_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfSumIso04_sumNeutralHadronEt") ) tree->Branch("mu_pfSumIso04_sumNeutralHadronEt", "std::vector<float>", &mu_pfSumIso04_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfSumIso04_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfSumIso04_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfSumIso04_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso04_sumPhotonEt") ) tree->Branch("mu_pfSumIso04_sumPhotonEt", "std::vector<float>", &mu_pfSumIso04_sumPhotonEt, buffersize);
   if( doWrite("mu_pfSumIso04_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfSumIso04_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfSumIso04_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso04_sumPUPt") ) tree->Branch("mu_pfSumIso04_sumPUPt", "std::vector<float>", &mu_pfSumIso04_sumPUPt, buffersize);

   if( doWrite("mu_isGlobalMuon") ) tree->Branch("mu_isGlobalMuon", "std::vector<int>", &mu_isGlobalMuon, buffersize);
   if( doWrite("mu_isTrackerMuon") ) tree->Branch("mu_isTrackerMuon", "std::vector<int>", &mu_isTrackerMuon, buffersize);
   if( doWrite("mu_isStandAloneMuon") ) tree->Branch("mu_isStandAloneMuon", "std::vector<int>", &mu_isStandAloneMuon, buffersize);
   if( doWrite("mu_isCaloMuon") ) tree->Branch("mu_isCaloMuon", "std::vector<int>", &mu_isCaloMuon, buffersize);
   if( doWrite("mu_isPFMuon") ) tree->Branch("mu_isPFMuon", "std::vector<int>", &mu_isPFMuon, buffersize);
   if( doWrite("mu_isRPCMuon") ) tree->Branch("mu_isRPCMuon", "std::vector<int>", &mu_isRPCMuon, buffersize);

   if( doWrite("mu_vx") ) tree->Branch("mu_vx", "std::vector<float>", &mu_vx, buffersize);
   if( doWrite("mu_vy") ) tree->Branch("mu_vy", "std::vector<float>", &mu_vy, buffersize);
   if( doWrite("mu_vz") ) tree->Branch("mu_vz", "std::vector<float>", &mu_vz, buffersize);

   if( doWrite("mu_numberOfMatches") ) tree->Branch("mu_numberOfMatches", "std::vector<int>", &mu_numberOfMatches, buffersize);
   if( doWrite("mu_numberOfMatchedStations") ) tree->Branch("mu_numberOfMatchedStations", "std::vector<int>", &mu_numberOfMatchedStations, buffersize);

   if( doWrite("mu_segmentCompatibility") ) tree->Branch("mu_segmentCompatibility", "std::vector<float>", &mu_segmentCompatibility, buffersize);
   if( doWrite("mu_caloCompatibility") ) tree->Branch("mu_caloCompatibility", "std::vector<float>", &mu_caloCompatibility, buffersize);

   if( doWrite("mu_combinedQuality_chi2LocalPosition") ) tree->Branch("mu_combinedQuality_chi2LocalPosition", "std::vector<float>", &mu_combinedQuality_chi2LocalPosition, buffersize);
   if( doWrite("mu_combinedQuality_trkKink") ) tree->Branch("mu_combinedQuality_trkKink", "std::vector<float>", &mu_combinedQuality_trkKink, buffersize);

   if( doWrite("mu_isLooseMuon") ) tree->Branch("mu_isLooseMuon", "std::vector<bool>", &mu_isLooseMuon, buffersize);
   if( doWrite("mu_isMediumMuon") ) tree->Branch("mu_isMediumMuon", "std::vector<bool>", &mu_isMediumMuon, buffersize);
   if( doWrite("mu_isTightMuon") ) tree->Branch("mu_isTightMuon", "std::vector<bool>", &mu_isTightMuon, buffersize);
   if( doWrite("mu_isSoftMuon") ) tree->Branch("mu_isSoftMuon", "std::vector<bool>", &mu_isSoftMuon, buffersize);
   if( doWrite("mu_isHighPtMuon") ) tree->Branch("mu_isHighPtMuon", "std::vector<bool>", &mu_isHighPtMuon, buffersize);

   if( doWrite("mu_isGoodMuon_AllGlobalMuons") ) tree->Branch("mu_isGoodMuon_AllGlobalMuons", "std::vector<bool>", &mu_isGoodMuon_AllGlobalMuons, buffersize);
   if( doWrite("mu_isGoodMuon_AllStandAloneMuons") ) tree->Branch("mu_isGoodMuon_AllStandAloneMuons", "std::vector<bool>", &mu_isGoodMuon_AllStandAloneMuons, buffersize);
   if( doWrite("mu_isGoodMuon_AllTrackerMuons") ) tree->Branch("mu_isGoodMuon_AllTrackerMuons", "std::vector<bool>", &mu_isGoodMuon_AllTrackerMuons, buffersize);
   if( doWrite("mu_isGoodMuon_TrackerMuonArbitrated") ) tree->Branch("mu_isGoodMuon_TrackerMuonArbitrated", "std::vector<bool>", &mu_isGoodMuon_TrackerMuonArbitrated, buffersize);
   if( doWrite("mu_isGoodMuon_AllArbitrated") ) tree->Branch("mu_isGoodMuon_AllArbitrated", "std::vector<bool>", &mu_isGoodMuon_AllArbitrated, buffersize);
   if( doWrite("mu_isGoodMuon_GlobalMuonPromptTight") ) tree->Branch("mu_isGoodMuon_GlobalMuonPromptTight", "std::vector<bool>", &mu_isGoodMuon_GlobalMuonPromptTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationTight") ) tree->Branch("mu_isGoodMuon_TMLastStationTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationTight, buffersize);
   if( doWrite("mu_isGoodMuon_TM2DCompatibilityLoose") ) tree->Branch("mu_isGoodMuon_TM2DCompatibilityLoose", "std::vector<bool>", &mu_isGoodMuon_TM2DCompatibilityLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TM2DCompatibilityTight") ) tree->Branch("mu_isGoodMuon_TM2DCompatibilityTight", "std::vector<bool>", &mu_isGoodMuon_TM2DCompatibilityTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationLoose") ) tree->Branch("mu_isGoodMuon_TMOneStationLoose", "std::vector<bool>", &mu_isGoodMuon_TMOneStationLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationTight") ) tree->Branch("mu_isGoodMuon_TMOneStationTight", "std::vector<bool>", &mu_isGoodMuon_TMOneStationTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedLowPtLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedLowPtLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedLowPtLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedLowPtTight") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedLowPtTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedLowPtTight, buffersize);
   if( doWrite("mu_isGoodMuon_GMTkChiCompatibility") ) tree->Branch("mu_isGoodMuon_GMTkChiCompatibility", "std::vector<bool>", &mu_isGoodMuon_GMTkChiCompatibility, buffersize);
   if( doWrite("mu_isGoodMuon_GMStaChiCompatibility") ) tree->Branch("mu_isGoodMuon_GMStaChiCompatibility", "std::vector<bool>", &mu_isGoodMuon_GMStaChiCompatibility, buffersize);
   if( doWrite("mu_isGoodMuon_GMTkKinkTight") ) tree->Branch("mu_isGoodMuon_GMTkKinkTight", "std::vector<bool>", &mu_isGoodMuon_GMTkKinkTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationAngLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationAngLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationAngLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationAngTight") ) tree->Branch("mu_isGoodMuon_TMLastStationAngTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationAngTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationAngLoose") ) tree->Branch("mu_isGoodMuon_TMOneStationAngLoose", "std::vector<bool>", &mu_isGoodMuon_TMOneStationAngLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationAngTight") ) tree->Branch("mu_isGoodMuon_TMOneStationAngTight", "std::vector<bool>", &mu_isGoodMuon_TMOneStationAngTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight, buffersize);

   if( doWrite("mu_calEnergy_em") ) tree->Branch("mu_calEnergy_em", "std::vector<float>", &mu_calEnergy_em, buffersize);
   if( doWrite("mu_calEnergy_had") ) tree->Branch("mu_calEnergy_had", "std::vector<float>", &mu_calEnergy_had, buffersize);
   if( doWrite("mu_calEnergy_ho") ) tree->Branch("mu_calEnergy_ho", "std::vector<float>", &mu_calEnergy_ho, buffersize);
   if( doWrite("mu_calEnergy_emS9") ) tree->Branch("mu_calEnergy_emS9", "std::vector<float>", &mu_calEnergy_emS9, buffersize);
   if( doWrite("mu_calEnergy_hadS9") ) tree->Branch("mu_calEnergy_hadS9", "std::vector<float>", &mu_calEnergy_hadS9, buffersize);
   if( doWrite("mu_calEnergy_hoS9") ) tree->Branch("mu_calEnergy_hoS9", "std::vector<float>", &mu_calEnergy_hoS9, buffersize);
   if( doWrite("mu_calEnergy_emS25") ) tree->Branch("mu_calEnergy_emS25", "std::vector<float>", &mu_calEnergy_emS25, buffersize);
   if( doWrite("mu_calEnergy_emMax") ) tree->Branch("mu_calEnergy_emMax", "std::vector<float>", &mu_calEnergy_emMax, buffersize);
   if( doWrite("mu_calEnergy_hadMax") ) tree->Branch("mu_calEnergy_hadMax", "std::vector<float>", &mu_calEnergy_hadMax, buffersize);
   if( doWrite("mu_calEnergy_ecal_time") ) tree->Branch("mu_calEnergy_ecal_time", "std::vector<float>", &mu_calEnergy_ecal_time, buffersize);
   if( doWrite("mu_calEnergy_hcal_time") ) tree->Branch("mu_calEnergy_hcal_time", "std::vector<float>", &mu_calEnergy_hcal_time, buffersize);
   if( doWrite("mu_calEnergy_ecal_rawId") ) tree->Branch("mu_calEnergy_ecal_rawId", "std::vector<int>", &mu_calEnergy_ecal_rawId, buffersize);
   if( doWrite("mu_calEnergy_hcal_rawId") ) tree->Branch("mu_calEnergy_hcal_rawId", "std::vector<int>", &mu_calEnergy_hcal_rawId, buffersize);

   if( doWrite("mu_isolationR03_trackerVetoPt") ) tree->Branch("mu_isolationR03_trackerVetoPt", "std::vector<float>", &mu_isolationR03_trackerVetoPt, buffersize);
   if( doWrite("mu_isolationR03_emVetoEt") ) tree->Branch("mu_isolationR03_emVetoEt", "std::vector<float>", &mu_isolationR03_emVetoEt, buffersize);
   if( doWrite("mu_isolationR03_hadVetoEt") ) tree->Branch("mu_isolationR03_hadVetoEt", "std::vector<float>", &mu_isolationR03_hadVetoEt, buffersize);
   if( doWrite("mu_isolationR03_hoVetoEt") ) tree->Branch("mu_isolationR03_hoVetoEt", "std::vector<float>", &mu_isolationR03_hoVetoEt, buffersize);
   if( doWrite("mu_isolationR03_sumPt") ) tree->Branch("mu_isolationR03_sumPt", "std::vector<float>", &mu_isolationR03_sumPt, buffersize);
   if( doWrite("mu_isolationR03_emEt") ) tree->Branch("mu_isolationR03_emEt", "std::vector<float>", &mu_isolationR03_emEt, buffersize);
   if( doWrite("mu_isolationR03_hadEt") ) tree->Branch("mu_isolationR03_hadEt", "std::vector<float>", &mu_isolationR03_hadEt, buffersize);
   if( doWrite("mu_isolationR03_hoEt") ) tree->Branch("mu_isolationR03_hoEt", "std::vector<float>", &mu_isolationR03_hoEt, buffersize);
   if( doWrite("mu_isolationR03_nTracks") ) tree->Branch("mu_isolationR03_nTracks", "std::vector<int>", &mu_isolationR03_nTracks, buffersize);
   if( doWrite("mu_isolationR03_nJets") ) tree->Branch("mu_isolationR03_nJets", "std::vector<int>", &mu_isolationR03_nJets, buffersize);

   if( doWrite("mu_isolationR05_trackerVetoPt") ) tree->Branch("mu_isolationR05_trackerVetoPt", "std::vector<float>", &mu_isolationR05_trackerVetoPt, buffersize);
   if( doWrite("mu_isolationR05_emVetoEt") ) tree->Branch("mu_isolationR05_emVetoEt", "std::vector<float>", &mu_isolationR05_emVetoEt, buffersize);
   if( doWrite("mu_isolationR05_hadVetoEt") ) tree->Branch("mu_isolationR05_hadVetoEt", "std::vector<float>", &mu_isolationR05_hadVetoEt, buffersize);
   if( doWrite("mu_isolationR05_hoVetoEt") ) tree->Branch("mu_isolationR05_hoVetoEt", "std::vector<float>", &mu_isolationR05_hoVetoEt, buffersize);
   if( doWrite("mu_isolationR05_sumPt") ) tree->Branch("mu_isolationR05_sumPt", "std::vector<float>", &mu_isolationR05_sumPt, buffersize);
   if( doWrite("mu_isolationR05_emEt") ) tree->Branch("mu_isolationR05_emEt", "std::vector<float>", &mu_isolationR05_emEt, buffersize);
   if( doWrite("mu_isolationR05_hadEt") ) tree->Branch("mu_isolationR05_hadEt", "std::vector<float>", &mu_isolationR05_hadEt, buffersize);
   if( doWrite("mu_isolationR05_hoEt") ) tree->Branch("mu_isolationR05_hoEt", "std::vector<float>", &mu_isolationR05_hoEt, buffersize);
   if( doWrite("mu_isolationR05_nTracks") ) tree->Branch("mu_isolationR05_nTracks", "std::vector<int>", &mu_isolationR05_nTracks, buffersize);
   if( doWrite("mu_isolationR05_nJets") ) tree->Branch("mu_isolationR05_nJets", "std::vector<int>", &mu_isolationR05_nJets, buffersize);

   if( doWrite("mu_hasGlobalTrack") ) tree->Branch("mu_hasGlobalTrack", "std::vector<int>", &mu_hasGlobalTrack, buffersize);
   if( doWrite("mu_globalTrack_d0") ) tree->Branch("mu_globalTrack_d0", "std::vector<float>", &mu_globalTrack_d0, buffersize);
   if( doWrite("mu_globalTrack_z0") ) tree->Branch("mu_globalTrack_z0", "std::vector<float>", &mu_globalTrack_z0, buffersize);
   if( doWrite("mu_globalTrack_d0Error") ) tree->Branch("mu_globalTrack_d0Error", "std::vector<float>", &mu_globalTrack_d0Error, buffersize);
   if( doWrite("mu_globalTrack_z0Error") ) tree->Branch("mu_globalTrack_z0Error", "std::vector<float>", &mu_globalTrack_z0Error, buffersize);
   if( doWrite("mu_globalTrack_PV_dxy") ) tree->Branch("mu_globalTrack_PV_dxy", "std::vector<float>", &mu_globalTrack_PV_dxy, buffersize);
   if( doWrite("mu_globalTrack_PV_dz") ) tree->Branch("mu_globalTrack_PV_dz", "std::vector<float>", &mu_globalTrack_PV_dz, buffersize);
   if( doWrite("mu_globalTrack_RP_dxy") ) tree->Branch("mu_globalTrack_RP_dxy", "std::vector<float>", &mu_globalTrack_RP_dxy, buffersize);
   if( doWrite("mu_globalTrack_RP_dz") ) tree->Branch("mu_globalTrack_RP_dz", "std::vector<float>", &mu_globalTrack_RP_dz, buffersize);
   if( doWrite("mu_globalTrack_BS_dxy") ) tree->Branch("mu_globalTrack_BS_dxy", "std::vector<float>", &mu_globalTrack_BS_dxy, buffersize);
   if( doWrite("mu_globalTrack_BS_dz") ) tree->Branch("mu_globalTrack_BS_dz", "std::vector<float>", &mu_globalTrack_BS_dz, buffersize);
   if( doWrite("mu_globalTrack_dxyError") ) tree->Branch("mu_globalTrack_dxyError", "std::vector<float>", &mu_globalTrack_dxyError, buffersize);
   if( doWrite("mu_globalTrack_dzError") ) tree->Branch("mu_globalTrack_dzError", "std::vector<float>", &mu_globalTrack_dzError, buffersize);
   if( doWrite("mu_globalTrack_normalizedChi2") ) tree->Branch("mu_globalTrack_normalizedChi2", "std::vector<float>", &mu_globalTrack_normalizedChi2, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidHits") ) tree->Branch("mu_globalTrack_numberOfValidHits", "std::vector<int>", &mu_globalTrack_numberOfValidHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidMuonHits") ) tree->Branch("mu_globalTrack_numberOfValidMuonHits", "std::vector<int>", &mu_globalTrack_numberOfValidMuonHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfLostHits") ) tree->Branch("mu_globalTrack_numberOfLostHits", "std::vector<int>", &mu_globalTrack_numberOfLostHits, buffersize);
   if( doWrite("mu_globalTrack_pt") ) tree->Branch("mu_globalTrack_pt", "std::vector<float>", &mu_globalTrack_pt, buffersize);
   if( doWrite("mu_globalTrack_eta") ) tree->Branch("mu_globalTrack_eta", "std::vector<float>", &mu_globalTrack_eta, buffersize);
   if( doWrite("mu_globalTrack_phi") ) tree->Branch("mu_globalTrack_phi", "std::vector<float>", &mu_globalTrack_phi, buffersize);
   if( doWrite("mu_globalTrack_ptError") ) tree->Branch("mu_globalTrack_ptError", "std::vector<float>", &mu_globalTrack_ptError, buffersize);
   if( doWrite("mu_globalTrack_etaError") ) tree->Branch("mu_globalTrack_etaError", "std::vector<float>", &mu_globalTrack_etaError, buffersize);
   if( doWrite("mu_globalTrack_phiError") ) tree->Branch("mu_globalTrack_phiError", "std::vector<float>", &mu_globalTrack_phiError, buffersize);
   if( doWrite("mu_globalTrack_vx") ) tree->Branch("mu_globalTrack_vx", "std::vector<float>", &mu_globalTrack_vx, buffersize);
   if( doWrite("mu_globalTrack_vy") ) tree->Branch("mu_globalTrack_vy", "std::vector<float>", &mu_globalTrack_vy, buffersize);
   if( doWrite("mu_globalTrack_vz") ) tree->Branch("mu_globalTrack_vz", "std::vector<float>", &mu_globalTrack_vz, buffersize);
   if( doWrite("mu_globalTrack_qoverp") ) tree->Branch("mu_globalTrack_qoverp", "std::vector<float>", &mu_globalTrack_qoverp, buffersize);
   if( doWrite("mu_globalTrack_qoverpError") ) tree->Branch("mu_globalTrack_qoverpError", "std::vector<float>", &mu_globalTrack_qoverpError, buffersize);
   if( doWrite("mu_globalTrack_charge") ) tree->Branch("mu_globalTrack_charge", "std::vector<int>", &mu_globalTrack_charge, buffersize);
   if( doWrite("mu_globalTrack_trackerLayersWithMeasurement") ) tree->Branch("mu_globalTrack_trackerLayersWithMeasurement", "std::vector<int>", &mu_globalTrack_trackerLayersWithMeasurement, buffersize);
   if( doWrite("mu_globalTrack_pixelLayersWithMeasurement") ) tree->Branch("mu_globalTrack_pixelLayersWithMeasurement", "std::vector<int>", &mu_globalTrack_pixelLayersWithMeasurement, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("mu_globalTrack_trackerLayersWithoutMeasurement") ) tree->Branch("mu_globalTrack_trackerLayersWithoutMeasurement", "std::vector<int>", &mu_globalTrack_trackerLayersWithoutMeasurement, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidPixelHits") ) tree->Branch("mu_globalTrack_numberOfValidPixelHits", "std::vector<int>", &mu_globalTrack_numberOfValidPixelHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfLostPixelHits") ) tree->Branch("mu_globalTrack_numberOfLostPixelHits", "std::vector<int>", &mu_globalTrack_numberOfLostPixelHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfInnerHits") ) tree->Branch("mu_globalTrack_numberOfInnerHits", "std::vector<int>", &mu_globalTrack_numberOfInnerHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfOuterHits") ) tree->Branch("mu_globalTrack_numberOfOuterHits", "std::vector<int>", &mu_globalTrack_numberOfOuterHits, buffersize);
   if( doWrite("mu_globalTrack_validFraction") ) tree->Branch("mu_globalTrack_validFraction", "std::vector<float>", &mu_globalTrack_validFraction, buffersize);

   if( doWrite("mu_bestTrackType") ) tree->Branch("mu_bestTrackType", "std::vector<int>", &mu_bestTrackType, buffersize);
   if( doWrite("mu_hasBestTrack") ) tree->Branch("mu_hasBestTrack", "std::vector<int>", &mu_hasBestTrack, buffersize);
   if( doWrite("mu_bestTrack_d0") ) tree->Branch("mu_bestTrack_d0", "std::vector<float>", &mu_bestTrack_d0, buffersize);
   if( doWrite("mu_bestTrack_z0") ) tree->Branch("mu_bestTrack_z0", "std::vector<float>", &mu_bestTrack_z0, buffersize);
   if( doWrite("mu_bestTrack_d0Error") ) tree->Branch("mu_bestTrack_d0Error", "std::vector<float>", &mu_bestTrack_d0Error, buffersize);
   if( doWrite("mu_bestTrack_z0Error") ) tree->Branch("mu_bestTrack_z0Error", "std::vector<float>", &mu_bestTrack_z0Error, buffersize);
   if( doWrite("mu_bestTrack_PV_dxy") ) tree->Branch("mu_bestTrack_PV_dxy", "std::vector<float>", &mu_bestTrack_PV_dxy, buffersize);
   if( doWrite("mu_bestTrack_PV_dz") ) tree->Branch("mu_bestTrack_PV_dz", "std::vector<float>", &mu_bestTrack_PV_dz, buffersize);
   if( doWrite("mu_bestTrack_RP_dxy") ) tree->Branch("mu_bestTrack_RP_dxy", "std::vector<float>", &mu_bestTrack_RP_dxy, buffersize);
   if( doWrite("mu_bestTrack_RP_dz") ) tree->Branch("mu_bestTrack_RP_dz", "std::vector<float>", &mu_bestTrack_RP_dz, buffersize);
   if( doWrite("mu_bestTrack_BS_dxy") ) tree->Branch("mu_bestTrack_BS_dxy", "std::vector<float>", &mu_bestTrack_BS_dxy, buffersize);
   if( doWrite("mu_bestTrack_BS_dz") ) tree->Branch("mu_bestTrack_BS_dz", "std::vector<float>", &mu_bestTrack_BS_dz, buffersize);
   if( doWrite("mu_bestTrack_dxyError") ) tree->Branch("mu_bestTrack_dxyError", "std::vector<float>", &mu_bestTrack_dxyError, buffersize);
   if( doWrite("mu_bestTrack_dzError") ) tree->Branch("mu_bestTrack_dzError", "std::vector<float>", &mu_bestTrack_dzError, buffersize);
   if( doWrite("mu_bestTrack_normalizedChi2") ) tree->Branch("mu_bestTrack_normalizedChi2", "std::vector<float>", &mu_bestTrack_normalizedChi2, buffersize);
   if( doWrite("mu_bestTrack_numberOfValidHits") ) tree->Branch("mu_bestTrack_numberOfValidHits", "std::vector<int>", &mu_bestTrack_numberOfValidHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfLostHits") ) tree->Branch("mu_bestTrack_numberOfLostHits", "std::vector<int>", &mu_bestTrack_numberOfLostHits, buffersize);
   if( doWrite("mu_bestTrack_pt") ) tree->Branch("mu_bestTrack_pt", "std::vector<float>", &mu_bestTrack_pt, buffersize);
   if( doWrite("mu_bestTrack_eta") ) tree->Branch("mu_bestTrack_eta", "std::vector<float>", &mu_bestTrack_eta, buffersize);
   if( doWrite("mu_bestTrack_phi") ) tree->Branch("mu_bestTrack_phi", "std::vector<float>", &mu_bestTrack_phi, buffersize);
   if( doWrite("mu_bestTrack_ptError") ) tree->Branch("mu_bestTrack_ptError", "std::vector<float>", &mu_bestTrack_ptError, buffersize);
   if( doWrite("mu_bestTrack_etaError") ) tree->Branch("mu_bestTrack_etaError", "std::vector<float>", &mu_bestTrack_etaError, buffersize);
   if( doWrite("mu_bestTrack_phiError") ) tree->Branch("mu_bestTrack_phiError", "std::vector<float>", &mu_bestTrack_phiError, buffersize);
   if( doWrite("mu_bestTrack_vx") ) tree->Branch("mu_bestTrack_vx", "std::vector<float>", &mu_bestTrack_vx, buffersize);
   if( doWrite("mu_bestTrack_vy") ) tree->Branch("mu_bestTrack_vy", "std::vector<float>", &mu_bestTrack_vy, buffersize);
   if( doWrite("mu_bestTrack_vz") ) tree->Branch("mu_bestTrack_vz", "std::vector<float>", &mu_bestTrack_vz, buffersize);
   if( doWrite("mu_bestTrack_qoverp") ) tree->Branch("mu_bestTrack_qoverp", "std::vector<float>", &mu_bestTrack_qoverp, buffersize);
   if( doWrite("mu_bestTrack_qoverpError") ) tree->Branch("mu_bestTrack_qoverpError", "std::vector<float>", &mu_bestTrack_qoverpError, buffersize);
   if( doWrite("mu_bestTrack_charge") ) tree->Branch("mu_bestTrack_charge", "std::vector<int>", &mu_bestTrack_charge, buffersize);
   if( doWrite("mu_bestTrack_trackerLayersWithMeasurement") ) tree->Branch("mu_bestTrack_trackerLayersWithMeasurement", "std::vector<int>", &mu_bestTrack_trackerLayersWithMeasurement, buffersize);
   if( doWrite("mu_bestTrack_pixelLayersWithMeasurement") ) tree->Branch("mu_bestTrack_pixelLayersWithMeasurement", "std::vector<int>", &mu_bestTrack_pixelLayersWithMeasurement, buffersize);
   if( doWrite("mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("mu_bestTrack_trackerLayersWithoutMeasurement") ) tree->Branch("mu_bestTrack_trackerLayersWithoutMeasurement", "std::vector<int>", &mu_bestTrack_trackerLayersWithoutMeasurement, buffersize);
   if( doWrite("mu_bestTrack_numberOfValidPixelHits") ) tree->Branch("mu_bestTrack_numberOfValidPixelHits", "std::vector<int>", &mu_bestTrack_numberOfValidPixelHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfLostPixelHits") ) tree->Branch("mu_bestTrack_numberOfLostPixelHits", "std::vector<int>", &mu_bestTrack_numberOfLostPixelHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfInnerHits") ) tree->Branch("mu_bestTrack_numberOfInnerHits", "std::vector<int>", &mu_bestTrack_numberOfInnerHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfOuterHits") ) tree->Branch("mu_bestTrack_numberOfOuterHits", "std::vector<int>", &mu_bestTrack_numberOfOuterHits, buffersize);
   if( doWrite("mu_bestTrack_validFraction") ) tree->Branch("mu_bestTrack_validFraction", "std::vector<float>", &mu_bestTrack_validFraction, buffersize);

   if( doWrite("mu_hasInnerTrack") ) tree->Branch("mu_hasInnerTrack", "std::vector<int>", &mu_hasInnerTrack, buffersize);
   if( doWrite("mu_innerTrack_d0") ) tree->Branch("mu_innerTrack_d0", "std::vector<float>", &mu_innerTrack_d0, buffersize);
   if( doWrite("mu_innerTrack_z0") ) tree->Branch("mu_innerTrack_z0", "std::vector<float>", &mu_innerTrack_z0, buffersize);
   if( doWrite("mu_innerTrack_d0Error") ) tree->Branch("mu_innerTrack_d0Error", "std::vector<float>", &mu_innerTrack_d0Error, buffersize);
   if( doWrite("mu_innerTrack_z0Error") ) tree->Branch("mu_innerTrack_z0Error", "std::vector<float>", &mu_innerTrack_z0Error, buffersize);
   if( doWrite("mu_innerTrack_PV_dxy") ) tree->Branch("mu_innerTrack_PV_dxy", "std::vector<float>", &mu_innerTrack_PV_dxy, buffersize);
   if( doWrite("mu_innerTrack_PV_dz") ) tree->Branch("mu_innerTrack_PV_dz", "std::vector<float>", &mu_innerTrack_PV_dz, buffersize);
   if( doWrite("mu_innerTrack_RP_dxy") ) tree->Branch("mu_innerTrack_RP_dxy", "std::vector<float>", &mu_innerTrack_RP_dxy, buffersize);
   if( doWrite("mu_innerTrack_RP_dz") ) tree->Branch("mu_innerTrack_RP_dz", "std::vector<float>", &mu_innerTrack_RP_dz, buffersize);
   if( doWrite("mu_innerTrack_BS_dxy") ) tree->Branch("mu_innerTrack_BS_dxy", "std::vector<float>", &mu_innerTrack_BS_dxy, buffersize);
   if( doWrite("mu_innerTrack_BS_dz") ) tree->Branch("mu_innerTrack_BS_dz", "std::vector<float>", &mu_innerTrack_BS_dz, buffersize);
   if( doWrite("mu_innerTrack_dxyError") ) tree->Branch("mu_innerTrack_dxyError", "std::vector<float>", &mu_innerTrack_dxyError, buffersize);
   if( doWrite("mu_innerTrack_dzError") ) tree->Branch("mu_innerTrack_dzError", "std::vector<float>", &mu_innerTrack_dzError, buffersize);
   if( doWrite("mu_innerTrack_normalizedChi2") ) tree->Branch("mu_innerTrack_normalizedChi2", "std::vector<float>", &mu_innerTrack_normalizedChi2, buffersize);
   if( doWrite("mu_innerTrack_numberOfValidHits") ) tree->Branch("mu_innerTrack_numberOfValidHits", "std::vector<int>", &mu_innerTrack_numberOfValidHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfLostHits") ) tree->Branch("mu_innerTrack_numberOfLostHits", "std::vector<int>", &mu_innerTrack_numberOfLostHits, buffersize);
   if( doWrite("mu_innerTrack_pt") ) tree->Branch("mu_innerTrack_pt", "std::vector<float>", &mu_innerTrack_pt, buffersize);
   if( doWrite("mu_innerTrack_eta") ) tree->Branch("mu_innerTrack_eta", "std::vector<float>", &mu_innerTrack_eta, buffersize);
   if( doWrite("mu_innerTrack_phi") ) tree->Branch("mu_innerTrack_phi", "std::vector<float>", &mu_innerTrack_phi, buffersize);
   if( doWrite("mu_innerTrack_ptError") ) tree->Branch("mu_innerTrack_ptError", "std::vector<float>", &mu_innerTrack_ptError, buffersize);
   if( doWrite("mu_innerTrack_etaError") ) tree->Branch("mu_innerTrack_etaError", "std::vector<float>", &mu_innerTrack_etaError, buffersize);
   if( doWrite("mu_innerTrack_phiError") ) tree->Branch("mu_innerTrack_phiError", "std::vector<float>", &mu_innerTrack_phiError, buffersize);
   if( doWrite("mu_innerTrack_vx") ) tree->Branch("mu_innerTrack_vx", "std::vector<float>", &mu_innerTrack_vx, buffersize);
   if( doWrite("mu_innerTrack_vy") ) tree->Branch("mu_innerTrack_vy", "std::vector<float>", &mu_innerTrack_vy, buffersize);
   if( doWrite("mu_innerTrack_vz") ) tree->Branch("mu_innerTrack_vz", "std::vector<float>", &mu_innerTrack_vz, buffersize);
   if( doWrite("mu_innerTrack_qoverp") ) tree->Branch("mu_innerTrack_qoverp", "std::vector<float>", &mu_innerTrack_qoverp, buffersize);
   if( doWrite("mu_innerTrack_qoverpError") ) tree->Branch("mu_innerTrack_qoverpError", "std::vector<float>", &mu_innerTrack_qoverpError, buffersize);
   if( doWrite("mu_innerTrack_charge") ) tree->Branch("mu_innerTrack_charge", "std::vector<int>", &mu_innerTrack_charge, buffersize);
   if( doWrite("mu_innerTrack_trackerLayersWithMeasurement") ) tree->Branch("mu_innerTrack_trackerLayersWithMeasurement", "std::vector<int>", &mu_innerTrack_trackerLayersWithMeasurement, buffersize);
   if( doWrite("mu_innerTrack_pixelLayersWithMeasurement") ) tree->Branch("mu_innerTrack_pixelLayersWithMeasurement", "std::vector<int>", &mu_innerTrack_pixelLayersWithMeasurement, buffersize);
   if( doWrite("mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("mu_innerTrack_trackerLayersWithoutMeasurement") ) tree->Branch("mu_innerTrack_trackerLayersWithoutMeasurement", "std::vector<int>", &mu_innerTrack_trackerLayersWithoutMeasurement, buffersize);
   if( doWrite("mu_innerTrack_numberOfValidPixelHits") ) tree->Branch("mu_innerTrack_numberOfValidPixelHits", "std::vector<int>", &mu_innerTrack_numberOfValidPixelHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfLostPixelHits") ) tree->Branch("mu_innerTrack_numberOfLostPixelHits", "std::vector<int>", &mu_innerTrack_numberOfLostPixelHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfInnerHits") ) tree->Branch("mu_innerTrack_numberOfInnerHits", "std::vector<int>", &mu_innerTrack_numberOfInnerHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfOuterHits") ) tree->Branch("mu_innerTrack_numberOfOuterHits", "std::vector<int>", &mu_innerTrack_numberOfOuterHits, buffersize);
   if( doWrite("mu_innerTrack_validFraction") ) tree->Branch("mu_innerTrack_validFraction", "std::vector<float>", &mu_innerTrack_validFraction, buffersize);

   if( doWrite("mu_type") ) tree->Branch("mu_type", "std::vector<int>", &mu_type, buffersize);

   if( doWrite("mu_hasMCMatch") ) tree->Branch("mu_hasMCMatch", "std::vector<int>", &mu_hasMCMatch, buffersize);
   if( doWrite("mu_gen_pt") ) tree->Branch("mu_gen_pt", "std::vector<float>", &mu_gen_pt, buffersize);
   if( doWrite("mu_gen_eta") ) tree->Branch("mu_gen_eta", "std::vector<float>", &mu_gen_eta, buffersize);
   if( doWrite("mu_gen_phi") ) tree->Branch("mu_gen_phi", "std::vector<float>", &mu_gen_phi, buffersize);
   if( doWrite("mu_gen_m") ) tree->Branch("mu_gen_m", "std::vector<float>", &mu_gen_m, buffersize);
   if( doWrite("mu_gen_status") ) tree->Branch("mu_gen_status", "std::vector<int>", &mu_gen_status, buffersize);
   if( doWrite("mu_gen_id") ) tree->Branch("mu_gen_id", "std::vector<int>", &mu_gen_id, buffersize);
   if( doWrite("mu_gen_charge") ) tree->Branch("mu_gen_charge", "std::vector<int>", &mu_gen_charge, buffersize);
   if( doWrite("mu_gen_dr") ) tree->Branch("mu_gen_dr", "std::vector<float>", &mu_gen_dr, buffersize);

   if( doWrite("mu_hasMCMatchPAT") ) tree->Branch("mu_hasMCMatchPAT", "std::vector<int>", &mu_hasMCMatchPAT, buffersize);
   if( doWrite("mu_genPAT_pt") ) tree->Branch("mu_genPAT_pt", "std::vector<float>", &mu_genPAT_pt, buffersize);
   if( doWrite("mu_genPAT_eta") ) tree->Branch("mu_genPAT_eta", "std::vector<float>", &mu_genPAT_eta, buffersize);
   if( doWrite("mu_genPAT_phi") ) tree->Branch("mu_genPAT_phi", "std::vector<float>", &mu_genPAT_phi, buffersize);
   if( doWrite("mu_genPAT_m") ) tree->Branch("mu_genPAT_m", "std::vector<float>", &mu_genPAT_m, buffersize);
   if( doWrite("mu_genPAT_status") ) tree->Branch("mu_genPAT_status", "std::vector<int>", &mu_genPAT_status, buffersize);
   if( doWrite("mu_genPAT_id") ) tree->Branch("mu_genPAT_id", "std::vector<int>", &mu_genPAT_id, buffersize);
   if( doWrite("mu_genPAT_charge") ) tree->Branch("mu_genPAT_charge", "std::vector<int>", &mu_genPAT_charge, buffersize);

   if( doWrite("tau_n") ) tree->Branch("tau_n", &tau_n, "tau_n/I", buffersize);
   if( doWrite("tau_pt") ) tree->Branch("tau_pt", "std::vector<float>", &tau_pt, buffersize);
   if( doWrite("tau_eta") ) tree->Branch("tau_eta", "std::vector<float>", &tau_eta, buffersize);
   if( doWrite("tau_phi") ) tree->Branch("tau_phi", "std::vector<float>", &tau_phi, buffersize);
   if( doWrite("tau_m") ) tree->Branch("tau_m", "std::vector<float>", &tau_m, buffersize);
   if( doWrite("tau_E") ) tree->Branch("tau_E", "std::vector<float>", &tau_E, buffersize);
   if( doWrite("tau_id") ) tree->Branch("tau_id", "std::vector<int>", &tau_id, buffersize);
   if( doWrite("tau_charge") ) tree->Branch("tau_charge", "std::vector<int>", &tau_charge, buffersize);

   if( doWrite("tau_hasLeadChargedHadrCand") ) tree->Branch("tau_hasLeadChargedHadrCand", "std::vector<bool>", &tau_hasLeadChargedHadrCand, buffersize);
   if( doWrite("tau_leadingTrackPt") ) tree->Branch("tau_leadingTrackPt", "std::vector<float>", &tau_leadingTrackPt, buffersize);
   if( doWrite("tau_leadingTrackCharge") ) tree->Branch("tau_leadingTrackCharge", "std::vector<int>", &tau_leadingTrackCharge, buffersize);
   if( doWrite("tau_leadingTrackDz") ) tree->Branch("tau_leadingTrackDz", "std::vector<float>", &tau_leadingTrackDz, buffersize);
   if( doWrite("tau_leadingTrackDxy") ) tree->Branch("tau_leadingTrackDxy", "std::vector<float>", &tau_leadingTrackDxy, buffersize);

   if( doWrite("tau_decayMode") ) tree->Branch("tau_decayMode", "std::vector<int>", &tau_decayMode, buffersize);
   if( doWrite("tau_decayModeFindingOldDMs") ) tree->Branch("tau_decayModeFindingOldDMs", "std::vector<float>", &tau_decayModeFindingOldDMs, buffersize);
   if( doWrite("tau_decayModeFindingNewDMs") ) tree->Branch("tau_decayModeFindingNewDMs", "std::vector<float>", &tau_decayModeFindingNewDMs, buffersize);

   if( doWrite("tau_puCorrPtSum") ) tree->Branch("tau_puCorrPtSum", "std::vector<float>", &tau_puCorrPtSum, buffersize);
   if( doWrite("tau_neutralIsoPtSum") ) tree->Branch("tau_neutralIsoPtSum", "std::vector<float>", &tau_neutralIsoPtSum, buffersize);
   if( doWrite("tau_chargedIsoPtSum") ) tree->Branch("tau_chargedIsoPtSum", "std::vector<float>", &tau_chargedIsoPtSum, buffersize);
   if( doWrite("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits") ) tree->Branch("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", "std::vector<float>", &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, buffersize);

   if( doWrite("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, buffersize);
   if( doWrite("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, buffersize);
   if( doWrite("tau_byTightCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits, buffersize);

   if( doWrite("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
   if( doWrite("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
   if( doWrite("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
   if( doWrite("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT, buffersize);

   if( doWrite("tau_againstMuonLoose3") ) tree->Branch("tau_againstMuonLoose3", "std::vector<float>", &tau_againstMuonLoose3, buffersize);
   if( doWrite("tau_againstMuonTight3") ) tree->Branch("tau_againstMuonTight3", "std::vector<float>", &tau_againstMuonTight3, buffersize);

   if( doWrite("tau_pfEssential_jet_pt") ) tree->Branch("tau_pfEssential_jet_pt", "std::vector<float>", &tau_pfEssential_jet_pt, buffersize);
   if( doWrite("tau_pfEssential_jet_eta") ) tree->Branch("tau_pfEssential_jet_eta", "std::vector<float>", &tau_pfEssential_jet_eta, buffersize);
   if( doWrite("tau_pfEssential_jet_phi") ) tree->Branch("tau_pfEssential_jet_phi", "std::vector<float>", &tau_pfEssential_jet_phi, buffersize);
   if( doWrite("tau_pfEssential_jet_m") ) tree->Branch("tau_pfEssential_jet_m", "std::vector<float>", &tau_pfEssential_jet_m, buffersize);

   if( doWrite("tau_pfEssential_jetCorr_pt") ) tree->Branch("tau_pfEssential_jetCorr_pt", "std::vector<float>", &tau_pfEssential_jetCorr_pt, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_eta") ) tree->Branch("tau_pfEssential_jetCorr_eta", "std::vector<float>", &tau_pfEssential_jetCorr_eta, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_phi") ) tree->Branch("tau_pfEssential_jetCorr_phi", "std::vector<float>", &tau_pfEssential_jetCorr_phi, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_m") ) tree->Branch("tau_pfEssential_jetCorr_m", "std::vector<float>", &tau_pfEssential_jetCorr_m, buffersize);

   if( doWrite("tau_pfEssential_hasSV") ) tree->Branch("tau_pfEssential_hasSV", "std::vector<bool>", &tau_pfEssential_hasSV, buffersize);
   if( doWrite("tau_pfEssential_sv_x") ) tree->Branch("tau_pfEssential_sv_x", "std::vector<float>", &tau_pfEssential_sv_x, buffersize);
   if( doWrite("tau_pfEssential_sv_y") ) tree->Branch("tau_pfEssential_sv_y", "std::vector<float>", &tau_pfEssential_sv_y, buffersize);
   if( doWrite("tau_pfEssential_sv_z") ) tree->Branch("tau_pfEssential_sv_z", "std::vector<float>", &tau_pfEssential_sv_z, buffersize);

   if( doWrite("tau_pfEssential_flightLengthSig") ) tree->Branch("tau_pfEssential_flightLengthSig", "std::vector<float>", &tau_pfEssential_flightLengthSig, buffersize);
   if( doWrite("tau_pfEssential_dxy") ) tree->Branch("tau_pfEssential_dxy", "std::vector<float>", &tau_pfEssential_dxy, buffersize);
   if( doWrite("tau_pfEssential_dxy_error") ) tree->Branch("tau_pfEssential_dxy_error", "std::vector<float>", &tau_pfEssential_dxy_error, buffersize);
   if( doWrite("tau_pfEssential_dxy_Sig") ) tree->Branch("tau_pfEssential_dxy_Sig", "std::vector<float>", &tau_pfEssential_dxy_Sig, buffersize);

   if( doWrite("jet_n") ) tree->Branch("jet_n", &jet_n, "jet_n/I", buffersize);
   if( doWrite("jet_pt") ) tree->Branch("jet_pt", "std::vector<float>", &jet_pt, buffersize);
   if( doWrite("jet_eta") ) tree->Branch("jet_eta", "std::vector<float>", &jet_eta, buffersize);
   if( doWrite("jet_phi") ) tree->Branch("jet_phi", "std::vector<float>", &jet_phi, buffersize);
   if( doWrite("jet_m") ) tree->Branch("jet_m", "std::vector<float>", &jet_m, buffersize);
   if( doWrite("jet_E") ) tree->Branch("jet_E", "std::vector<float>", &jet_E, buffersize);

   if( doWrite("jet_ntrk") ) tree->Branch("jet_ntrk", "std::vector<int>", &jet_ntrk, buffersize);

   if( doWrite("jet_JBP") ) tree->Branch("jet_JBP", "std::vector<float>", &jet_JBP, buffersize);
   if( doWrite("jet_JP") ) tree->Branch("jet_JP", "std::vector<float>", &jet_JP, buffersize);
   if( doWrite("jet_TCHP") ) tree->Branch("jet_TCHP", "std::vector<float>", &jet_TCHP, buffersize);
   if( doWrite("jet_TCHE") ) tree->Branch("jet_TCHE", "std::vector<float>", &jet_TCHE, buffersize);
   if( doWrite("jet_SSVHE") ) tree->Branch("jet_SSVHE", "std::vector<float>", &jet_SSVHE, buffersize);
   if( doWrite("jet_SSVHP") ) tree->Branch("jet_SSVHP", "std::vector<float>", &jet_SSVHP, buffersize);
   if( doWrite("jet_CMVA") ) tree->Branch("jet_CMVA", "std::vector<float>", &jet_CMVA, buffersize);
   if( doWrite("jet_CSVv2") ) tree->Branch("jet_CSVv2", "std::vector<float>", &jet_CSVv2, buffersize);
   if( doWrite("jet_cMVAv2") ) tree->Branch("jet_cMVAv2", "std::vector<float>", &jet_cMVAv2, buffersize);
   if( doWrite("jet_CharmCvsL") ) tree->Branch("jet_CharmCvsL", "std::vector<float>", &jet_CharmCvsL, buffersize);
   if( doWrite("jet_CharmCvsB") ) tree->Branch("jet_CharmCvsB", "std::vector<float>", &jet_CharmCvsB, buffersize);
   if( doWrite("jet_partonFlavour") ) tree->Branch("jet_partonFlavour", "std::vector<int>", &jet_partonFlavour, buffersize);
   if( doWrite("jet_hadronFlavour") ) tree->Branch("jet_hadronFlavour", "std::vector<int>", &jet_hadronFlavour, buffersize);

   if( doWrite("jet_neutralHadronEnergy") ) tree->Branch("jet_neutralHadronEnergy", "std::vector<float>", &jet_neutralHadronEnergy, buffersize);
   if( doWrite("jet_neutralEmEnergy") ) tree->Branch("jet_neutralEmEnergy", "std::vector<float>", &jet_neutralEmEnergy, buffersize);
   if( doWrite("jet_chargedHadronEnergy") ) tree->Branch("jet_chargedHadronEnergy", "std::vector<float>", &jet_chargedHadronEnergy, buffersize);
   if( doWrite("jet_chargedEmEnergy") ) tree->Branch("jet_chargedEmEnergy", "std::vector<float>", &jet_chargedEmEnergy, buffersize);
   if( doWrite("jet_electronEnergy") ) tree->Branch("jet_electronEnergy", "std::vector<float>", &jet_electronEnergy, buffersize);
   if( doWrite("jet_muonEnergy") ) tree->Branch("jet_muonEnergy", "std::vector<float>", &jet_muonEnergy, buffersize);
   if( doWrite("jet_photonEnergy") ) tree->Branch("jet_photonEnergy", "std::vector<float>", &jet_photonEnergy, buffersize);

   if( doWrite("jet_charge") ) tree->Branch("jet_charge", "std::vector<int>", &jet_charge, buffersize);
   if( doWrite("jet_chargeVec") ) tree->Branch("jet_chargeVec", "std::vector<int>", &jet_chargeVec, buffersize);
   if( doWrite("jet_chargedMultiplicity") ) tree->Branch("jet_chargedMultiplicity", "std::vector<int>", &jet_chargedMultiplicity, buffersize);
   if( doWrite("jet_neutralMultiplicity") ) tree->Branch("jet_neutralMultiplicity", "std::vector<int>", &jet_neutralMultiplicity, buffersize);
   if( doWrite("jet_chargedHadronMultiplicity") ) tree->Branch("jet_chargedHadronMultiplicity", "std::vector<int>", &jet_chargedHadronMultiplicity, buffersize);

   if( doWrite("jet_jetArea") ) tree->Branch("jet_jetArea", "std::vector<float>", &jet_jetArea, buffersize);

   if( doWrite("jet_jecFactorUncorrected") ) tree->Branch("jet_jecFactorUncorrected", "std::vector<float>", &jet_jecFactorUncorrected, buffersize);
   if( doWrite("jet_jecFactorL1FastJet") ) tree->Branch("jet_jecFactorL1FastJet", "std::vector<float>", &jet_jecFactorL1FastJet, buffersize);
   if( doWrite("jet_jecFactorL2Relative") ) tree->Branch("jet_jecFactorL2Relative", "std::vector<float>", &jet_jecFactorL2Relative, buffersize);
   if( doWrite("jet_jecFactorL3Absolute") ) tree->Branch("jet_jecFactorL3Absolute", "std::vector<float>", &jet_jecFactorL3Absolute, buffersize);

   if( doWrite("jet_neutralHadronEnergyFraction") ) tree->Branch("jet_neutralHadronEnergyFraction", "std::vector<float>", &jet_neutralHadronEnergyFraction, buffersize);
   if( doWrite("jet_neutralEmEnergyFraction") ) tree->Branch("jet_neutralEmEnergyFraction", "std::vector<float>", &jet_neutralEmEnergyFraction, buffersize);
   if( doWrite("jet_chargedHadronEnergyFraction") ) tree->Branch("jet_chargedHadronEnergyFraction", "std::vector<float>", &jet_chargedHadronEnergyFraction, buffersize);
   if( doWrite("jet_muonEnergyFraction") ) tree->Branch("jet_muonEnergyFraction", "std::vector<float>", &jet_muonEnergyFraction, buffersize);
   if( doWrite("jet_chargedEmEnergyFraction") ) tree->Branch("jet_chargedEmEnergyFraction", "std::vector<float>", &jet_chargedEmEnergyFraction, buffersize);

   if( doWrite("jet_Unc") ) tree->Branch("jet_Unc", "std::vector<float>", &jet_Unc, buffersize);

   if( doWrite("jet_pileupJetId") ) tree->Branch("jet_pileupJetId", "std::vector<float>", &jet_pileupJetId, buffersize);

   if( doWrite("jet_looseJetID") ) tree->Branch("jet_looseJetID", "std::vector<bool>", &jet_looseJetID, buffersize);
   if( doWrite("jet_tightJetID") ) tree->Branch("jet_tightJetID", "std::vector<bool>", &jet_tightJetID, buffersize);

   if( doWrite("jet_hasGenJet") ) tree->Branch("jet_hasGenJet", "std::vector<bool>", &jet_hasGenJet, buffersize);
   if( doWrite("jet_genJet_pt") ) tree->Branch("jet_genJet_pt", "std::vector<float>", &jet_genJet_pt, buffersize);
   if( doWrite("jet_genJet_eta") ) tree->Branch("jet_genJet_eta", "std::vector<float>", &jet_genJet_eta, buffersize);
   if( doWrite("jet_genJet_phi") ) tree->Branch("jet_genJet_phi", "std::vector<float>", &jet_genJet_phi, buffersize);
   if( doWrite("jet_genJet_m") ) tree->Branch("jet_genJet_m", "std::vector<float>", &jet_genJet_m, buffersize);
   if( doWrite("jet_genJet_E") ) tree->Branch("jet_genJet_E", "std::vector<float>", &jet_genJet_E, buffersize);
   if( doWrite("jet_genJet_status") ) tree->Branch("jet_genJet_status", "std::vector<int>", &jet_genJet_status, buffersize);
   if( doWrite("jet_genJet_id") ) tree->Branch("jet_genJet_id", "std::vector<int>", &jet_genJet_id, buffersize);

   if( doWrite("jet_hasGenParton") ) tree->Branch("jet_hasGenParton", "std::vector<bool>", &jet_hasGenParton, buffersize);
   if( doWrite("jet_genParton_pt") ) tree->Branch("jet_genParton_pt", "std::vector<float>", &jet_genParton_pt, buffersize);
   if( doWrite("jet_genParton_eta") ) tree->Branch("jet_genParton_eta", "std::vector<float>", &jet_genParton_eta, buffersize);
   if( doWrite("jet_genParton_phi") ) tree->Branch("jet_genParton_phi", "std::vector<float>", &jet_genParton_phi, buffersize);
   if( doWrite("jet_genParton_m") ) tree->Branch("jet_genParton_m", "std::vector<float>", &jet_genParton_m, buffersize);
   if( doWrite("jet_genParton_E") ) tree->Branch("jet_genParton_E", "std::vector<float>", &jet_genParton_E, buffersize);
   if( doWrite("jet_genParton_status") ) tree->Branch("jet_genParton_status", "std::vector<int>", &jet_genParton_status, buffersize);
   if( doWrite("jet_genParton_id") ) tree->Branch("jet_genParton_id", "std::vector<int>", &jet_genParton_id, buffersize);

   //------------------------
   //  genJet collection
   //------------------------

   if( doWrite("genJet_n") ) tree->Branch("genJet_n", &genJet_n, "genJet_n/I", buffersize);
   if( doWrite("genJet_pt") ) tree->Branch("genJet_pt", "std::vector<float>", &genJet_pt, buffersize);
   if( doWrite("genJet_eta") ) tree->Branch("genJet_eta", "std::vector<float>", &genJet_eta, buffersize);
   if( doWrite("genJet_phi") ) tree->Branch("genJet_phi", "std::vector<float>", &genJet_phi, buffersize);
   if( doWrite("genJet_m") ) tree->Branch("genJet_m", "std::vector<float>", &genJet_m, buffersize);
   if( doWrite("genJet_E") ) tree->Branch("genJet_E", "std::vector<float>", &genJet_E, buffersize);
   if( doWrite("genJet_emEnergy") ) tree->Branch("genJet_emEnergy", "std::vector<float>", &genJet_emEnergy, buffersize);
   if( doWrite("genJet_hadEnergy") ) tree->Branch("genJet_hadEnergy", "std::vector<float>", &genJet_hadEnergy, buffersize);
   if( doWrite("genJet_invisibleEnergy") ) tree->Branch("genJet_invisibleEnergy", "std::vector<float>", &genJet_invisibleEnergy, buffersize);
   if( doWrite("genJet_auxiliaryEnergy") ) tree->Branch("genJet_auxiliaryEnergy", "std::vector<float>", &genJet_auxiliaryEnergy, buffersize);
   if( doWrite("genJet_flavour") ) tree->Branch("genJet_flavour", "std::vector<int>", &genJet_flavour, buffersize);

   if( doWrite("gen_PVz") ) tree->Branch("gen_PVz", &gen_PVz, "gen_PVz/F", buffersize);

   if( doWrite("gen_all") )
     {
	tree->Branch("gen_n", &gen_n, "gen_n/I", buffersize);
	tree->Branch("gen_pt", "std::vector<float>", &gen_pt, buffersize);
	tree->Branch("gen_eta", "std::vector<float>", &gen_eta, buffersize);
	tree->Branch("gen_phi", "std::vector<float>", &gen_phi, buffersize);
	tree->Branch("gen_m", "std::vector<float>", &gen_m, buffersize);
	tree->Branch("gen_E", "std::vector<float>", &gen_E, buffersize);
	tree->Branch("gen_status", "std::vector<int>", &gen_status, buffersize);
	tree->Branch("gen_id", "std::vector<int>", &gen_id, buffersize);
	tree->Branch("gen_charge", "std::vector<int>", &gen_charge, buffersize);
	tree->Branch("gen_index", "std::vector<int>", &gen_index, buffersize);
	tree->Branch("gen_mother_index", "std::vector<int>", &gen_mother_index, buffersize);
	tree->Branch("gen_daughter_n", "std::vector<int>", &gen_daughter_n, buffersize);
	tree->Branch("gen_daughter_index", "std::vector<std::vector<int> >", &gen_daughter_index, buffersize);
     }
   
   if( doWrite("genTTX_id") ) tree->Branch("genTTX_id", &genTTX_id, "genTTX_id/I", buffersize);
}

bool FlatTree::doWrite(const std::string& name)
{
   std::map<std::string,bool>::iterator it = conf.find(name);
   if( it != conf.end() )
     return it->second;
   return 0;
}
