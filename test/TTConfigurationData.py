import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

process = Framework.create(True, eras.Run2_25ns, '74X_dataRun2_v2', cms.PSet(
    tt = cms.PSet(
        type = cms.string('tt_analyzer'),
        prefix = cms.string('tt_'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            electronPtCut = cms.untracked.double(20),
            electronEtaCut = cms.untracked.double(2.5),
            electronVetoIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-veto'),
            electronLooseIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-loose'),
            electronMediumIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-medium'),
            electronTightIDName = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-tight'),
            
            muonPtCut = cms.untracked.double(20),
            muonEtaCut = cms.untracked.double(2.4),
            muonBaseIsoCut = cms.untracked.double(.20), # Loose cut recommended for dilepton analysis

            jetPtCut = cms.untracked.double(30),
            jetEtaCut = cms.untracked.double(2.5),
            #jetPUID = cms.untracked.double(-9999999),
            jetDRleptonCut = cms.untracked.double(0.3),
            #jetID = cms.untracked.string('tight'), # not tightLeptonVeto since DeltaR(l,j) cut should be enough
            jetCSVv2Name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            jetCSVv2L = cms.untracked.double(0.605),
            jetCSVv2M = cms.untracked.double(0.89),
            jetCSVv2T = cms.untracked.double(0.97),

            hltDRCut = cms.untracked.double(0.05), # DeltaR cut for trigger matching
            hltDPtCut = cms.untracked.double(0.1), #Delta(Pt)/Pt cut for trigger matching
            ),
        categories_parameters = cms.PSet(
            MllCutSF = cms.untracked.double(20),
            MllCutDF = cms.untracked.double(0),
            MllZVetoCutLow = cms.untracked.double(86),
            MllZVetoCutHigh = cms.untracked.double(116),
            HLTDoubleMuon = cms.untracked.vstring('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2'),
            HLTDoubleEG = cms.untracked.vstring('HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3'),
            HLTMuonEG = cms.untracked.vstring('HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3', 'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3'),
            ),
        )
    ), redoJEC=True
    )

Framework.schedule(process, ['tt']) 

process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/292/00000/E20B7746-8245-E511-B593-02163E0135AD.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
