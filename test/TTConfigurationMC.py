import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

process = Framework.create(False, eras.Run2_25ns, '74X_mcRun2_asymptotic_v2', cms.PSet(
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
            muonLooseIsoCut = cms.untracked.double(.20), # Loose cut recommended for dilepton analysis
            muonTightIsoCut = cms.untracked.double(.12),

            jetPtCut = cms.untracked.double(30),
            jetEtaCut = cms.untracked.double(2.5),
            #jetPUID = cms.untracked.double(-9999999),
            jetDRleptonCut = cms.untracked.double(0.3),
            jetID = cms.untracked.string('tight'), # not tightLeptonVeto since DeltaR(l,j) cut should be enough
            jetCSVv2Name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            jetCSVv2L = cms.untracked.double(0.605),
            jetCSVv2M = cms.untracked.double(0.89),
            jetCSVv2T = cms.untracked.double(0.97),

            hltDRCut = cms.untracked.double(0.05), # DeltaR cut for trigger matching
            hltDPtCut = cms.untracked.double(0.1), #Delta(Pt)/Pt cut for trigger matching
            ),
        categories_parameters = cms.PSet(
            MllCutSF = cms.untracked.double(20),
            MllCutDF = cms.untracked.double(20),
            MllZVetoCutLow = cms.untracked.double(76),
            MllZVetoCutHigh = cms.untracked.double(116),
            HLTDoubleMuon = cms.untracked.vstring('HLT_Mu17_TrkIsoVVL_(Tk)?Mu8_TrkIsoVVL_DZ_v.*'),
            HLTDoubleEG = cms.untracked.vstring('HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.*'),
            HLTMuonEG = cms.untracked.vstring('HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.*', 'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.*'),
            ),
        )
    ), redoJEC=False
    )

del process.framework.producers.fat_jets

Framework.schedule(process, ['tt'])

#process.source.firstEvent = cms.untracked.uint32(13083444)
#process.source.firstLuminosityBlock = cms.untracked.uint32(52386)

# Tricky gen event from /store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/00000/0014DC94-DC5C-E511-82FB-7845C4FC39F5.root
# First one is g g -> t tbar with one W -> bbar c
# Second is b bar -> t tbar semi-leptonic
#process.source.eventsToProcess = cms.untracked.VEventRange(
        #'1:52386:13083444',
        #'1:34020:8496854'
        #)

# Other tricky gen events, with lots of ISR
# From file:/nfs/scratch/fynu/swertz/CMSSW_7_4_15/src/cp3_llbb/TTAnalysis/test/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_miniAODv2_oneFile.root
#process.source.eventsToProcess = cms.untracked.VEventRange(
        #'1:321521:80300260',
        #'1:357590:89308562',
        #'1:387992:96901374'
        #)
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source.fileNames = cms.untracked.vstring(
        'file:/nfs/scratch/fynu/swertz/CMSSW_7_4_15/src/cp3_llbb/TTAnalysis/test/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_miniAODv2_oneFile.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
