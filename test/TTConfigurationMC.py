import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

process = Framework.create(False, eras.Run2_25ns, '74X_mcRun2_asymptotic_v2', cms.PSet(
    bTagsLoose = cms.PSet(
        type = cms.string('btags_analyzer'),
        prefix = cms.string('btags_CSVv2_loose'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            discr_name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            discr_cut = cms.untracked.double(0.605),
            eta_cut = cms.untracked.double(2.4),
            pt_cut = cms.untracked.double(30)
            )
        ),

    bTagsMedium = cms.PSet(
        type = cms.string('btags_analyzer'),
        prefix = cms.string('btags_CSVv2_medium'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            discr_name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            discr_cut = cms.untracked.double(0.89),
            eta_cut = cms.untracked.double(2.4),
            pt_cut = cms.untracked.double(30)
            )
        ),

    bTagsTight = cms.PSet(
        type = cms.string('btags_analyzer'),
        prefix = cms.string('btags_CSVv2_tight'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            discr_name = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
            discr_cut = cms.untracked.double(0.97),
            eta_cut = cms.untracked.double(2.4),
            pt_cut = cms.untracked.double(30)
            )
        ),

    tt = cms.PSet(
        type = cms.string('tt_analyzer'),
        prefix = cms.string('tt_'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            #electronIsoCut = cms.untracked.double(0.11), # Iso cut already applied in cutBasedID
            electronPtCut = cms.untracked.double(20),
            electronEtaCut = cms.untracked.double(2.4),
            electron_veto_wp_name = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-veto'),
            electron_medium_wp_name = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-medium'),
            electron_loose_wp_name = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-loose'),
            electron_tight_wp_name = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-tight'),
            electron_selectedID_name = cms.untracked.string('cutBasedElectronID-Spring15-25ns-V1-standalone-tight'),
            
            muonIsoCut = cms.untracked.double(.20), # Loose cut recommended for dilepton analysis
            muonPtCut = cms.untracked.double(20),
            muonEtaCut = cms.untracked.double(2.4),
            muon_selectedID_wp = cms.untracked.string('tight'),
            
            jetPtCut = cms.untracked.double(30),
            jetEtaCut = cms.untracked.double(2.5),
            ),
        categories_parameters = cms.PSet(
            mll_cut = cms.untracked.double(20),
            mll_ZVetoCut_low = cms.untracked.double(86),
            mll_ZVetoCut_high = cms.untracked.double(116)
            ),
        )
    ), redoJEC=True
    )

Framework.schedule(process, ['bTagsLoose', 'bTagsMedium', 'bTagsTight', 'tt'])

process.source.fileNames = cms.untracked.vstring(
        'file:///home/fynu/sbrochet/storage/MINIAODSIM/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9_reduced.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
