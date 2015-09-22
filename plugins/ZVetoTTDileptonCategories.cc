#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>

#include <cp3_llbb/TTAnalysis/interface/TTAnalyzer.h>
#include <cp3_llbb/TTAnalysis/interface/ZVetoTTDileptonCategories.h>

// ***** ***** *****
// Dilepton Z veto category
// ***** ***** *****

bool ZVetoTTDileptonCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const { return true; }

//void ZVetoTTDileptonCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
//void ZVetoTTDileptonCategory::register_cuts(CutManager& manager) {
//    manager.new_cut("ll_mass_Zveto", "116 > mll > 86");
//}

bool ZVetoTTDileptonCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    // El-El
    for(unsigned int iele1 = 0; iele1 < tt.selectedElectrons.size(); iele1++){
        for(unsigned int iele2 = iele1+1; iele2 < tt.selectedElectrons.size(); iele2++){
            if( electrons.charge[ tt.selectedElectrons[iele1] ] * electrons.charge[ tt.selectedElectrons[iele2] ] < 0 ){
                LorentzVector dilep = electrons.p4[ tt.selectedElectrons[iele1] ] + electrons.p4[ tt.selectedElectrons[iele2] ];
                if( dilep.M() > m_mll_cut_low && dilep.M() < m_mll_cut_high)
                    return false;
            }
        }
    }

    // Mu-Mu
    for(unsigned int imu1 = 0; imu1 < tt.selectedMuons.size(); imu1++){
        for(unsigned int imu2 = imu1+1; imu2 < tt.selectedMuons.size(); imu2++){
            if( muons.charge[ tt.selectedMuons[imu1] ] * muons.charge[ tt.selectedMuons[imu2] ] < 0 ){
                LorentzVector dilep = muons.p4[ tt.selectedMuons[imu1] ] + muons.p4[ tt.selectedMuons[imu2] ];
                if( dilep.M() > m_mll_cut_low && dilep.M() < m_mll_cut_high)
                    return false;
            }
        }
    }

    return true;
}

