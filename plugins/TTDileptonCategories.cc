#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>

#include <cp3_llbb/TTAnalysis/interface/TTDileptonCategories.h>
#include <cp3_llbb/TTAnalysis/interface/TTAnalyzer.h>

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool TTMuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return muons.p4.size() >= 2 ;
};

bool TTMuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const TTAnalyzer& tt_analyzer = analyzers.get<TTAnalyzer>("tt");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    int selMu1 = tt_analyzer.selectedLeadingMuMu.first;
    int selMu2 = tt_analyzer.selectedLeadingMuMu.second;

    if( selMu1 < 0 )
        return false;

    if( muons.p4[selMu1].Pt() < electrons.p4[ tt_analyzer.selectedLeadingElMu.first ].Pt() )
        return false;

    if( muons.p4[selMu1].Pt() < electrons.p4[ tt_analyzer.selectedLeadingElEl.first ].Pt() )
        return false;

    if( muons.p4[selMu2].Pt() < electrons.p4[ tt_analyzer.selectedLeadingMuEl.second ].Pt() )
        return false;

    return true;
};

void TTMuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
};

void TTMuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const TTAnalyzer& tt_analyzer = analyzers.get<TTAnalyzer>("tt");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");

    int selMu1 = tt_analyzer.selectedLeadingMuMu.first;
    int selMu2 = tt_analyzer.selectedLeadingMuMu.second;

    if( selMu1 >= 0 ){
        if( (muons.p4[selMu1] + muons.p4[selMu2]).M() > m_mll_cut)
            manager.pass_cut("ll_mass");
    }
}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****
bool TTMuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return muons.p4.size() >= 1 && electrons.p4.size() >= 1;
};

bool TTMuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const TTAnalyzer& tt_analyzer = analyzers.get<TTAnalyzer>("tt");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    int selMu = tt_analyzer.selectedLeadingMuEl.first;
    int selEl = tt_analyzer.selectedLeadingMuEl.second;

    if( selMu < 0 )
        return false;

    if( muons.p4[selMu].Pt() < electrons.p4[ tt_analyzer.selectedLeadingElMu.first ].Pt() )
        return false;

    if( muons.p4[selMu].Pt() < electrons.p4[ tt_analyzer.selectedLeadingElEl.first ].Pt() )
        return false;

    if( electrons.p4[selEl].Pt() < muons.p4[ tt_analyzer.selectedLeadingMuMu.second ].Pt() )
        return false;

    return true;
};

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool TTElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return muons.p4.size() >= 1 && electrons.p4.size() >= 1;
};

bool TTElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const TTAnalyzer& tt_analyzer = analyzers.get<TTAnalyzer>("tt");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    int selEl = tt_analyzer.selectedLeadingElMu.first;
    int selMu = tt_analyzer.selectedLeadingElMu.second;

    if( selEl < 0 )
        return false;

    if( electrons.p4[selEl].Pt() < muons.p4[ tt_analyzer.selectedLeadingMuEl.first ].Pt() )
        return false;

    if( electrons.p4[selEl].Pt() < muons.p4[ tt_analyzer.selectedLeadingMuMu.first ].Pt() )
        return false;

    if( muons.p4[selMu].Pt() < electrons.p4[ tt_analyzer.selectedLeadingElEl.second ].Pt() )
        return false;

    return true;
};

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool TTElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return electrons.p4.size() >= 2 ;
};

bool TTElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const TTAnalyzer& tt_analyzer = analyzers.get<TTAnalyzer>("tt");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    int selEl1 = tt_analyzer.selectedLeadingElEl.first;
    int selEl2 = tt_analyzer.selectedLeadingElEl.second;

    if( selEl1 < 0 )
        return false;

    if( electrons.p4[selEl1].Pt() < muons.p4[ tt_analyzer.selectedLeadingMuEl.first ].Pt() )
        return false;

    if( electrons.p4[selEl1].Pt() < muons.p4[ tt_analyzer.selectedLeadingMuMu.first ].Pt() )
        return false;

    if( electrons.p4[selEl2].Pt() < muons.p4[ tt_analyzer.selectedLeadingElMu.second ].Pt() )
        return false;

    return true;
};

void TTElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
};

void TTElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const TTAnalyzer& tt_analyzer = analyzers.get<TTAnalyzer>("tt");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    int selEl1 = tt_analyzer.selectedLeadingElEl.first;
    int selEl2 = tt_analyzer.selectedLeadingElEl.second;

    if( selEl1 >= 0 ){
        if( (electrons.p4[selEl1] + electrons.p4[selEl2]).M() > m_mll_cut)
            manager.pass_cut("ll_mass");
    }
}
