#include <cp3_llbb/TTAnalysis/interface/TTAnalyzer.h>

#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

void TTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers) {


    //std::cout << " Starting an analysis " << std::endl;

    const JetsProducer& jets = producers.get<JetsProducer>("jets");
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");

    BRANCH(bb, LorentzVector);
    BRANCH(ee, LorentzVector);
    BRANCH(mm, LorentzVector);

    bb.SetPxPyPzE(0.,0.,0.,0.);
    ee.SetPxPyPzE(0.,0.,0.,0.);
    mm.SetPxPyPzE(0.,0.,0.,0.);

    if (jets.p4.size() > 1) {
        bb = jets.p4.at(0)+jets.p4.at(1);
    }

    if (electrons.p4.size() > 1) {
        ee = electrons.p4.at(0)+electrons.p4.at(1);
    }
    
    if (muons.p4.size() > 1) {
        mm = muons.p4.at(0)+muons.p4.at(1);
    }

}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, TTAnalyzer, "tt_analyzer");
