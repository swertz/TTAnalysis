#ifndef TTANALYZER_H
#define TTANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>

#include <cp3_llbb/TTAnalysis/interface/Categories.h>

class TTAnalyzer: public Framework::Analyzer {
    public:
        TTAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config) {

        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&) override;

        virtual void registerCategories(CategoryManager& manager) {
            manager.new_category("mumu", "Category with leading leptons as two muons", &mumu_category);
            manager.new_category("elel", "Category with leading leptons as two electrons", &elel_category);
            manager.new_category("muel", "Category with leading leptons as muon, electron", &muel_category);
            manager.new_category("elmu", "Category with leading leptons as electron, muon", &elmu_category);
        }

    private:
        MuMuCategory mumu_category;
        ElElCategory elel_category;
        MuElCategory muel_category;
        ElMuCategory elmu_category;

};

#endif
