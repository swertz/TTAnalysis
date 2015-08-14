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
            manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons");
            /*manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons");
            manager.new_category<MuElCategory>("muel", "Category with leading leptons as muon, electron");
            manager.new_category<ElMuCategory>("elmu", "Category with leading leptons as electron, muon");*/
        }

    private:

};

#endif
