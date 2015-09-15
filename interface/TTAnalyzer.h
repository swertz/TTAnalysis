#ifndef TTANALYZER_H
#define TTANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>

class TTAnalyzer: public Framework::Analyzer {
    public:
        TTAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config) {

        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;

    private:

};

#endif
