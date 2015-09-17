#ifndef TTANALYZER_H
#define TTANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>

class TTAnalyzer: public Framework::Analyzer {
    public:
        TTAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config),
            m_electronIsoCut( config.getUntrackedParameter<double>("electronIsoCut") ),
            m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut") ),
            m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut") ),
            m_muonIsoCut( config.getUntrackedParameter<double>("muonIsoCut") ),
            m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut") ),
            m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut") )
        {
            m_electron_loose_wp_name = config.getUntrackedParameter<std::string>("electrons_loose_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-loose");
            m_electron_tight_wp_name = config.getUntrackedParameter<std::string>("electrons_tight_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-tight");
        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;

        BRANCH(isolatedElectrons, std::vector<int>);
        BRANCH(isolatedMuons, std::vector<int>);
        BRANCH(tightElectrons, std::vector<int>);
        BRANCH(tightMuons, std::vector<int>);
        BRANCH(looseElectrons, std::vector<int>);
        BRANCH(looseMuons, std::vector<int>);
        BRANCH(selectedElectrons, std::vector<int>);
        BRANCH(selectedMuons, std::vector<int>);

        BRANCH(selectedLeadingMuMu, std::pair<int, int>);
        BRANCH(selectedLeadingMuEl, std::pair<int, int>);
        BRANCH(selectedLeadingElMu, std::pair<int, int>);
        BRANCH(selectedLeadingElEl, std::pair<int, int>);

    private:

        const float m_electronIsoCut, m_electronEtaCut, m_electronPtCut;
        const float m_muonIsoCut, m_muonEtaCut, m_muonPtCut;

        std::string m_electron_loose_wp_name;
        std::string m_electron_tight_wp_name;

};

#endif
