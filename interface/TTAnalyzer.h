#ifndef TTANALYZER_H
#define TTANALYZER_H

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/Analyzer.h>

class TTAnalyzer: public Framework::Analyzer {
    public:
        TTAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config),
            
            m_electronIsoCut( config.getUntrackedParameter<double>("electronIsoCut", 10.) ),
            m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut") ),
            m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut") ),
            m_electron_veto_wp_name( config.getUntrackedParameter<std::string>("electron_veto_wp_name") ),
            m_electron_loose_wp_name( config.getUntrackedParameter<std::string>("electron_loose_wp_name") ),
            m_electron_medium_wp_name( config.getUntrackedParameter<std::string>("electron_medium_wp_name") ),
            m_electron_tight_wp_name( config.getUntrackedParameter<std::string>("electron_tight_wp_name") ),
            m_electron_selectedID_name( config.getUntrackedParameter<std::string>("electron_selectedID_name") ),
            
            m_muonIsoCut( config.getUntrackedParameter<double>("muonIsoCut") ),
            m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut") ),
            m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut") ),
            m_muon_selectedID_wp( config.getUntrackedParameter<std::string>("muon_selectedID_wp") ),
            
            m_jetEtaCut( config.getUntrackedParameter<double>("jetEtaCut") ),
            m_jetPtCut( config.getUntrackedParameter<double>("jetPtCut") )
        {
        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;

        BRANCH(selectedElectrons, std::vector<uint8_t>);
        BRANCH(isolatedElectrons, std::vector<uint8_t>);
        BRANCH(tightElectrons, std::vector<uint8_t>);
        BRANCH(looseElectrons, std::vector<uint8_t>);
        BRANCH(mediumElectrons, std::vector<uint8_t>);
        BRANCH(vetoElectrons, std::vector<uint8_t>);
        
        BRANCH(selectedMuons, std::vector<uint8_t>);
        BRANCH(isolatedMuons, std::vector<uint8_t>);
        BRANCH(tightMuons, std::vector<uint8_t>);
        BRANCH(mediumMuons, std::vector<uint8_t>);
        BRANCH(looseMuons, std::vector<uint8_t>);

        BRANCH(selectedLeadingMuMu, std::pair<int, int>);
        BRANCH(selectedLeadingMuEl, std::pair<int, int>);
        BRANCH(selectedLeadingElMu, std::pair<int, int>);
        BRANCH(selectedLeadingElEl, std::pair<int, int>);
        
        BRANCH(selectedJets, std::vector<uint8_t>);

    private:

        const float m_electronIsoCut, m_electronEtaCut, m_electronPtCut;
        const std::string m_electron_veto_wp_name;
        const std::string m_electron_loose_wp_name;
        const std::string m_electron_medium_wp_name;
        const std::string m_electron_tight_wp_name;
        const std::string m_electron_selectedID_name;

        const float m_muonIsoCut, m_muonEtaCut, m_muonPtCut;
        const std::string m_muon_selectedID_wp;

        const float m_jetEtaCut, m_jetPtCut;

        static inline bool muonIDAccessor(const MuonsProducer& muons, const uint8_t index, const std::string muonID){
            if(muonID == "loose")
                return muons.isLoose[index];
            
            if(muonID == "medium")
                return muons.isMedium[index];
            
            if(muonID == "tight")
                return muons.isTight[index];
            
            throw edm::Exception(edm::errors::NotFound, "Unknown muonID passed to analyzer");
        }
};

#endif
