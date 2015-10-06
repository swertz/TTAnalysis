#ifndef TTANALYZER_H
#define TTANALYZER_H

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/Analyzer.h>

struct Lepton {
  LorentzVector p4;
  uint8_t idx;
  bool isMu;
  bool isEl;
}

struct DiLepton {
  LorentzVector p4;
  std::pair<int, int> idxs;
  bool isElEl;
  bool isElMu;
  bool isMuEl;
  bool isMuMu;
}

struct DiJet {
  LorentzVector p4;
  std::pair<uint8_t, uint8_t> idxs;
}

class TTAnalyzer: public Framework::Analyzer {
    public:
        TTAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config),
            
            m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut") ),
            m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut") ),
            m_electronVetoIDName( config.getUntrackedParameter<std::string>("electronVetoIDName") ),
            m_electronLooseIDName( config.getUntrackedParameter<std::string>("electronLooseIDName") ),
            m_electronMediumIDName( config.getUntrackedParameter<std::string>("electronMediumIDName") ),
            m_electronTightIDName( config.getUntrackedParameter<std::string>("electronTightIDName") ),
            m_electronSelectedIDName( config.getUntrackedParameter<std::string>("electronSelectedIDName") ),
            
            m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut") ),
            m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut") ),
            m_muonBaseIsoCut( config.getUntrackedParameter<double>("muonBaseIsoCut") ),
            m_muonSelectedIsoCut( config.getUntrackedParameter<double>("muonSelectedIsoCut") ),
            m_muonSelectedID( config.getUntrackedParameter<std::string>("muonSelectedID") ),
            
            m_MllBaseCutSF( config.getUntrackedParameter<double>("MllBaseCutSF", 0.) ),
            m_MllBaseCutDF( config.getUntrackedParameter<double>("MllBaseCutDF", 0.) ),
            
            m_jetPtCut( config.getUntrackedParameter<double>("jetPtCut") )
            m_jetEtaCut( config.getUntrackedParameter<double>("jetEtaCut") ),
            m_jetPUID( config.getUntrackedParameter<double>("jetPUID", -1000) ),
            m_jetID( config.getUntrackedParameter<std::string>("jetID") ),
            
            m_hltEtaCut( config.getUntrackedParameter<double>("hltEtaCut") ),
            m_hltPtCut( config.getUntrackedParameter<double>("hltPtCut", 999999) ),
        {
        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;

        BRANCH(selectedElectrons, std::vector<uint8_t>); // index points to electrons array
        BRANCH(isolatedElectrons, std::vector<uint8_t>);
        BRANCH(tightElectrons, std::vector<uint8_t>);
        BRANCH(looseElectrons, std::vector<uint8_t>);
        BRANCH(mediumElectrons, std::vector<uint8_t>);
        BRANCH(vetoElectrons, std::vector<uint8_t>);
        
        BRANCH(selectedMuons, std::vector<uint8_t>); // index points to muons array
        BRANCH(tightMuons, std::vector<uint8_t>);
        BRANCH(mediumMuons, std::vector<uint8_t>);
        BRANCH(looseMuons, std::vector<uint8_t>);

        BRANCH(leadingSelectedMuMu, std::pair<int, int>); // indices point to electrons/muons arrays
        BRANCH(leadingSelectedMuEl, std::pair<int, int>);
        BRANCH(leadingSelectedElMu, std::pair<int, int>);
        BRANCH(leadingSelectedElEl, std::pair<int, int>);
        
        BRANCH(selectedJets, std::vector<uint8_t>);
        BRANCH(selectedJets_diBJetExcluded, std::vector<uint8_t>);
        BRANCH(selectedBJets_CSVv2L, std::vector<uint8_t>);
        BRANCH(selectedBJets_CSVv2M, std::vector<uint8_t>);
        BRANCH(selectedBJets_CSVv2T, std::vector<uint8_t>);

        BRANCH(lepton_p4, std::vector<LorentzVector>);
        BRANCH(lepton_idx, std::vector<uint8_t>); // index points to electrons/muons arrays
        BRANCH(lepton_isMu, std::vector<bool>);
        BRANCH(lepton_isEl, std::vector<bool>);
        
        BRANCH(diLepton_p4, LorentzVector);
        BRANCH(diLepton_idxs, std::pair<uint8_t, uint8_t>); // index points to lepton array
        BRANCH(diLepton_isMuMu, bool);
        BRANCH(diLepton_isMuEl, bool);
        BRANCH(diLepton_isElMu, bool);
        BRANCH(diLepton_isElEl, bool);
        BRANCH(diLepton_DR, float);
        BRANCH(diLepton_DPhi, float);
        BRANCH(diLepton_DEta, float);

        BRANCH(diJet_p4, LorentzVector);
        BRANCH(diJet_idxs, std::pair<uint8_t, uint8_t>);
        BRANCH(diJet_DR, float);
        BRANCH(diJet_DPhi, float);
        BRANCH(diJet_DEta, float);

        BRANCH(diBJet_p4, LorentzVector);
        BRANCH(diBJet_idxs, std::pair<uint8_t, uint8_t>);
        BRANCH(diBJet_DR, float);
        BRANCH(diBJet_DPhi, float);
        BRANCH(diBJet_DEta, float);

        BRANCH(ll_jj_p4, LorentzVector);
        BRANCH(ll_jj_DR, float);
        BRANCH(ll_jj_DPhi, float);
        BRANCH(ll_jj_DEta, float);
        BRANCH(lj_minDR, float);

        BRANCH(ll_bb_p4, LorentzVector);
        BRANCH(ll_bb_DR, float);
        BRANCH(ll_bb_DPhi, float);
        BRANCH(ll_bb_DEta, float);
        BRANCH(lb_minDR, float);

        BRANCH(ll_jj_pfMET_p4, LorentzVector);
        BRANCH(ll_bb_pfMET_p4, LorentzVector);

    private:

        const float m_electronPtCut, m_electronEtaCut;
        const std::string m_electronVetoIDName;
        const std::string m_electronLooseIDName;
        const std::string m_electronMediumIDName;
        const std::string m_electronTightIDName;
        const std::string m_electronSelectedIDName;

        const float m_muonPtCut, m_muonEtaCut, m_muonBaseIsoCut, m_muonSelectedIsoCut;
        const std::string m_muonSelectedID;

        const float m_MllBaseCutSF, m_MllBaseCutDF;

        const float m_jetPtCut, m_jetEtaCut, m_jetPUID;
        const std::string m_jetID;
        const float m_bJetCSVv2L, m_bJetCSVv2M, m_bJetCSVv2T;

        const float m_hltEtaCut, m_hltPtCut;

        std::vector<Lepton> m_leptons;
        DiLepton m_diLepton;
        DiJet m_diJet, m_diBJet_PtChosen, m_diBJet_DiscrChosen;
        
        static inline bool muonIDAccessor(const MuonsProducer& muons, const uint8_t index, const std::string& muonID){
            if(index >= muons.p4.size())
              throw edm::Exception(edm::errors::StdException, "Invalid muon index passed to ID accessor");

            if(muonID == "loose")
                return muons.isLoose[index];
            
            if(muonID == "medium")
                return muons.isMedium[index];
            
            if(muonID == "tight")
                return muons.isTight[index];
            
            throw edm::Exception(edm::errors::NotFound, "Unknown muonID passed to analyzer");
        }
        
        static inline bool jetIDAccessor(const JetsProducer& jets, const uint8_t index, const std::string& jetID){
            if(index >= jets.p4.size())
              throw edm::Exception(edm::errors::StdException, "Invalid jet index passed to ID accessor");
            
            if(jetID == "loose")
                return jets.passLooseID[index];
            
            if(jetID == "tight")
                return jets.passTightID[index];
            
            if(jetID == "tightLeptonVeto")
                return jets.passTightLeptonVetoID[index];
            
            throw edm::Exception(edm::errors::NotFound, "Unknown jetID passed to analyzer");
        }
};

#endif
