#ifndef TTANALYZER_H
#define TTANALYZER_H

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/Analyzer.h>

#define myLorentzVector ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>

template <typename T1, typename T2>
struct ResetterT<std::map<T1, T2>>: Resetter {
    public:
        ResetterT(std::map<T1, T2>& data)
            : m_data(data) {
            }

        virtual void reset() {
            m_data.clear();
        }

    private:
        std::map<T1, T2>& m_data;
};

float DeltaEta(const myLorentzVector &v1, const myLorentzVector &v2){
  return abs(v1.Eta() - v2.Eta());
}

class jetBTagDiscriminantSorter {
  
  public:

    jetBTagDiscriminantSorter(const JetsProducer& jets, const std::string& taggerName): m_jetsProducer(jets), m_taggerName(taggerName) {}
    bool operator()(uint8_t idxJet1, uint8_t idxJet2){
      return m_jetsProducer.getBTagDiscriminant(idxJet1, m_taggerName) > m_jetsProducer.getBTagDiscriminant(idxJet2, m_taggerName);
    }

  private:

    const JetsProducer& m_jetsProducer;
    const std::string m_taggerName;
};

struct Lepton {
  myLorentzVector p4;
  uint8_t idx; // stores index to electron/muon arrays
  bool isMu;
  bool isEl;
};

struct DiLepton {
  myLorentzVector p4;
  std::pair<int, int> idxs; // stores indices to electron/muon arrays
  std::pair<int, int> lidxs; // stores indices to Lepton array
  bool isElEl;
  bool isElMu;
  bool isMuEl;
  bool isMuMu;
};

struct DiJet {
  myLorentzVector p4;
  std::pair<int, int> idxs;
};

class TTAnalyzer: public Framework::Analyzer {
    public:
        TTAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config),
            
            m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut", 20) ),
            m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut", 2.5) ),
            m_electronVetoIDName( config.getUntrackedParameter<std::string>("electronVetoIDName") ),
            m_electronLooseIDName( config.getUntrackedParameter<std::string>("electronLooseIDName") ),
            m_electronMediumIDName( config.getUntrackedParameter<std::string>("electronMediumIDName") ),
            m_electronTightIDName( config.getUntrackedParameter<std::string>("electronTightIDName") ),
            m_electronSelectedIDName( config.getUntrackedParameter<std::string>("electronSelectedIDName") ),
            
            m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut", 20) ),
            m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut", 2.4) ),
            m_muonBaseIsoCut( config.getUntrackedParameter<double>("muonBaseIsoCut", 0.2) ),
            m_muonSelectedIsoCut( config.getUntrackedParameter<double>("muonSelectedIsoCut", 0.2) ),
            m_muonSelectedID( config.getUntrackedParameter<std::string>("muonSelectedID") ),
            
            m_MllBaseCutSF( config.getUntrackedParameter<double>("MllBaseCutSF", 0.) ),
            m_MllBaseCutDF( config.getUntrackedParameter<double>("MllBaseCutDF", 0.) ),
            
            m_jetPtCut( config.getUntrackedParameter<double>("jetPtCut", 30) ),
            m_jetEtaCut( config.getUntrackedParameter<double>("jetEtaCut", 2.5) ),
            m_jetPUID( config.getUntrackedParameter<double>("jetPUID", -1000) ),
            m_jetDRleptonCut( config.getUntrackedParameter<double>("jetDRleptonCut", 0.3) ),
            m_jetID( config.getUntrackedParameter<std::string>("jetID") ),
            m_jetCSVv2Name( config.getUntrackedParameter<std::string>("jetCSVv2Name") ),
            m_jetCSVv2L( config.getUntrackedParameter<double>("jetCSVv2L", 0.605) ),
            m_jetCSVv2M( config.getUntrackedParameter<double>("jetCSVv2M", 0.89) ),
            m_jetCSVv2T( config.getUntrackedParameter<double>("jetCSVv2T", 0.97) ),
            
            m_hltDRCut( config.getUntrackedParameter<double>("hltDrCut", 100) ),
            m_hltPtCut( config.getUntrackedParameter<double>("hltPtCut", 999999) )
        {
        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;

        BRANCH(selectedElectrons, std::vector<uint8_t>); // index points to electrons array
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

        BRANCH(lepton_p4, std::vector<myLorentzVector>);
        BRANCH(lepton_idx, std::vector<uint8_t>); // index points to electrons/muons arrays
        BRANCH(lepton_isMu, std::vector<bool>);
        BRANCH(lepton_isEl, std::vector<bool>);
        
        BRANCH(diLepton_p4, myLorentzVector);
        BRANCH(diLepton_idxs, std::pair<int, int>); // index points to electrons/muons arrays
        BRANCH(diLepton_lidxs, std::pair<int, int>); // index points to lepton array
        BRANCH(diLepton_isMuMu, bool);
        BRANCH(diLepton_isMuEl, bool);
        BRANCH(diLepton_isElMu, bool);
        BRANCH(diLepton_isElEl, bool);
        BRANCH(diLepton_DR, float);
        BRANCH(diLepton_DPhi, float);
        BRANCH(diLepton_DEta, float);

        BRANCH(selectedJets, std::vector<uint8_t>); // all following indices point to jets array
        BRANCH(selectedJets_tightID, std::vector<uint8_t>);
        BRANCH(selectedJets_tightID_DRcut, std::vector<uint8_t>);
        BRANCH(selectedJets_diBJetExcluded, std::map<std::string, std::vector<uint8_t>>); // see note below on b-jet indices
        BRANCH(selectedBJets_CSVv2L, std::vector<uint8_t>);
        BRANCH(selectedBJets_CSVv2M, std::vector<uint8_t>);
        BRANCH(selectedBJets_CSVv2T, std::vector<uint8_t>);
        BRANCH(selectedNonBJets_CSVv2L, std::vector<uint8_t>);
        BRANCH(selectedNonBJets_CSVv2M, std::vector<uint8_t>);
        BRANCH(selectedNonBJets_CSVv2T, std::vector<uint8_t>);

        BRANCH(diJet_p4, myLorentzVector);
        BRANCH(diJet_idxs, std::pair<int, int>); // indices point to jets array
        BRANCH(diJet_DR, float);
        BRANCH(diJet_DPhi, float);
        BRANCH(diJet_DEta, float);

        // For b-jet, the string index corresponds to "WPs-ORDER" where WPs is "LL", "MM" or "TT",
        // and ORDER is "Pt" for Pt-ordered b-jet choice and "CSVv2" for CSVv2-ordered b-jet choice.
        
        BRANCH(diBJet_p4, std::map<std::string, myLorentzVector>);
        BRANCH(diBJet_idxs, std::map<std::string, std::pair<int, int>>); // indices point to jets array
        BRANCH(diBJet_DR, std::map<std::string, float>);
        BRANCH(diBJet_DPhi, std::map<std::string, float>);
        BRANCH(diBJet_DEta, std::map<std::string, float>);

        BRANCH(ll_jj_p4, myLorentzVector);
        BRANCH(ll_jj_DR, float);
        BRANCH(ll_jj_DPhi, float);
        BRANCH(ll_jj_DEta, float);
        BRANCH(lj_minDR, float);
        BRANCH(lj_noDRcut_minDR, float);

        BRANCH(ll_bb_p4, std::map<std::string, myLorentzVector>);
        BRANCH(ll_bb_DR, std::map<std::string, float>);
        BRANCH(ll_bb_DPhi, std::map<std::string, float>);
        BRANCH(ll_bb_DEta, std::map<std::string, float>);
        BRANCH(lb_minDR, std::map<std::string, float>);

        BRANCH(ll_jj_pfMET_p4, myLorentzVector);
        BRANCH(ll_pfMET_DPhi, float);
        BRANCH(jj_pfMET_DPhi, float);
        BRANCH(ll_jj_pfMET_DPhi, float);
        
        BRANCH(ll_jj_noHFMET_p4, myLorentzVector);
        BRANCH(ll_noHFMET_DPhi, float);
        BRANCH(jj_noHFMET_DPhi, float);
        BRANCH(ll_jj_noHFMET_DPhi, float);
        
        BRANCH(ll_bb_pfMET_p4, std::map<std::string, myLorentzVector>);
        BRANCH(bb_pfMET_DPhi, std::map<std::string, float>);
        BRANCH(ll_bb_pfMET_DPhi, std::map<std::string, float>);
        
        BRANCH(ll_bb_noHFMET_p4, std::map<std::string, myLorentzVector>);
        BRANCH(bb_noHFMET_DPhi, std::map<std::string, float>);
        BRANCH(ll_bb_noHFMET_DPhi, std::map<std::string, float>);

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

        const float m_jetPtCut, m_jetEtaCut, m_jetPUID, m_jetDRleptonCut;
        const std::string m_jetID, m_jetCSVv2Name;
        const float m_jetCSVv2L, m_jetCSVv2M, m_jetCSVv2T;

        const float m_hltDRCut, m_hltPtCut;

        std::vector<Lepton> m_leptons;
        DiLepton m_diLepton;
        DiJet m_diJet;
        std::map<std::string, DiJet> m_diBJet_PtChosen, m_diBJet_CSVv2Chosen;
        
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
