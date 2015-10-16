#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/TTAnalysis/interface/TTDileptonCategories.h>
#include <cp3_llbb/TTAnalysis/interface/TTAnalyzer.h>

using namespace TTAnalysis;

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
  return electrons.p4.size() >= 2;
}

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton MuMu among all ID pairs is found, keep event in this category
  bool oneDiLepton(false);

  for(const LepLepID::LepLepID& id: LepLepID::it) {
    if(tt.diLeptons_LepIDs[id].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id][0]];
      if(m_diLepton.isElEl){
        oneDiLepton = true;
        break;
      }
    }
  }

  return oneDiLepton;
}

void ElElCategory::register_cuts(CutManager& manager) {
  
  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
    manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
    manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
    manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
    manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
    manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);
  }

}

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    if(tt.diLeptons_LepIDs[id.first].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id.first][0]];
      
      if(m_diLepton.isElEl) {
        manager.pass_cut(baseStrCategory + postFix);

        if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
          // We have fired a trigger. Now, check that it is actually a DoubleEG trigger
          if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::DoubleEG) )
            manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
        }
        
        if(m_diLepton.p4.M() > m_MllCutSF)
          manager.pass_cut(baseStrMllCut + postFix);
        
        if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
          manager.pass_cut(baseStrMllZVetoCut + postFix);
        
        if(m_diLepton.isOS)
          manager.pass_cut(baseStrDiLeptonIsOS + postFix);
      }
    
    } else if(tt.diLeptons_LepIDs[id.first].size() >= 2) {
      manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
    }
  }

}

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
  return electrons.p4.size() >= 1 && muons.p4.size() >= 1;
}

bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton ElMu among all ID pairs is found, keep event in this category
  bool oneDiLepton(false);

  for(const LepLepID::LepLepID& id: LepLepID::it) {
    if(tt.diLeptons_LepIDs[id].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id][0]];
      if(m_diLepton.isElMu){
        oneDiLepton = true;
        break;
      }
    }
  }

  return oneDiLepton;
}

void ElMuCategory::register_cuts(CutManager& manager) {
  
  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
    manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
    manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
    manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
    manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
    manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);
  }

}

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    if(tt.diLeptons_LepIDs[id.first].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id.first][0]];
      
      if(m_diLepton.isElEl) {
        manager.pass_cut(baseStrCategory + postFix);

        if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
          // We have fired a trigger. Now, check that it is actually a MuonEG trigger
          if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::MuonEG) )
            manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
        }
        
        if(m_diLepton.p4.M() > m_MllCutDF)
          manager.pass_cut(baseStrMllCut + postFix);
        
        if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
          manager.pass_cut(baseStrMllZVetoCut + postFix);
        
        if(m_diLepton.isOS)
          manager.pass_cut(baseStrDiLeptonIsOS + postFix);
      }
    
    } else if(tt.diLeptons_LepIDs[id.first].size() >= 2) {
      manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
    }
  }

}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****
bool MuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
  return electrons.p4.size() >= 1 && muons.p4.size() >= 1;
}

bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton MuEl among all ID pairs is found, keep event in this category
  bool oneDiLepton(false);

  for(const LepLepID::LepLepID& id: LepLepID::it) {
    if(tt.diLeptons_LepIDs[id].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id][0]];
      if(m_diLepton.isMuEl){
        oneDiLepton = true;
        break;
      }
    }
  }

  return oneDiLepton;
}

void MuElCategory::register_cuts(CutManager& manager) {
  
  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
    manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
    manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
    manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
    manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
    manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);
  }

}

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    if(tt.diLeptons_LepIDs[id.first].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id.first][0]];
      
      if(m_diLepton.isElEl) {
        manager.pass_cut(baseStrCategory + postFix);

        if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
          // We have fired a trigger. Now, check that it is actually a MuonEG trigger
          if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::DoubleEG) )
            manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
        }
        
        if(m_diLepton.p4.M() > m_MllCutDF)
          manager.pass_cut(baseStrMllCut + postFix);
        
        if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
          manager.pass_cut(baseStrMllZVetoCut + postFix);
        
        if(m_diLepton.isOS)
          manager.pass_cut(baseStrDiLeptonIsOS + postFix);
      }
    
    } else if(tt.diLeptons_LepIDs[id.first].size() >= 2) {
      manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
    }
  }

}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
  return muons.p4.size() >= 2;
}

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton MuMu among all ID pairs is found, keep event in this category
  bool oneDiLepton(false);

  for(const LepLepID::LepLepID& id: LepLepID::it) {
    if(tt.diLeptons_LepIDs[id].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id][0]];
      if(m_diLepton.isMuMu){
        oneDiLepton = true;
        break;
      }
    }
  }

  return oneDiLepton;
}

void MuMuCategory::register_cuts(CutManager& manager) {
  
  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
    manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
    manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
    manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
    manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
    manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);
  }

}

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(auto const& id: LepLepID::map){
    std::string postFix("_");
    postFix += id.second;
    
    if(tt.diLeptons_LepIDs[id.first].size() >= 1) {
      const DiLepton& m_diLepton = tt.diLeptons[tt.diLeptons_LepIDs[id.first][0]];
      
      if(m_diLepton.isElEl) {
        manager.pass_cut(baseStrCategory + postFix);

        if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
          // We have fired a trigger. Now, check that it is actually a DoubleMuon trigger
          if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::DoubleMuon) )
            manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
        }
        
        if(m_diLepton.p4.M() > m_MllCutSF)
          manager.pass_cut(baseStrMllCut + postFix);
        
        if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
          manager.pass_cut(baseStrMllZVetoCut + postFix);
        
        if(m_diLepton.isOS)
          manager.pass_cut(baseStrDiLeptonIsOS + postFix);
      }
    
    } else if(tt.diLeptons_LepIDs[id.first].size() >= 2) {
      manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
    }
  }

}

