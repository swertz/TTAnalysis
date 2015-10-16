#pragma once

#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/TTAnalysis/interface/Types.h>

namespace TTAnalysis{

class DileptonCategory: public Category {
  public:
    virtual void configure(const edm::ParameterSet& conf) override {
      m_MllCutSF = conf.getUntrackedParameter<double>("MllCutSF", 20);
      m_MllCutDF = conf.getUntrackedParameter<double>("MllCutDF", 0);
      m_MllZVetoCutLow = conf.getUntrackedParameter<double>("MllZVetoCutLow", 86);
      m_MllZVetoCutHigh = conf.getUntrackedParameter<double>("MllZVetoCutHigh", 116);
      m_HLTDoubleMuon = conf.getUntrackedParameter<std::vector<std::string>>("HLTDoubleMuon");
      m_HLTDoubleEG = conf.getUntrackedParameter<std::vector<std::string>>("HLTDoubleEG");
      m_HLTMuonEG = conf.getUntrackedParameter<std::vector<std::string>>("HLTMuonEG");
    }

    DileptonCategory():
      baseStrCategory("Category"),
      baseStrExtraDiLeptonVeto("ExtraDiLeptonVeto"),
      baseStrDiLeptonTriggerMatch("DiLeptonTriggerMatch"),
      baseStrMllCut("Mll"),
      baseStrMllZVetoCut("MllZVeto"),
      baseStrDiLeptonIsOS("DiLeptonIsOS")
      {}

  protected:
    float m_MllCutSF, m_MllCutDF, m_MllZVetoCutLow, m_MllZVetoCutHigh;

    std::vector<std::string> m_HLTDoubleMuon;
    std::vector<std::string> m_HLTDoubleEG;
    std::vector<std::string> m_HLTMuonEG;

    std::string baseStrCategory;
    std::string baseStrExtraDiLeptonVeto;
    std::string baseStrDiLeptonTriggerMatch;
    std::string baseStrMllCut;
    std::string baseStrMllZVetoCut;
    std::string baseStrDiLeptonIsOS;

    enum class HLT { DoubleMuon, DoubleEG, MuonEG };

    // Check that the hlt objects at indices hltIdx1, hltIdx2 have fired one and the same of the trigger paths
    // in the group specified by pathGroup.
    bool checkHLT(const HLTProducer& hlt, uint8_t hltIdx1, uint8_t hltIdx2, HLT pathGroup) const {
      const std::vector<std::string>* chosenPathGroup(nullptr);

      switch(pathGroup){
        
        case HLT::DoubleMuon:
          chosenPathGroup = &m_HLTDoubleMuon;
          break;

        case HLT::DoubleEG:
          chosenPathGroup = &m_HLTDoubleEG;
          break;

        case HLT::MuonEG:
          chosenPathGroup = &m_HLTMuonEG;
          break;

        default:
          break;
      
      }

      if( !chosenPathGroup || hltIdx1 >= hlt.object_paths.size() || hltIdx2 >= hlt.object_paths.size() )
        return false;

      std::vector<std::string> matchedTriggersObj1, matchedTriggersObj2, commonMatchedTriggers;

      std::set_intersection( (*chosenPathGroup).begin(), (*chosenPathGroup).end(), hlt.object_paths[hltIdx1].begin(), hlt.object_paths[hltIdx1].end(), std::back_inserter(matchedTriggersObj1));
      std::set_intersection( (*chosenPathGroup).begin(), (*chosenPathGroup).end(), hlt.object_paths[hltIdx2].begin(), hlt.object_paths[hltIdx2].end(), std::back_inserter(matchedTriggersObj2));

      std::set_intersection( matchedTriggersObj1.begin(), matchedTriggersObj1.end(), matchedTriggersObj2.begin(), matchedTriggersObj2.end(), std::back_inserter(commonMatchedTriggers));

      return commonMatchedTriggers.size() > 0;
    }
      
};

class ElElCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class ElMuCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

 class MuElCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class MuMuCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

}

