#pragma once

#include <cp3_llbb/TTAnalysis/interface/Types.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

namespace TTAnalysis {

  float DeltaEta(const myLorentzVector &v1, const myLorentzVector &v2){
    return abs(v1.Eta() - v2.Eta());
  }
  
  // Used by std::sort to sort jets according to decreasing b-tagging discriminant value
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

}
