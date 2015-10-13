#pragma once

#include <utility>
#include <vector>
#include <map>
#include <limits>
#include <array>

#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>

#include <cp3_llbb/TreeWrapper/interface/Resetter.h>

#define BRANCHNOCLEAR(NAME, ...) __VA_ARGS__& NAME = tree[#NAME].write<__VA_ARGS__>(false)

// Needed because of gcc bug when using typedef and std::map
#define myLorentzVector ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>

namespace TTAnalysis{

  namespace LepID{
    enum LepID{ L, M, T, V, Count };
    // Ugly way to allow iterating over all items in the enumeration ( for(const LepID::LepID& id: LepID::it) )
    const std::array<LepID, Count> it = { L, M, T, V };
  }
  
  namespace LepLepID{
    enum LepLepID{ LL, LM, ML, LT, TL, MM, MT, TM, TT, Count };
    const std::array<LepLepID, Count> it = { LL, LM, ML, LT, TL, MM, MT, TM, TT };
  }
  
  namespace JetID{
    enum JetID{ L, T, TLV, Count };
    const std::array<JetID, Count> it = { L, T, TLV };
  }
  
  namespace BWP{
    enum BWP{ L, M, T, Count };
    const std::array<BWP, Count> it = { L, M, T };
  }
  
  namespace BBWP{
    enum BBWP{ LL, LM, ML, LT, TL, MM, MT, TM, TT, Count };
    const std::array<BBWP, Count> it = { LL, LM, ML, LT, TL, MM, MT, TM, TT };
  }

  struct BaseObject {
    BaseObject(myLorentzVector p4): p4(p4) {}
    BaseObject() {}

    myLorentzVector p4;
  };

  struct Lepton: BaseObject {
    Lepton(myLorentzVector p4, uint8_t idx, bool isEl, bool isMu, bool isLoose = false, bool isMedium = false, bool isTight = false, bool isVeto = false):
      BaseObject(p4), idx(idx), charge(charge), isEl(isEl), isMu(isMu)
      {
        lepID.push_back(isLoose);
        lepID.push_back(isMedium);
        lepID.push_back(isTight);
        if(isEl)
          lepID.push_back(isVeto);
        else
          lepID.push_back(false);
      }
    Lepton() {}
    
    uint8_t idx; // stores index to electron/muon arrays
    uint8_t charge;
    bool isEl;
    bool isMu;
    std::vector<bool> lepID; // lepton IDs: loose-medium-tight(-veto)
  };
  
  struct DiLepton: BaseObject {
    DiLepton():
      lepIDs(LepLepID::Count, false)
      {}
    
    std::pair<int, int> idxs; // stores indices to electron/muon arrays
    std::pair<int, int> lidxs; // stores indices to Lepton array
    bool isElEl, isElMu, isMuEl, isMuMu;
    bool isOS; // opposite sign
    bool isSF; // same flavour
    std::vector<bool> lepIDs;
    float DR;
    float DEta;
    float DPhi;
  };
  
  struct DiJet: BaseObject {
    DiJet():
      minDRjl_lepIDs(LepLepID::Count, std::numeric_limits<float>::max()),
      CSVv2_WPs(BBWP::Count, false)
      {}
    
    std::pair<int, int> idxs; // stores indices to jets array
    //std::pair<int, int> jidxs; // stores indices to TTAnalysis::Jet array (NOT implemented... necessary??)
    float minDRjl_lepIDs;
    std::vector<bool> CSVv2_WPs;
    float DR;
    float DEta;
    float DPhi;
  };

  struct DiLepDiJet: BaseObject {
    DiLepDiJet(DiLepton& diLepton, const int diLepIdx, DiJet& diJet, const int diJetIdx):
      BaseObject(diLepton.p4 + diJet.p4),
      diLepton(diLepton),
      diLepIdx(diLepIdx),
      diJet(diJet),
      diJetIdx(diJetIxd),
      DR_ll_jj( ROOT::Math::VectorUtil::DeltaR(diLepton.p4, diJet.p4) ),
      DEta_ll_jj( ROOT::Math::DeltaEta(diLepton.p4, diJet.p4) ),
      DPhi_ll_jj( ROOT::Math::VectorUtil::DeltaPhi(diLepton.p4, diJet.p4) )
      {}

    const DiLepton& diLepton;
    const int diLepIdx;
    const DiJet& diJet;
    const int diJetIdx;

    float DR_ll_jj, DEta_ll_jj, DPhi_ll_jj;
    
    float minDRjl, maxDRjl;
    float minDEtajl, maxDEtajl;
    float minDPhijl, maxDPhijl;
  };

  struct DiLepDiJetMet: DiLepDiJet {
    DiLepDiJetMet(DiLepDiJet& diLepDiJet, uint8_t diLepDiJetIdx, LorentzVector& MetP4, bool hasNoHFMet = false):
      DiLepDiJet(diLepDiJet.diLepton, diLepDiJet.diLepIdx, diLepDiJet.diJet, diLepDiJet.diJetIdx),
      diLepDiJetIdx(diLepDiJetIdx),
      minDRjl(DiLepDiJet.minDRjl),
      maxDRjl(DiLepDiJet.maxDRjl),
      minDEtajl(DiLepDiJet.minDEtajl),
      maxDEtajl(DiLepDiJet.maxDEtajl),
      minDPhijl(DiLepDiJet.minDPhijl),
      maxDPhijl(DiLepDiJet.maxDPhijl),
      hasNoHFMet(hasNoHFMet)
    {
      p4 += MetP4;

      DR_ll_Met = ROOT::Math::VectorUtil::DeltaR(diLepton.p4, MetP4);
      DR_jj_Met = ROOT::Math::VectorUtil::DeltaR(diJet.p4, MetP4);
      
      DEta_ll_Met = DeltaEta(diLepton.p4, MetP4);
      DEta_jj_Met = DeltaEta(diJet.p4, MetP4);
      
      DPhi_ll_Met = ROOT::Math::VectorUtil::DeltaPhi(diLepton.p4, MetP4);
      DPhi_jj_Met = ROOT::Math::VectorUtil::DeltaPhi(diJet.p4, MetP4);
      
      DR_lljj_Met = ROOT::Math::VectorUtil::DeltaR(diLepton.p4 + diJet.p4, MetP4);
      DEta_lljj_Met = ROOT::Math::DeltaEta(diLepton.p4 + diJet.p4, MetP4);
      DPhi_lljj_Met = ROOT::Math::VectorUtil::DeltaPhi(diLepton.p4 + diJet.p4, MetP4);
    }

    uint8_t diLepDiJetIdx:
    bool hasNoHFMet;

    float DR_ll_Met, DR_jj_Met;
    float DEta_ll_Met, DEta_jj_Met;
    float DPhi_ll_Met, DPhi_jj_Met;

    float DR_lljj_Met;
    float DEta_lljj_Met;
    float DPhi_lljj_Met;

    float minDR_l_Met, minDR_j_Met;
    float maxDR_l_Met, maxDR_j_Met;
    float minDEta_l_Met, minDEta_j_Met;
    float maxDEta_l_Met, maxDEta_j_Met;
    float minDPhi_l_Met, minDPhi_j_Met;
    float maxDPhi_l_Met, maxDPhi_j_Met;
  };

}
