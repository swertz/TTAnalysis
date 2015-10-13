#pragma once

#include <utility>
#include <vector>
#include <map>
#include <limits>
#include <array>

#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>

#include <cp3_llbb/TreeWrapper/interface/Resetter.h>

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
    DiLepDiJet(DiLepton diLepton, const int diLepIdx, DiJet diJet, const int diJetIdx):
      diLepton(diLepton),
      diLepIdx(diLepIdx),
      diJet(diJet),
      diJetIdx(diJetIxd)
      {}

    const DiLepton& diLepton;
    const int diLepIdx;
    const DiJet& diJet;
    const int diJetIdx;

    float minDRjl, maxDRjl;
    float minDEtajl, maxDEtajl;
    float minDPhijl, maxDPhijl;
}
