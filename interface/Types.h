#pragma once

#include <utility>
#include <vector>
#include <map>

#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>

#include <cp3_llbb/TreeWrapper/interface/Resetter.h>

// Needed because of gcc bug when using typedef and std::map
#define myLorentzVector ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>

namespace TTAnalysis{

  enum class lepID{ L, M, T, V };
  enum class jetID{ L, T, TLV };
  enum class btWP{ L, M, T };

  struct BaseObject {
    BaseObject(myLorentzVector p4): p4(p4) {}
    BaseObject() {}

    myLorentzVector p4;

    virtual void clear(){
      p4.SetCoordinates(0,0,0,0);
    }  
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
      }
    Lepton() {}
    
    uint8_t idx; // stores index to electron/muon arrays
    uint8_t charge;
    bool isEl;
    bool isMu;
    std::vector<bool> lepID; // lepton IDs: loose-medium-tight(-veto)

    /*virtual void clear(){
      BaseObject::clear();
      idx = -1;
      charge = 0;
      isEl = false;
      isMu = false;
      lepID.clear();
    }*/
  };
  
  struct DiLepton: BaseObject {
    /*DiLepton(myLorentzVector p4, std::pair<int, int> idxs, std::pair<int, int> lidxs, bool isElEl, bool isElMu, bool isMuEl, bool isMuMu, bool isOS, bool isSF, float DR, float DEta, float DPhi):
      BaseObject(p4), idxs(idxs), lidxs(lidxs), isElEl(isElEl), isElMu(isElMu), isMuEl(isMuEl), isMuMu(isMuMu), isOS(isOS), isSF(isSF), DR(DR), DEta(DEta), DPhi(DPhi)
      {}
    DiLepton() {}*/
    
    std::pair<int, int> idxs; // stores indices to electron/muon arrays
    std::pair<int, int> lidxs; // stores indices to Lepton array
    bool isElEl, isElMu, isMuEl, isMuMu;
    bool isOS; // opposite sign
    bool isSF; // same flavour
    bool lepID_LL, lepID_LM, lepID_ML, lepID_LT, lepID_TL, lepID_MM, lepID_MT, lepID_TM, lepID_TT;
    float DR;
    float DEta;
    float DPhi;
    
    /*virtual void clear(){
      BaseObject::clear();
      idxs.first = -1;
      idxs.second = -1;
      lidxs.first = -1;
      lidxs.second = -1;
      isElEl = false;
      isElMu = false;
      isMuEl = false;
      isMuMu = false;
      isOS = false;
      isSF = false;
      DR = -1;
      DEta = -1;
      DPhi = -1;
    }*/
  };
  
  struct DiJet: BaseObject {
    DiJet(myLorentzVector p4, std::pair<int, int> idxs):
      BaseObject(p4), idxs(idxs)
      {}
    DiJet() {}
    
    std::pair<int, int> idxs; // stores indices to jets array
    std::pair<int, int> jidxs; // stores indices to TTAnalysis::Jet array
    float DR;
    float DEta;
    float DPhi;
    
    /*virtual void clear(){
      BaseObject::clear();
      idxs.first = -1;
      idxs.second = -1;
    }*/
  };

}

// Needed for TreeWrapper to handle the objects correctly

/*template <>
struct ResetterT<TTAnalysis::BaseObject>: Resetter {
    public:
        ResetterT(TTAnalysis::BaseObject& data): mdata(data) {}

        virtual void reset() {
            mdata.clear();
        }

    private:
        TTAnalysis::BaseObject& mdata;
};

template <typename T1, typename T2>
struct ResetterT<std::map<T1, T2>>: Resetter {
    public:
        ResetterT(std::map<T1, T2>& data): mdata(data) {}

        virtual void reset() {
            mdata.clear();
        }

    private:
        std::map<T1, T2>& mdata;
};
*/