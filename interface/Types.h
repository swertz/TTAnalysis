#pragma once

#include <utility>
#include <vector>
#include <map>

//#include <Math/Vector4D.h>
#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>

#include <cp3_llbb/TreeWrapper/interface/Resetter.h>

// Needed because of gcc bug when using typedef and std::map
#define myLorentzVector ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>

namespace TTAnalysis{

  struct BaseObject {
    BaseObject(myLorentzVector p4): p4(p4) {}
    BaseObject() {}

    myLorentzVector p4;

    void clear(void){
      p4.SetCoordinates(0,0,0,0);
    }
  };

  struct Lepton: BaseObject {
    Lepton(myLorentzVector p4, uint8_t idx, bool isMu, bool isEl):
      BaseObject(p4), idx(idx), isMu(isMu), isEl(isEl)
      {}
    Lepton() {}
    
    uint8_t idx; // stores index to electron/muon arrays
    bool isMu;
    bool isEl;

    void clear(void){
      BaseObject::clear();
      idx = 0;
      isMu = false;
      isEl = false;
    }
  };
  
  struct DiLepton: BaseObject {
    DiLepton(myLorentzVector p4, std::pair<int, int> idxs, std::pair<int, int> lidxs, bool isElEl, bool isElMu, bool isMuEl, bool isMuMu):
      BaseObject(p4), idxs(idxs), lidxs(lidxs), isElEl(isElEl), isElMu(isElMu), isMuEl(isMuEl), isMuMu(isMuMu)
      {}
    DiLepton() {}
    
    std::pair<int, int> idxs; // stores indices to electron/muon arrays
    std::pair<int, int> lidxs; // stores indices to Lepton array
    bool isElEl;
    bool isElMu;
    bool isMuEl;
    bool isMuMu;
    
    void clear(void){
      BaseObject::clear();
      idxs.first = -1;
      idxs.second = -1;
      lidxs.first = -1;
      lidxs.second = -1;
      isElEl = false;
      isElMu = false;
      isMuEl = false;
      isMuMu = false;
    }
  };
  
  struct DiJet: BaseObject {
    DiJet(myLorentzVector p4, std::pair<int, int> idxs):
      BaseObject(p4), idxs(idxs)
      {}
    DiJet() {}
    
    std::pair<int, int> idxs;
    
    void clear(void){
      BaseObject::clear();
      idxs.first = -1;
      idxs.second = -1;
    }
  };

}

// Needed for TreeWrapper to handle the objects correctly

template <>
struct ResetterT<TTAnalysis::BaseObject>: Resetter {
    public:
        ResetterT(TTAnalysis::BaseObject& data)
            : mdata(data) {
            }

        virtual void reset() {
            mdata.clear();
        }

    private:
        TTAnalysis::BaseObject& mdata;
};

template <typename T1, typename T2>
struct ResetterT<std::map<T1, T2>>: Resetter {
    public:
        ResetterT(std::map<T1, T2>& data)
            : mdata(data) {
            }

        virtual void reset() {
            mdata.clear();
        }

    private:
        std::map<T1, T2>& mdata;
};

