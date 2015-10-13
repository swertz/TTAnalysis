#include <cp3_llbb/TTAnalysis/interface/Types.h>
#include <cp3_llbb/TTAnalysis/interface/Tools.h>
#include <cp3_llbb/TTAnalysis/interface/TTAnalyzer.h>
#include <cp3_llbb/TTAnalysis/interface/NoZTTDileptonCategories.h>
#include <cp3_llbb/TTAnalysis/interface/TTDileptonCategories.h>

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>

// To access VectorUtil::DeltaR() more easily
using namespace ROOT::Math;

using namespace TTAnalysis;

void TTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const CategoryManager& categories) {

    ///////////////////////////
    //       ELECTRONS       //
    ///////////////////////////

    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    for(uint8_t ielectron = 0; ielectron < electrons.p4.size(); ielectron++){
      if( electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut ){
        
        leptons.push_back( 
          Lepton(
            electron.p4[ielectron], 
            ielectron, 
            electron.charge[ielectron], 
            true, false,
            electrons.ids[ielectron][m_electronLooseIDName],
            electrons.ids[ielectron][m_electronMediumIDName],
            electrons.ids[ielectron][m_electronTightIDName],
            electrons.ids[ielectron][m_electronVetoIDName]
          )
        );
        
        if( electrons.ids[ielectron][m_electronVetoIDName] )
          vetoElectrons.push_back(ielectron);
        if( electrons.ids[ielectron][m_electronLooseIDName] )
          looseElectrons.push_back(ielectron);
        if( electrons.ids[ielectron][m_electronMediumIDName] )
          mediumElectrons.push_back(ielectron);
        if( electrons.ids[ielectron][m_electronTightIDName] )
          tightElectrons.push_back(ielectron);
      }
    }

    ///////////////////////////
    //       MUONS           //
    ///////////////////////////
    
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");

    for(uint8_t imuon = 0; imuon < muons.p4.size(); imuon++){
      if(muons.relativeIsoR04_withEA[imuon] < m_muonBaseIsoCut && muons.p4[imuon].Pt() > m_muonPtCut && abs(muons.p4[imuon].Eta()) < m_muonEtaCut ){
        
        leptons.push_back(
          Lepton(
            muon.p4[imuon], 
            imuon,
            muon.charge[imuon], 
            false, true,
            muons.isLoose[imuon],
            muons.isMedium[imuon],
            muons.isTight[imuon]
          )
        );

        if( muons.isLoose[imuon] )
          looseMuons.push_back(imuon);
        if( muons.isMedium[imuon] )
          mediumMuons.push_back(imuon);
        if( muons.isTight[imuon] )
          tightMuons.push_back(imuon);
      }
    }

    // Sort the leptons vector according to Pt
    std::sort(leptons.begin(), leptons.end(), [](const Lepton& a, const Lepton &b){ return b.p4.Pt() > a.p4.Pt(); });

    ///////////////////////////
    //       DILEPTONS       //
    ///////////////////////////

    for(uint8_t i1 = 0; i1 < leptons.size(); i1++){
      for(uint8_t i2 = i1 + 1; i2 < leptons.size(); i2++){
        const Lepton& l1 = leptons[i1];
        const Lepton& l2 = leptons[i2];

        DiLepton m_diLepton();

        m_diLepton.p4 = l1.p4 + l2.p4; 
        m_diLepton.idxs = std::make_pair(l1.idx; l2.idx); 
        m_diLepton.lidxs = std::make_pair(i1; i2); 
        m_diLepton.isElEl = l1.isEl && l2.isEl;
        m_diLepton.isElMu = l1.isEl && l2.isMu;
        m_diLepton.isMuEl = l1.isMu && l2.isEl;
        m_diLepton.isMuMu = l1.isMu && l2.isMu;
        m_diLepton.isOS = l1.charge != l2.charge;
        m_diLepton.isSF = m_diLepton.isElEl || m_diLepton.isMuMu;
        m_diLepton.lepIDs[LepLepID::LL] = l1.lepID[LepID::L] && l2.lepID[LepID::L];
        m_diLepton.lepIDs[LepLepID::LM] = l1.lepID[LepID::L] && l2.lepID[LepID::M];
        m_diLepton.lepIDs[LepLepID::ML] = l1.lepID[LepID::M] && l2.lepID[LepID::L];
        m_diLepton.lepIDs[LepLepID::LT] = l1.lepID[LepID::L] && l2.lepID[LepID::T];
        m_diLepton.lepIDs[LepLepID::TL] = l1.lepID[LepID::T] && l2.lepID[LepID::L];
        m_diLepton.lepIDs[LepLepID::MM] = l1.lepID[LepID::M] && l2.lepID[LepID::M];
        m_diLepton.lepIDs[LepLepID::MT] = l1.lepID[LepID::M] && l2.lepID[LepID::T];
        m_diLepton.lepIDs[LepLepID::TM] = l1.lepID[LepID::T] && l2.lepID[LepID::M];
        m_diLepton.lepIDs[LepLepID::TT] = l1.lepID[LepID::T] && l2.lepID[LepID::T];
        m_diLepton.DR = VectorUtil::DeltaR(l1.p4; l2.p4);
        m_diLepton.DEta = TTAnalysis::DeltaEta(l1.p4; l2.p4);
        m_diLepton.DPhi = VectorUtil::DeltaPhi(l1.p4; l2.p4);

        diLeptons.push_back(m_diLepton);
      }
    }

    for(uint8_t i = 0; i < diLeptons.size(); i++){
      const diLepton& m_diLepton = diLeptons[i];
      
      for(const LepLepID::LepLepID& id: LepLepID::it){
        if(m_diLepton.lepIDs[id])
          diLeptons_LepIDs[id].push_back(i);
      }
      
    }

    ///////////////////////////
    //       JETS            //
    ///////////////////////////

    const JetsProducer& jets = producers.get<JetsProducer>("jets");

    // First find the jets passing kinematic cuts and tight jetID

    for(uint8_t ijet = 0; ijet < jets.p4.size(); ijet++){
        
      // Save the jets that pass the kinematic cuts
      if( abs(jets.p4[ijet].Eta()) < m_jetEtaCut && jets.p4[ijet].Pt() > m_jetPtCut){
        selectedJets.push_back(ijet);
          
        // Save the jets that pass the kinematic cuts and tight jetID
        if( abs(jets.p4[ijet].Eta()) < m_jetEtaCut && jets.p4[ijet].Pt() > m_jetPtCut && jets.passTightID[ijet]){
          selectedJets_tightID.push_back(ijet);
          
          // Save the jets that pass the kinematic cuts and tight jetID and DR(l,j)>cut using selected leptons, for each DiLepton ID pair
          for(const LepLepID::LepLepID& id: LepLepID::it){
            bool passDRcut(true);
            for(const DiLepton& m_diLepton: diLeptons){
              if( m_diLepton.lepIDs[id] && min(VectorUtil::DeltaR(jets.p4[ijet], leptons[m_diLepton.lidxs.first].p4), VectorUtil::DeltaR(jets.p4[ijet], leptons[m_diLepton.lidxs.second].p4)) < m_jetDRleptonCut ){
                passDRcut[id] = false;
                break;
              }
            }
            if(passDRcut)
              selectedJets_tightID_DRcut[id].push_back(ijet);
          }
        } // end selected + tightID
      } // end selected
    } // end for

    // Next, construct DiJets out of these (not accounting for minDRjl here)

    for(uint8_t j1 = 0; j1 < selectedJets_tightID.size(); j1++){
      for(uint8_t j2 = 0; j2 < selectedJets_tightID.size(); j2++){
        const uint8_t jet1 = selectedJets_tightID[j1];
        const uint8_t jet2 = selectedJets_tightID[j2];

        m_diJet = DiJet(); 
        m_diJet.p4 = jets.p4[jet1] + jets.p4[jet2];
        m_diJet.idxs = std::make_pair(jet1, jet2);
        
        m_diJet.DR = VectorUtil::DeltaR(jets.p4[jet1], jets.p4[jet2]);
        m_diJet.DPhi = VectorUtil::DeltaPhi(jets.p4[jet1], jets.p4[jet2]);
        m_diJet.DEta = DeltaEta(jets.p4[jet1], jets.p4[jet2]);
        
        m_diJet.CSVv2_WPs[BBWP::LL] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2L && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2L; 
        m_diJet.CSVv2_WPs[BBWP::LM] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2L && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2M; 
        m_diJet.CSVv2_WPs[BBWP::ML] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2M && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2L; 
        m_diJet.CSVv2_WPs[BBWP::LT] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2L && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2T; 
        m_diJet.CSVv2_WPs[BBWP::TL] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2T && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2L; 
        m_diJet.CSVv2_WPs[BBWP::MM] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2M && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2M; 
        m_diJet.CSVv2_WPs[BBWP::MT] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2M && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2T; 
        m_diJet.CSVv2_WPs[BBWP::TM] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2T && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2M; 
        m_diJet.CSVv2_WPs[BBWP::TT] = jets.getBTagDiscriminant(jet1, m_jetCSVv2Name) > m_jetCSVv2T && jets.getBTagDiscriminant(jet2, m_jetCSVv2Name) > m_jetCSVv2T;
        
        for(const auto& m_diLep: diLeptons){
          const float minDR = std::min( { VectorUtil::DeltaR(jets.p4[ijet1], leptons[m_diLep.lidxs.first].p4), VectorUtil::DeltaR(jets.p4[ijet1], leptons[m_diLep.lidxs.second].p4), VectorUtil::DeltaR(jets.p4[ijet2], leptons[m_diLep.lidxs.first].p4), VectorUtil::DeltaR(jets.p4[ijet2], leptons[m_diLep.lidxs.second].p4) } );
          for(const LepLepID::LepLepID& id: LepLepID::it){
            if( minDR < m_diJet.minDRjl_lepIDs[id] && m_diLep.lepIDs[id] )
              m_diJet.minDRjl_lepIDs[id] = minDR;
          }

        }

        diJets.push_back(m_diJet); 
      }
    }

    // Save the DiJets which have minDRjl>cut, for each leptonID pair
    for(uint8_t i = 0; i < diJets.size(); i++){
      for(const LepLepID::LepLepID& id: LepLepID::it){
        if(diJets.minDRjl_lepIDs[id] > m_jetDRleptonCut)
          diJets_DRCut[id].push_back(i);
      }
    }

    ///////////////////////////
    //       B-JETS          //
    ///////////////////////////
  
    // Save b-jets one by one, this time taking into account minDRjl for different DiLepton ID pairs
    // Do it Pt- or CSVv2-ordered
    for(const LepLepID::LepLepID& id: LepLepID::it){
      for(uint8_t i = 0; i < selectedJets_tightID_DRcut[id].size(); i++){
        uint8_t ijet = selectedJets_tightID_DRcut[id][i];
        if(jets.getBTagDiscriminant(i, m_jetCSVv2Name) > m_jetCSVv2L)
          selectedBJets_DRCut_BWPs_PtOrdered[id][BWP::L].push_back(ijet);
        if(jets.getBTagDiscriminant(i, m_jetCSVv2Name) > m_jetCSVv2M)
          selectedBJets_DRCut_BWPs_PtOrdered[id][BWP::M].push_back(ijet);
        if(jets.getBTagDiscriminant(i, m_jetCSVv2Name) > m_jetCSVv2T)
          selectedBJets_DRCut_BWPs_PtOrdered[id][BWP::T].push_back(ijet);
      }
    }
    selectedBJets_DRCut_BWPs_CSVv2Ordered = selectedBJets_DRCut_PtOrdered;
    for(const LepLepID::LepLepID& id: LepLepID::it){
      for(const BWP::BWP& wp: BWP::it){
        std::sort(selectedBJets_DRCut_BWPs_CSVv2Ordered[id][wp].begin(), selectedBJets_DRCut_BWPs_CSVv2Ordered[id][wp].end(), jetBTagDiscriminantSorter(jets, m_jetCSVv2Name));
      }
    }

    // Save di-b-jets, for each CSVv2 working point pair, taking into account minDRjl for each different DiLepton ID pair
    // Do it Pt- or CSVv2-ordered
    for(uint8_t i = 0; i < diJets.size(); i++){
      const Dijet& m_diJet = diJets[i];
      for(const LepLepID::LepLepID& id: LepLepID::it){
        for(const BBWP::BBWP& wp: BBWP::it){
          if(m_diJet.CSVv2_WPs[wp] && m_diJet.minDRjl_lepIDs[id] > m_jetDRleptonCut)
            diBJets_DRCut_BBWPs_PtOrdered[id][wp].push_back(i);
        }
      }
    }
    diBJets_DRCut_BBWPs_CSVv2Ordered = diBJets_DRCut_BBWPs_PtOrdered;
    for(const LepLepID::LepLepID& id: LepLepID::it){
      for(const BBWP::BBWP& wp: BBWP::it)
        std::sort(diBJets_DRCut_BBWPs_CSVv2Ordered[id][wp].begin(), diBJets_DRCut_BBWPs_CSVv2Ordered[id][wp].end(), diJetBTagDiscriminantSorter(jets, m_jetCSVv2Name, diJets));
    }


    /*// Save the selected jets that pass (or not) the different CSVv2 working points
    for(uint8_t i = 0; i < selectedJets_tightID_DRcut.size(); i++){
      uint8_t ijet = selectedJets_tightID_DRcut[i];

      float csvv2 = jets.getBTagDiscriminant(ijet, m_jetCSVv2Name);

      if(csvv2 >= m_jetCSVv2L){
        selectedBJets_CSVv2L.push_back(ijet);
        
        if(csvv2 >= m_jetCSVv2M){
          selectedBJets_CSVv2M.push_back(ijet);
          
          if(csvv2 >= m_jetCSVv2T){
            selectedBJets_CSVv2T.push_back(ijet);
          }else{
            selectedNonBJets_CSVv2T.push_back(ijet);
          } // end Tight
        
        }else{
          selectedNonBJets_CSVv2M.push_back(ijet);
        } // end Medium
      
      }else{
        selectedNonBJets_CSVv2L.push_back(ijet);
      } // end Loose
    
    } // end for

    // Choose pairs of L b-jets, from Pt- or CSVv2-ordered selected jets

    std::vector<uint8_t> selectedBJets_CSVv2L_CSVv2Ordered = selectedBJets_CSVv2L;
    std::sort(selectedBJets_CSVv2L_CSVv2Ordered.begin(), selectedBJets_CSVv2L_CSVv2Ordered.end(), jetBTagDiscriminantSorter(jets, m_jetCSVv2Name));
    if(selectedBJets_CSVv2L.size() >= 2){
      uint8_t jet1 = selectedBJets_CSVv2L[0];
      uint8_t jet2 = selectedBJets_CSVv2L[1];
      m_diBJet_PtChosen["LL"] = { jets.p4[jet1] + jets.p4[jet2], std::make_pair(jet1, jet2) };
      
      jet1 = selectedBJets_CSVv2L_CSVv2Ordered[0];
      jet2 = selectedBJets_CSVv2L_CSVv2Ordered[1];
      m_diBJet_CSVv2Chosen["LL"] = { jets.p4[jet1] + jets.p4[jet2], std::make_pair(jet1, jet2) };
    }else{
      m_diBJet_PtChosen["LL"] = { myLorentzVector(), std::make_pair(-1, -1) };
      m_diBJet_CSVv2Chosen["LL"] = { myLorentzVector(), std::make_pair(-1, -1) };
    }

    // Choose pairs of M b-jets, from Pt- or CSVv2-ordered selected jets

    std::vector<uint8_t> selectedBJets_CSVv2M_CSVv2Ordered = selectedBJets_CSVv2M;
    std::sort(selectedBJets_CSVv2M_CSVv2Ordered.begin(), selectedBJets_CSVv2M_CSVv2Ordered.end(), jetBTagDiscriminantSorter(jets, m_jetCSVv2Name));
    if(selectedBJets_CSVv2M.size() >= 2){
      uint8_t jet1 = selectedBJets_CSVv2M[0];
      uint8_t jet2 = selectedBJets_CSVv2M[1];
      m_diBJet_PtChosen["MM"] = { jets.p4[jet1] + jets.p4[jet2], std::make_pair(jet1, jet2) };
      
      jet1 = selectedBJets_CSVv2M_CSVv2Ordered[0];
      jet2 = selectedBJets_CSVv2M_CSVv2Ordered[1];
      m_diBJet_CSVv2Chosen["MM"] = { jets.p4[jet1] + jets.p4[jet2], std::make_pair(jet1, jet2) };
    }else{
      m_diBJet_PtChosen["MM"] = { myLorentzVector(), std::make_pair(-1, -1) };
      m_diBJet_CSVv2Chosen["MM"] = { myLorentzVector(), std::make_pair(-1, -1) };
    }

    // Choose pairs of T b-jets, from Pt- or CSVv2-ordered selected jets

    std::vector<uint8_t> selectedBJets_CSVv2T_CSVv2Ordered = selectedBJets_CSVv2T;
    std::sort(selectedBJets_CSVv2T_CSVv2Ordered.begin(), selectedBJets_CSVv2T_CSVv2Ordered.end(), jetBTagDiscriminantSorter(jets, m_jetCSVv2Name));
    if(selectedBJets_CSVv2T.size() >= 2){
      uint8_t jet1 = selectedBJets_CSVv2T[0];
      uint8_t jet2 = selectedBJets_CSVv2T[1];
      m_diBJet_PtChosen["TT"] = { jets.p4[jet1] + jets.p4[jet2], std::make_pair(jet1, jet2) };
      
      jet1 = selectedBJets_CSVv2T_CSVv2Ordered[0];
      jet2 = selectedBJets_CSVv2T_CSVv2Ordered[1];
      m_diBJet_CSVv2Chosen["TT"] = { jets.p4[jet1] + jets.p4[jet2], std::make_pair(jet1, jet2) };
    }else{
      m_diBJet_PtChosen["TT"] = { myLorentzVector(), std::make_pair(-1, -1) };
      m_diBJet_CSVv2Chosen["TT"] = { myLorentzVector(), std::make_pair(-1, -1) };
    }

    // Save the interesting quantitites to disk, for all diBJet working points and both ordering schemes.
    // Always take care that undefined quantities can still be accessed => it will be the categories' role
    // to determine whether a quantity is well-defined or not.

    for(const std::string &wp: {"LL", "MM", "TT"}){
      // Pt-ordered
      diBJet_p4[wp + "-Pt"] = m_diBJet_PtChosen[wp].p4;
      diBJet_idxs[wp + "-Pt"] = std::make_pair(m_diBJet_PtChosen[wp].idxs.first, m_diBJet_PtChosen[wp].idxs.second);
      if(m_diBJet_PtChosen[wp].idxs.first >=0 ){
        diBJet_DR[wp + "-Pt"] = VectorUtil::DeltaR(jets.p4[m_diBJet_PtChosen[wp].idxs.first], jets.p4[m_diBJet_PtChosen[wp].idxs.second]);
        diBJet_DPhi[wp + "-Pt"] = VectorUtil::DeltaPhi(jets.p4[m_diBJet_PtChosen[wp].idxs.first], jets.p4[m_diBJet_PtChosen[wp].idxs.second]);
        diBJet_DEta[wp + "-Pt"] = DeltaEta(jets.p4[m_diBJet_PtChosen[wp].idxs.first], jets.p4[m_diBJet_PtChosen[wp].idxs.second]);
      }else{
        diBJet_DR[wp + "-Pt"] = 0;
        diBJet_DPhi[wp + "-Pt"] = 0;
        diBJet_DEta[wp + "-Pt"] = 0;
      }
      selectedJets_diBJetExcluded[wp + "-Pt"] = std::vector<uint8_t>();
      for(uint8_t ijet = 0; ijet < selectedJets_tightID_DRcut.size(); ijet++){
        if(m_diBJet_PtChosen[wp].idxs.first != ijet && m_diBJet_PtChosen[wp].idxs.second != ijet)
          selectedJets_diBJetExcluded[wp + "-Pt"].push_back(ijet);
      }

      // CSVv2-ordered
      diBJet_p4[wp + "-CSVv2"] = m_diBJet_CSVv2Chosen[wp].p4;
      diBJet_idxs[wp + "-CSVv2"] = std::make_pair(m_diBJet_CSVv2Chosen[wp].idxs.first, m_diBJet_CSVv2Chosen[wp].idxs.second);
      if(m_diBJet_CSVv2Chosen[wp].idxs.first >=0 ){
        diBJet_DR[wp + "-CSVv2"] = VectorUtil::DeltaR(jets.p4[m_diBJet_CSVv2Chosen[wp].idxs.first], jets.p4[m_diBJet_CSVv2Chosen[wp].idxs.second]);
        diBJet_DPhi[wp + "-CSVv2"] = VectorUtil::DeltaPhi(jets.p4[m_diBJet_CSVv2Chosen[wp].idxs.first], jets.p4[m_diBJet_CSVv2Chosen[wp].idxs.second]);
        diBJet_DEta[wp + "-CSVv2"] = DeltaEta(jets.p4[m_diBJet_CSVv2Chosen[wp].idxs.first], jets.p4[m_diBJet_CSVv2Chosen[wp].idxs.second]);
      }else{
        diBJet_DR[wp + "-CSVv2"] = 0;
        diBJet_DPhi[wp + "-CSVv2"] = 0;
        diBJet_DEta[wp + "-CSVv2"] = 0;
      }
      selectedJets_diBJetExcluded[wp + "-CSVv2"] = std::vector<uint8_t>();
      for(uint8_t ijet = 0; ijet < selectedJets_tightID_DRcut.size(); ijet++){
        if(m_diBJet_CSVv2Chosen[wp].idxs.first != ijet && m_diBJet_CSVv2Chosen[wp].idxs.second != ijet)
          selectedJets_diBJetExcluded[wp + "-CSVv2"].push_back(ijet);
      }

    }*/

    ///////////////////////////
    //    EVENT VARIABLES    //
    ///////////////////////////
    
    // leptons-(b-)jets

    int diLepDiJetCounter(0);

    for(uint8_t dilep = 0; dilep < diLeptons.size(); dilep++){
      const DiLepton& m_diLepton = diLeptons[dilep];
      
      for(uint8_t dijet = 0; dijet < diJets.size(); dijet++){
        const DiJet& m_diJet =  diJets[dijet];
        
        DiLepDiJet m_diLepDiJet(m_diLepton, dilep, m_diJet, dijet);

        m_diLepDiJet.minDRjl = min( {
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.second].p4),
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.second].p4)
            } );
        m_diLepDiJet.maxDRjl = max( {
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.second].p4),
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.second].p4)
            } );
        m_diLepDiJet.minDEtajl = min( {
            DeltaEta(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.first].p4),
            DeltaEta(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.second].p4),
            DeltaEta(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.first].p4),
            DeltaEta(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.second].p4)
            } );
        m_diLepDiJet.maxDEtajl = max( {
            DeltaEta(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.first].p4),
            DeltaEta(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.second].p4),
            DeltaEta(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.first].p4),
            DeltaEta(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.second].p4)
            } );
        m_diLepDiJet.minDPhijl = min( {
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.second].p4),
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.second].p4)
            } );
        m_diLepDiJet.maxDPhijl = max( {
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, jets[m_diJet.idxs.second].p4),
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.first].p4),
            VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, jets[m_diJet.idxs.second].p4)
            } );

        diLepDiJets.push_back(m_diLepDiJet);

        for(const LepLepID::LepLepID& id: LepLepID::it){
          if(m_diLepton.lepIDs[id] && m_diLepDiJet.minDRjl > m_jetDRleptonCut){
            diLepDiJets_DRCut[id].push_back(diLepDiJetCounter);
            for(const BBWP::BBWP& wp: BBWP::it){
              if(m_diBJet.CSVv2_WPs[wp])
                diLepDiBJets_DRCut_BBWPs_PtOrdered[id][wp].push_back(diLepDiJetCounter);
            }
          }
        }

        diLepDiJetCounter++;
      }
    }
        
   
    /*ll_jj_p4 = diLepton_p4 + diJet_p4;
    ll_jj_DR = VectorUtil::DeltaR(diLepton_p4, diJet_p4);
    ll_jj_DPhi = VectorUtil::DeltaPhi(diLepton_p4, diJet_p4);
    ll_jj_DEta = DeltaEta(diLepton_p4, diJet_p4);

    lj_minDR = 10;
    lj_noDRcut_minDR = 10;
    if(m_diLepton.idxs.first >= 0){
      for(uint8_t ijet = 0; ijet < selectedJets_tightID_DRcut.size(); ijet++){
        float dr = std::min(
            VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.first].p4, jets.p4[ijet]),
            VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.second].p4, jets.p4[ijet])
            );
        if(dr < lj_minDR)
          lj_minDR = dr;
      }
      for(uint8_t ijet = 0; ijet < selectedJets_tightID.size(); ijet++){
        float dr = std::min(
            VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.first].p4, jets.p4[ijet]),
            VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.second].p4, jets.p4[ijet])
            );
        if(dr < lj_noDRcut_minDR)
          lj_noDRcut_minDR = dr;
      }
    }
    
    // leptons-bjets

    for(const std::string &wp: {"LL", "MM", "TT"}){
      for(const std::string &order: {"-Pt", "-CSVv2"}){
        ll_bb_p4[wp+order] = diLepton_p4 + diBJet_p4[wp+order];
        ll_bb_DR[wp+order] = VectorUtil::DeltaR(diLepton_p4, diBJet_p4[wp+order]);
        ll_bb_DPhi[wp+order] = VectorUtil::DeltaPhi(diLepton_p4, diBJet_p4[wp+order]);
        ll_bb_DEta[wp+order] = DeltaEta(diLepton_p4, diBJet_p4[wp+order]);

        lb_minDR[wp+order] = 10;
        if(m_diLepton.idxs.first >= 0){
          float dr = std::min( {
              VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.first].p4, jets.p4[diBJet_idxs[wp+order].first]),
              VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.second].p4, jets.p4[diBJet_idxs[wp+order].first]),
              VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.first].p4, jets.p4[diBJet_idxs[wp+order].second]),
              VectorUtil::DeltaR(m_leptons[m_diLepton.lidxs.second].p4, jets.p4[diBJet_idxs[wp+order].second])
              } );
          if(dr < lb_minDR[wp+order])
            lb_minDR[wp+order] = dr;
        }
      }
    }*/

    // leptons-(b-)jets-MET

    const METProducer &met = producers.get<METProducer>("met");
    const METProducer &noHFmet = producers.get<METProducer>("nohf_met");

    ll_jj_pfMET_p4 = diLepton_p4 + diJet_p4 + met.p4;
    ll_pfMET_DPhi = VectorUtil::DeltaPhi(diLepton_p4, met.p4);
    jj_pfMET_DPhi = VectorUtil::DeltaPhi(diJet_p4, met.p4);
    ll_jj_pfMET_DPhi = VectorUtil::DeltaPhi(diLepton_p4 + diJet_p4, met.p4);
    
    ll_jj_noHFMET_p4 = diLepton_p4 + diJet_p4 + noHFmet.p4;
    ll_noHFMET_DPhi = VectorUtil::DeltaPhi(diLepton_p4, noHFmet.p4);
    jj_noHFMET_DPhi = VectorUtil::DeltaPhi(diJet_p4, noHFmet.p4);
    ll_jj_noHFMET_DPhi = VectorUtil::DeltaPhi(diLepton_p4 + diJet_p4, noHFmet.p4);

    for(const std::string &wp: {"LL", "MM", "TT"}){
      for(const std::string &order: {"-Pt", "-CSVv2"}){
        ll_bb_p4[wp+order] = diLepton_p4 + diBJet_p4[wp+order];
    
        ll_bb_pfMET_p4[wp+order] = diLepton_p4 + diBJet_p4[wp+order] + met.p4;
        bb_pfMET_DPhi[wp+order] = VectorUtil::DeltaPhi(diBJet_p4[wp+order], met.p4);
        ll_bb_pfMET_DPhi[wp+order] = VectorUtil::DeltaPhi(diLepton_p4 + diBJet_p4[wp+order], met.p4);
        
        ll_bb_noHFMET_p4[wp+order] = diLepton_p4 + diBJet_p4[wp+order] + noHFmet.p4;
        bb_noHFMET_DPhi[wp+order] = VectorUtil::DeltaPhi(diBJet_p4[wp+order], noHFmet.p4);
        ll_bb_noHFMET_DPhi[wp+order] = VectorUtil::DeltaPhi(diLepton_p4 + diBJet_p4[wp+order], noHFmet.p4);
      }
    }

    ///////////////////////////
    //       TRIGGER         //
    ///////////////////////////

    ///////////////////////////
    //       GEN INFO        //
    ///////////////////////////

}

void TTAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
  manager.new_category<TTMuMuCategory>("mumu", "Category with leading leptons as two muons", config);
  manager.new_category<TTElElCategory>("elel", "Category with leading leptons as two electrons", config);
  manager.new_category<TTMuElCategory>("muel", "Category with leading leptons as muon, electron", config);
  manager.new_category<TTElMuCategory>("elmu", "Category with leading leptons as electron, muon", config);
  
  manager.new_category<NoZTTMuMuCategory>("noZmumu", "Category with leading leptons as two muons, excluding the Z peak", config);
  manager.new_category<NoZTTElElCategory>("noZelel", "Category with leading leptons as two electrons, excluding the Z peak", config);
}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, TTAnalyzer, "tt_analyzer");
