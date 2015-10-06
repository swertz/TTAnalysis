#include <cp3_llbb/TTAnalysis/interface/TTAnalyzer.h>
#include <cp3_llbb/TTAnalysis/interface/NoZTTDileptonCategories.h>
#include <cp3_llbb/TTAnalysis/interface/TTDileptonCategories.h>

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

void TTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const CategoryManager& categories) {

    ///////////////////////////
    //       ELECTRONS       //
    ///////////////////////////

    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

    for(unsigned int ielectron = 0; ielectron < electrons.p4.size(); ielectron++){
        if( electrons.ids[ielectron][m_electronVetoIDName] && electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut )
            vetoElectrons.push_back(ielectron);
        if( electrons.ids[ielectron][m_electronLooseIDName] && electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut )
            looseElectrons.push_back(ielectron);
        if( electrons.ids[ielectron][m_electronMediumIDName] && electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut )
            mediumElectrons.push_back(ielectron);
        if( electrons.ids[ielectron][m_electronTightIDName] && electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut )
            tightElectrons.push_back(ielectron);

        if( electrons.ids[ielectron][m_electronSelectedIDName] && electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut ){
            selectedElectrons.push_back(ielectron);
            m_leptons.push_back( { electrons.p4[ielectron], ielectron, true, false } );
        }
    }

    ///////////////////////////
    //       MUONS           //
    ///////////////////////////
    
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");

    for(unsigned int imuon = 0; imuon < muons.p4.size(); imuon++){
        if(muons.relativeIsoR04_withEA[imuon] < m_muonBaseIsoCut &&  muons.isLoose[imuon] && muons.p4[imuon].Pt() > m_muonPtCut && abs(muons.p4[imuon].Eta()) < m_muonEtaCut )
            looseMuons.push_back(imuon);
        if(muons.relativeIsoR04_withEA[imuon] < m_muonBaseIsoCut &&  muons.isMedium[imuon] && muons.p4[imuon].Pt() > m_muonPtCut && abs(muons.p4[imuon].Eta()) < m_muonEtaCut )
            mediumMuons.push_back(imuon);
        if(muons.relativeIsoR04_withEA[imuon] < m_muonBaseIsoCut &&  muons.isTight[imuon] && muons.p4[imuon].Pt() > m_muonPtCut && abs(muons.p4[imuon].Eta()) < m_muonEtaCut )
            tightMuons.push_back(imuon);

        if(muonIDAccessor(muons, imuon, m_muon_selectedID) && muons.relativeIsoR04_withEA[imuon] < m_muonSelectedIsoCut && muons.p4[imuon].Pt() > m_muonPtCut && abs(muons.p4[imuon].Eta()) < m_muonEtaCut ){
            selectedMuons.push_back(imuon);
            m_leptons.push_back( { muons.p4[imuon], imuon, false, true } );
        }
    }

    // Sort the m_leptons vector according to Pt and write the content to disk
    std::sort(m_leptons.begin(), m_leptons.end(), [](Lepton a, Lepton b){ return b.p4.Pt() > a.p4.Pt(); });
    for(const Lepton lepton&: m_leptons){
      lepton_p4.push_back(lepton.p4);
      lepton_idx.push_back(lepton.idx);
      lepton_isMu.push_back(lepton.isMu);
      lepton_isEl.push_back(lepton.isEl);
    }

    ///////////////////////////
    //       DILEPTONS       //
    ///////////////////////////

    // First find the highest-Pt opposite-charge pairs for each couple of flavours, on selected objets only
    // Already apply a cut on Mll
    
    leadingSelectedElEl = std::make_pair(-1, -1);
    leadingSelectedElMu = std::make_pair(-1, -1);
    leadingSelectedMuEl = std::make_pair(-1, -1);
    leadingSelectedMuMu = std::make_pair(-1, -1);

    float tempSumPtElEl = 0, tempSumPtElMu = 0, tempSumPtMuEl = 0, tempSumPtMuMu = 0; 

    // El-El
    for(unsigned int iele1 = 0; iele1 < selectedElectrons.size(); iele1++){
        for(unsigned int iele2 = iele1+1; iele2 < selectedElectrons.size(); iele2++){
            if( electrons.charge[ selectedElectrons[iele1] ] * electrons.charge[ selectedElectrons[iele2] ] < 0 && (electrons.p4[ selectedElectrons[iele1] ] + electrons.p4[ selectedElectrons[iele2] ]).M() > m_MllBaseCutSF ){
                leadingSelectedElEl = std::make_pair(iele1, iele2);
                tempSumPtElEl = electrons.p4[ selectedElectrons[iele1] ].Pt() + electrons.p4[ selectedElectrons[iele2] ].Pt(); 
                break;
            }
        }
        if(leadingSelectedElEl.first >= 0)
            break;
    }

    // El-Mu
    for(unsigned int iele = 0; iele < selectedElectrons.size(); iele++){
        for(unsigned int imu = 0; imu < selectedMuons.size(); imu++){
            if( electrons.p4[ selectedElectrons[iele] ].Pt() > muons.p4[ selectedMuons[imu] ].Pt() && electrons.charge[ selectedElectrons[iele] ] * muons.charge[ selectedMuons[imu] ] < 0 && (electrons.p4[ selectedElectrons[iele] ] + muons.p4[ selectedMuons[imu] ]).M() > m_MllBaseCutDF){
                leadingSelectedElMu = std::make_pair(iele, imu);
                tempSumPtElMu = electrons.p4[ selectedElectrons[iele] ].Pt() + muons.p4[ selectedMuons[imu] ].Pt(); 
                break;
            }
        }
        if(leadingSelectedElMu.first >= 0)
            break;
    }
    
    // Mu-El
    for(unsigned int imu = 0; imu < selectedMuons.size(); imu++){
        for(unsigned int iele = 0; iele < selectedElectrons.size(); iele++){
            if( electrons.p4[ selectedElectrons[iele] ].Pt() < muons.p4[ selectedMuons[imu] ].Pt() && electrons.charge[ selectedElectrons[iele] ] * muons.charge[ selectedMuons[imu] ] < 0 && (electrons.p4[ selectedElectrons[iele] ] + muons.p4[ selectedMuons[imu] ]).M() > m_MllBaseCutDF){
                leadingSelectedMuEl = std::make_pair(imu, iele);
                tempSumPtMuEl = electrons.p4[ selectedElectrons[iele] ].Pt() + muons.p4[ selectedMuons[imu] ].Pt(); 
                break;
            }
        }
        if(leadingSelectedMuEl.first >= 0)
            break;
    }

    // Mu-Mu
    for(unsigned int imu1 = 0; imu1 < selectedMuons.size(); imu1++){
        for(unsigned int imu2 = imu1+1; imu2 < selectedMuons.size(); imu2++){
            if( muons.charge[ selectedMuons[imu1] ] * muons.charge[ selectedMuons[imu2] ] < 0 && (muons.p4[ selectedMuons[imu1] ] + muons.p4[ selectedMuons[imu2] ]).M() > m_MllBaseCutSF){
                leadingSelectedMuMu = std::make_pair(imu1, imu2);
                tempSumPtMuMu = muons.p4[ selectedMuons[imu1] ].Pt() + muons.p4[ selectedMuons[imu2] ].Pt(); 
                break;
            }
        }
        if(leadingSelectedMuMu.first >= 0)
            break;
    }

    // Next, we define the selected lepton pair as the highest Sum(Pt) pair

    if( (tempSumPtElEl > 0) && (tempSumPtElEl > tempSumPtElMu) && (tempSumPtElEl > tempSumPtMuEl) && (tempSumPtElEl > tempSumPtMuMu) ){
      m_diLepton = { 
        electrons.p4[leadingSelectedElEl.first] + electrons.p4[leadingSelectedElEl.second], 
        std::make_pair<int, int>(leadingSelectedElEl.first, leadingSelectedElEl.second),
        true, false, false, false
        };
    }else if( (tempSumPtElMu > 0) && (tempSumPtElMu > tempSumPtElEl) && (tempSumPtElMu > tempSumPtMuEl) && (tempSumPtElMu > tempSumPtMuMu) ){
      m_diLepton = { 
        electrons.p4[leadingSelectedElMu.first] + muons.p4[leadingSelectedElMu.second], 
        std::make_pair<int, int>(leadingSelectedElMu.first, leadingSelectedElMu.second),
        false, true, false, false
        };
    }else if( (tempSumPtMuEl > 0) && (tempSumPtMuEl > tempSumPtElEl) && (tempSumPtMuEl > tempSumPtElMu) && (tempSumPtMuEl > tempSumPtMuMu) ){
      m_diLepton = { 
        muons.p4[leadingSelectedMuEl.first] + electrons.p4[leadingSelectedMuEl.second], 
        std::make_pair<int, int>(leadingSelectedMuEl.first, leadingSelectedMuEl.second),
        false, false, true, false
        };
    }else if( (tempSumPtMuMu > 0) && (tempSumPtMuMu > tempSumPtElEl) && (tempSumPtMuMu > tempSumPtElMu) && (tempSumPtMuMu > tempSumPtMuEl) ){
      m_diLepton = { 
        muons.p4[leadingSelectedMuMu.first] + muons.p4[leadingSelectedMuMu.second], 
        std::make_pair<int, int>(leadingSelectedMuMu.first, leadingSelectedMuMu.second),
        false, false, false, true
        };
    }else{
      m_diLepton = { LorentzVector(0,0,0,0), std::make_pair<int, int>(-1, -1), false, false, false, false };
    }

    


    ///////////////////////////
    //       JETS            //
    ///////////////////////////

    // Save the jets that pass the cuts

    const JetsProducer& jets = producers.get<JetsProducer>("jets");

    for(unsigned int ijet = 0; ijet < jets.p4.size(); ijet++){
        if( abs(jets.p4[ijet].Eta()) < m_jetEtaCut && jets.p4[ijet].Pt() > m_jetPtCut)
            selectedJets.push_back(ijet);
    }

    ///////////////////////////
    //       B-JETS          //
    ///////////////////////////

    ///////////////////////////
    //    EVENT VARIABLES    //
    ///////////////////////////

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
