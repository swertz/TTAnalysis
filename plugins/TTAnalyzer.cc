#include <cp3_llbb/TTAnalysis/interface/TTAnalyzer.h>
#include <cp3_llbb/TTAnalysis/interface/NoZTTDileptonCategories.h>
#include <cp3_llbb/TTAnalysis/interface/TTDileptonCategories.h>

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>

void TTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const CategoryManager& categories) {

    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");

    for(unsigned int ielectron = 0; ielectron < electrons.p4.size(); ielectron++){
        if( electrons.ids[ielectron][m_electron_loose_wp_name] )
            looseElectrons.push_back(ielectron);
        if( electrons.ids[ielectron][m_electron_tight_wp_name] )
            tightElectrons.push_back(ielectron);
        if( electrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut )
            isolatedElectrons.push_back(ielectron);

        if( electrons.ids[ielectron][m_electron_tight_wp_name] && electrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut && electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut )
          selectedElectrons.push_back(ielectron);
    }

    for(unsigned int imuon = 0; imuon < muons.p4.size(); imuon++){
        if( muons.isLoose[imuon] )
            looseMuons.push_back(imuon);
        if( muons.isTight[imuon] )
            tightMuons.push_back(imuon);
        if( muons.relativeIsoR04_withEA[imuon] < m_muonIsoCut )
            isolatedMuons.push_back(imuon);

        if( muons.isTight[imuon] && muons.relativeIsoR04_withEA[imuon] < m_muonIsoCut && muons.p4[imuon].Pt() > m_muonPtCut && abs(muons.p4[imuon].Eta()) < m_muonEtaCut )
          selectedMuons.push_back(imuon);
    }

    // Find the highest-Pt opposite-charge pairs for each couple of flavours, on selected objets only
    selectedLeadingElEl = std::make_pair(-1, -1);
    selectedLeadingElMu = std::make_pair(-1, -1);
    selectedLeadingMuEl = std::make_pair(-1, -1);
    selectedLeadingMuMu = std::make_pair(-1, -1);

    // El-El
    for(unsigned int iele1 = 0; iele1 < selectedElectrons.size(); iele1++){
        for(unsigned int iele2 = iele1+1; iele2 < selectedElectrons.size(); iele2++){
            if( electrons.charge[ selectedElectrons[iele1] ] * electrons.charge[ selectedElectrons[iele2] ] < 0 ){
                selectedLeadingElEl = std::make_pair(iele1, iele2);
                break;
            }
        }
        if(selectedLeadingElEl.first >= 0)
            break;
    }

    // El-Mu
    for(unsigned int iele = 0; iele < selectedElectrons.size(); iele++){
        for(unsigned int imu = 0; imu < selectedMuons.size(); imu++){
            if( electrons.p4[ selectedElectrons[iele] ].Pt() > muons.p4[ selectedMuons[imu] ].Pt() && electrons.charge[ selectedElectrons[iele] ] * muons.charge[ selectedMuons[imu] ] < 0 ){
                selectedLeadingElMu = std::make_pair(iele, imu);
                break;
            }
        }
        if(selectedLeadingElMu.first >= 0)
            break;
    }
    
    // Mu-El
    for(unsigned int imu = 0; imu < selectedMuons.size(); imu++){
        for(unsigned int iele = 0; iele < selectedElectrons.size(); iele++){
            if( electrons.p4[ selectedElectrons[iele] ].Pt() < muons.p4[ selectedMuons[imu] ].Pt() && electrons.charge[ selectedElectrons[iele] ] * muons.charge[ selectedMuons[imu] ] < 0 ){
                selectedLeadingMuEl = std::make_pair(imu, iele);
                break;
            }
        }
        if(selectedLeadingMuEl.first >= 0)
            break;
    }

    // Mu-Mu
    for(unsigned int imu1 = 0; imu1 < selectedMuons.size(); imu1++){
        for(unsigned int imu2 = imu1+1; imu2 < selectedMuons.size(); imu2++){
            if( muons.charge[ selectedMuons[imu1] ] * muons.charge[ selectedMuons[imu2] ] < 0 ){
                selectedLeadingMuMu = std::make_pair(imu1, imu2);
                break;
            }
        }
        if(selectedLeadingMuMu.first >= 0)
            break;
    }

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
