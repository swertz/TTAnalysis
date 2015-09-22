#ifndef ZVETOTTDILEPTONCATEGORIES_H
#define ZVETOTTDILEPTONCATEGORIES_H

#include <cp3_llbb/Framework/interface/Category.h>

class ZVetoTTDileptonCategory: public Category {
    public:
        virtual void configure(const edm::ParameterSet& conf) {
            m_mll_cut_low = conf.getUntrackedParameter<double>("mll_cut_low", 86);
            m_mll_cut_high = conf.getUntrackedParameter<double>("mll_cut_high", 116);
        }
        virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
        virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
        //virtual void register_cuts(CutManager& manager) override;
        //virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;

    protected:
        float m_mll_cut_low, m_mll_cut_high;
};

#endif
