#pragma once
#ifndef SIREN_CharmMesonDecay3Body_H
#define SIREN_CharmMesonDecay3Body_H

// CharmMesonDecay3Body -- Pythia-style 3-body phase-space decay for D mesons
//
// Sister class to CharmMesonDecay (the legacy 2-body-cascade implementation).
// Both inherit from Decay and share the same decay-width machinery
// (DifferentialDecayWidth, TotalDecayWidthForFinalState, FinalStateProbability,
// and the signature catalog). The only behavioural difference is in
// SampleFinalState, where this class generates the final-state kinematics by
// sampling 3-body phase space (following Pythia's ParticleDecays::threeBody)
// with V-A matrix element reweighting. It also mixes D -> K l nu with
// D -> K*(892) l nu channels per event, in line with PDG branching ratios.

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <math.h>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>


#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"

#include "SIREN/interactions/Decay.h"
#include "SIREN/utilities/Interpolator.h"


namespace siren {
namespace interactions {

class CharmMesonDecay3Body : public Decay {
friend cereal::access;
private:
    const std::set<siren::dataclasses::Particle::ParticleType> primary_types = {siren::dataclasses::Particle::ParticleType::D0, siren::dataclasses::Particle::ParticleType::DPlus};
    siren::utilities::Interpolator1D<double> inverseCdf; // for dGamma (form-factor model; not used by FinalStateProbability)
    // Shared closure helpers: SampleFinalState's density and FinalStateProbability
    // both build on these, so Sample == Density by construction.
    static double KStarMass();
    double SampledQ2Density(double mD, double mK, double ml, double q2, bool apply_va) const;
    double SampledQ2Normalization(double mD, double mK, double ml, bool apply_va) const;
    // Per-component normalization cache (not serialized; keyed by mass set).
    mutable std::map<long, double> norm_cache;
public:
    CharmMesonDecay3Body();
    CharmMesonDecay3Body(siren::dataclasses::Particle::ParticleType primary);
    virtual bool equal(Decay const & other) const override;
    static double particleMass(siren::dataclasses::ParticleType particle);
    // Analytic angle-average of the accepted V-A weight (q^2 density factor).
    // Public so the closure/regression test can check it against a numeric
    // quadrature oracle; pure function of the decay masses and m23.
    double VAWeightAngleAverage(double mD, double mK, double ml, double m23) const;
    double TotalDecayWidth(dataclasses::InteractionRecord const &) const override;
    double TotalDecayWidth(siren::dataclasses::Particle::ParticleType primary) const override;
    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const override;
    double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const override;
    double DifferentialDecayWidth(std::vector<double> constants, double Q2, double mD, double mK) const;
    void SampleFinalStateHadronic(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::Particle::ParticleType primary) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    std::vector<double> FormFactorFromRecord(dataclasses::CrossSectionDistributionRecord const & record) const;
    void computeDiffGammaCDF(std::vector<double> constants, double mD, double mK);

public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("CharmMesonDecay3Body only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            // Default-constructible with a fixed const primary_types, so cereal
            // default-constructs then calls this load(). Read PrimaryTypes into a
            // temporary to consume the stream symmetrically with save(); a
            // load_and_construct would be bypassed for a default-constructible
            // type, leaving the body unread and corrupting following archive data.
            std::set<siren::dataclasses::Particle::ParticleType> _primary_types;
            archive(::cereal::make_nvp("PrimaryTypes", _primary_types));
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("CharmMesonDecay3Body only supports version <= 0!");
        }
    }

};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::CharmMesonDecay3Body, 0);
CEREAL_REGISTER_TYPE(siren::interactions::CharmMesonDecay3Body);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::interactions::CharmMesonDecay3Body);

#endif // SIREN_CharmMesonDecay3Body_H
