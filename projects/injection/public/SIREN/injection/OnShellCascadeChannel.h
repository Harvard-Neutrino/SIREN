#pragma once
#include <array>
#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

namespace siren { namespace injection {

// Prompt P -> spectator + R, R -> d1 + d2, with equal pair masses.
// The normalized internal proposal is (1 + kappa*c*c)/(2*(1+kappa/3)).
// This is a proposal parameter; physical amplitudes remain model-owned.
// Final-state order is spectator, d1, d2. Complete configurations are rotated,
// and both daughter-targeting densities enter the returned mixture density.
class OnShellCascadeChannel : public PhaseSpaceChannel {
public:
    OnShellCascadeChannel(std::shared_ptr<siren::geometry::Geometry const> target,
        double pair_mass, double kappa=0,
        std::array<double,3> orientation_weights={0.1,0.45,0.45},
        double first_daughter_probability=0.5, double volume=-1.0);
    void Sample(std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        siren::dataclasses::InteractionRecord &) const override;
    double Density(std::shared_ptr<siren::detector::DetectorModel const>,
        siren::dataclasses::InteractionRecord const &) const override;
    double InternalDensity(siren::dataclasses::InteractionRecord const &) const;
    std::string Name() const override { return "OnShellCascade"; }
    PhaseSpaceTopology Topology() const override { return PhaseSpaceTopology::Decay3Body; }
    PhaseSpaceMeasure Measure() const override { return PhaseSpaceMeasure::OnShellCascade(pair_mass_); }
    template<class Archive> void save(Archive & ar, std::uint32_t version) const {
        if (version != 1) throw std::runtime_error("Legacy OnShellCascadeChannel density rejected; regenerate the archive");
        ar(cereal::make_nvp("Target",target_), cereal::make_nvp("PairMass",pair_mass_),
           cereal::make_nvp("Kappa",kappa_), cereal::make_nvp("OrientationWeights",weights_),
           cereal::make_nvp("FirstDaughterProbability",rho_), cereal::make_nvp("TargetVolume",target_volume_),
           cereal::virtual_base_class<PhaseSpaceChannel>(this));
    }
    template<class Archive> void load(Archive & ar, std::uint32_t version) {
        if (version != 1) throw std::runtime_error("Legacy OnShellCascadeChannel density rejected; regenerate the archive");
        ar(cereal::make_nvp("Target",target_), cereal::make_nvp("PairMass",pair_mass_),
           cereal::make_nvp("Kappa",kappa_), cereal::make_nvp("OrientationWeights",weights_),
           cereal::make_nvp("FirstDaughterProbability",rho_), cereal::make_nvp("TargetVolume",target_volume_),
           cereal::virtual_base_class<PhaseSpaceChannel>(this));
        Validate(true);
        BuildComponents();
    }
private:
    friend class cereal::access;
    OnShellCascadeChannel()=default;
    void Validate(bool archived) const;
    // Envelope and volume proposals, built once from the checked volume.
    void BuildComponents();
    std::shared_ptr<siren::geometry::Geometry const> target_;
    double pair_mass_=0, kappa_=0, rho_=0.5;
    std::array<double,3> weights_{0.1,0.45,0.45};
    double target_volume_=-1;
    std::shared_ptr<PhaseSpaceChannel const> envelope_, directed_;
};
}}
CEREAL_CLASS_VERSION(siren::injection::OnShellCascadeChannel,1);
CEREAL_REGISTER_TYPE(siren::injection::OnShellCascadeChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::injection::PhaseSpaceChannel,siren::injection::OnShellCascadeChannel);
