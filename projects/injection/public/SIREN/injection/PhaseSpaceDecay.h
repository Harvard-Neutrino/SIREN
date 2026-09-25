#pragma once
#include "SIREN/interactions/Decay.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

namespace siren { namespace injection {

// A native decay whose PHYSICAL normalized final-state law is explicitly given
// by a channel. An importance proposal belongs on the injection process, not
// here. Widths and final-state masses are supplied by the owning physics model.
class PhaseSpaceDecay : public siren::interactions::Decay {
public:
    PhaseSpaceDecay(siren::dataclasses::InteractionSignature signature,
        std::vector<double> masses, double partial_width, double total_width,
        std::shared_ptr<PhaseSpaceChannel> physical_channel);
    bool equal(siren::interactions::Decay const &) const override;
    double ParentDecayWidth(siren::dataclasses::InteractionRecord const &) const override;
    double TotalDecayWidthAllFinalStates(siren::dataclasses::InteractionRecord const &) const override;
    double TotalDecayWidth(siren::dataclasses::ParticleType) const override;
    double TotalDecayWidth(siren::dataclasses::InteractionRecord const &) const override;
    double DifferentialDecayWidth(siren::dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(siren::dataclasses::CrossSectionDistributionRecord &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override;
    std::vector<double> SecondaryMasses(std::vector<siren::dataclasses::ParticleType> const &) const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(
        siren::dataclasses::ParticleType) const override;
    double FinalStateProbability(siren::dataclasses::InteractionRecord const &) const override;
    std::vector<std::string> DensityVariables() const override;
    PhaseSpaceTopology Topology() const override;
    PhaseSpaceMeasure Measure() const override;
    PhaseSpaceTopology TopologyForSignature(siren::dataclasses::InteractionSignature const & s) const override {
        if (s != signature_) throw std::invalid_argument("Unknown PhaseSpaceDecay signature");
        return Topology();
    }
    PhaseSpaceMeasure MeasureForSignature(siren::dataclasses::InteractionSignature const & s) const override {
        if (s != signature_) throw std::invalid_argument("Unknown PhaseSpaceDecay signature");
        return Measure();
    }
    std::shared_ptr<PhaseSpaceChannel> GetPhysicalChannel() const { return physical_channel_; }
    template<class Archive> void save(Archive & ar, std::uint32_t version) const {
        if (version != 1) throw std::runtime_error("Legacy PhaseSpaceDecay branching contract rejected; regenerate the archive");
        ar(cereal::make_nvp("Signature", signature_), cereal::make_nvp("Masses", masses_),
           cereal::make_nvp("PartialWidth", partial_width_), cereal::make_nvp("TotalWidth", total_width_),
           cereal::make_nvp("PhysicalChannel", physical_channel_),
           cereal::virtual_base_class<siren::interactions::Decay>(this));
    }
    template<class Archive> void load(Archive & ar, std::uint32_t version) {
        if (version != 1) throw std::runtime_error("Legacy PhaseSpaceDecay branching contract rejected; regenerate the archive");
        ar(cereal::make_nvp("Signature", signature_), cereal::make_nvp("Masses", masses_),
           cereal::make_nvp("PartialWidth", partial_width_), cereal::make_nvp("TotalWidth", total_width_),
           cereal::make_nvp("PhysicalChannel", physical_channel_),
           cereal::virtual_base_class<siren::interactions::Decay>(this));
        Validate();
    }
private:
    friend class cereal::access;
    PhaseSpaceDecay() = default;
    void Validate() const;
    siren::dataclasses::InteractionSignature signature_;
    std::vector<double> masses_;
    double partial_width_ = 0, total_width_ = 0;
    std::shared_ptr<PhaseSpaceChannel> physical_channel_;
};
}}
CEREAL_CLASS_VERSION(siren::injection::PhaseSpaceDecay, 1);
CEREAL_REGISTER_TYPE(siren::injection::PhaseSpaceDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::injection::PhaseSpaceDecay);
