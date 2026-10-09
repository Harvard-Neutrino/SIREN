#pragma once
#ifndef SIREN_PhaseSpaceChannel_H
#define SIREN_PhaseSpaceChannel_H

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include "SIREN/dataclasses/PhaseSpaceConvention.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace injection {

// The convention types live in dataclasses.
using PhaseSpaceTopology = siren::dataclasses::PhaseSpaceTopology;
using PhaseSpaceMeasure = siren::dataclasses::PhaseSpaceMeasure;
using PhaseSpaceConvention = siren::dataclasses::PhaseSpaceConvention;
using siren::dataclasses::PhaseSpaceTopologyName;
using siren::dataclasses::PhaseSpaceMeasureName;
using siren::dataclasses::MeasureConvertibilityGroup;
using siren::dataclasses::PhaseSpaceDensityConvertible;

// A proposal must evaluate its density at points drawn by every peer channel.
class PhaseSpaceChannel {
public:
    virtual ~PhaseSpaceChannel() = default;

    template<class Archive>
    void save(Archive &, std::uint32_t const version) const {
        if (version != 0) {
            throw std::runtime_error(
                "PhaseSpaceChannel only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive &, std::uint32_t const version) {
        if (version != 0) {
            throw std::runtime_error(
                "PhaseSpaceChannel only supports version <= 0!");
        }
    }

    virtual void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const = 0;

    virtual double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const = 0;

    virtual std::string Name() const = 0;

    // Structural shape of the final state.
    virtual PhaseSpaceTopology Topology() const = 0;

    // What variables the density is differential in.
    virtual PhaseSpaceMeasure Measure() const = 0;

};

double ConvertDensity(
    double density,
    PhaseSpaceMeasure const & from,
    PhaseSpaceMeasure const & to,
    PhaseSpaceTopology topology,
    siren::dataclasses::InteractionRecord const & record);

// Mixture density is sum_i weights[i] * channels[i]->Density(...).
struct MultiChannelPhaseSpace {
    std::vector<std::shared_ptr<PhaseSpaceChannel>> channels;
    std::vector<double> weights;  // alpha_i, must sum to 1

    // Skips the fatal channel-compatibility checks. Stored in archives.
    bool allow_incompatible_ = false;

    // Leaves channels and weights empty; assign them, then call Normalize().
    MultiChannelPhaseSpace() = default;

    explicit MultiChannelPhaseSpace(
        std::vector<std::shared_ptr<PhaseSpaceChannel>> channels,
        std::vector<double> weights = {},
        bool allow_incompatible = false);

    // Normalize weights in place to sum 1.  Empty weights -> uniform 1/N.
    // Throws ConfigurationError on channels/weights length mismatch or a
    // non-positive weight sum.
    void Normalize();

    // Severity-tagged compatibility diagnostic.  Fatal entries block a mixture
    // (measures not convertible); Info entries report a supported auto-conversion.
    struct ChannelDiagnostic {
        enum class Severity { Info, Fatal };
        Severity severity;
        std::string message;
    };

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if (version == 0) {
            archive(::cereal::make_nvp("Channels", channels));
            archive(::cereal::make_nvp("Weights", weights));
            archive(::cereal::make_nvp(
                "AllowIncompatible", allow_incompatible_));
        } else {
            throw std::runtime_error(
                "MultiChannelPhaseSpace only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if (version != 0) {
            throw std::runtime_error(
                "MultiChannelPhaseSpace only supports version <= 0!");
        }

        std::vector<std::shared_ptr<PhaseSpaceChannel>> loaded_channels;
        std::vector<double> loaded_weights;
        bool loaded_allow_incompatible = false;
        archive(::cereal::make_nvp("Channels", loaded_channels));
        archive(::cereal::make_nvp("Weights", loaded_weights));
        archive(::cereal::make_nvp(
            "AllowIncompatible", loaded_allow_incompatible));

        MultiChannelPhaseSpace loaded;
        loaded.channels = std::move(loaded_channels);
        loaded.weights = std::move(loaded_weights);
        loaded.allow_incompatible_ = loaded_allow_incompatible;
        loaded.RequireNormalizedWeights("load");
        loaded.ThrowOnIncompatibility();
        *this = std::move(loaded);
    }

    // Sample from the multi-channel mixture.
    // Returns the index of the channel that was used.
    int Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const;

    double DensityIn(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record,
        PhaseSpaceMeasure const & measure
    ) const;

    // DensityIn with the full convention: additionally checks that the
    // mixture's topology matches, throwing MeasureCompatibilityError when it
    // does not.
    double DensityIn(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record,
        PhaseSpaceConvention const & convention
    ) const;

    std::vector<double> DensityBreakdown(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const;

    // Return the common topology. Throws if channels disagree.
    PhaseSpaceTopology CommonTopology() const;

    PhaseSpaceMeasure CommonMeasure() const;

    // CommonTopology() and CommonMeasure() as one convention.
    PhaseSpaceConvention CommonConvention() const;

    // Validate topology and measure compatibility, returning every diagnostic
    // (Fatal and Info) with a severity tag.  Empty if all checks pass.
    std::vector<ChannelDiagnostic> ValidateChannelsDetailed() const;

    std::vector<std::string> ValidateChannelDensities(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord template_record,
        int samples_per_channel = 100
    ) const;

private:
    mutable bool convention_cache_valid_ = false;
    mutable std::size_t convention_fingerprint_ = 0;
    mutable PhaseSpaceTopology cached_common_topology_ =
        PhaseSpaceTopology::Unspecified;
    mutable PhaseSpaceMeasure cached_common_measure_ =
        PhaseSpaceMeasure::Unspecified();
    mutable std::vector<PhaseSpaceMeasure> cached_channel_measures_;
    mutable std::string cached_topology_error_;
    mutable std::vector<ChannelDiagnostic> cached_compatibility_diagnostics_;

    std::size_t ConventionFingerprint() const;
    void EnsureConventionCache() const;

    void ThrowOnIncompatibility() const;
    // Throws ConfigurationError unless the weights are valid and sum to one.
    void RequireNormalizedWeights(char const * where) const;

    double ComputeContributions(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record,
        std::vector<double> * weighted,
        std::vector<double> * bare
    ) const;

};

class NestedMixtureChannel : public PhaseSpaceChannel {
public:
    std::shared_ptr<MultiChannelPhaseSpace> mixture;
    std::string label = "NestedMixture";

    NestedMixtureChannel() = default;
    explicit NestedMixtureChannel(
        std::shared_ptr<MultiChannelPhaseSpace> mixture_)
        : mixture(mixture_) {}

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if (version == 0) {
            archive(::cereal::make_nvp(
                "PhaseSpaceChannel",
                ::cereal::virtual_base_class<PhaseSpaceChannel>(this)));
            archive(::cereal::make_nvp("Mixture", mixture));
            archive(::cereal::make_nvp("Label", label));
        } else {
            throw std::runtime_error(
                "NestedMixtureChannel only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if (version == 0) {
            archive(::cereal::make_nvp(
                "PhaseSpaceChannel",
                ::cereal::virtual_base_class<PhaseSpaceChannel>(this)));
            archive(::cereal::make_nvp("Mixture", mixture));
            archive(::cereal::make_nvp("Label", label));
        } else {
            throw std::runtime_error(
                "NestedMixtureChannel only supports version <= 0!");
        }
    }

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const override;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const override;

    std::string Name() const override;
    PhaseSpaceTopology Topology() const override;
    PhaseSpaceMeasure Measure() const override;

};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(siren::injection::PhaseSpaceChannel, 0);
CEREAL_CLASS_VERSION(siren::injection::MultiChannelPhaseSpace, 0);
CEREAL_CLASS_VERSION(siren::injection::NestedMixtureChannel, 0);
CEREAL_REGISTER_TYPE(siren::injection::NestedMixtureChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::NestedMixtureChannel);

#endif // SIREN_PhaseSpaceChannel_H
