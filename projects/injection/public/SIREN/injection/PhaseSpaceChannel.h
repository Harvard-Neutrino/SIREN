#pragma once
#ifndef SIREN_PhaseSpaceChannel_H
#define SIREN_PhaseSpaceChannel_H

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cereal/access.hpp>
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
class MultiChannelPhaseSpace {
public:
    explicit MultiChannelPhaseSpace(
        std::vector<std::shared_ptr<PhaseSpaceChannel>> channels,
        std::vector<double> weights = {});

    std::vector<std::shared_ptr<PhaseSpaceChannel>> const & Channels() const {
        return channels_;
    }
    std::vector<double> const & Weights() const { return weights_; }

    // Normalize the supplied weights, or use uniform weights for an empty
    // vector. A failed update leaves the existing weights intact.
    void SetWeights(std::vector<double> weights);

    // Severity-tagged compatibility diagnostic.  Fatal entries block a mixture
    // (measures not convertible); Info entries report a supported auto-conversion.
    struct ChannelDiagnostic {
        enum class Severity { Info, Fatal };
        Severity severity;
        std::string message;
    };

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

    // Accumulate w^2 g_i/g in the common measure. Nested mixtures share the
    // outer density denominator.
    void Accumulate(std::shared_ptr<siren::detector::DetectorModel const> detector,
                    siren::dataclasses::InteractionRecord const & record,
                    double weight, bool recurse = true);
    void UpdateWeights(std::string const & rule = "alpha_sqrt_W",
                       double damping = 0.5, double min_weight = 0.001,
                       bool recurse = true);
    void ResetAccumulators(bool recurse = true);
    std::uint64_t TuningCount() const { return tuning_count_; }

    // Return the convention shared by, or convertible from, all channels.
    // Throws MeasureCompatibilityError if no common convention exists.
    PhaseSpaceTopology CommonTopology() const;

    PhaseSpaceMeasure CommonMeasure() const;

    // CommonTopology() and CommonMeasure() as one convention.
    PhaseSpaceConvention CommonConvention() const;

    // Validate topology and measure compatibility, returning every diagnostic
    // (Fatal and Info) with a severity tag. Info describes auto-conversions.
    std::vector<ChannelDiagnostic> ValidateChannelsDetailed() const;

private:
    friend class cereal::access;
    // Cereal creates an empty object before loading and validating its fields.
    MultiChannelPhaseSpace() = default;

    std::vector<std::shared_ptr<PhaseSpaceChannel>> channels_;
    std::vector<double> weights_;

    std::vector<double> tuning_log_sums_;
    std::vector<double> tuning_weights_;
    std::vector<std::shared_ptr<PhaseSpaceChannel>> tuning_channels_;
    std::uint64_t tuning_count_ = 0;
    void CheckTuningState() const;

    struct ConventionEvaluation {
        PhaseSpaceConvention common;
        std::vector<PhaseSpaceMeasure> measures;
        std::vector<ChannelDiagnostic> diagnostics;
    };

    ConventionEvaluation EvaluateConventions() const;
    ConventionEvaluation RequireCompatibleConvention() const;
    void ValidateArchivedWeights() const;

    double ComputeDensity(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record,
        ConventionEvaluation const & conventions,
        std::vector<double> * weighted = nullptr,
        std::vector<double> * bare = nullptr
    ) const;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if (version == 0) {
            if (tuning_count_ != 0) {
                throw std::runtime_error("Finish or reset mixture tuning before serialization");
            }
            archive(::cereal::make_nvp("Channels", channels_));
            archive(::cereal::make_nvp("Weights", weights_));
            // Retain the version-0 wire field for existing readers.
            bool legacy_allow_incompatible = false;
            archive(::cereal::make_nvp(
                "AllowIncompatible", legacy_allow_incompatible));
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
        bool legacy_allow_incompatible = false;
        archive(::cereal::make_nvp("Channels", loaded_channels));
        archive(::cereal::make_nvp("Weights", loaded_weights));
        archive(::cereal::make_nvp(
            "AllowIncompatible", legacy_allow_incompatible));

        MultiChannelPhaseSpace loaded;
        loaded.channels_ = std::move(loaded_channels);
        loaded.weights_ = std::move(loaded_weights);
        // Preserve exact weights for deterministic continuation. The legacy
        // compatibility override never changes validation on load.
        loaded.ValidateArchivedWeights();
        loaded.RequireCompatibleConvention();
        *this = std::move(loaded);
    }
};

// Probe peer densities at points sampled by each channel.
std::vector<std::string> ValidateChannelDensities(
    MultiChannelPhaseSpace const & mixture,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord template_record,
    int samples_per_channel = 100);

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
