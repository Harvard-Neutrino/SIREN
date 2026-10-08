#pragma once
#ifndef SIREN_PrimaryExternalDistribution_H
#define SIREN_PrimaryExternalDistribution_H

#include <memory>                                        // for shared_ptr
#include <array>                                         // for array
#include <set>                                           // for set
#include <string>                                        // for string
#include <tuple>                                         // for tuple
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>

#include "SIREN/distributions/Distributions.h"
#include "SIREN/distributions/DistributionVariable.h"
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }
namespace cereal { class access; }

namespace siren {
namespace distributions {

// Table-driven primary distribution: each row fixes the primary's kinematics
// (and optionally its vertex, initial position, time, and extra interaction
// parameters).
//
// Rows are sampled UNIFORMLY unless explicit sampling weights are supplied:
// the supplied parent list is taken to be the intended injection ensemble
// (a dk2nu production's importance reweighting deliberately oversamples the
// phase space it cares about, and injection respects that choice).
//
// A column named "weight" declares per-row PHYSICAL weights: row i represents
// weight_i physical primaries per unit exposure (for dk2nu tables,
// nimpwt/POT). The weights encode how to return from the supplied ensemble to
// the physical one: on the physical side of a process the distribution
// reports the physical row density through PhysicalDensity (never cancelled
// against the injection side, whose GenerationProbability stays the sampling
// density) and the weight total through the PhysicallyNormalizedDistribution
// interface. One instance shared between the injection and physical sides is
// therefore exact for uniform sampling and for any explicit sampling bias.
// Tables without a "weight" column keep the legacy semantics: flat physical
// density, no normalization, injection/physical cancellation as usual.
//
// SEGMENT MODE (explicit opt-in via SetSegmentColumn / the segment_column
// constructor argument): the named column (metres) turns each row from a
// point primary into a straight track segment starting at x0/y0/z0 along the
// momentum direction (Geant4 photon steps, for example). The table must then
// carry x0/y0/z0 and px/py/pz and must not carry x/y/z: the interaction
// vertex is sampled uniformly along the segment, this distribution declares
// InteractionVertex and therefore owns the injection bounds, and
// InjectionBounds returns the segment end points so the weighter integrates
// the interaction probability and the normalized position density along the
// segment. GenerationProbability carries the extra factor 1/length (the
// sampled longitudinal density) and is evaluated on THIS table's support:
// zero for a record whose initial position, direction or vertex is not on
// the row's segment, so pooled proposals over overlapping segments combine
// correctly. DensityVariables then also reports
// "PrimaryPositionLongitudinal" while PhysicalDensityVariables does not,
// because the physical longitudinal density comes from the weighter's
// normalized position factor. The owning process must use
// VertexWeightingMode::ExternalBounds(); the Injector rejects any other mode,
// and rejects ExternalBounds() for point tables. A column merely named
// "length" is ordinary metadata unless opted in, and archives written before
// segment mode existed (versions 0 and 1) always load as point tables.
//
// ROW LAYOUT: records cache the sampled row index, and every table asked to
// evaluate a record reads that index into its own rows, with row densities
// relative to its own uniform-over-rows measure. Tables pooled in one Weighter
// or paired across its injection/physical sides must therefore share a row
// layout (RowLayoutMismatch): the same primaries at the same indices. The
// Weighter rejects other combinations at configuration.
class PrimaryExternalDistribution : virtual public VertexPositionDistribution,
                                    virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
protected:
    PrimaryExternalDistribution() {};
    void LoadInputFile(std::string const & _filename);
private:
    std::string filename;
    std::vector<std::vector<double>> input_data;
    std::vector<std::string> keys;
    std::vector<double> sampling_weights_;
    double sampling_weights_sum_ = 0;
    std::vector<double> sampling_cdf_;
    // Physical row weights from the "weight" column (empty without one).
    // Derived state: rebuilt from input_data on construction, refiltering,
    // and deserialization; never archived.
    std::vector<double> physical_weights_;
    double physical_weights_sum_ = 0;
    bool init_pos_set = false;
    bool vertex_set = false;
    bool mom_set = false;
    // Segment mode: the opted-in column name is archived (version 2); the
    // flag and index are derived from it and the keys on construction and
    // deserialization.
    std::string segment_column_;
    bool length_set = false;
    size_t length_index_ = 0;
    double emin = 0;
    std::set<DistributionVariable> set_variables_;
    mutable std::array<double, 3> _cached_position = {0.0, 0.0, 0.0};
    void DeriveColumnFlags();
    void ComputeSetVariables();
    void BuildSamplingCDF();
    void DerivePhysicalRowWeights();
    void ValidateSegmentLengths() const;
    double SegmentLength(siren::dataclasses::InteractionRecord const & record) const;
    double RowSamplingDensity(siren::dataclasses::InteractionRecord const & record) const;
    bool OnOwnSegment(siren::dataclasses::InteractionRecord const & record, double & length) const;
public:
    PrimaryExternalDistribution(std::string _filename);
    PrimaryExternalDistribution(std::string _filename, double emin);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data, double emin);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data, std::vector<double> _sampling_weights);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data, std::vector<double> _sampling_weights, double emin);
    PrimaryExternalDistribution(PrimaryExternalDistribution const & other) = default;
    size_t GetPhysicalNumEvents() const;
    void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual double PhysicalDensity(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual bool PhysicalDensityDiffers() const override;
    virtual bool ProvidesExternalBounds() const override;
    virtual std::vector<std::string> PhysicalDensityVariables() const override;
    // Opt into segment mode using the named length column (metres); an empty
    // name returns to point semantics. Re-derives declared variables and
    // validates the table.
    void SetSegmentColumn(std::string const & column);
    std::string const & GetSegmentColumn() const { return segment_column_; }
    // Empty when the two tables may evaluate each other's records (same row
    // count and the same primary at every row index, ignoring the "weight"
    // column and either table's segment length column); otherwise a
    // description of the first difference. Weighter requires this of every
    // external table it pools or pairs across the injection/physical sides.
    std::string RowLayoutMismatch(PrimaryExternalDistribution const & other) const;
    virtual std::set<DistributionVariable> SetVariables() const override;
    virtual std::set<DistributionVariable> RequiredVariables() const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
        } else if(version == 1) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            archive(::cereal::make_nvp("SamplingWeights", sampling_weights_));
        } else if(version == 2) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            archive(::cereal::make_nvp("SamplingWeights", sampling_weights_));
            archive(::cereal::make_nvp("SegmentColumn", segment_column_));
        } else {
            throw std::runtime_error("PrimaryExternalDistribution only supports version <= 2!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            segment_column_.clear();
            ComputeSetVariables();
            DerivePhysicalRowWeights();
        } else if(version == 1) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            archive(::cereal::make_nvp("SamplingWeights", sampling_weights_));
            // Segment mode did not exist: a column named "length" is metadata.
            segment_column_.clear();
            ComputeSetVariables();
            DerivePhysicalRowWeights();
            BuildSamplingCDF();
        } else if(version == 2) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            archive(::cereal::make_nvp("SamplingWeights", sampling_weights_));
            archive(::cereal::make_nvp("SegmentColumn", segment_column_));
            ComputeSetVariables();
            DerivePhysicalRowWeights();
            BuildSamplingCDF();
            ValidateSegmentLengths();
        } else {
            throw std::runtime_error("PrimaryExternalDistribution only supports version <= 2!");
        }
    }
private:
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PrimaryExternalDistribution, 2);
CEREAL_REGISTER_TYPE(siren::distributions::PrimaryExternalDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::PrimaryExternalDistribution);

#endif // SIREN_PrimaryExternalDistribution_H
