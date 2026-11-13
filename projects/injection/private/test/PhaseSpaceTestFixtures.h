#pragma once

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/utilities/Random.h"

#include <memory>
#include <string>

namespace siren::injection::test {

using siren::dataclasses::InteractionRecord;

inline InteractionRecord TwoBodyDecayRecord() {
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.5};
    record.secondary_momenta = {
        {11.0, 12.0, 13.0, 14.0},
        {21.0, 22.0, 23.0, 24.0}};
    return record;
}

class ConstantChannel final : public PhaseSpaceChannel {
public:
    explicit ConstantChannel(
        double density = 1.0,
        PhaseSpaceTopology topology = PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure measure = PhaseSpaceMeasure::SolidAngleRest())
        : density_(density), topology_(topology), measure_(measure) {}

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord &) const override {}

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const &) const override
    {
        return density_;
    }

    std::string Name() const override { return "Constant"; }
    PhaseSpaceTopology Topology() const override { return topology_; }
    PhaseSpaceMeasure Measure() const override { return measure_; }
    void SetConvention(PhaseSpaceTopology topology, PhaseSpaceMeasure measure) {
        topology_ = topology;
        measure_ = measure;
    }

private:
    double density_;
    PhaseSpaceTopology topology_;
    PhaseSpaceMeasure measure_;
};

} // namespace siren::injection::test
