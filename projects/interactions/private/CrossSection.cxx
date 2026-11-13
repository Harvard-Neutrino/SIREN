#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleMasses.h"
#include "SIREN/dataclasses/PhaseSpaceConvention.h"

namespace siren {
namespace interactions {

CrossSection::CrossSection() {}

void CrossSection::SampleFinalState(dataclasses::InteractionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> rand) const {
    siren::dataclasses::CrossSectionDistributionRecord csdr(record);
    this->SampleFinalState(csdr, rand);
    csdr.Finalize(record);
}

double CrossSection::SampleInteractionTime(siren::dataclasses::CrossSectionDistributionRecord const & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // Identity default: keep the flight-time value already on the record.
    return record.GetInteractionTime();
}

double CrossSection::TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const {
    std::vector<siren::dataclasses::InteractionSignature> signatures = this->GetPossibleSignaturesFromParents(record.signature.primary_type, record.signature.target_type);
    siren::dataclasses::InteractionRecord fake_record = record;
    double total_cross_section = 0;
    for(auto signature : signatures) {
        fake_record.signature = signature;
        total_cross_section += this->TotalCrossSection(fake_record);
    }
    return total_cross_section;
}

std::vector<double> CrossSection::SecondaryMasses(dataclasses::InteractionRecord const & record) const {
    return SecondaryMasses(record.signature.secondary_types);
}

std::vector<double> CrossSection::SecondaryMasses(std::vector<dataclasses::ParticleType> const & secondary_types) const {
    std::vector<double> masses;
    masses.reserve(secondary_types.size());
    for(auto const & type : secondary_types) {
        masses.push_back(siren::dataclasses::GetParticleMass(type));
    }
    return masses;
}

std::vector<double> CrossSection::SecondaryHelicities(dataclasses::InteractionRecord const & record) const {
    return std::vector<double>(record.signature.secondary_types.size(), 0.0);
}

bool CrossSection::operator==(CrossSection const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

siren::dataclasses::PhaseSpaceTopology CrossSection::Topology() const {
    using T = siren::dataclasses::PhaseSpaceTopology;
    auto signatures = GetPossibleSignatures();
    if (signatures.empty()) return T::Unspecified;
    size_t n = signatures.front().secondary_types.size();
    for (auto const & sig : signatures) {
        if (sig.secondary_types.size() != n) return T::Unspecified;
    }
    if (n == 2) return T::Scatter2to2;
    if (n == 3) return T::Scatter2to3;
    return T::Unspecified;
}

siren::dataclasses::PhaseSpaceMeasure CrossSection::Measure() const {
    return siren::dataclasses::PhaseSpaceMeasure::Unspecified();
}

siren::dataclasses::PhaseSpaceTopology CrossSection::TopologyForSignature(
    siren::dataclasses::InteractionSignature const & signature) const
{
    using T = siren::dataclasses::PhaseSpaceTopology;
    size_t n = signature.secondary_types.size();
    if (n == 2) return T::Scatter2to2;
    if (n == 3) return T::Scatter2to3;
    return T::Unspecified;
}

siren::dataclasses::PhaseSpaceMeasure CrossSection::MeasureForSignature(
    siren::dataclasses::InteractionSignature const &) const
{
    return Measure();
}

} // namespace interactions
} // namespace siren
