#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleMasses.h"
#include "SIREN/dataclasses/PhaseSpaceConvention.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <iostream>

namespace siren {
namespace interactions {

namespace {

// Infer the scattering measure a model's density is differential in from its
// DensityVariables() strings. Explicit azimuth variables select the *Phi
// joint measures; a (Q2, y) pair is the two-dimensional MandelstamQ2Y, not
// the one-dimensional MandelstamQ2.
siren::dataclasses::PhaseSpaceMeasure InferMeasureFromDensityVariables(
    std::vector<std::string> const & variables)
{
    using M = siren::dataclasses::PhaseSpaceMeasure;

    auto lower_match = [&variables](std::string const & needle) {
        for (auto const & v : variables) {
            std::string lv = v;
            std::transform(lv.begin(), lv.end(), lv.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lv.find(needle) != std::string::npos) return true;
        }
        return false;
    };
    auto lower_eq = [&variables](std::string const & needle) {
        for (auto const & v : variables) {
            std::string lv = v;
            std::transform(lv.begin(), lv.end(), lv.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lv == needle) return true;
        }
        return false;
    };

    bool has_bjorken_x =
        lower_match("bjorken x") ||
        lower_match("bjorken_x") ||
        lower_eq("x");
    bool has_bjorken_y =
        lower_match("bjorken y") ||
        lower_match("bjorken_y") ||
        lower_eq("y");
    bool has_q2 =
        lower_match("q^2") ||
        lower_match("mandelstam") ||
        lower_eq("q2") ||
        lower_eq("t") ||
        lower_eq("-t");
    bool has_phi =
        lower_match("azimuth") ||
        lower_eq("phi");

    if (has_bjorken_x && has_bjorken_y)
        return has_phi ? M::BjorkenXYPhi() : M::BjorkenXY();
    if (has_q2 && has_bjorken_y)
        return has_phi ? M::MandelstamQ2YPhi() : M::MandelstamQ2Y();
    if (has_q2)
        return has_phi ? M::MandelstamQ2Phi() : M::MandelstamQ2();
    if (has_bjorken_y) {
        // A y-only density for a fixed-mass 2->2 final state is differential
        // in dy, not in the two-dimensional DIS measure dx dy.
        return has_phi ? M::FixedMassYPhi() : M::FixedMassY();
    }
    return M::Unspecified();
}

} // anonymous namespace

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
    using M = siren::dataclasses::PhaseSpaceMeasure;
    auto variables = DensityVariables();

    // Warn at most once per process that this base fallback is being used.
    static std::atomic<bool> warned{false};
    auto warn_once = [](std::string const & detail) {
        if (!warned.exchange(true)) {
            std::cerr << "Warning: CrossSection subclass does not override Measure(); "
                      << detail << " Override Measure() to silence this warning."
                      << std::endl;
        }
    };

    // A malformed signature set (empty, or of mixed secondary count) has no well-defined
    // topology or measure; fall back to Unspecified as Topology() does.
    auto signatures = GetPossibleSignatures();
    if (signatures.empty()) {
        warn_once("auto-detected Unspecified (no signatures defined).");
        return M::Unspecified();
    }
    size_t n = signatures.front().secondary_types.size();
    for (auto const & sig : signatures) {
        if (sig.secondary_types.size() != n) {
            warn_once("auto-detected Unspecified (mixed secondary counts).");
            return M::Unspecified();
        }
    }

    M result = InferMeasureFromDensityVariables(variables);

    warn_once("auto-detected " +
              siren::dataclasses::PhaseSpaceMeasureName(result) +
              " from DensityVariables().");
    return result;
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
    siren::dataclasses::InteractionSignature const & signature) const
{
    using M = siren::dataclasses::PhaseSpaceMeasure;
    size_t n = signature.secondary_types.size();
    if (n != 2 && n != 3) return M::Unspecified();

    return InferMeasureFromDensityVariables(DensityVariables());
}

} // namespace interactions
} // namespace siren
