#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleMasses.h"
#include "SIREN/dataclasses/PhaseSpaceConvention.h"

#include <algorithm>
#include <cctype>
#include <iostream>

namespace siren {
namespace interactions {

CrossSection::CrossSection() {}

void CrossSection::SampleFinalState(dataclasses::InteractionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> rand) const {
    siren::dataclasses::CrossSectionDistributionRecord csdr(record);
    this->SampleFinalState(csdr, rand);
    csdr.Finalize(record);
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
    return siren::dataclasses::MeasureFromConvention(Convention());
}

siren::dataclasses::PhaseSpaceConvention CrossSection::Convention() const {
    using C = siren::dataclasses::PhaseSpaceConvention;
    auto variables = DensityVariables();

    auto lower_match = [](std::vector<std::string> const & vars,
                          std::string const & needle) {
        for (auto const & v : vars) {
            std::string lv = v;
            std::transform(lv.begin(), lv.end(), lv.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lv.find(needle) != std::string::npos) return true;
        }
        return false;
    };
    auto lower_eq = [](std::vector<std::string> const & vars,
                       std::string const & needle) {
        for (auto const & v : vars) {
            std::string lv = v;
            std::transform(lv.begin(), lv.end(), lv.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lv == needle) return true;
        }
        return false;
    };

    bool has_bjorken_x =
        lower_match(variables, "bjorken x") ||
        lower_match(variables, "bjorken_x") ||
        lower_eq(variables, "x");
    bool has_bjorken_y =
        lower_match(variables, "bjorken y") ||
        lower_match(variables, "bjorken_y") ||
        lower_eq(variables, "y");
    bool has_q2 =
        lower_match(variables, "q^2") ||
        lower_match(variables, "mandelstam") ||
        lower_eq(variables, "q2") ||
        lower_eq(variables, "t") ||
        lower_eq(variables, "-t");

    C result = C::Custom;
    if (has_bjorken_x && has_bjorken_y) {
        result = C::BjorkenXY;
    } else if (has_q2) {
        result = C::MandelstamST;
    } else if (has_bjorken_y) {
        result = C::BjorkenXY;
    }

    std::cerr << "Warning: CrossSection subclass does not override Convention(); "
              << "auto-detected "
              << siren::dataclasses::PhaseSpaceConventionName(result)
              << " from DensityVariables(). Override Convention() to silence "
              << "this warning." << std::endl;
    return result;
}

} // namespace interactions
} // namespace siren
