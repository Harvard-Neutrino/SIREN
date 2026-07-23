#include "SIREN/interactions/TrivialCrossSection.h"

#include <tuple>                                              // for tie
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <algorithm>                                          // for find
#include <stdexcept>                                          // for runtime...

#include "SIREN/interactions/CrossSection.h"         // for CrossSe...
#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle
#include "SIREN/utilities/Random.h"                  // for SIREN_r...

namespace siren {
namespace interactions {

TrivialCrossSection::TrivialCrossSection(double cross_section_cm2,
        std::vector<siren::dataclasses::ParticleType> primary_types,
        std::vector<siren::dataclasses::ParticleType> target_types)
    : targets_(std::move(target_types)) {
    if(not (cross_section_cm2 > 0))
        throw std::runtime_error("TrivialCrossSection: the cross section must be positive");
    if(primary_types.empty())
        throw std::runtime_error("TrivialCrossSection: need at least one primary type");
    if(targets_.empty())
        throw std::runtime_error("TrivialCrossSection: need at least one target type");
    for(auto primary : primary_types) {
        energies_[primary] = {0.0, 1e9};
        cross_sections_[primary] = {cross_section_cm2, cross_section_cm2};
    }
}

TrivialCrossSection::TrivialCrossSection(std::map<siren::dataclasses::ParticleType,
            std::pair<std::vector<double>, std::vector<double>>> tables,
        std::vector<siren::dataclasses::ParticleType> target_types)
    : targets_(std::move(target_types)) {
    if(tables.empty())
        throw std::runtime_error("TrivialCrossSection: need at least one primary table");
    if(targets_.empty())
        throw std::runtime_error("TrivialCrossSection: need at least one target type");
    for(auto & entry : tables) {
        std::vector<double> & energies = entry.second.first;
        std::vector<double> & cross_sections = entry.second.second;
        if(energies.empty())
            throw std::runtime_error("TrivialCrossSection: a table must have at least one point");
        if(energies.size() != cross_sections.size())
            throw std::runtime_error("TrivialCrossSection: energies and cross sections must have the same length");
        for(size_t i = 0; i + 1 < energies.size(); ++i) {
            if(not (energies[i] < energies[i + 1]))
                throw std::runtime_error("TrivialCrossSection: energies must be strictly ascending");
        }
        for(double xs : cross_sections) {
            if(xs < 0)
                throw std::runtime_error("TrivialCrossSection: cross sections must be non-negative");
        }
        energies_[entry.first] = std::move(energies);
        cross_sections_[entry.first] = std::move(cross_sections);
    }
}

bool TrivialCrossSection::equal(CrossSection const & other) const {
    const TrivialCrossSection* x = dynamic_cast<const TrivialCrossSection*>(&other);

    if(!x)
        return false;
    else
        return std::tie(energies_, cross_sections_, targets_)
            == std::tie(x->energies_, x->cross_sections_, x->targets_);
}

double TrivialCrossSection::TotalCrossSection(siren::dataclasses::ParticleType primary, double energy) const {
    auto it = energies_.find(primary);
    if(it == energies_.end())
        return 0.0;
    std::vector<double> const & energies = it->second;
    std::vector<double> const & cross_sections = cross_sections_.at(primary);
    if(energy < energies.front())
        return 0.0;
    if(energy >= energies.back())
        return cross_sections.back();
    size_t i = std::distance(energies.begin(),
            std::upper_bound(energies.begin(), energies.end(), energy)) - 1;
    double t = (energy - energies[i]) / (energies[i + 1] - energies[i]);
    return cross_sections[i] + t * (cross_sections[i + 1] - cross_sections[i]);
}

double TrivialCrossSection::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    return TotalCrossSection(interaction.signature.primary_type, interaction.primary_momentum[0]);
}

double TrivialCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    siren::dataclasses::InteractionSignature const & signature = interaction.signature;
    bool matches = energies_.count(signature.primary_type) > 0
        and std::find(targets_.begin(), targets_.end(), signature.target_type) != targets_.end()
        and signature.secondary_types.size() == 2
        and signature.secondary_types[0] == signature.primary_type
        and signature.secondary_types[1] == signature.target_type;
    if(not matches)
        return 0.0;
    return TotalCrossSection(interaction);
}

double TrivialCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    auto it = energies_.find(interaction.signature.primary_type);
    if(it == energies_.end())
        return 0.0;
    return it->second.front();
}

void TrivialCrossSection::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random>) const {
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();

    siren::dataclasses::SecondaryParticleRecord & neutrino = secondaries[0];
    neutrino.SetFourMomentum({record.primary_momentum[0], record.primary_momentum[1],
                              record.primary_momentum[2], record.primary_momentum[3]});
    neutrino.SetMass(record.primary_mass);
    neutrino.SetHelicity(record.primary_helicity);

    siren::dataclasses::SecondaryParticleRecord & target = secondaries[1];
    target.SetFourMomentum({record.target_mass, 0.0, 0.0, 0.0});
    target.SetMass(record.target_mass);
    target.SetHelicity(record.target_helicity);
}

double TrivialCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<siren::dataclasses::ParticleType> TrivialCrossSection::GetPossibleTargets() const {
    return targets_;
}

std::vector<siren::dataclasses::ParticleType> TrivialCrossSection::GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const {
    if(energies_.count(primary_type) == 0)
        return {};
    return targets_;
}

std::vector<siren::dataclasses::ParticleType> TrivialCrossSection::GetPossiblePrimaries() const {
    std::vector<siren::dataclasses::ParticleType> primaries;
    primaries.reserve(energies_.size());
    for(auto const & entry : energies_)
        primaries.push_back(entry.first);
    return primaries;
}

std::vector<dataclasses::InteractionSignature> TrivialCrossSection::GetPossibleSignatures() const {
    std::vector<siren::dataclasses::InteractionSignature> signatures;
    for(auto const & entry : energies_) {
        for(auto target : targets_) {
            siren::dataclasses::InteractionSignature signature;
            signature.primary_type = entry.first;
            signature.target_type = target;
            signature.secondary_types = {entry.first, target};
            signatures.push_back(signature);
        }
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> TrivialCrossSection::GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const {
    std::vector<siren::dataclasses::InteractionSignature> signatures;
    if(energies_.count(primary_type) == 0)
        return signatures;
    if(std::find(targets_.begin(), targets_.end(), target_type) == targets_.end())
        return signatures;
    siren::dataclasses::InteractionSignature signature;
    signature.primary_type = primary_type;
    signature.target_type = target_type;
    signature.secondary_types = {primary_type, target_type};
    signatures.push_back(signature);
    return signatures;
}

std::vector<std::string> TrivialCrossSection::DensityVariables() const {
    return std::vector<std::string>{};
}

} // namespace interactions
} // namespace siren
