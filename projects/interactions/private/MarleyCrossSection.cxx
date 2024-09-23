#include "SIREN/interactions/MarleyCrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "marley/Generator.hh"  // MARLEY generator
#include "marley/JSON.hh"  // For handling MARLEY JSON configuration


namespace siren {
namespace interactions {

//MarleyCrossSection::MarleyCrossSection() {}
MarleyCrossSection::MarleyCrossSection(const std::string& marley_config) {
    InitializeMarley(marley_config);
}

void MarleyCrossSection::InitializeMarley(const std::string& marley_config) {
    marley::JSON marley_json = marley::JSON::load(marley_config);
    marley::JSONConfig config(marley_json);
    marley_generator_ = config.create_generator();
}

void CrossSection::SampleFinalState(dataclasses::InteractionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> rand) const {
    siren::dataclasses::CrossSectionDistributionRecord csdr(record);
    this->SampleFinalState(csdr, rand);
    csdr.Finalize(record);
}

double MarleyCrossSection::TotalCrossSection(siren::dataclasses::InteractionRecord const & record) const {

    int pdg_a = record.GetPrimaryPDG();  // PDG code of the primary particle
    double KEa = record.GetPrimaryKineticEnergy();  // Kinetic energy of the primary particle
    int pdg_atom = record.GetTargetPDG();  // PDG code of the target atom

    return marley_generator_.total_xs(pdg_a, KEa, pdg_atom);
}

//double CrossSection::TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const {
//    std::vector<siren::dataclasses::InteractionSignature> signatures = this->GetPossibleSignaturesFromParents(record.signature.primary_type, record.signature.target_type);
//    siren::dataclasses::InteractionRecord fake_record = record;
//    double total_cross_section = 0;
//    for(auto signature : signatures) {
//        fake_record.signature = signature;
//        total_cross_section += this->TotalCrossSection(fake_record);
//    }
//    return total_cross_section;
//}


bool MarleyCrossSection::operator==(CrossSection const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

} // namespace interactions
} // namespace siren

