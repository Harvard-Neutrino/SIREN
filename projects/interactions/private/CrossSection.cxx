#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"

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

bool CrossSection::operator==(CrossSection const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

} // namespace interactions
} // namespace siren
