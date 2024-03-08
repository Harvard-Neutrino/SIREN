#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"

namespace SI {
namespace interactions {

CrossSection::CrossSection() {}

void CrossSection::SampleFinalState(dataclasses::InteractionRecord & record, std::shared_ptr<SI::utilities::LI_random> rand) const {
    SI::dataclasses::CrossSectionDistributionRecord csdr(record);
    this->SampleFinalState(csdr, rand);
    csdr.Finalize(record);
}

double CrossSection::TotalCrossSectionAllFinalStates(SI::dataclasses::InteractionRecord const & record) const {
    std::vector<SI::dataclasses::InteractionSignature> signatures = this->GetPossibleSignaturesFromParents(record.signature.primary_type, record.signature.target_type);
    SI::dataclasses::InteractionRecord fake_record = record;
    double total_cross_section = 0;
    for(auto signature : signatures) {
        fake_record.signature = signature;
        total_cross_section += this->TotalCrossSection(record);
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
} // namespace SI

