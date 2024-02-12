#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"

namespace LI {
namespace interactions {

CrossSection::CrossSection() {}

void CrossSection::SampleFinalState(dataclasses::InteractionRecord & record, std::shared_ptr<LI::utilities::LI_random> rand) const {
    LI::dataclasses::CrossSectionDistributionRecord csdr(record);
    this->SampleFinalState(csdr, rand);
    csdr.Finalize(record);
}

double CrossSection::TotalCrossSectionAllFinalStates(LI::dataclasses::InteractionRecord const & record) const {
    std::vector<LI::dataclasses::InteractionSignature> signatures = this->GetPossibleSignaturesFromParents(record.signature.primary_type, record.signature.target_type);
    LI::dataclasses::InteractionRecord fake_record = record;
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
} // namespace LI

