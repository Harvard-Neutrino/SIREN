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

bool CrossSection::operator==(CrossSection const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

} // namespace interactions
} // namespace LI

