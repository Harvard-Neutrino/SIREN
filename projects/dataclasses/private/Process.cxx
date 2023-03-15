#include "LeptonInjector/dataclasses/Process.h"

namespace LI {
namespace dataclasses {

bool Process::operator==(Process const & other) const {
    return std::tie(
        primary_type,
        cross_sections)
        ==
        std::tie(
        other.primary_type,
        other.cross_sections);
}

} // namespace dataclasses
} // namespace LI
