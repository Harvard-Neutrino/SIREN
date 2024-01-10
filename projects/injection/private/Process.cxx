#include "LeptonInjector/injection/Process.h"

#include <tuple>

namespace LI {
namespace injection {

bool Process::operator==(Process const & other) const {
    return std::tie(
        primary_type,
        cross_sections)
        ==
        std::tie(
        other.primary_type,
        other.cross_sections);
}

bool Process::MatchesHead(std::shared_ptr<Process> const & other) const {
    return std::tie(
        primary_type,
        cross_sections)
        ==
        std::tie(
        other->primary_type,
        other->cross_sections);
}

} // namespace injection
} // namespace LI
