#pragma once
#ifndef SIREN_StringMapulation_H
#define SIREN_StringMapulation_H

#include <string>

namespace siren {
namespace utilities {

constexpr char const * tab = "  ";

std::string add_prefix(std::string const & input, std::string const & prefix);

} // namespace utilities
} // namespace siren

#endif // SIREN_StringMapulation_H
