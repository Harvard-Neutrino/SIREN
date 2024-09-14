#pragma once
#ifndef SIREN_StringMapulation_H
#define SIREN_StringMapulation_H

#include <string>

namespace siren {
namespace utilities {

constexpr char const * tab = "  ";

std::string add_prefix(std::string const & input, std::string const & prefix);
std::string indent(std::string const & input, size_t n_indent = 1);

} // namespace utilities
} // namespace siren

#endif // SIREN_StringMapulation_H
