#include <string>
#include <vector>
#include <sstream>

#include "SIREN/utilities/StringManipulation.h"

namespace siren {
namespace utilities {

std::string add_prefix(std::string const & input, std::string const & prefix) {
    std::istringstream iss(input);
    std::vector<std::string> lines;
    std::string line;
    ssize_t last_non_empty_line = -1;
    size_t line_number = 0;

    // Read each line and track the last non-empty line
    while (std::getline(iss, line)) {
        lines.push_back(line);
        if (!line.empty()) {
            last_non_empty_line = line_number;
        }
        line_number++;
    }

    std::ostringstream oss;

    // Add prefix to each line up to the last non-empty line
    if (last_non_empty_line >= 0) {
        for (size_t i = 0; i <= static_cast<size_t>(last_non_empty_line); ++i) {
            oss << prefix << lines[i] << '\n';
        }
        // Ensure the string ends with an empty newline
        oss << '\n';
    }

    return oss.str();
}

} // namespace utilities
} // namespace siren
