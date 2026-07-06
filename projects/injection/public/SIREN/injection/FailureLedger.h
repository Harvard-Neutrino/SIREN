#ifndef SIREN_FailureLedger_H
#define SIREN_FailureLedger_H

#include <map>
#include <string>
#include <cstdint>

#include "SIREN/utilities/Errors.h"

namespace siren {
namespace injection {

// Aggregated failure accounting keyed by (depth, parent pdg, FailureReason).
// Storage and per-record cost are O(distinct reasons), not O(failed events):
// repeated identical failures increment a count and keep a single exemplar
// message. Transient runtime diagnostics only; never serialized.
struct FailureLedger {
    struct Key {
        int depth;
        int parent_pdg;
        siren::utilities::FailureReason reason;
        bool operator<(Key const & other) const {
            if(depth != other.depth) return depth < other.depth;
            if(parent_pdg != other.parent_pdg) return parent_pdg < other.parent_pdg;
            return reason < other.reason;
        }
    };
    struct Entry {
        uint64_t count = 0;
        std::string exemplar;
    };

    std::map<Key, Entry> entries;

    void Record(int depth, int parent_pdg, siren::utilities::FailureReason reason,
                std::string const & message) {
        Entry & entry = entries[Key{depth, parent_pdg, reason}];
        entry.count += 1;
        if(entry.exemplar.empty()) {
            entry.exemplar = message;
        }
    }

    void Clear() {
        entries.clear();
    }
};

} // namespace injection
} // namespace siren

#endif // SIREN_FailureLedger_H
