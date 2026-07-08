#include "SIREN/dataclasses/InteractionTree.h"

#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <ostream>
#include <unordered_map>

namespace siren {
namespace dataclasses {

int InteractionTreeDatum::depth(InteractionTree const & owner) const {
    int depth = 0;
    std::uint32_t p = parent_index;
    // An acyclic parent chain visits at most tree.size() nodes; exceeding that
    // means the indices form a cycle (a corrupt archive), so stop rather than
    // spin forever.
    while(p != kNoParent) {
        ++depth;
        if(static_cast<std::size_t>(depth) > owner.tree.size())
            throw std::runtime_error("InteractionTreeDatum::depth: parent_index chain exceeds tree size");
        p = owner.tree.at(p)->parent_index;
    }
    return depth;
}

bool InteractionTreeDatum::operator==(InteractionTreeDatum const & other) const {
    return node_id == other.node_id
        and parent_index == other.parent_index
        and daughter_indices == other.daughter_indices
        and record == other.record;
}

bool InteractionTreeHeader::operator==(InteractionTreeHeader const & other) const {
    return event_number == other.event_number
        and weights == other.weights
        and provenance == other.provenance;
}

bool InteractionTree::operator==(InteractionTree const & other) const {
    if(not (header == other.header)) return false;
    if(tree.size() != other.tree.size()) return false;
    for(std::size_t i=0; i<tree.size(); ++i) {
        bool const lhs = bool(tree[i]);
        bool const rhs = bool(other.tree[i]);
        if(lhs != rhs) return false;
        if(lhs and not (*tree[i] == *other.tree[i])) return false;
    }
    return true;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(std::shared_ptr<InteractionTreeDatum> datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    datum->node_id = static_cast<std::uint32_t>(tree.size());
    datum->daughter_indices.clear();
    if (parent) {
        datum->parent_index = parent->node_id;
        parent->daughter_indices.push_back(datum->node_id);
    } else {
        datum->parent_index = kNoParent;
    }
    tree.push_back(datum);
    // Enforced invariant: node_id is exactly the position in tree.
    assert(tree.back()->node_id == tree.size() - 1);
    return datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionTreeDatum& datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> _datum = std::make_shared<InteractionTreeDatum>(datum);
    return add_entry(_datum, parent);
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionRecord& record,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    return add_entry(datum, parent);
}

void InteractionTree::rebuild_indices_from_legacy() {
    // tree is already populated in parent-before-daughter order by the legacy
    // writer. Map each node's shared_ptr identity to its position, then derive
    // node_id/parent_index/daughter_indices from the legacy pointer edges.
    std::unordered_map<InteractionTreeDatum const *, std::uint32_t> index_of;
    index_of.reserve(tree.size());
    for(std::size_t i=0; i<tree.size(); ++i) {
        index_of[tree[i].get()] = static_cast<std::uint32_t>(i);
    }
    for(std::size_t i=0; i<tree.size(); ++i) {
        std::shared_ptr<InteractionTreeDatum> & datum = tree[i];
        datum->node_id = static_cast<std::uint32_t>(i);
        if(datum->legacy_parent) {
            auto it = index_of.find(datum->legacy_parent.get());
            datum->parent_index = (it != index_of.end()) ? it->second : kNoParent;
        } else {
            datum->parent_index = kNoParent;
        }
        datum->daughter_indices.clear();
        for(auto const & d : datum->legacy_daughters) {
            auto it = index_of.find(d.get());
            if(it != index_of.end()) datum->daughter_indices.push_back(it->second);
        }
    }
    // Break the legacy shared_ptr cycle so nothing leaks into the running process.
    for(std::shared_ptr<InteractionTreeDatum> & datum : tree) {
        datum->legacy_parent = nullptr;
        datum->legacy_daughters.clear();
    }
    validate_indices();
}

void InteractionTree::validate_indices() const {
    std::uint32_t const n = static_cast<std::uint32_t>(tree.size());
    for(std::uint32_t i = 0; i < n; ++i) {
        std::shared_ptr<InteractionTreeDatum> const & datum = tree[i];
        if(not datum)
            throw std::runtime_error("InteractionTree: null datum in tree vector");
        if(datum->node_id != i)
            throw std::runtime_error("InteractionTree: node_id does not equal its position");
        if(datum->parent_index != kNoParent) {
            if(datum->parent_index >= n)
                throw std::runtime_error("InteractionTree: parent_index out of range");
            if(not tree[datum->parent_index])
                throw std::runtime_error("InteractionTree: parent_index points to null datum");
        }
        for(std::uint32_t d : datum->daughter_indices) {
            if(d >= n)
                throw std::runtime_error("InteractionTree: daughter_index out of range");
            if(not tree[d])
                throw std::runtime_error("InteractionTree: daughter_index points to null datum");
        }
    }
}

void SaveInteractionTrees(std::vector<std::shared_ptr<InteractionTree>>& trees, std::string const & filename) {
    std::string const path = filename + ".siren_events";
    std::ofstream os(path, std::ios::binary);
    if(!os.is_open()) {
        throw std::runtime_error(
            "Cannot open '" + path + "' for writing. Does the output "
            "directory exist?");
    }
    ::cereal::BinaryOutputArchive archive(os);
    archive(trees);
}

std::vector<std::shared_ptr<InteractionTree>> LoadInteractionTrees(std::string const & filename) {
    std::string const path = filename + ".siren_events";
    std::ifstream is(path, std::ios::binary);
    if(!is.is_open()) {
        throw std::runtime_error(
            "Cannot open '" + path + "' for reading. Does the file exist?");
    }
    ::cereal::BinaryInputArchive archive(is);
    std::vector<std::shared_ptr<InteractionTree>> trees;
    archive(trees);
    return trees;
}

} // namespace dataclasses
} // namespace siren

std::string to_str(siren::dataclasses::InteractionTreeDatum const & datum) {
    using siren::dataclasses::kNoParent;
    std::stringstream ss;
    ss << "[ InteractionTreeDatum (" << &datum << ")\n";
    ss << "  node_id: " << datum.node_id << '\n';
    if(datum.parent_index == kNoParent) ss << "  parent_index: (root)\n";
    else ss << "  parent_index: " << datum.parent_index << '\n';
    ss << "  daughter_indices: [";
    for(std::size_t i=0; i<datum.daughter_indices.size(); ++i) {
        ss << (i ? ", " : "") << datum.daughter_indices[i];
    }
    ss << "]\n";
    ss << "  record: " << to_str(datum.record) << '\n';
    ss << ']';
    return ss.str();
}

std::string to_str(siren::dataclasses::InteractionTreeHeader const & header) {
    std::stringstream ss;
    ss << "[ InteractionTreeHeader (" << &header << ")\n";
    ss << "  event_number: " << header.event_number << '\n';
    ss << "  weights: [";
    for(std::size_t i=0; i<header.weights.size(); ++i) ss << (i ? ", " : "") << header.weights[i];
    ss << "]\n";
    ss << "  provenance: {";
    bool first = true;
    for(auto const & kv : header.provenance) {
        ss << (first ? "" : ", ") << kv.first << ": " << kv.second;
        first = false;
    }
    ss << "}\n]";
    return ss.str();
}

std::string to_str(siren::dataclasses::InteractionTree const & tree) {
    std::stringstream ss;
    ss << "[ InteractionTree (" << &tree << ")\n";
    ss << "  " << to_str(tree.header) << '\n';
    ss << "  nodes: " << tree.tree.size() << '\n';
    ss << ']';
    return ss.str();
}

std::string to_repr(siren::dataclasses::InteractionTreeDatum const & datum) {
    using siren::dataclasses::kNoParent;
    std::stringstream ss;
    ss << "InteractionTreeDatum(node_id=" << datum.node_id << ", parent_index=";
    if(datum.parent_index == kNoParent) ss << "None"; else ss << datum.parent_index;
    ss << ", daughters=" << datum.daughter_indices.size() << ")";
    return ss.str();
}

std::string to_repr(siren::dataclasses::InteractionTreeHeader const & header) {
    std::stringstream ss;
    ss << "InteractionTreeHeader(event_number=" << header.event_number
       << ", weights=" << header.weights.size()
       << ", provenance=" << header.provenance.size() << ")";
    return ss.str();
}

std::string to_repr(siren::dataclasses::InteractionTree const & tree) {
    std::stringstream ss;
    ss << "InteractionTree(nodes=" << tree.tree.size()
       << ", event_number=" << tree.header.event_number << ")";
    return ss.str();
}

std::ostream & operator<<(std::ostream & os, siren::dataclasses::InteractionTreeDatum const & datum) {
    os << to_str(datum); return os;
}
std::ostream & operator<<(std::ostream & os, siren::dataclasses::InteractionTreeHeader const & header) {
    os << to_str(header); return os;
}
std::ostream & operator<<(std::ostream & os, siren::dataclasses::InteractionTree const & tree) {
    os << to_str(tree); return os;
}
