#include "SIREN/dataclasses/InteractionTree.h"

#include <memory>

namespace siren {
namespace dataclasses {

int InteractionTreeDatum::depth() const {
    int depth = 0;
    if(parent==NULL) return depth;
    std::shared_ptr<InteractionTreeDatum> test = std::make_shared<InteractionTreeDatum>(*parent);
    while(true) {
        ++depth;
        if(test->parent==NULL) return depth;
        test = std::make_shared<InteractionTreeDatum>(*(test->parent));
    }
    return -1;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(std::shared_ptr<InteractionTreeDatum> datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    if (parent) {
        datum->parent = parent;
        parent->daughters.push_back(datum);
    }
    tree.push_back(datum);
    return datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionTreeDatum& datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> _datum = std::make_shared<InteractionTreeDatum>(datum);
    if (parent) {
        _datum->parent = parent;
        parent->daughters.push_back(_datum);
    }
    tree.push_back(_datum);
    return _datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionRecord& record,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    if (parent) {
        datum->parent = parent;
        parent->daughters.push_back(datum);
    }
    tree.push_back(datum);
    return datum;
}

void SaveInteractionTrees(std::vector<std::shared_ptr<InteractionTree>>& trees, std::string const & filename) {
    std::ofstream os(filename+".siren_events", std::ios::binary);
    ::cereal::BinaryOutputArchive archive(os);
    archive(trees);
}

std::vector<std::shared_ptr<InteractionTree>> LoadInteractionTrees(std::string const & filename) {
    std::ifstream is(filename+".siren_events", std::ios::binary);
    ::cereal::BinaryInputArchive archive(is);
    std::vector<std::shared_ptr<InteractionTree>> trees;
    archive(trees);
    return trees;
}

} // namespace dataclasses
} // namespace siren
