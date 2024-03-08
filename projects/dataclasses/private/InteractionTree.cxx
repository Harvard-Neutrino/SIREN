#include "SIREN/dataclasses/InteractionTree.h"

#include <memory>

namespace SI {
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
    tree.insert(datum);
    return datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionTreeDatum& datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> _datum = std::make_shared<InteractionTreeDatum>(datum);
    if (parent) {
        _datum->parent = parent;
        parent->daughters.push_back(_datum);
    }
    tree.insert(_datum);
    return _datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionRecord& record,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    if (parent) {
        datum->parent = parent;
        parent->daughters.push_back(datum);
    }
    tree.insert(datum);
    return datum;
}

} // namespace dataclasses
} // namespace SI
