#include "LeptonInjector/dataclasses/InteractionTree.h"

namespace LI {
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

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionRecord& record,
                                                                 std::shared_ptr<InteractionTreeDatum> parent) {
  std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
  datum->parent = parent;
  if (parent) {
    parent->daughters.push_back(datum);
  }
  tree.insert(datum);
  return datum;
}

} // namespace dataclasses
} // namespace LI
