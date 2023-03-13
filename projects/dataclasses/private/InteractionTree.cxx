#include "LeptonInjector/dataclasses/InteractionTree.h"

namespace LI {
namespace dataclasses {
  
int InteractionTreeDatum::depth() const {
  int depth = 0;
  std::shared_ptr<InteractionTreeDatum> test = parent;
  while(true) {
    if(!test) return depth;
    test = test->parent;
    ++depth;
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
