#include "LeptonInjector/dataclasses/InteractionTree.h"

namespace LI {
namespace dataclasses {
  
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
