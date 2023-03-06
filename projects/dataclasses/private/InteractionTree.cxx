#include "LeptonInjector/dataclasses/InteractionTree.h"

namespace LI {
namespace dataclasses {
  
void InteractionTree::add_entry(dataclasses::InteractionRecord& record,
                                std::shared_ptr<dataclassess:InteractionTreeDatum> parent = NULL) {
  std::shared_ptr<dataclasses::InteractionTreeDatum> datum = std::make_shared<dataclasses::InteractionTreeDatum>(record);
  datum->parent = parent;
  if (parent) {
    parent->daughters.append(datum);
  }
}

} // namespace dataclasses
} // namespace LI
