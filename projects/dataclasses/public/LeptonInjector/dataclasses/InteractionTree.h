#pragma once
#ifndef LI_InteractionTree_H
#define LI_InteractionTree_H

#endif // LI_InteractionTree_H

#include "LeptonInjector/dataclasses/InteractionRecord.h"


namespace LI {
namespace dataclasses {

struct InteractionTreeDatum {
  InteractionTreeDatum(dataclasses::InteractionRecord record) : record(record) {}
  dataclasses::InteractionRecord record;
  std::shared_ptr<dataclasses::InteractionTreeDatum> parent;
  std::vector<std::shared_ptr<dataclasses::InteractionTreeDatum>> daughters;
};

struct InteractionTree {
  std::set<std::shared_ptr<dataclasses::InteractionTreeDatum>> tree;
  void add_entry(dataclasses::InteractionRecord& record,
                 std::shared_ptr<dataclasses::InteractionTreeDatum> parent);
};

} // namespace dataclasses
} // namespace LI
