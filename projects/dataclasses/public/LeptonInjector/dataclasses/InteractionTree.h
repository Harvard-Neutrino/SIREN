#pragma once
#ifndef LI_InteractionTree_H
#define LI_InteractionTree_H

#endif // LI_InteractionTree_H

#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>


namespace LI {
namespace dataclasses {

struct InteractionTreeDatum {
  InteractionTreeDatum(dataclasses::InteractionRecord& record) : record(record) {}
  dataclasses::InteractionRecord record;
  std::shared_ptr<dataclasses::InteractionTreeDatum> parent;
  std::vector<std::shared_ptr<dataclasses::InteractionTreeDatum>> daughters;
  int depth() const;
  template<class Archive>
  void serialize(Archive & archive, std::uint32_t const version) {
      if(version == 0) {
          archive(::cereal::make_nvp("Record", record));
          archive(::cereal::make_nvp("Parent", parent));
          archive(::cereal::make_nvp("Daughters", daughters));
      } else {
          throw std::runtime_error("InteractionTreeDatum only supports version <= 0!");
      }
  };
};

struct InteractionTree {
  std::set<std::shared_ptr<dataclasses::InteractionTreeDatum>> tree;
  std::shared_ptr<InteractionTreeDatum> add_entry(dataclasses::InteractionRecord& record,
                                                  std::shared_ptr<dataclasses::InteractionTreeDatum> parent = NULL);
  template<class Archive>
  void serialize(Archive & archive, std::uint32_t const version) {
      if(version == 0) {
          archive(::cereal::make_nvp("Tree", tree));
      } else {
          throw std::runtime_error("InteractionTree only supports version <= 0!");
      }
  };
};

} // namespace dataclasses
} // namespace LI
