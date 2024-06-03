#pragma once
#ifndef SIREN_InteractionTree_H
#define SIREN_InteractionTree_H

#include "SIREN/dataclasses/InteractionRecord.h"

#include <set>                                             // for set
#include <memory>                                          // for shared_ptr
#include <vector>                                          // for vector
#include <cstdint>                                         // for uint32_t
#include <stddef.h>                                        // for NULL
#include <stdexcept>                                       // for runtime_error
#include <fstream>                                         // for if/ofstream

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace siren {
namespace dataclasses {

struct InteractionTreeDatum {
  InteractionTreeDatum() {}
  InteractionTreeDatum(dataclasses::InteractionRecord& record) : record(record) {}
  dataclasses::InteractionRecord record;
  std::shared_ptr<dataclasses::InteractionTreeDatum> parent = NULL;
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
  std::vector<std::shared_ptr<dataclasses::InteractionTreeDatum>> tree;
  std::shared_ptr<InteractionTreeDatum> add_entry(std::shared_ptr<dataclasses::InteractionTreeDatum> datum,
                                                  std::shared_ptr<dataclasses::InteractionTreeDatum> parent = NULL);
  std::shared_ptr<InteractionTreeDatum> add_entry(dataclasses::InteractionTreeDatum& datum,
                                                  std::shared_ptr<dataclasses::InteractionTreeDatum> parent = NULL);
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

void SaveInteractionTrees(std::vector<std::shared_ptr<InteractionTree>>& trees, std::string const & filename);
std::vector<std::shared_ptr<InteractionTree>> LoadInteractionTrees(std::string const & filename);

} // namespace dataclasses
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::InteractionTreeDatum, 0);
CEREAL_CLASS_VERSION(siren::dataclasses::InteractionTree, 0);

#endif // SIREN_InteractionTree_H

