#pragma once
#ifndef SIREN_InteractionTree_H
#define SIREN_InteractionTree_H

#include "SIREN/dataclasses/InteractionRecord.h"

#include <map>                                             // for map
#include <set>                                             // for set
#include <string>                                          // for string
#include <limits>                                          // for numeric_limits
#include <memory>                                          // for shared_ptr
#include <vector>                                          // for vector
#include <cstdint>                                         // for uint32_t
#include <iosfwd>                                          // for ostream
#include <stddef.h>                                        // for NULL
#include <stdexcept>                                       // for runtime_error
#include <fstream>                                         // for if/ofstream

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace siren { namespace dataclasses { struct InteractionTreeDatum; } }
namespace siren { namespace dataclasses { struct InteractionTreeHeader; } }
namespace siren { namespace dataclasses { struct InteractionTree; } }

std::ostream & operator<<(std::ostream & os, siren::dataclasses::InteractionTreeDatum const & datum);
std::ostream & operator<<(std::ostream & os, siren::dataclasses::InteractionTreeHeader const & header);
std::ostream & operator<<(std::ostream & os, siren::dataclasses::InteractionTree const & tree);

std::string to_str(siren::dataclasses::InteractionTreeDatum const & datum);
std::string to_str(siren::dataclasses::InteractionTreeHeader const & header);
std::string to_str(siren::dataclasses::InteractionTree const & tree);

std::string to_repr(siren::dataclasses::InteractionTreeDatum const & datum);
std::string to_repr(siren::dataclasses::InteractionTreeHeader const & header);
std::string to_repr(siren::dataclasses::InteractionTree const & tree);

namespace siren {
namespace dataclasses {

// Sentinel value for InteractionTreeDatum::parent_index / node_id meaning "no
// parent" (i.e. a root/primary interaction). Chosen as the max uint32_t so that
// a valid index is never confused with "unset".
static constexpr std::uint32_t kNoParent = std::numeric_limits<std::uint32_t>::max();

struct InteractionTreeDatum {
  InteractionTreeDatum() {}
  InteractionTreeDatum(dataclasses::InteractionRecord& record) : record(record) {}

  dataclasses::InteractionRecord record;

  // Stable identity == position of this datum in the owning InteractionTree::tree
  // vector. Assigned by InteractionTree::add_entry and never changed afterwards.
  std::uint32_t node_id = kNoParent;
  // Index of the parent interaction within the owning tree, or kNoParent for a
  // root. An integer back-reference cannot form an ownership cycle, and depth()
  // can walk it without copying ancestors.
  std::uint32_t parent_index = kNoParent;
  // Indices of the daughter interactions within the owning tree.
  std::vector<std::uint32_t> daughter_indices;

  bool is_root() const { return parent_index == kNoParent; }

  // Depth of this datum below the root, resolving parent_index against its owning
  // tree. O(depth) integer pointer-chase, zero record copies. There is
  // deliberately NO zero-argument depth(): a datum holds no back-pointer to its
  // tree (that would recreate the ownership cycle), so depth requires the owner.
  int depth(InteractionTree const & owner) const;

  bool operator==(InteractionTreeDatum const & other) const;

  // Legacy shared_ptr edges. Populated ONLY while loading a version-0 archive so
  // that InteractionTree::load can rebuild the integer indices; cleared
  // immediately afterwards and never written by the version-1 save path.
  std::shared_ptr<dataclasses::InteractionTreeDatum> legacy_parent = nullptr;
  std::vector<std::shared_ptr<dataclasses::InteractionTreeDatum>> legacy_daughters;

  template<class Archive>
  void save(Archive & archive, std::uint32_t const version) const {
      if(version == 1) {
          archive(::cereal::make_nvp("Record", record));
          archive(::cereal::make_nvp("NodeId", node_id));
          archive(::cereal::make_nvp("ParentIndex", parent_index));
          archive(::cereal::make_nvp("DaughterIndices", daughter_indices));
      } else {
          throw std::runtime_error("InteractionTreeDatum save only supports version 1!");
      }
  }
  template<class Archive>
  void load(Archive & archive, std::uint32_t const version) {
      if(version == 1) {
          archive(::cereal::make_nvp("Record", record));
          archive(::cereal::make_nvp("NodeId", node_id));
          archive(::cereal::make_nvp("ParentIndex", parent_index));
          archive(::cereal::make_nvp("DaughterIndices", daughter_indices));
      } else if(version == 0) {
          // Legacy layout: Record, Parent (shared_ptr), Daughters (vector<shared_ptr>).
          // cereal pointer-tracking resolves the aliased nodes; indices are
          // reconstructed by InteractionTree::load once every node exists.
          archive(::cereal::make_nvp("Record", record));
          archive(::cereal::make_nvp("Parent", legacy_parent));
          archive(::cereal::make_nvp("Daughters", legacy_daughters));
      } else {
          throw std::runtime_error("InteractionTreeDatum only supports version <= 1!");
      }
  }
};

// Event-level metadata for a tree. Maps onto HepMC3 GenEvent (event_number,
// weights) and GenRunInfo (provenance) at export time.
struct InteractionTreeHeader {
  std::uint64_t event_number = 0;
  std::vector<double> weights;                    // per-injector; names carried in provenance
  std::map<std::string, std::string> provenance;  // generator name/version, run tag, etc.

  bool operator==(InteractionTreeHeader const & other) const;

  template<class Archive>
  void serialize(Archive & archive, std::uint32_t const version) {
      if(version == 0) {
          archive(::cereal::make_nvp("EventNumber", event_number));
          archive(::cereal::make_nvp("Weights", weights));
          archive(::cereal::make_nvp("Provenance", provenance));
      } else {
          throw std::runtime_error("InteractionTreeHeader only supports version 0!");
      }
  }
};

struct InteractionTree {
  // Flat vector of all nodes, parents strictly before their daughters. The
  // enforced invariant is that every datum's node_id equals its position here.
  std::vector<std::shared_ptr<dataclasses::InteractionTreeDatum>> tree;
  InteractionTreeHeader header;

  std::shared_ptr<InteractionTreeDatum> add_entry(std::shared_ptr<dataclasses::InteractionTreeDatum> datum,
                                                  std::shared_ptr<dataclasses::InteractionTreeDatum> parent = nullptr);
  std::shared_ptr<InteractionTreeDatum> add_entry(dataclasses::InteractionTreeDatum& datum,
                                                  std::shared_ptr<dataclasses::InteractionTreeDatum> parent = nullptr);
  std::shared_ptr<InteractionTreeDatum> add_entry(dataclasses::InteractionRecord& record,
                                                  std::shared_ptr<dataclasses::InteractionTreeDatum> parent = nullptr);

  std::shared_ptr<InteractionTreeDatum> const & at(std::uint32_t node_id) const { return tree.at(node_id); }
  int depth(std::uint32_t node_id) const { return tree.at(node_id)->depth(*this); }

  bool operator==(InteractionTree const & other) const;

  template<class Archive>
  void save(Archive & archive, std::uint32_t const version) const {
      if(version == 1) {
          archive(::cereal::make_nvp("Tree", tree));
          archive(::cereal::make_nvp("Header", header));
      } else {
          throw std::runtime_error("InteractionTree save only supports version 1!");
      }
  }
  template<class Archive>
  void load(Archive & archive, std::uint32_t const version) {
      if(version == 1) {
          archive(::cereal::make_nvp("Tree", tree));
          archive(::cereal::make_nvp("Header", header));
      } else if(version == 0) {
          archive(::cereal::make_nvp("Tree", tree));
          rebuild_indices_from_legacy();
      } else {
          throw std::runtime_error("InteractionTree only supports version <= 1!");
      }
  }

private:
  // Reconstruct node_id/parent_index/daughter_indices from the legacy shared_ptr
  // edges after a version-0 load, then clear the legacy pointers so no reference
  // cycle survives into the running process.
  void rebuild_indices_from_legacy();
};

void SaveInteractionTrees(std::vector<std::shared_ptr<InteractionTree>>& trees, std::string const & filename);
std::vector<std::shared_ptr<InteractionTree>> LoadInteractionTrees(std::string const & filename);

} // namespace dataclasses
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::InteractionTreeDatum, 1);
CEREAL_CLASS_VERSION(siren::dataclasses::InteractionTreeHeader, 0);
CEREAL_CLASS_VERSION(siren::dataclasses::InteractionTree, 1);

#endif // SIREN_InteractionTree_H
