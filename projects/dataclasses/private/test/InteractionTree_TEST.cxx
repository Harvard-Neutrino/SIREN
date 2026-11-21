
#include <math.h>
#include <cmath>
#include <random>
#include <sstream>
#include <iostream>
#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/dataclasses/InteractionTree.h"

using namespace siren::dataclasses;

// Legacy (version-0) mirror of the pointer-edge container layout: parent/daughters are
// shared_ptr edges. Used only to write a genuine version-0 ".siren_events"
// archive so the real (version-1) LoadInteractionTrees can be exercised on its
// backward-compatibility branch (rebuild_indices_from_legacy).
namespace legacy_v0 {
struct LegacyDatum {
    siren::dataclasses::InteractionRecord record;
    std::shared_ptr<LegacyDatum> parent = nullptr;
    std::vector<std::shared_ptr<LegacyDatum>> daughters;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Record", record));
            archive(::cereal::make_nvp("Parent", parent));
            archive(::cereal::make_nvp("Daughters", daughters));
        } else {
            throw std::runtime_error("legacy_v0::LegacyDatum only writes version 0!");
        }
    }
};
struct LegacyTree {
    std::vector<std::shared_ptr<LegacyDatum>> tree;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Tree", tree));
        } else {
            throw std::runtime_error("legacy_v0::LegacyTree only writes version 0!");
        }
    }
};
} // namespace legacy_v0

CEREAL_CLASS_VERSION(legacy_v0::LegacyDatum, 0);
CEREAL_CLASS_VERSION(legacy_v0::LegacyTree, 0);

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

// Helper: build an InteractionRecord with a fixed EPlus/EMinus signature so the
// tests below can populate node records without repeating the boilerplate.
InteractionRecord MakeRecord() {
    InteractionRecord record;
    record.signature.primary_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.target_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);
    return record;
}

TEST(DatumConstructor, InteractionRecord)
{
    InteractionRecord record;
    record.signature.primary_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.target_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);

    InteractionTreeDatum datum(record);
    EXPECT_EQ(datum.record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum.record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum.record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum.record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
}

TEST(DatumConstructor, DefaultsAreRoot)
{
    // A freshly constructed datum has no owning tree yet; node_id and
    // parent_index both default to the kNoParent sentinel and it reads as root.
    InteractionTreeDatum datum;
    EXPECT_EQ(datum.node_id, kNoParent);
    EXPECT_EQ(datum.parent_index, kNoParent);
    EXPECT_TRUE(datum.daughter_indices.empty());
    EXPECT_TRUE(datum.is_root());
}

TEST(TreeConstructor, Default)
{
    InteractionTree tree;
    EXPECT_EQ(tree.tree.size(), 0);
}

TEST(TreeAddEntry, Record)
{
    InteractionTree tree;
    InteractionRecord record = MakeRecord();

    std::shared_ptr<InteractionTreeDatum> datum = tree.add_entry(record);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
    // The first entry is a root at position 0.
    EXPECT_EQ(datum->node_id, 0u);
    EXPECT_EQ(datum->parent_index, kNoParent);
    EXPECT_TRUE(datum->is_root());
}

TEST(TreeAddEntry, DatumReference)
{
    InteractionTree tree;
    InteractionRecord record = MakeRecord();

    InteractionTreeDatum datum = InteractionTreeDatum(record);
    std::shared_ptr<InteractionTreeDatum> datum2 = tree.add_entry(datum);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum2->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
    EXPECT_EQ(datum2->node_id, 0u);
}

TEST(TreeAddEntry, DatumPointer)
{
    InteractionTree tree;
    InteractionRecord record = MakeRecord();

    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    std::shared_ptr<InteractionTreeDatum> datum2 = tree.add_entry(datum);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum2->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
    // add_entry with a shared_ptr stores that very object and assigns its node_id.
    EXPECT_EQ(datum2, datum);
    EXPECT_EQ(datum->node_id, 0u);
}

TEST(TreeAddEntry, Parent)
{
    InteractionTree tree;
    InteractionRecord record = MakeRecord();

    auto datum = tree.add_entry(record);
    auto datum2 = tree.add_entry(record, datum);
    EXPECT_EQ(tree.tree.size(), 2);
    EXPECT_EQ(datum2->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);

    // New integer-edge API: the child points back at the parent by node_id, and
    // the parent lists the child among its daughter_indices.
    EXPECT_EQ(datum2->parent_index, datum->node_id);
    EXPECT_EQ(datum2->node_id, 1u);
    EXPECT_FALSE(datum2->is_root());
    EXPECT_TRUE(datum->is_root());
    EXPECT_EQ(datum->daughter_indices.size(), 1u);
    EXPECT_EQ(datum->daughter_indices[0], datum2->node_id);
}

TEST(TreeInvariant, NodeIdEqualsPosition)
{
    // Build a small multi-node tree and assert every datum's node_id matches its
    // position in the flat tree vector.
    InteractionTree tree;
    InteractionRecord record = MakeRecord();

    auto root = tree.add_entry(record);
    auto childA = tree.add_entry(record, root);
    auto childB = tree.add_entry(record, root);
    auto grandchild = tree.add_entry(record, childA);
    (void)childB;
    (void)grandchild;

    EXPECT_EQ(tree.tree.size(), 4u);
    for(std::uint32_t k=0; k<tree.tree.size(); ++k) {
        EXPECT_EQ(tree.tree[k]->node_id, k);
    }
}

TEST(TreeStructure, ParentIndexAndDaughters)
{
    // Root + two children + one grandchild under the first child.
    InteractionTree tree;
    InteractionRecord record = MakeRecord();

    auto root = tree.add_entry(record);
    auto childA = tree.add_entry(record, root);
    auto childB = tree.add_entry(record, root);
    auto grandchild = tree.add_entry(record, childA);

    // Positions / node_ids.
    EXPECT_EQ(root->node_id, 0u);
    EXPECT_EQ(childA->node_id, 1u);
    EXPECT_EQ(childB->node_id, 2u);
    EXPECT_EQ(grandchild->node_id, 3u);

    // Root: no parent, is_root, both children recorded as daughters.
    EXPECT_EQ(root->parent_index, kNoParent);
    EXPECT_TRUE(root->is_root());
    ASSERT_EQ(root->daughter_indices.size(), 2u);
    EXPECT_EQ(root->daughter_indices[0], childA->node_id);
    EXPECT_EQ(root->daughter_indices[1], childB->node_id);

    // Children point back at the root.
    EXPECT_EQ(childA->parent_index, root->node_id);
    EXPECT_FALSE(childA->is_root());
    EXPECT_EQ(childB->parent_index, root->node_id);
    EXPECT_FALSE(childB->is_root());

    // childA owns the grandchild; childB has none.
    ASSERT_EQ(childA->daughter_indices.size(), 1u);
    EXPECT_EQ(childA->daughter_indices[0], grandchild->node_id);
    EXPECT_TRUE(childB->daughter_indices.empty());

    // Grandchild points back at childA.
    EXPECT_EQ(grandchild->parent_index, childA->node_id);
    EXPECT_FALSE(grandchild->is_root());
    EXPECT_TRUE(grandchild->daughter_indices.empty());
}

TEST(TreeDepth, DatumAndTreeAgree)
{
    // Root(0) -> childA(1) -> grandchild(3); childB(2) is a sibling of childA.
    InteractionTree tree;
    InteractionRecord record = MakeRecord();

    auto root = tree.add_entry(record);
    auto childA = tree.add_entry(record, root);
    auto childB = tree.add_entry(record, root);
    auto grandchild = tree.add_entry(record, childA);

    // Depth via the datum (resolving against its owner).
    EXPECT_EQ(root->depth(tree), 0);
    EXPECT_EQ(childA->depth(tree), 1);
    EXPECT_EQ(childB->depth(tree), 1);
    EXPECT_EQ(grandchild->depth(tree), 2);

    // Depth via the tree convenience wrapper keyed by node_id.
    EXPECT_EQ(tree.depth(root->node_id), 0);
    EXPECT_EQ(tree.depth(childA->node_id), 1);
    EXPECT_EQ(tree.depth(childB->node_id), 1);
    EXPECT_EQ(tree.depth(grandchild->node_id), 2);

    // at() returns the same handles that add_entry produced.
    EXPECT_EQ(tree.at(grandchild->node_id), grandchild);
}

// Build a fixed 2-level tree (root + two children) for equality tests.
static InteractionTree BuildSampleTree() {
    InteractionTree tree;
    InteractionRecord record = MakeRecord();
    auto root = tree.add_entry(record);
    tree.add_entry(record, root);
    tree.add_entry(record, root);
    tree.header.event_number = 42;
    tree.header.weights.push_back(1.5);
    tree.header.provenance["generator"] = "siren-test";
    return tree;
}

TEST(DatumEquality, IdenticalAndMutations)
{
    InteractionTree a = BuildSampleTree();
    InteractionTree b = BuildSampleTree();

    // Independently built identical data compare equal.
    ASSERT_EQ(a.tree.size(), b.tree.size());
    EXPECT_TRUE(*a.tree[1] == *b.tree[1]);

    // Mutating node_id breaks equality.
    {
        InteractionTreeDatum mutated = *b.tree[1];
        mutated.node_id = 99u;
        EXPECT_FALSE(*a.tree[1] == mutated);
    }
    // Mutating parent_index breaks equality.
    {
        InteractionTreeDatum mutated = *b.tree[1];
        mutated.parent_index = kNoParent;
        EXPECT_FALSE(*a.tree[1] == mutated);
    }
    // Mutating daughter_indices breaks equality.
    {
        InteractionTreeDatum mutated = *b.tree[0];
        mutated.daughter_indices.push_back(123u);
        EXPECT_FALSE(*a.tree[0] == mutated);
    }
    // Mutating the record breaks equality.
    {
        InteractionTreeDatum mutated = *b.tree[1];
        mutated.record.signature.primary_type = siren::dataclasses::ParticleType::EMinus;
        EXPECT_FALSE(*a.tree[1] == mutated);
    }
}

TEST(HeaderEquality, IdenticalAndMutations)
{
    InteractionTreeHeader h1;
    h1.event_number = 7;
    h1.weights = {1.0, 2.0};
    h1.provenance["gen"] = "v1";

    InteractionTreeHeader h2 = h1;
    EXPECT_TRUE(h1 == h2);

    h2.event_number = 8;
    EXPECT_FALSE(h1 == h2);

    h2 = h1;
    h2.weights.push_back(3.0);
    EXPECT_FALSE(h1 == h2);

    h2 = h1;
    h2.provenance["gen"] = "v2";
    EXPECT_FALSE(h1 == h2);
}

TEST(TreeEquality, IdenticalAndMutations)
{
    InteractionTree a = BuildSampleTree();
    InteractionTree b = BuildSampleTree();

    // Two independently-built identical trees compare equal.
    EXPECT_TRUE(a == b);

    // Mutating one node's node_id makes them unequal.
    {
        InteractionTree c = BuildSampleTree();
        c.tree[1]->node_id = 99u;
        EXPECT_FALSE(a == c);
    }
    // Mutating a node's parent_index makes them unequal.
    {
        InteractionTree c = BuildSampleTree();
        c.tree[1]->parent_index = kNoParent;
        EXPECT_FALSE(a == c);
    }
    // Mutating a node's daughter_indices makes them unequal.
    {
        InteractionTree c = BuildSampleTree();
        c.tree[0]->daughter_indices.push_back(123u);
        EXPECT_FALSE(a == c);
    }
    // Mutating a node's record makes them unequal.
    {
        InteractionTree c = BuildSampleTree();
        c.tree[1]->record.signature.primary_type = siren::dataclasses::ParticleType::EMinus;
        EXPECT_FALSE(a == c);
    }
    // Mutating the header makes them unequal.
    {
        InteractionTree c = BuildSampleTree();
        c.header.event_number = 43;
        EXPECT_FALSE(a == c);
    }
}

TEST(TreeSerialization, JSONRoundTripV1)
{
    // Populate a header (event_number, one weight, one provenance entry) and a
    // 2-level tree, then round-trip through a cereal JSON archive.
    InteractionTree original;
    InteractionRecord record = MakeRecord();
    auto root = original.add_entry(record);
    original.add_entry(record, root);
    original.add_entry(record, root);
    original.header.event_number = 12345;
    original.header.weights.push_back(0.75);
    original.header.provenance["generator"] = "siren-serialization-test";

    std::stringstream ss;
    {
        cereal::JSONOutputArchive oarchive(ss);
        oarchive(original);
    }

    InteractionTree loaded;
    {
        cereal::JSONInputArchive iarchive(ss);
        iarchive(loaded);
    }

    // Records, integer edges, and header all survive the round-trip.
    EXPECT_TRUE(loaded == original);
    EXPECT_EQ(loaded.tree.size(), original.tree.size());
    EXPECT_EQ(loaded.header.event_number, original.header.event_number);
    ASSERT_EQ(loaded.header.weights.size(), 1u);
    EXPECT_EQ(loaded.header.weights[0], 0.75);
    EXPECT_EQ(loaded.header.provenance.at("generator"), "siren-serialization-test");

    // The reconstructed node_id == position invariant still holds after load.
    for(std::uint32_t k=0; k<loaded.tree.size(); ++k) {
        EXPECT_EQ(loaded.tree[k]->node_id, k);
    }
}

TEST(LegacyLoad, V0FileRebuildsIndices)
{
    // Build a version-0 tree with shared_ptr edges:
    //   root -> {child0, child1};  child0 -> {grandchild}
    using legacy_v0::LegacyDatum;
    using legacy_v0::LegacyTree;
    auto root = std::make_shared<LegacyDatum>(); root->record = MakeRecord();
    auto child0 = std::make_shared<LegacyDatum>(); child0->record = MakeRecord();
    auto child1 = std::make_shared<LegacyDatum>(); child1->record = MakeRecord();
    auto grandchild = std::make_shared<LegacyDatum>(); grandchild->record = MakeRecord();
    child0->parent = root; child1->parent = root; grandchild->parent = child0;
    root->daughters = {child0, child1};
    child0->daughters = {grandchild};

    auto ltree = std::make_shared<LegacyTree>();
    ltree->tree = {root, child0, child1, grandchild}; // parent-before-daughter order
    std::vector<std::shared_ptr<LegacyTree>> ltrees = {ltree};

    // Write it exactly like SaveInteractionTrees does (binary, ".siren_events").
    std::string base = std::string(::testing::TempDir()) + "siren_v0_fixture";
    {
        std::ofstream os(base + ".siren_events", std::ios::binary);
        ::cereal::BinaryOutputArchive archive(os);
        archive(ltrees);
    }

    // Load with the real version-1 loader; it must take the version-0 branch.
    std::vector<std::shared_ptr<InteractionTree>> loaded = LoadInteractionTrees(base);
    std::remove((base + ".siren_events").c_str());

    ASSERT_EQ(loaded.size(), 1u);
    InteractionTree & t = *loaded[0];
    ASSERT_EQ(t.tree.size(), 4u);

    // node_id == position invariant reconstructed
    for(std::uint32_t k=0; k<t.tree.size(); ++k) EXPECT_EQ(t.tree[k]->node_id, k);

    // parent_index topology reconstructed from the legacy pointer edges
    EXPECT_TRUE(t.tree[0]->is_root());
    EXPECT_EQ(t.tree[0]->parent_index, kNoParent);
    EXPECT_EQ(t.tree[1]->parent_index, 0u);
    EXPECT_EQ(t.tree[2]->parent_index, 0u);
    EXPECT_EQ(t.tree[3]->parent_index, 1u);

    // daughter_indices reconstructed
    ASSERT_EQ(t.tree[0]->daughter_indices.size(), 2u);
    EXPECT_EQ(t.tree[0]->daughter_indices[0], 1u);
    EXPECT_EQ(t.tree[0]->daughter_indices[1], 2u);
    ASSERT_EQ(t.tree[1]->daughter_indices.size(), 1u);
    EXPECT_EQ(t.tree[1]->daughter_indices[0], 3u);
    EXPECT_TRUE(t.tree[2]->daughter_indices.empty());
    EXPECT_TRUE(t.tree[3]->daughter_indices.empty());

    // depth resolves through the rebuilt indices
    EXPECT_EQ(t.tree[3]->depth(t), 2);
    EXPECT_EQ(t.depth(3), 2);

    // legacy shared_ptr edges cleared: no reference cycle survives the load
    for(auto const & d : t.tree) {
        EXPECT_EQ(d->legacy_parent, nullptr);
        EXPECT_TRUE(d->legacy_daughters.empty());
    }

    // default header for an upgraded file
    EXPECT_EQ(t.header.event_number, 0u);
    EXPECT_TRUE(t.header.weights.empty());
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
