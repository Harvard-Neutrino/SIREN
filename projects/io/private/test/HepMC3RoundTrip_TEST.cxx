#include <cstdio>
#include <string>
#include <vector>
#include <memory>
#include <gtest/gtest.h>

#include "SIREN/io/HepMC3Writer.h"
#include "SIREN/io/HepMC3Reader.h"
#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleID.h"
#include "SIREN/dataclasses/ParticleType.h"

using namespace siren::dataclasses;

namespace {

// Two-level tree:
//   node0  primary interaction  NuMu + EPlus -> NuMu(1:1) + EMinus(1:2)   (has target)
//   node1  decay of the NuMu(1:1)            NuMu       -> EMinus + NuMu   (no target)
InteractionTree BuildTree() {
    InteractionRecord root;
    root.signature.primary_type = ParticleType::NuMu;
    root.signature.target_type = ParticleType::EPlus;
    root.signature.secondary_types = {ParticleType::NuMu, ParticleType::EMinus};
    root.primary_momentum = {10.0, 1.0, 2.0, 9.5};
    root.primary_mass = 0.0;
    root.primary_helicity = -1.0;
    root.primary_initial_position = {0.5, -0.25, 1.5};
    root.primary_initial_time = 2.0;
    root.target_mass = 0.000511;
    root.target_helicity = 0.5;
    root.interaction_vertex = {1.0, 2.0, 3.0};
    root.interaction_time = 5.0;
    root.secondary_ids = {ParticleID(1, 1), ParticleID(1, 2)};
    root.secondary_momenta = {{6.0, 0.5, 1.0, 5.8}, {4.0, 0.5, 1.0, 3.7}};
    root.secondary_masses = {0.0, 0.000511};
    root.secondary_helicities = {-1.0, 1.0};
    root.secondary_times = {5.0, 5.0};
    root.interaction_parameters = {{"Q2", 0.7}, {"bjorken_x", 0.3}};

    InteractionRecord dec;
    dec.signature.primary_type = ParticleType::NuMu;
    dec.signature.target_type = ParticleType::Decay;   // decay vertex -> no target
    dec.signature.secondary_types = {ParticleType::EMinus, ParticleType::NuMu};
    dec.primary_id = ParticleID(1, 1);                 // == root secondary 0 (shared particle)
    dec.primary_momentum = {6.0, 0.5, 1.0, 5.8};
    dec.primary_mass = 0.0;
    dec.primary_helicity = -1.0;
    // A daughter's initial position is its parent's interaction vertex (as
    // CreateSecondaryRecord sets in real output); non-zero so the round trip is
    // not vacuously satisfied and the reweighting-critical field is exercised.
    dec.primary_initial_position = {1.0, 2.0, 3.0};
    dec.primary_initial_time = 5.0;
    dec.interaction_vertex = {2.0, 3.0, 4.0};
    dec.interaction_time = 6.5;
    dec.secondary_ids = {ParticleID(1, 3), ParticleID(1, 4)};
    dec.secondary_momenta = {{3.0, 0.2, 0.5, 2.9}, {3.0, 0.3, 0.5, 2.9}};
    dec.secondary_masses = {0.000511, 0.0};
    dec.secondary_helicities = {1.0, -1.0};
    dec.secondary_times = {6.5, 6.5};
    dec.interaction_parameters = {{"width", 1.25e-3}};

    InteractionTree tree;
    std::shared_ptr<InteractionTreeDatum> r = tree.add_entry(root);
    tree.add_entry(dec, r);
    tree.header.event_number = 11;
    tree.header.weights = {2.5};
    return tree;
}

void ExpectVecClose(std::array<double, 4> const & a, std::array<double, 4> const & b, double tol) {
    for(int i = 0; i < 4; ++i) EXPECT_NEAR(a[i], b[i], tol) << "component " << i;
}
void ExpectVecClose(std::array<double, 3> const & a, std::array<double, 3> const & b, double tol) {
    for(int i = 0; i < 3; ++i) EXPECT_NEAR(a[i], b[i], tol) << "component " << i;
}

// Compare the reweighting-critical fields (ids are intentionally not stored).
void ExpectRecordClose(InteractionRecord const & orig, InteractionRecord const & got) {
    double const tol = 1e-6;
    EXPECT_EQ(orig.signature.primary_type, got.signature.primary_type);
    EXPECT_EQ(orig.signature.target_type, got.signature.target_type);
    ASSERT_EQ(orig.signature.secondary_types.size(), got.signature.secondary_types.size());
    for(size_t i = 0; i < orig.signature.secondary_types.size(); ++i)
        EXPECT_EQ(orig.signature.secondary_types[i], got.signature.secondary_types[i]);

    ExpectVecClose(orig.primary_momentum, got.primary_momentum, tol);
    EXPECT_NEAR(orig.primary_mass, got.primary_mass, tol);
    EXPECT_NEAR(orig.primary_helicity, got.primary_helicity, tol);
    ExpectVecClose(orig.primary_initial_position, got.primary_initial_position, tol);
    EXPECT_NEAR(orig.primary_initial_time, got.primary_initial_time, tol);
    EXPECT_NEAR(orig.target_mass, got.target_mass, tol);
    EXPECT_NEAR(orig.target_helicity, got.target_helicity, tol);
    ExpectVecClose(orig.interaction_vertex, got.interaction_vertex, tol);
    EXPECT_NEAR(orig.interaction_time, got.interaction_time, tol);

    ASSERT_EQ(orig.secondary_momenta.size(), got.secondary_momenta.size());
    for(size_t i = 0; i < orig.secondary_momenta.size(); ++i)
        ExpectVecClose(orig.secondary_momenta[i], got.secondary_momenta[i], tol);
    ASSERT_EQ(orig.secondary_masses.size(), got.secondary_masses.size());
    for(size_t i = 0; i < orig.secondary_masses.size(); ++i)
        EXPECT_NEAR(orig.secondary_masses[i], got.secondary_masses[i], tol);
    ASSERT_EQ(orig.secondary_helicities.size(), got.secondary_helicities.size());
    for(size_t i = 0; i < orig.secondary_helicities.size(); ++i)
        EXPECT_NEAR(orig.secondary_helicities[i], got.secondary_helicities[i], tol);
    ASSERT_EQ(orig.secondary_times.size(), got.secondary_times.size());
    for(size_t i = 0; i < orig.secondary_times.size(); ++i)
        EXPECT_NEAR(orig.secondary_times[i], got.secondary_times[i], tol);

    ASSERT_EQ(orig.interaction_parameters.size(), got.interaction_parameters.size());
    for(auto const & kv : orig.interaction_parameters) {
        auto it = got.interaction_parameters.find(kv.first);
        ASSERT_NE(it, got.interaction_parameters.end()) << "missing key " << kv.first;
        EXPECT_NEAR(kv.second, it->second, tol) << "key " << kv.first;
    }
}

} // namespace

TEST(HepMC3RoundTrip, RecordsAndTopology) {
    InteractionTree tree = BuildTree();

    std::string const path = std::string(::testing::TempDir()) + "siren_hepmc3_roundtrip.hepmc3";
    std::vector<std::shared_ptr<InteractionTree>> trees = {std::make_shared<InteractionTree>(tree)};
    siren::io::SaveInteractionTreesAsHepMC3(trees, path);

    std::vector<std::shared_ptr<InteractionTree>> loaded =
        siren::io::LoadInteractionTreesFromHepMC3(path);
    std::remove(path.c_str());

    ASSERT_EQ(loaded.size(), 1u);
    InteractionTree const & lt = *loaded[0];

    // header
    EXPECT_EQ(lt.header.event_number, 11u);
    ASSERT_FALSE(lt.header.weights.empty());
    EXPECT_NEAR(lt.header.weights[0], 2.5, 1e-9);

    // topology: two nodes, node0 root, node1 child of node0
    ASSERT_EQ(lt.tree.size(), 2u);
    EXPECT_TRUE(lt.tree[0]->is_root());
    EXPECT_EQ(lt.tree[0]->node_id, 0u);
    EXPECT_FALSE(lt.tree[1]->is_root());
    EXPECT_EQ(lt.tree[1]->parent_index, 0u);
    ASSERT_EQ(lt.tree[0]->daughter_indices.size(), 1u);
    EXPECT_EQ(lt.tree[0]->daughter_indices[0], 1u);

    // records match the reweight-critical fields
    ExpectRecordClose(tree.tree[0]->record, lt.tree[0]->record);
    ExpectRecordClose(tree.tree[1]->record, lt.tree[1]->record);
}

namespace {

InteractionRecord MkRec(ParticleType prim, ParticleType tgt, ParticleID prim_id,
                        std::vector<ParticleType> const & sec_types,
                        std::vector<ParticleID> const & sec_ids) {
    InteractionRecord r;
    r.signature.primary_type = prim;
    r.signature.target_type = tgt;
    r.signature.secondary_types = sec_types;
    r.primary_id = prim_id;
    r.primary_momentum = {5.0, 1.0, 1.0, 4.7};
    r.interaction_vertex = {0.1, 0.2, 0.3};
    r.interaction_time = 1.0;
    r.secondary_ids = sec_ids;
    for(size_t i = 0; i < sec_types.size(); ++i) {
        r.secondary_momenta.push_back({2.0, 0.1, 0.1, 1.9});
        r.secondary_masses.push_back(0.0);
        r.secondary_helicities.push_back(0.0);
        r.secondary_times.push_back(1.0);
    }
    return r;
}

} // namespace

// A 3-level tree with two daughters at the root, a grandchild, and leaf
// secondaries that do not spawn vertices -- exercises the reader's
// order-independent topology reconstruction.
TEST(HepMC3RoundTrip, DeepTreeTopology) {
    using PT = ParticleType;
    // node0 root -> A(1:1) [->node1], B(1:2) [->node2], leaf(1:3)
    InteractionRecord n0 = MkRec(PT::NuMu, PT::EPlus, ParticleID(),
        {PT::NuMu, PT::NuMu, PT::EMinus}, {ParticleID(1,1), ParticleID(1,2), ParticleID(1,3)});
    // node1 = decay of A(1:1) -> leaf(1:4), C(1:5) [->node3]
    InteractionRecord n1 = MkRec(PT::NuMu, PT::Decay, ParticleID(1,1),
        {PT::EMinus, PT::NuMu}, {ParticleID(1,4), ParticleID(1,5)});
    // node2 = decay of B(1:2) -> leaf(1:6)
    InteractionRecord n2 = MkRec(PT::NuMu, PT::Decay, ParticleID(1,2),
        {PT::EPlus}, {ParticleID(1,6)});
    // node3 = decay of C(1:5) -> leaf(1:7)
    InteractionRecord n3 = MkRec(PT::NuMu, PT::Decay, ParticleID(1,5),
        {PT::NuMu}, {ParticleID(1,7)});

    InteractionTree tree;
    auto d0 = tree.add_entry(n0);
    auto d1 = tree.add_entry(n1, d0);
    tree.add_entry(n2, d0);
    tree.add_entry(n3, d1);
    tree.header.event_number = 99;

    std::string const path = std::string(::testing::TempDir()) + "siren_hepmc3_deep.hepmc3";
    std::vector<std::shared_ptr<InteractionTree>> trees = {std::make_shared<InteractionTree>(tree)};
    siren::io::SaveInteractionTreesAsHepMC3(trees, path);
    std::vector<std::shared_ptr<InteractionTree>> loaded =
        siren::io::LoadInteractionTreesFromHepMC3(path);
    std::remove(path.c_str());

    ASSERT_EQ(loaded.size(), 1u);
    InteractionTree const & lt = *loaded[0];
    ASSERT_EQ(lt.tree.size(), 4u);

    // topology: 0 root; 1,2 children of 0; 3 child of 1
    EXPECT_TRUE(lt.tree[0]->is_root());
    EXPECT_EQ(lt.tree[1]->parent_index, 0u);
    EXPECT_EQ(lt.tree[2]->parent_index, 0u);
    EXPECT_EQ(lt.tree[3]->parent_index, 1u);
    EXPECT_EQ(lt.tree[0]->daughter_indices.size(), 2u);
    EXPECT_EQ(lt.tree[1]->daughter_indices.size(), 1u);
    EXPECT_TRUE(lt.tree[2]->daughter_indices.empty());
    EXPECT_TRUE(lt.tree[3]->daughter_indices.empty());

    // secondary counts survive (3 at root incl. the leaf, 2, 1, 1)
    EXPECT_EQ(lt.tree[0]->record.signature.secondary_types.size(), 3u);
    EXPECT_EQ(lt.tree[1]->record.signature.secondary_types.size(), 2u);
    EXPECT_EQ(lt.tree[2]->record.signature.secondary_types.size(), 1u);
    EXPECT_EQ(lt.tree[3]->record.signature.secondary_types.size(), 1u);
}

// A non-root interaction that has a real target must round-trip through the
// status-22 (non-primary target) encoding, while the root keeps the single
// status-20 target. Guards the writer 20->22 split and the matching reader
// change that treats 20 and 22 as targets.
TEST(HepMC3RoundTrip, NonRootTargetSurvives) {
    using PT = ParticleType;
    // root: NuMu + EPlus -> NuMu(1:1)   [root target -> status 20]
    InteractionRecord root = MkRec(PT::NuMu, PT::EPlus, ParticleID(),
        {PT::NuMu}, {ParticleID(1, 1)});
    // child: NuMu(1:1) + HNucleus -> EMinus(1:2)   [non-root target -> status 22]
    InteractionRecord child = MkRec(PT::NuMu, PT::HNucleus, ParticleID(1, 1),
        {PT::EMinus}, {ParticleID(1, 2)});

    InteractionTree tree;
    auto d0 = tree.add_entry(root);
    tree.add_entry(child, d0);
    tree.header.event_number = 5;

    std::string const path = std::string(::testing::TempDir()) + "siren_hepmc3_nonroot_target.hepmc3";
    std::vector<std::shared_ptr<InteractionTree>> trees = {std::make_shared<InteractionTree>(tree)};
    siren::io::SaveInteractionTreesAsHepMC3(trees, path);
    std::vector<std::shared_ptr<InteractionTree>> loaded =
        siren::io::LoadInteractionTreesFromHepMC3(path);
    std::remove(path.c_str());

    ASSERT_EQ(loaded.size(), 1u);
    InteractionTree const & lt = *loaded[0];
    ASSERT_EQ(lt.tree.size(), 2u);
    EXPECT_TRUE(lt.tree[0]->is_root());
    EXPECT_EQ(lt.tree[0]->record.signature.target_type, PT::EPlus);      // root target survives
    EXPECT_FALSE(lt.tree[1]->is_root());
    EXPECT_EQ(lt.tree[1]->parent_index, 0u);
    EXPECT_EQ(lt.tree[1]->record.signature.target_type, PT::HNucleus);   // non-root target survives
}

int main(int argc, char ** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
