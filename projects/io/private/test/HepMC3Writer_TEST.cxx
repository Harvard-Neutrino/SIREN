#include <cstdio>
#include <string>
#include <vector>
#include <memory>
#include <gtest/gtest.h>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/Attribute.h>
#include <HepMC3/ReaderAscii.h>

#include "SIREN/io/HepMC3Writer.h"
#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleID.h"
#include "SIREN/dataclasses/ParticleType.h"

using namespace siren::dataclasses;

namespace {

// A two-level tree: a primary interaction (NuMu on EPlus -> NuMu + EMinus) whose
// first secondary (NuMu, id 1:1) is itself the primary of a daughter vertex.
InteractionTree BuildTree() {
    InteractionRecord root;
    root.signature.primary_type = ParticleType::NuMu;
    root.signature.target_type = ParticleType::EPlus;
    root.signature.secondary_types = {ParticleType::NuMu, ParticleType::EMinus};
    root.primary_momentum = {10.0, 0.0, 0.0, 10.0}; // [E, px, py, pz]
    root.primary_mass = 0.0;
    root.primary_helicity = -1.0;
    root.primary_initial_position = {0.0, 0.0, 0.0};
    root.primary_initial_time = 0.0;
    root.target_mass = 0.000511;
    root.target_helicity = 0.5;
    root.interaction_vertex = {1.0, 2.0, 3.0}; // internal meters
    root.interaction_time = 5.0;
    root.secondary_ids = {ParticleID(1, 1), ParticleID(1, 2)};
    root.secondary_momenta = {{6.0, 0.0, 0.0, 6.0}, {4.0, 0.0, 0.0, 4.0}};
    root.secondary_masses = {0.0, 0.000511};
    root.secondary_helicities = {-1.0, 1.0};
    root.secondary_times = {5.0, 5.0};
    root.interaction_parameters = {{"foo", 1.5}, {"bar", -2.25}};

    InteractionRecord dau;
    dau.signature.primary_type = ParticleType::NuMu;
    dau.signature.target_type = ParticleType::EPlus;
    dau.signature.secondary_types = {ParticleType::NuMu};
    dau.primary_id = ParticleID(1, 1); // == root.secondary_ids[0] -> shared particle
    dau.primary_momentum = {6.0, 0.0, 0.0, 6.0};
    dau.primary_mass = 0.0;
    dau.primary_helicity = -1.0;
    dau.target_mass = 0.000511;
    dau.interaction_vertex = {2.0, 3.0, 4.0};
    dau.interaction_time = 6.0;
    dau.secondary_ids = {ParticleID(1, 3)};
    dau.secondary_momenta = {{6.0, 0.0, 0.0, 6.0}};
    dau.secondary_masses = {0.0};
    dau.secondary_helicities = {-1.0};
    dau.secondary_times = {6.0};

    InteractionTree tree;
    std::shared_ptr<InteractionTreeDatum> r = tree.add_entry(root);
    tree.add_entry(dau, r);
    tree.header.event_number = 42;
    tree.header.weights = {2.0};
    return tree;
}

} // namespace

TEST(HepMC3Writer, RoundTripStructure) {
    InteractionTree tree = BuildTree();

    std::string const path = std::string(::testing::TempDir()) + "siren_hepmc3_writer_test.hepmc3";

    std::vector<std::shared_ptr<InteractionTree>> trees = {
        std::make_shared<InteractionTree>(tree)};
    siren::io::SaveInteractionTreesAsHepMC3(trees, path);

    HepMC3::ReaderAscii reader(path);
    ASSERT_FALSE(reader.failed());

    int nevents = 0;
    HepMC3::GenEvent last;
    while(true) {
        HepMC3::GenEvent evt;
        reader.read_event(evt);
        if(reader.failed()) break;
        last = evt;
        ++nevents;
    }
    reader.close();
    std::remove(path.c_str());

    EXPECT_EQ(nevents, 1);
    EXPECT_EQ(last.event_number(), 42);

    ASSERT_FALSE(last.weights().empty());
    EXPECT_DOUBLE_EQ(last.weights()[0], 2.0);

    // One vertex per interaction record.
    EXPECT_EQ(last.vertices().size(), 2u);
    // beam + root target + shared NuMu + EMinus + daughter target + daughter NuMu.
    EXPECT_GE(last.particles().size(), 5u);

    // Exactly one status-4 beam particle carrying the primary's pid and energy.
    int beams = 0;
    for(auto const & p : last.particles()) {
        if(p->status() == 4) {
            ++beams;
            EXPECT_EQ(p->pid(), static_cast<int>(ParticleType::NuMu));
            EXPECT_NEAR(p->momentum().e(), 10.0, 1e-9);
        }
    }
    EXPECT_EQ(beams, 1);

    // Primary vertex position converted internal-meters -> CM ({1,2,3} -> {100,200,300}).
    bool found_primary_vertex = false;
    for(auto const & v : last.vertices()) {
        if(v->status() == 1) {
            found_primary_vertex = true;
            EXPECT_NEAR(v->position().x(), 100.0, 1e-6);
            EXPECT_NEAR(v->position().y(), 200.0, 1e-6);
            EXPECT_NEAR(v->position().z(), 300.0, 1e-6);
        }
    }
    EXPECT_TRUE(found_primary_vertex);

    // interaction_parameters survive as per-key DoubleAttributes on a vertex.
    bool found_foo = false;
    for(auto const & v : last.vertices()) {
        auto a = v->attribute<HepMC3::DoubleAttribute>("siren.param.foo");
        if(a) {
            found_foo = true;
            EXPECT_DOUBLE_EQ(a->value(), 1.5);
        }
    }
    EXPECT_TRUE(found_foo);
}

int main(int argc, char ** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
