#include <cstdio>
#include <cstdint>
#include <array>
#include <map>
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
#include "SIREN/utilities/Constants.h"

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

// The vertex 4-position must be Minkowski self-consistent: the time slot carries
// c*t re-expressed in the SAME CM length unit used for the spatial slots (divisor
// Constants::cm). This full-precision GenVertex position() t-slot is the primary
// subject and is asserted against a number computed directly from Constants -- if
// the writer ever changes the divisor on one slot but not the others, the 4-vector
// stops being Minkowski and this test fails.
//
// The per-event lab_pos time component (4th entry, in seconds) is a KNOWN LOSSY
// slot: it is emitted through HepMC3::VectorDoubleAttribute, whose to_string uses
// std::to_string(double) == "%f" with only 6 fractional digits, so any lab time
// below ~5e-7 s (i.e. every realistic sub-microsecond neutrino lab time) is
// truncated to 0.000000 on write. This test pins that truncated behavior so the
// finding is visible; the full-precision time is available from the vertex t-slot
// above (t_seconds = position().t() * cm / c / second).
TEST(HepMC3Writer, VertexTimeSlotScaleConsistency) {
    namespace C = siren::utilities::Constants;

    // A single-vertex tree with a known vertex position and interaction time.
    double const t = 5.0;                        // internal time units
    std::array<double, 3> const x = {1.0, 2.0, 3.0}; // internal meters
    InteractionRecord r;
    r.signature.primary_type = ParticleType::NuMu;
    r.signature.target_type = ParticleType::EPlus;
    r.signature.secondary_types = {ParticleType::NuMu};
    r.primary_momentum = {10.0, 0.0, 0.0, 10.0};
    r.primary_mass = 0.0;
    r.primary_helicity = -1.0;
    r.primary_initial_position = {0.0, 0.0, 0.0};
    r.primary_initial_time = 0.0;
    r.target_mass = 0.000511;
    r.target_helicity = 0.5;
    r.interaction_vertex = x;
    r.interaction_time = t;
    r.secondary_ids = {ParticleID(1, 1)};
    r.secondary_momenta = {{10.0, 0.0, 0.0, 10.0}};
    r.secondary_masses = {0.0};
    r.secondary_helicities = {-1.0};
    r.secondary_times = {t};

    InteractionTree tree;
    tree.add_entry(r);
    tree.header.event_number = 7;
    tree.header.weights = {1.0};

    std::string const path =
        std::string(::testing::TempDir()) + "siren_hepmc3_writer_tslot.hepmc3";
    std::vector<std::shared_ptr<InteractionTree>> trees = {
        std::make_shared<InteractionTree>(tree)};
    siren::io::SaveInteractionTreesAsHepMC3(trees, path);

    HepMC3::ReaderAscii reader(path);
    ASSERT_FALSE(reader.failed());
    HepMC3::GenEvent evt;
    reader.read_event(evt);
    ASSERT_FALSE(reader.failed());
    reader.close();

    ASSERT_EQ(evt.vertices().size(), 1u);
    HepMC3::FourVector const pos = evt.vertices().front()->position();

    // Spatial slots: internal meters divided by cm.
    double const spatial_scale = 1.0 / C::cm;
    EXPECT_NEAR(pos.x(), x[0] * spatial_scale, 1e-9);
    EXPECT_NEAR(pos.y(), x[1] * spatial_scale, 1e-9);
    EXPECT_NEAR(pos.z(), x[2] * spatial_scale, 1e-9);

    // Time slot: c*t divided by the SAME cm divisor (Minkowski self-consistency).
    double const expected_t_slot = (C::c * t) / C::cm;
    EXPECT_NEAR(pos.t(), expected_t_slot, 1e-9);
    // Cross-check the numeric anchor: c*5/cm = 0.299792458*5/1e-2.
    EXPECT_NEAR(expected_t_slot, 149.89622899999998, 1e-9);
    // Dividing the time slot by (c/cm) recovers the internal time exactly, which is
    // only true when the slot used the same length divisor as the spatial slots.
    EXPECT_NEAR(pos.t() / (C::c / C::cm), t, 1e-12);

    // The per-event lab_pos vector: spatial slots are the CM vertex position (same
    // cm divisor as above) and survive at full precision because they are O(100).
    auto lab = evt.attribute<HepMC3::VectorDoubleAttribute>("lab_pos");
    ASSERT_TRUE(static_cast<bool>(lab));
    std::vector<double> const lab_pos = lab->value();
    ASSERT_EQ(lab_pos.size(), 4u);
    EXPECT_NEAR(lab_pos[0], x[0] * spatial_scale, 1e-9);
    EXPECT_NEAR(lab_pos[1], x[1] * spatial_scale, 1e-9);
    EXPECT_NEAR(lab_pos[2], x[2] * spatial_scale, 1e-9);

    // The intended lab time is t / Constants::second == 5e-9 s. But the
    // VectorDoubleAttribute %f formatting (6 fractional digits) truncates it to
    // exactly 0 on write. Pin the truncated value: this is the documented finding,
    // not a passing round trip. The recoverable, lossless lab time lives in the
    // vertex t-slot, cross-checked here to be exact.
    double const intended_lab_time = t / C::second; // 5e-9 s
    EXPECT_GT(intended_lab_time, 0.0);
    EXPECT_LT(intended_lab_time, 5e-7); // below the %f truncation floor
    EXPECT_EQ(lab_pos[3], 0.0);         // lossy: the intended 5e-9 s is gone
    // The same physical time, taken from the full-precision Minkowski vertex slot,
    // is preserved to machine precision (this is the value consumers should use).
    double const lab_time_from_vertex = pos.t() * C::cm / C::c / C::second;
    EXPECT_NEAR(lab_time_from_vertex, intended_lab_time, intended_lab_time * 1e-12);

    std::remove(path.c_str());
}

namespace {

// A single-vertex tree with the given primary type (distinct primaries -> distinct
// root signatures -> distinct process ids).
InteractionTree BuildSimpleTree(ParticleType primary, std::uint64_t event_number) {
    InteractionRecord r;
    r.signature.primary_type = primary;
    r.signature.target_type = ParticleType::EPlus;
    r.signature.secondary_types = {ParticleType::NuMu};
    r.primary_momentum = {10.0, 0.0, 0.0, 10.0};
    r.primary_mass = 0.0;
    r.primary_helicity = -1.0;
    r.primary_initial_position = {0.0, 0.0, 0.0};
    r.primary_initial_time = 0.0;
    r.target_mass = 0.000511;
    r.target_helicity = 0.5;
    r.interaction_vertex = {1.0, 2.0, 3.0};
    r.interaction_time = 5.0;
    r.secondary_ids = {ParticleID(1, 1)};
    r.secondary_momenta = {{10.0, 0.0, 0.0, 10.0}};
    r.secondary_masses = {0.0};
    r.secondary_helicities = {-1.0};
    r.secondary_times = {5.0};

    InteractionTree tree;
    tree.add_entry(r);
    tree.header.event_number = event_number;
    tree.header.weights = {1.0};
    return tree;
}

// Map each event's beam pid to its signal_process_id by reading a written file.
std::map<int, int> BeamPidToProcessId(std::string const & path) {
    std::map<int, int> result;
    HepMC3::ReaderAscii reader(path);
    EXPECT_FALSE(reader.failed());
    while(true) {
        HepMC3::GenEvent evt;
        reader.read_event(evt);
        if(reader.failed()) break;
        auto pid_attr = evt.attribute<HepMC3::IntAttribute>("signal_process_id");
        int beam_pid = 0;
        for(auto const & p : evt.particles())
            if(p->status() == 4) { beam_pid = p->pid(); break; }
        if(pid_attr) result[beam_pid] = pid_attr->value();
    }
    reader.close();
    return result;
}

} // namespace

// Process ids are assigned by a deterministic signature sort, not tree-encounter
// order, so the same physics gets the same id no matter what order the trees are
// written in. Writing {A,B} and {B,A} must yield the same beam-pid -> id map.
TEST(HepMC3Writer, DeterministicProcessIds) {
    InteractionTree a = BuildSimpleTree(ParticleType::NuMu, 1);
    InteractionTree b = BuildSimpleTree(ParticleType::NuEBar, 2);

    std::string const p1 = std::string(::testing::TempDir()) + "siren_hepmc3_detid_ab.hepmc3";
    std::string const p2 = std::string(::testing::TempDir()) + "siren_hepmc3_detid_ba.hepmc3";

    std::vector<std::shared_ptr<InteractionTree>> ab = {
        std::make_shared<InteractionTree>(a), std::make_shared<InteractionTree>(b)};
    std::vector<std::shared_ptr<InteractionTree>> ba = {
        std::make_shared<InteractionTree>(b), std::make_shared<InteractionTree>(a)};
    siren::io::SaveInteractionTreesAsHepMC3(ab, p1);
    siren::io::SaveInteractionTreesAsHepMC3(ba, p2);

    std::map<int, int> const m1 = BeamPidToProcessId(p1);
    std::map<int, int> const m2 = BeamPidToProcessId(p2);
    std::remove(p1.c_str());
    std::remove(p2.c_str());

    ASSERT_EQ(m1.size(), 2u);
    EXPECT_EQ(m1, m2); // identical id assignment regardless of write order

    // Ids come from the "Other" band and are numbered in signature order: NuEBar
    // (pdg -12) sorts before NuMu (pdg 14), so NuEBar -> 700 and NuMu -> 701.
    EXPECT_EQ(m1.at(static_cast<int>(ParticleType::NuEBar)), 700);
    EXPECT_EQ(m1.at(static_cast<int>(ParticleType::NuMu)), 701);
}

int main(int argc, char ** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
