// projects/dataclasses/private/test/InteractionRecord_TEST.cxx
//
// Full test suite for InteractionRecord and related record classes.
// This file is intended to REPLACE the current InteractionRecord_TEST.cxx
// (it includes main()).

#include <array>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/ParticleID.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/utilities/Constants.h"

using namespace siren::dataclasses;

namespace {

template <size_t N>
void ExpectArrayNear(std::array<double, N> const& a,
                     std::array<double, N> const& b,
                     double tol = 1e-12) {
    for (size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(a.at(i), b.at(i), tol) << "index i=" << i;
    }
}

InteractionRecord MakePopulatedInteractionRecord() {
    InteractionRecord r;
    r.signature.primary_type = ParticleType::EMinus;
    r.signature.target_type = ParticleType::PPlus;
    r.signature.secondary_types = {ParticleType::EMinus, ParticleType::EMinus};

    r.primary_id = ParticleID(1, 1);
    r.primary_initial_position = {1.0, 2.0, 3.0};
    r.primary_initial_time = 0.5;
    r.primary_mass = 8.0;
    r.primary_momentum = {10.0, 6.0, 0.0, 0.0};
    r.primary_helicity = 0.5;

    r.target_id = ParticleID(2, 2);
    r.target_mass = 1.0;
    r.target_helicity = -0.5;

    r.interaction_vertex = {4.0, 2.0, 3.0};
    r.interaction_time = 17.0;

    r.secondary_ids = {ParticleID(3, 3), ParticleID(4, 4)};
    r.secondary_masses = {0.5, 0.25};
    r.secondary_momenta = {
        std::array<double, 4>{1.0, 0.0, 1.0, 0.0},
        std::array<double, 4>{2.0, 0.0, 0.0, 2.0},
    };
    r.secondary_helicities = {1.0, -1.0};
    r.secondary_times = {17.0, 21.0};

    r.interaction_parameters = {{"alpha", 1.0}, {"beta", 2.0}};
    return r;
}

}  // namespace

///////////////////////////////////////////////
// InteractionRecord: comparison + ordering
///////////////////////////////////////////////

TEST(InteractionRecord_Comparison, Equality_CoversAllFields) {
    InteractionRecord A;
    InteractionRecord B;
    EXPECT_TRUE(A == B);

    // signature
    A.signature.primary_type = ParticleType::EMinus;
    EXPECT_FALSE(A == B);
    B.signature.primary_type = ParticleType::EMinus;
    EXPECT_TRUE(A == B);

    A.signature.target_type = ParticleType::PPlus;
    EXPECT_FALSE(A == B);
    B.signature.target_type = ParticleType::PPlus;
    EXPECT_TRUE(A == B);

    A.signature.secondary_types.push_back(ParticleType::EMinus);
    EXPECT_FALSE(A == B);
    B.signature.secondary_types.push_back(ParticleType::EMinus);
    EXPECT_TRUE(A == B);

    // primary
    A.primary_id = ParticleID(1, 1);
    EXPECT_FALSE(A == B);
    B.primary_id = ParticleID(1, 1);
    EXPECT_TRUE(A == B);

    A.primary_initial_position = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A == B);
    B.primary_initial_position = {1.0, 2.0, 3.0};
    EXPECT_TRUE(A == B);

    A.primary_initial_time = 1.0;
    EXPECT_FALSE(A == B);
    B.primary_initial_time = 1.0;
    EXPECT_TRUE(A == B);

    A.primary_mass = 1.0;
    EXPECT_FALSE(A == B);
    B.primary_mass = 1.0;
    EXPECT_TRUE(A == B);

    A.primary_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_FALSE(A == B);
    B.primary_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_TRUE(A == B);

    A.primary_helicity = 1.0;
    EXPECT_FALSE(A == B);
    B.primary_helicity = 1.0;
    EXPECT_TRUE(A == B);

    // target id is part of equality
    A.target_id = ParticleID(9, 9);
    EXPECT_FALSE(A == B);
    B.target_id = ParticleID(9, 9);
    EXPECT_TRUE(A == B);

    A.target_mass = 1.0;
    EXPECT_FALSE(A == B);
    B.target_mass = 1.0;
    EXPECT_TRUE(A == B);

    A.target_helicity = 1.0;
    EXPECT_FALSE(A == B);
    B.target_helicity = 1.0;
    EXPECT_TRUE(A == B);

    // vertex
    A.interaction_vertex = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A == B);
    B.interaction_vertex = {1.0, 2.0, 3.0};
    EXPECT_TRUE(A == B);

    A.interaction_time = 1.0;
    EXPECT_FALSE(A == B);
    B.interaction_time = 1.0;
    EXPECT_TRUE(A == B);

    // secondaries: ids + masses + momenta + helicities
    A.secondary_ids.push_back(ParticleID(7, 7));
    EXPECT_FALSE(A == B);
    B.secondary_ids.push_back(ParticleID(7, 7));
    EXPECT_TRUE(A == B);

    A.secondary_masses.push_back(1.0);
    EXPECT_FALSE(A == B);
    B.secondary_masses.push_back(1.0);
    EXPECT_TRUE(A == B);

    A.secondary_momenta.push_back({1.0, 2.0, 3.0, 4.0});
    EXPECT_FALSE(A == B);
    B.secondary_momenta.push_back({1.0, 2.0, 3.0, 4.0});
    EXPECT_TRUE(A == B);

    A.secondary_helicities.push_back(1.0);
    EXPECT_FALSE(A == B);
    B.secondary_helicities.push_back(1.0);
    EXPECT_TRUE(A == B);

    A.secondary_times.push_back(1.0);
    EXPECT_FALSE(A == B);
    B.secondary_times.push_back(1.0);
    EXPECT_TRUE(A == B);

    // interaction parameters
    A.interaction_parameters["test"] = 1.0;
    EXPECT_FALSE(A == B);
    B.interaction_parameters["test"] = 1.0;
    EXPECT_TRUE(A == B);

    A.interaction_parameters["test"] = 2.0;
    EXPECT_FALSE(A == B);
    B.interaction_parameters["test"] = 2.0;
    EXPECT_TRUE(A == B);
}

TEST(InteractionRecord_Comparison, LessThan_IsStrictWeakOrdering_Smoke) {
    InteractionRecord A;
    InteractionRecord B;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.signature.primary_type = ParticleType::EMinus;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.signature.primary_type = ParticleType::EMinus;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    // Target ID participates too
    A.target_id = ParticleID(1, 1);
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.target_id = ParticleID(1, 1);
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);
}

///////////////////////////////////////////////
// InteractionRecord: serialization
///////////////////////////////////////////////

TEST(InteractionRecord_Serialization, Json_RoundTrip) {
    InteractionRecord in = MakePopulatedInteractionRecord();

    std::stringstream ss;
    {
        cereal::JSONOutputArchive oar(ss);
        oar(in);
    }

    InteractionRecord out;
    {
        cereal::JSONInputArchive iar(ss);
        iar(out);
    }

    EXPECT_TRUE(in == out);
}

TEST(InteractionRecord_Serialization, Binary_RoundTrip) {
    InteractionRecord in = MakePopulatedInteractionRecord();

    std::stringstream ss;
    {
        cereal::BinaryOutputArchive oar(ss);
        oar(in);
    }

    InteractionRecord out;
    {
        cereal::BinaryInputArchive iar(ss);
        iar(out);
    }

    EXPECT_TRUE(in == out);
}

TEST(InteractionRecord_Serialization, SaveUnsupportedVersionThrows) {
    InteractionRecord r;

    std::stringstream ss;
    cereal::JSONOutputArchive oar(ss);

    EXPECT_THROW(r.save(oar, 2), std::runtime_error);
}

TEST(InteractionRecord_Serialization, LoadUnsupportedVersionThrows) {
    InteractionRecord r;

    std::stringstream ss;
    ss << "{}";
    cereal::JSONInputArchive iar(ss);

    EXPECT_THROW(r.load(iar, 2), std::runtime_error);
}

TEST(InteractionRecord_Serialization, Version0_LoadsWithDefaultTimes) {
    InteractionRecord in = MakePopulatedInteractionRecord();

    // Write a version 0 archive (pre-time format) directly.
    std::stringstream ss;
    {
        cereal::JSONOutputArchive oar(ss);
        in.save(oar, 0);
    }

    // The version 0 format must not contain any of the new time fields.
    std::string v0_json = ss.str();
    EXPECT_EQ(v0_json.find("PrimaryInitialTime"), std::string::npos);
    EXPECT_EQ(v0_json.find("InteractionTime"), std::string::npos);
    EXPECT_EQ(v0_json.find("SecondaryTimes"), std::string::npos);

    InteractionRecord out;
    out.primary_initial_time = 123.0;
    out.interaction_time = 456.0;
    out.secondary_times = {789.0};
    {
        cereal::JSONInputArchive iar(ss);
        out.load(iar, 0);
    }

    // All non-time fields survive; times revert to their defaults, with
    // secondary_times sized to match the other secondary vectors.
    EXPECT_EQ(out.signature.primary_type, in.signature.primary_type);
    EXPECT_EQ(out.primary_id, in.primary_id);
    ExpectArrayNear<3>(out.primary_initial_position, in.primary_initial_position);
    ExpectArrayNear<3>(out.interaction_vertex, in.interaction_vertex);
    EXPECT_EQ(out.secondary_ids, in.secondary_ids);
    EXPECT_EQ(out.interaction_parameters, in.interaction_parameters);
    EXPECT_DOUBLE_EQ(out.primary_initial_time, 0.0);
    EXPECT_DOUBLE_EQ(out.interaction_time, 0.0);
    ASSERT_EQ(out.secondary_times.size(), in.secondary_momenta.size());
    EXPECT_DOUBLE_EQ(out.secondary_times.at(0), 0.0);
    EXPECT_DOUBLE_EQ(out.secondary_times.at(1), 0.0);
}

TEST(InteractionRecord_Strings, ToStrAndToReprAreNonEmpty) {
    InteractionRecord r = MakePopulatedInteractionRecord();

    std::string s = to_str(r);
    std::string rep = to_repr(r);

    EXPECT_FALSE(s.empty());
    EXPECT_FALSE(rep.empty());

    // Basic sanity that key information is present
    EXPECT_NE(s.find("InteractionRecord"), std::string::npos);
    EXPECT_NE(rep.find("InteractionRecord("), std::string::npos);

    std::stringstream os;
    os << r;
    EXPECT_FALSE(os.str().empty());
}

///////////////////////////////////////////////
// PrimaryDistributionRecord: set/get + derived updates
///////////////////////////////////////////////

TEST(PrimaryDistributionRecord, ConstructionAndInvariants) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    EXPECT_EQ(pr.GetType(), ParticleType::EMinus);
    EXPECT_TRUE(static_cast<bool>(pr.GetID()));

    // Helicity getter never throws; default is 0.
    EXPECT_DOUBLE_EQ(pr.GetHelicity(), 0.0);
}

TEST(PrimaryDistributionRecord, SetParticleRejectsMismatchedIDOrType) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    Particle p;
    p.id = ParticleID(123, 456);
    p.type = ParticleType::EMinus;
    p.mass = 1.0;
    p.momentum = {2.0, 1.0, 0.0, 0.0};
    p.position = {0.0, 0.0, 0.0};
    p.length = 1.0;
    p.helicity = 0.0;

    EXPECT_THROW(pr.SetParticle(p), std::runtime_error);

    // Now match ID but mismatch type
    p.id = pr.GetID();
    p.type = ParticleType::PPlus;
    EXPECT_THROW(pr.SetParticle(p), std::runtime_error);
}

TEST(PrimaryDistributionRecord, SetAndGetFourMomentumRoundTrip) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    pr.SetFourMomentum({10.0, 6.0, 0.0, 0.0});
    auto fm = pr.GetFourMomentum();
    ExpectArrayNear<4>(fm, {10.0, 6.0, 0.0, 0.0});
    ExpectArrayNear<3>(pr.GetThreeMomentum(), {6.0, 0.0, 0.0});
    EXPECT_DOUBLE_EQ(pr.GetEnergy(), 10.0);
}

TEST(PrimaryDistributionRecord, UpdateMassFromEnergyAndMomentum) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    pr.SetEnergy(10.0);
    pr.SetThreeMomentum({6.0, 0.0, 0.0});

    EXPECT_NEAR(pr.GetMass(), 8.0, 1e-12);
}

TEST(PrimaryDistributionRecord, UpdateEnergyFromMassAndMomentum) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    pr.SetMass(8.0);
    pr.SetThreeMomentum({6.0, 0.0, 0.0});

    EXPECT_NEAR(pr.GetEnergy(), 10.0, 1e-12);
}

TEST(PrimaryDistributionRecord, UpdateMomentumFromEnergyMassDirection) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    pr.SetEnergy(10.0);
    pr.SetMass(8.0);
    pr.SetDirection({1.0, 0.0, 0.0});

    ExpectArrayNear<3>(pr.GetThreeMomentum(), {6.0, 0.0, 0.0});
}

TEST(PrimaryDistributionRecord, UpdateDirectionFromMomentum) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    pr.SetThreeMomentum({0.0, 3.0, 4.0});
    auto d = pr.GetDirection();

    // magnitude = 5 -> direction (0, 0.6, 0.8)
    EXPECT_NEAR(d.at(0), 0.0, 1e-12);
    EXPECT_NEAR(d.at(1), 0.6, 1e-12);
    EXPECT_NEAR(d.at(2), 0.8, 1e-12);
}

TEST(PrimaryDistributionRecord, UpdateDirectionAndLengthFromPositions) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    pr.SetInitialPosition({1.0, 2.0, 3.0});
    pr.SetInteractionVertex({2.0, 4.0, 5.0});

    auto d = pr.GetDirection();  // normalized (1,2,2) / 3
    EXPECT_NEAR(d.at(0), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(d.at(1), 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(d.at(2), 2.0 / 3.0, 1e-12);

    EXPECT_NEAR(pr.GetLength(), 3.0, 1e-12);
}

TEST(PrimaryDistributionRecord, ClosestApproachClosure_FromPCAAndDistances) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    // Choose geometry:
    // dir = +x, PCA at (0,2,3), initial at (1,2,3), vertex at (4,2,3)
    pr.SetDirection({1.0, 0.0, 0.0});
    pr.SetPointOfClosestApproach({0.0, 2.0, 3.0});
    pr.SetInitialDistanceFromClosestApproach(1.0);
    pr.SetVertexDistanceFromClosestApproach(4.0);

    ExpectArrayNear<3>(pr.GetInitialPosition(), {1.0, 2.0, 3.0});
    ExpectArrayNear<3>(pr.GetInteractionVertex(), {4.0, 2.0, 3.0});
    EXPECT_NEAR(pr.GetLength(), 3.0, 1e-12);
}

TEST(PrimaryDistributionRecord, ClosestApproach_ComputesDistancesFromPositionAndDirection) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    pr.SetDirection({1.0, 0.0, 0.0});
    pr.SetInitialPosition({1.0, 2.0, 3.0});
    pr.SetInteractionVertex({4.0, 2.0, 3.0});

    // In this implementation, distance-from-closest-approach is the signed projection onto dir to the plane
    // through the origin orthogonal to dir (i.e., PCA w.r.t. the origin line).
    EXPECT_NEAR(pr.GetInitialDistanceFromClosestApproach(), 1.0, 1e-12);
    EXPECT_NEAR(pr.GetVertexDistanceFromClosestApproach(), 4.0, 1e-12);

    // PCA point should be (0,2,3)
    ExpectArrayNear<3>(pr.GetPointOfClosestApproach(), {0.0, 2.0, 3.0});
}

TEST(PrimaryDistributionRecord, ThrowsWhenInsufficientInformation) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    EXPECT_THROW(pr.GetMass(), std::runtime_error);
    EXPECT_THROW(pr.GetEnergy(), std::runtime_error);
    EXPECT_THROW(pr.GetDirection(), std::runtime_error);
    EXPECT_THROW(pr.GetLength(), std::runtime_error);
    EXPECT_THROW(pr.GetInitialPosition(), std::runtime_error);
    EXPECT_THROW(pr.GetInteractionVertex(), std::runtime_error);
    EXPECT_THROW(pr.GetPointOfClosestApproach(), std::runtime_error);
}

TEST(PrimaryDistributionRecord, FinalizeAvailable_FillsWhatItCan) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);
    pr.SetMass(1.0);  // only mass is set

    InteractionRecord out;
    pr.FinalizeAvailable(out);

    EXPECT_EQ(out.signature.primary_type, ParticleType::EMinus);
    EXPECT_EQ(out.primary_id, pr.GetID());

    // Mass should have been written; reminders should stay default (no throw)
    EXPECT_DOUBLE_EQ(out.primary_mass, 1.0);
    ExpectArrayNear<4>(out.primary_momentum, {0.0, 0.0, 0.0, 0.0});
    ExpectArrayNear<3>(out.primary_initial_position, {0.0, 0.0, 0.0});
    ExpectArrayNear<3>(out.interaction_vertex, {0.0, 0.0, 0.0});
}

TEST(PrimaryDistributionRecord, Finalize_RequiresDerivableFields) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    InteractionRecord out;
    // Finalize calls GetInteractionVertex(), GetInitialPosition(), GetMass(), GetFourMomentum()
    EXPECT_THROW(pr.Finalize(out), std::runtime_error);

    // Now provide enough to fully finalize:
    pr.SetMass(8.0);
    pr.SetEnergy(10.0);
    pr.SetDirection({1.0, 0.0, 0.0});
    pr.SetLength(3.0);
    pr.SetInteractionVertex({4.0, 2.0, 3.0});

    EXPECT_NO_THROW(pr.Finalize(out));
    EXPECT_EQ(out.signature.primary_type, ParticleType::EMinus);
    EXPECT_EQ(out.primary_id, pr.GetID());
    ExpectArrayNear<3>(out.interaction_vertex, {4.0, 2.0, 3.0});
    ExpectArrayNear<3>(out.primary_initial_position, {1.0, 2.0, 3.0});
    EXPECT_NEAR(out.primary_mass, 8.0, 1e-12);
    ExpectArrayNear<4>(out.primary_momentum, {10.0, 6.0, 0.0, 0.0});
}

TEST(PrimaryDistributionRecord, GetParticle_UsesDefaultsWhenUnset) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);

    Particle p = pr.GetParticle();
    EXPECT_EQ(p.id, pr.GetID());
    EXPECT_EQ(p.type, pr.GetType());

    // Nothing set => fallbacks
    EXPECT_DOUBLE_EQ(p.mass, 0.0);
    ExpectArrayNear<4>(p.momentum, {0.0, 0.0, 0.0, 0.0});
    ExpectArrayNear<3>(p.position, {0.0, 0.0, 0.0});
    EXPECT_DOUBLE_EQ(p.length, 0.0);
    EXPECT_DOUBLE_EQ(p.helicity, 0.0);
}

///////////////////////////////////////////////
// SecondaryParticleRecord
///////////////////////////////////////////////

TEST(SecondaryParticleRecord, ConstructionReflectsParentRecordAndSignature) {
    InteractionRecord parent;
    parent.signature.primary_type = ParticleType::EMinus;
    parent.signature.target_type = ParticleType::PPlus;
    parent.signature.secondary_types = {ParticleType::EMinus};

    parent.interaction_vertex = {9.0, 8.0, 7.0};

    SecondaryParticleRecord sr(parent, 0);

    EXPECT_EQ(sr.GetType(), ParticleType::EMinus);
    ExpectArrayNear<3>(sr.GetInitialPosition(), {9.0, 8.0, 7.0});
    EXPECT_TRUE(static_cast<bool>(sr.GetID()));
}

TEST(SecondaryParticleRecord, SetParticleRejectsMismatchedIDOrType) {
    InteractionRecord parent;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {0.0, 0.0, 0.0};

    SecondaryParticleRecord sr(parent, 0);

    Particle p;
    p.id = ParticleID(123, 456);
    p.type = ParticleType::EMinus;
    p.mass = 1.0;
    p.momentum = {2.0, 1.0, 0.0, 0.0};
    p.position = {0.0, 0.0, 0.0};
    p.length = 0.0;
    p.helicity = 0.0;

    EXPECT_THROW(sr.SetParticle(p), std::runtime_error);

    p.id = sr.GetID();
    p.type = ParticleType::PPlus;
    EXPECT_THROW(sr.SetParticle(p), std::runtime_error);
}

TEST(SecondaryParticleRecord, UpdateMassFromEnergyAndMomentum) {
    InteractionRecord parent;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {0.0, 0.0, 0.0};

    SecondaryParticleRecord sr(parent, 0);

    sr.SetEnergy(10.0);
    sr.SetThreeMomentum({6.0, 0.0, 0.0});

    EXPECT_NEAR(sr.GetMass(), 8.0, 1e-12);
}

TEST(SecondaryParticleRecord, KineticEnergyConvention_IsEnergyMinusMass) {
    InteractionRecord parent;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {0.0, 0.0, 0.0};

    SecondaryParticleRecord sr(parent, 0);

    sr.SetMass(8.0);
    sr.SetEnergy(10.0);

    EXPECT_NEAR(sr.GetKineticEnergy(), 2.0, 1e-12);
}

TEST(SecondaryParticleRecord, FinalizeWritesIntoPreSizedVectors) {
    InteractionRecord parent;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {0.0, 0.0, 0.0};
    parent.interaction_time = 7.0;

    // SecondaryParticleRecord::Finalize uses .at(), so parent must be sized.
    parent.secondary_ids.resize(1);
    parent.secondary_masses.resize(1);
    parent.secondary_momenta.resize(1);
    parent.secondary_helicities.resize(1);
    parent.secondary_times.resize(1);

    SecondaryParticleRecord sr(parent, 0);
    sr.SetMass(8.0);
    sr.SetFourMomentum({10.0, 6.0, 0.0, 0.0});
    sr.SetHelicity(1.0);

    EXPECT_NO_THROW(sr.Finalize(parent));

    EXPECT_EQ(parent.secondary_ids.at(0), sr.GetID());
    EXPECT_NEAR(parent.secondary_masses.at(0), 8.0, 1e-12);
    ExpectArrayNear<4>(parent.secondary_momenta.at(0), {10.0, 6.0, 0.0, 0.0});
    EXPECT_NEAR(parent.secondary_helicities.at(0), 1.0, 1e-12);
    // Default production time is the vertex time of the record being written.
    EXPECT_NEAR(parent.secondary_times.at(0), 7.0, 1e-12);

    // An explicit per-particle time overrides the default.
    sr.SetTime(11.0);
    EXPECT_NEAR(sr.GetTime(), 11.0, 1e-12);
    sr.Finalize(parent);
    EXPECT_NEAR(parent.secondary_times.at(0), 11.0, 1e-12);
}

///////////////////////////////////////////////
// CrossSectionDistributionRecord
///////////////////////////////////////////////

TEST(CrossSectionDistributionRecord, ConstructionBindsToInteractionRecordAndCreatesSecondaries) {
    InteractionRecord in;
    in.signature.primary_type = ParticleType::EMinus;
    in.signature.target_type = ParticleType::PPlus;
    in.signature.secondary_types = {ParticleType::EMinus, ParticleType::EMinus};

    in.primary_id = ParticleID(1, 1);
    in.primary_initial_position = {1.0, 2.0, 3.0};
    in.primary_mass = 8.0;
    in.primary_momentum = {10.0, 6.0, 0.0, 0.0};
    in.primary_helicity = 0.5;
    in.interaction_vertex = {4.0, 2.0, 3.0};

    // Leave target_id unset to force generation.
    EXPECT_FALSE(static_cast<bool>(in.target_id));

    CrossSectionDistributionRecord csr(in);

    EXPECT_EQ(csr.GetPrimaryID(), in.primary_id);
    EXPECT_EQ(csr.GetPrimaryType(), in.signature.primary_type);
    ExpectArrayNear<3>(csr.GetInteractionVertex(), in.interaction_vertex);

    EXPECT_TRUE(static_cast<bool>(csr.GetTargetID()));
    EXPECT_EQ(csr.GetTargetType(), in.signature.target_type);

    EXPECT_EQ(csr.GetSecondaryParticleRecords().size(), in.signature.secondary_types.size());
    EXPECT_EQ(csr.GetSecondaryParticleRecord(0).GetType(), ParticleType::EMinus);
}

TEST(CrossSectionDistributionRecord, FinalizePopulatesTargetParamsAndSecondaries) {
    InteractionRecord in;
    in.signature.primary_type = ParticleType::EMinus;
    in.signature.target_type = ParticleType::PPlus;
    in.signature.secondary_types = {ParticleType::EMinus};

    in.primary_id = ParticleID(1, 1);
    in.primary_initial_position = {1.0, 2.0, 3.0};
    in.primary_mass = 8.0;
    in.primary_momentum = {10.0, 6.0, 0.0, 0.0};
    in.primary_helicity = 0.5;
    in.interaction_vertex = {4.0, 2.0, 3.0};

    CrossSectionDistributionRecord csr(in);

    csr.SetTargetMass(2.0);
    csr.SetTargetHelicity(-1.0);
    csr.SetInteractionParameter("alpha", 1.0);

    auto& s0 = csr.GetSecondaryParticleRecord(0);
    s0.SetMass(0.5);
    s0.SetFourMomentum({1.0, 0.0, 1.0, 0.0});
    s0.SetHelicity(1.0);

    InteractionRecord out = in;  // start from input; finalize overwrites/sets relevant fields
    csr.Finalize(out);

    EXPECT_TRUE(static_cast<bool>(out.target_id));
    EXPECT_NEAR(out.target_mass, 2.0, 1e-12);
    EXPECT_NEAR(out.target_helicity, -1.0, 1e-12);
    EXPECT_EQ(out.interaction_parameters.at("alpha"), 1.0);

    ASSERT_EQ(out.secondary_ids.size(), 1u);
    ASSERT_EQ(out.secondary_masses.size(), 1u);
    ASSERT_EQ(out.secondary_momenta.size(), 1u);
    ASSERT_EQ(out.secondary_helicities.size(), 1u);

    EXPECT_TRUE(static_cast<bool>(out.secondary_ids.at(0)));
    EXPECT_NEAR(out.secondary_masses.at(0), 0.5, 1e-12);
    ExpectArrayNear<4>(out.secondary_momenta.at(0), {1.0, 0.0, 1.0, 0.0});
    EXPECT_NEAR(out.secondary_helicities.at(0), 1.0, 1e-12);
}

TEST(CrossSectionDistributionRecord_Strings, ToStrAndToReprAreNonEmpty) {
    InteractionRecord in = MakePopulatedInteractionRecord();
    CrossSectionDistributionRecord csr(in);

    std::string s = to_str(csr);
    std::string rep = to_repr(csr);

    EXPECT_FALSE(s.empty());
    EXPECT_FALSE(rep.empty());

    std::stringstream os;
    os << csr;
    EXPECT_FALSE(os.str().empty());
}

///////////////////////////////////////////////
// SecondaryDistributionRecord
///////////////////////////////////////////////

TEST(SecondaryDistributionRecord, GetLengthThrowsUntilSet) {
    InteractionRecord r;
    r.signature.primary_type = ParticleType::EMinus;
    r.primary_momentum = {10.0, 6.0, 8.0, 0.0};
    r.primary_initial_position = {1.0, 2.0, 3.0};

    SecondaryDistributionRecord sdr(r);

    EXPECT_THROW(sdr.GetLength(), std::runtime_error);
    sdr.SetLength(5.0);
    EXPECT_NEAR(sdr.GetLength(), 5.0, 1e-12);
}

TEST(SecondaryDistributionRecord, ConstructionFromRecordEnsuresPrimaryIDIsSet) {
    InteractionRecord r;
    r.signature.primary_type = ParticleType::EMinus;
    r.primary_momentum = {10.0, 6.0, 8.0, 0.0};
    r.primary_initial_position = {1.0, 2.0, 3.0};

    EXPECT_FALSE(static_cast<bool>(r.primary_id));

    SecondaryDistributionRecord sdr(r);

    // ctor mutates r.primary_id if unset
    EXPECT_TRUE(static_cast<bool>(r.primary_id));
    EXPECT_EQ(sdr.id, r.primary_id);
}

TEST(SecondaryDistributionRecord, FinalizeComputesInteractionVertexFromLengthAndDirection) {
    InteractionRecord r;
    r.signature.primary_type = ParticleType::EMinus;

    // energy!=0 triggers direction computation; spatial magnitude = 10 -> dir = (0.6,0.8,0)
    r.primary_momentum = {10.0, 6.0, 8.0, 0.0};
    r.primary_initial_position = {1.0, 2.0, 3.0};
    r.primary_mass = 8.0;
    r.primary_helicity = 0.0;

    SecondaryDistributionRecord sdr(r);
    sdr.SetLength(5.0);

    InteractionRecord out;
    sdr.Finalize(out);

    ExpectArrayNear<3>(out.primary_initial_position, {1.0, 2.0, 3.0});
    // vertex = initial + 5 * (0.6,0.8,0) = (4,6,3)
    ExpectArrayNear<3>(out.interaction_vertex, {4.0, 6.0, 3.0});
}

TEST(SecondaryDistributionRecord, ConstructFromParentSecondary_UsesSecondaryAsPrimary) {
    InteractionRecord parent;
    parent.signature.primary_type = ParticleType::EMinus;
    parent.signature.target_type = ParticleType::PPlus;
    parent.signature.secondary_types = {ParticleType::EMinus};

    parent.interaction_vertex = {10.0, 0.0, 0.0};

    parent.secondary_ids = {ParticleID(42, 7)};
    parent.secondary_masses = {0.5};
    parent.secondary_momenta = {std::array<double, 4>{5.0, 3.0, 4.0, 0.0}};  // spatial mag=5 -> dir=(0.6,0.8,0)
    parent.secondary_helicities = {1.0};

    SecondaryDistributionRecord sdr(parent, 0);
    sdr.SetLength(10.0);

    InteractionRecord out;
    sdr.Finalize(out);

    EXPECT_EQ(out.signature.primary_type, ParticleType::EMinus);
    EXPECT_EQ(out.primary_id, ParticleID(42, 7));
    ExpectArrayNear<3>(out.primary_initial_position, parent.interaction_vertex);
    // vertex = (10,0,0) + 10*(0.6,0.8,0) = (16,8,0)
    ExpectArrayNear<3>(out.interaction_vertex, {16.0, 8.0, 0.0});
}

///////////////////////////////////////////////
// to_str/to_repr/operator<< for record classes
///////////////////////////////////////////////

TEST(RecordStrings, PrimaryAndSecondaryReprsAreNonEmpty) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);
    pr.SetMass(8.0);
    pr.SetFourMomentum({10.0, 6.0, 0.0, 0.0});
    pr.SetInitialPosition({1.0, 2.0, 3.0});
    pr.SetInteractionVertex({4.0, 2.0, 3.0});

    std::string pr_s = to_str(pr);
    std::string pr_r = to_repr(pr);
    EXPECT_FALSE(pr_s.empty());
    EXPECT_FALSE(pr_r.empty());

    std::stringstream pr_os;
    pr_os << pr;
    EXPECT_FALSE(pr_os.str().empty());

    InteractionRecord parent;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {0.0, 0.0, 0.0};

    SecondaryParticleRecord sr(parent, 0);
    sr.SetMass(8.0);
    sr.SetFourMomentum({10.0, 6.0, 0.0, 0.0});

    std::string sr_s = to_str(sr);
    std::string sr_r = to_repr(sr);
    EXPECT_FALSE(sr_s.empty());
    EXPECT_FALSE(sr_r.empty());

    std::stringstream sr_os;
    sr_os << sr;
    EXPECT_FALSE(sr_os.str().empty());
}

///////////////////////////////////////////////
// Time propagation
///////////////////////////////////////////////

TEST(TimePropagation, PrimaryRecordComputesFlightTime) {
    using siren::utilities::Constants::c;

    PrimaryDistributionRecord pr(ParticleType::EMinus);
    pr.SetMass(8.0);
    pr.SetEnergy(10.0);  // |p| = 6 -> beta = 0.6
    pr.SetDirection({1.0, 0.0, 0.0});
    pr.SetInitialPosition({1.0, 2.0, 3.0});
    pr.SetInteractionVertex({4.0, 2.0, 3.0});  // distance = 3 (internal meters)

    double const flight = 3.0 / (0.6 * c);

    // Automatic rule: initial time defaults to 0.
    EXPECT_DOUBLE_EQ(pr.GetInitialTime(), 0.0);
    EXPECT_NEAR(pr.GetInteractionTime(), flight, 1e-9);

    InteractionRecord out;
    pr.Finalize(out);
    EXPECT_DOUBLE_EQ(out.primary_initial_time, 0.0);
    EXPECT_NEAR(out.interaction_time, flight, 1e-9);
}

TEST(TimePropagation, PrimaryRecordInitialTimeShiftsInteractionTime) {
    using siren::utilities::Constants::c;

    PrimaryDistributionRecord pr(ParticleType::EMinus);
    pr.SetMass(8.0);
    pr.SetEnergy(10.0);
    pr.SetDirection({1.0, 0.0, 0.0});
    pr.SetInitialPosition({1.0, 2.0, 3.0});
    pr.SetInteractionVertex({4.0, 2.0, 3.0});
    pr.SetInitialTime(2.5);

    double const flight = 3.0 / (0.6 * c);

    InteractionRecord out;
    pr.Finalize(out);
    EXPECT_DOUBLE_EQ(out.primary_initial_time, 2.5);
    EXPECT_NEAR(out.interaction_time, 2.5 + flight, 1e-9);
}

TEST(TimePropagation, PrimaryRecordInteractionTimeOverride) {
    PrimaryDistributionRecord pr(ParticleType::EMinus);
    pr.SetMass(8.0);
    pr.SetEnergy(10.0);
    pr.SetDirection({1.0, 0.0, 0.0});
    pr.SetInitialPosition({1.0, 2.0, 3.0});
    pr.SetInteractionVertex({4.0, 2.0, 3.0});
    pr.SetInteractionTime(42.0);

    InteractionRecord out;
    pr.Finalize(out);
    EXPECT_DOUBLE_EQ(out.interaction_time, 42.0);
}

TEST(TimePropagation, PrimaryRecordAtRestHasZeroFlightTime) {
    PrimaryDistributionRecord pr(ParticleType::MuMinus);
    pr.SetMass(0.105);
    pr.SetThreeMomentum({0.0, 0.0, 0.0});
    pr.SetInitialPosition({1.0, 2.0, 3.0});
    pr.SetInteractionVertex({1.0, 2.0, 3.0});
    pr.SetInitialTime(4.0);

    InteractionRecord out;
    pr.Finalize(out);
    EXPECT_DOUBLE_EQ(out.interaction_time, 4.0);
}

TEST(TimePropagation, CrossSectionRecordFillsSecondaryTimesFromVertexTime) {
    InteractionRecord in = MakePopulatedInteractionRecord();
    in.interaction_time = 17.0;

    CrossSectionDistributionRecord csr(in);
    EXPECT_DOUBLE_EQ(csr.GetInteractionTime(), 17.0);

    auto & s0 = csr.GetSecondaryParticleRecord(0);
    s0.SetMass(0.5);
    s0.SetFourMomentum({1.0, 0.0, 1.0, 0.0});
    s0.SetHelicity(1.0);
    auto & s1 = csr.GetSecondaryParticleRecord(1);
    s1.SetMass(0.25);
    s1.SetFourMomentum({2.0, 0.0, 0.0, 2.0});
    s1.SetHelicity(-1.0);

    InteractionRecord out = in;
    csr.Finalize(out);
    ASSERT_EQ(out.secondary_times.size(), 2u);
    EXPECT_DOUBLE_EQ(out.secondary_times.at(0), 17.0);
    EXPECT_DOUBLE_EQ(out.secondary_times.at(1), 17.0);
}

TEST(TimePropagation, CrossSectionRecordVertexTimeOverridePropagatesToSecondaries) {
    InteractionRecord in = MakePopulatedInteractionRecord();
    in.interaction_time = 17.0;

    CrossSectionDistributionRecord csr(in);
    csr.SetInteractionTime(20.0);
    EXPECT_DOUBLE_EQ(csr.GetInteractionTime(), 20.0);
    // The secondaries' default production times follow the override.
    EXPECT_DOUBLE_EQ(csr.GetSecondaryParticleRecord(0).GetTime(), 20.0);

    auto & s0 = csr.GetSecondaryParticleRecord(0);
    s0.SetMass(0.5);
    s0.SetFourMomentum({1.0, 0.0, 1.0, 0.0});
    s0.SetHelicity(1.0);
    auto & s1 = csr.GetSecondaryParticleRecord(1);
    s1.SetMass(0.25);
    s1.SetFourMomentum({2.0, 0.0, 0.0, 2.0});
    s1.SetHelicity(-1.0);
    // A per-particle time takes precedence over the vertex time.
    s1.SetTime(25.0);

    InteractionRecord out = in;
    csr.Finalize(out);
    EXPECT_DOUBLE_EQ(out.interaction_time, 20.0);
    ASSERT_EQ(out.secondary_times.size(), 2u);
    EXPECT_DOUBLE_EQ(out.secondary_times.at(0), 20.0);
    EXPECT_DOUBLE_EQ(out.secondary_times.at(1), 25.0);
}

TEST(TimePropagation, SecondaryRecordThreadsProductionTimeFromParent) {
    using siren::utilities::Constants::c;

    InteractionRecord parent;
    parent.signature.primary_type = ParticleType::EMinus;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {10.0, 0.0, 0.0};
    parent.interaction_time = 17.0;
    parent.secondary_ids = {ParticleID(42, 7)};
    parent.secondary_masses = {0.0};
    parent.secondary_momenta = {std::array<double, 4>{5.0, 3.0, 4.0, 0.0}};  // beta = 1
    parent.secondary_helicities = {1.0};
    parent.secondary_times = {21.0};

    SecondaryDistributionRecord sdr(parent, 0);
    EXPECT_DOUBLE_EQ(sdr.initial_time, 21.0);

    // Interaction time requires the length.
    EXPECT_THROW(sdr.GetInteractionTime(), std::runtime_error);
    sdr.SetLength(10.0);
    EXPECT_NEAR(sdr.GetInteractionTime(), 21.0 + 10.0 / c, 1e-9);

    InteractionRecord out;
    sdr.Finalize(out);
    EXPECT_DOUBLE_EQ(out.primary_initial_time, 21.0);
    EXPECT_NEAR(out.interaction_time, 21.0 + 10.0 / c, 1e-9);
}

TEST(TimePropagation, SecondaryRecordFallsBackToParentVertexTime) {
    InteractionRecord parent;
    parent.signature.primary_type = ParticleType::EMinus;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {10.0, 0.0, 0.0};
    parent.interaction_time = 17.0;
    parent.secondary_ids = {ParticleID(42, 7)};
    parent.secondary_masses = {0.0};
    parent.secondary_momenta = {std::array<double, 4>{5.0, 3.0, 4.0, 0.0}};
    parent.secondary_helicities = {1.0};
    // No secondary_times, as for records loaded from version 0 archives.

    SecondaryDistributionRecord sdr(parent, 0);
    EXPECT_DOUBLE_EQ(sdr.initial_time, 17.0);
}

TEST(TimePropagation, SecondaryRecordInteractionTimeOverride) {
    InteractionRecord parent;
    parent.signature.primary_type = ParticleType::EMinus;
    parent.signature.secondary_types = {ParticleType::EMinus};
    parent.interaction_vertex = {10.0, 0.0, 0.0};
    parent.interaction_time = 17.0;
    parent.secondary_ids = {ParticleID(42, 7)};
    parent.secondary_masses = {0.0};
    parent.secondary_momenta = {std::array<double, 4>{5.0, 3.0, 4.0, 0.0}};
    parent.secondary_helicities = {1.0};
    parent.secondary_times = {21.0};

    SecondaryDistributionRecord sdr(parent, 0);
    sdr.SetLength(10.0);
    sdr.SetInteractionTime(99.0);

    InteractionRecord out;
    sdr.Finalize(out);
    EXPECT_DOUBLE_EQ(out.interaction_time, 99.0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

