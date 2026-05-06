#include <cmath>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>

#include <gtest/gtest.h>

#include "SIREN/interactions/HNLDISFromSpline.h"
#include "SIREN/interactions/HNLDipoleDISFromSpline.h"
#include "SIREN/interactions/DISFromSpline.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/ParticleMasses.h"
#include "SIREN/utilities/Constants.h"

using namespace siren::interactions;
using namespace siren::dataclasses;
using namespace siren::utilities;

// ---------------------------------------------------------------------------
// Constants and ParticleMasses
// ---------------------------------------------------------------------------

TEST(Constants, MassesPositive) {
    EXPECT_GT(Constants::protonMass, 0);
    EXPECT_GT(Constants::neutronMass, 0);
    EXPECT_GT(Constants::electronMass, 0);
    EXPECT_GT(Constants::muonMass, 0);
    EXPECT_GT(Constants::tauMass, 0);
    EXPECT_GT(Constants::wMass, 0);
    EXPECT_GT(Constants::zMass, 0);
    EXPECT_GT(Constants::Pi0Mass, 0);
}

TEST(Constants, MassOrdering) {
    EXPECT_LT(Constants::electronMass, Constants::muonMass);
    EXPECT_LT(Constants::muonMass, Constants::tauMass);
    EXPECT_LT(Constants::protonMass, Constants::neutronMass);
    EXPECT_LT(Constants::Pi0Mass, Constants::protonMass);
}

TEST(Constants, CKMUnitarity) {
    // First row of CKM should satisfy |Vud|^2 + |Vus|^2 + |Vub|^2 ~ 1
    double row1 = Constants::Vud * Constants::Vud
                + Constants::Vus * Constants::Vus
                + Constants::Vub * Constants::Vub;
    EXPECT_NEAR(row1, 1.0, 0.01);

    // Second row
    double row2 = Constants::Vcd * Constants::Vcd
                + Constants::Vcs * Constants::Vcs
                + Constants::Vcb * Constants::Vcb;
    EXPECT_NEAR(row2, 1.0, 0.01);
}

TEST(ParticleMasses, NeutrinosMassless) {
    EXPECT_EQ(GetParticleMass(ParticleType::NuE), 0.0);
    EXPECT_EQ(GetParticleMass(ParticleType::NuMu), 0.0);
    EXPECT_EQ(GetParticleMass(ParticleType::NuTau), 0.0);
    EXPECT_EQ(GetParticleMass(ParticleType::NuEBar), 0.0);
}

TEST(ParticleMasses, ChargedLeptons) {
    EXPECT_DOUBLE_EQ(GetParticleMass(ParticleType::EMinus), Constants::electronMass);
    EXPECT_DOUBLE_EQ(GetParticleMass(ParticleType::MuMinus), Constants::muonMass);
    EXPECT_DOUBLE_EQ(GetParticleMass(ParticleType::TauMinus), Constants::tauMass);
}

TEST(ParticleMasses, AntiparticleSymmetry) {
    EXPECT_DOUBLE_EQ(GetParticleMass(ParticleType::MuMinus),
                     GetParticleMass(ParticleType::MuPlus));
    EXPECT_DOUBLE_EQ(GetParticleMass(ParticleType::PiPlus),
                     GetParticleMass(ParticleType::PiMinus));
    EXPECT_DOUBLE_EQ(GetParticleMass(ParticleType::KPlus),
                     GetParticleMass(ParticleType::KMinus));
    EXPECT_DOUBLE_EQ(GetParticleMass(ParticleType::WPlus),
                     GetParticleMass(ParticleType::WMinus));
}

TEST(ParticleMasses, UnknownReturnsZero) {
    EXPECT_EQ(GetParticleMass(ParticleType::unknown), 0.0);
}

// ---------------------------------------------------------------------------
// HNLDipoleDISFromSpline construction and validation
// ---------------------------------------------------------------------------

TEST(HNLDipoleDISFromSpline, DefaultConstructor) {
    HNLDipoleDISFromSpline xs;
    EXPECT_EQ(xs.GetPossiblePrimaries().size(), 0u);
}

TEST(HNLDipoleDISFromSpline, SetUnitsWrongCouplingSize) {
    HNLDipoleDISFromSpline xs;
    // Default-constructed has empty coupling vector; SetUnits rejects size != 3
    EXPECT_THROW(xs.SetUnits("invGeV"), std::runtime_error);
    EXPECT_THROW(xs.SetUnits("cm"), std::runtime_error);
    EXPECT_THROW(xs.SetUnits("m"), std::runtime_error);
}

TEST(HNLDipoleDISFromSpline, SetUnitsInvalidThrows) {
    // Default-constructed object has empty dipole_coupling_
    // SetUnits should throw because coupling size != 3
    HNLDipoleDISFromSpline xs;
    EXPECT_THROW(xs.SetUnits("cm"), std::runtime_error);
}

// ---------------------------------------------------------------------------
// HNLDISFromSpline construction and validation
// ---------------------------------------------------------------------------

TEST(HNLDISFromSpline, DefaultConstructor) {
    HNLDISFromSpline xs;
    EXPECT_EQ(xs.GetPossiblePrimaries().size(), 0u);
}

TEST(HNLDISFromSpline, SetUnitsInvalidThrows) {
    HNLDISFromSpline xs;
    EXPECT_THROW(xs.SetUnits("m"), std::runtime_error);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
