#include <gtest/gtest.h>
#include <limits>
#include "SIREN/math/FourMomentum.h"
#include "SIREN/dataclasses/RecordKinematics.h"

using siren::math::FourMomentum;
using siren::math::Vector3D;

TEST(FourMomentum, RawComponentsAndSignedTransfers) {
    auto raw = FourMomentum::FromComponents({-2., 3., 4., 0.});
    EXPECT_EQ(raw.Components(), (std::array<double, 4>{-2., 3., 4., 0.}));
    EXPECT_DOUBLE_EQ(raw.MassSquared(), -21.);
    EXPECT_FALSE(raw.HasKnownMass());
    EXPECT_THROW(raw.Mass(), std::domain_error);
    auto frame = FourMomentum::FromMass(2., {0., 0., 5.});
    EXPECT_THROW(raw.BoostToLab(frame), std::domain_error);
    EXPECT_THROW(raw.BoostToRest(frame), std::domain_error);
    auto a = FourMomentum::FromMass(1., {0., 0., 3.});
    auto b = FourMomentum::FromMass(1., {0., 0., -3.});
    auto transfer = a - b;
    EXPECT_DOUBLE_EQ(transfer.MassSquared(), -36.);
    EXPECT_DOUBLE_EQ((a + b).Energy(), a.Energy() + b.Energy());
    EXPECT_NEAR(FourMomentum::InvariantMass(a, b), 2.*std::sqrt(10.), 1e-14);
    EXPECT_THROW(FourMomentum::InvariantMass(raw, b), std::domain_error);
}

TEST(FourMomentum, KnownMassSurvivesCancellationAndBoost) {
    auto frame = FourMomentum::FromMass(2., {0., 0., 2e9});
    auto rest = FourMomentum::FromMass(0.001, {0., 0., 0.});
    auto lab = rest.BoostToLab(frame);
    EXPECT_TRUE(lab.HasKnownMass());
    EXPECT_DOUBLE_EQ(lab.Mass(), 0.001);
    EXPECT_NEAR(lab.Components()[3], 1e6, 1e-9);
    EXPECT_NEAR(lab.Energy(), frame.Energy()*0.0005, 1e-9);
    EXPECT_DOUBLE_EQ(lab.MassSquared(), 0.);
    EXPECT_DOUBLE_EQ(FourMomentum::FromComponents(lab.Components()).Mass(), 0.);
}

TEST(FourMomentum, RestMasslessAndBothEnergySigns) {
    auto frame = FourMomentum::FromMass(2., {0., 0., 0.});
    for (bool negative : {false, true}) {
        auto photon = FourMomentum::FromMass(0., {3., 4., 0.}, negative);
        EXPECT_EQ(photon.BoostToLab(frame).Components(), photon.Components());
        EXPECT_EQ(photon.BoostToRest(frame).Components(), photon.Components());
        EXPECT_DOUBLE_EQ(photon.Mass(), 0.);
        auto moving = FourMomentum::FromMass(2., {0., 0., 3.});
        auto lab = photon.BoostToLab(moving);
        EXPECT_NEAR(lab.Energy(), std::sqrt(13.)/2.*photon.Energy(), 2e-14);
        EXPECT_NEAR(lab.Components()[3], 1.5*photon.Energy(), 2e-14);
        auto back = lab.BoostToRest(moving);
        for (size_t i=0; i<4; ++i)
            EXPECT_NEAR(back.Components()[i], photon.Components()[i], 2e-14);
    }
}

TEST(FourMomentum, BoostRejectsInvalidFramesAndInputs) {
    auto p = FourMomentum::FromMass(1., {0., 0., 0.});
    for (auto frame : {FourMomentum::FromMass(0., {0., 0., 1.}),
                       FourMomentum::FromMass(1., {0., 0., 1.}, true),
                       FourMomentum::FromComponents({1., 2., 0., 0.})})
        EXPECT_THROW(p.BoostToLab(frame), std::domain_error);
    EXPECT_THROW(FourMomentum::FromMass(-1., {0., 0., 0.}), std::domain_error);
    EXPECT_THROW(FourMomentum::FromMass(1., {INFINITY, 0., 0.}), std::domain_error);
    EXPECT_THROW(FourMomentum::FromMass(INFINITY, {0., 0., 0.}), std::domain_error);
    EXPECT_THROW(FourMomentum::FromComponents({NAN, 0., 0., 0.}).Mass(), std::domain_error);
    EXPECT_THROW(FourMomentum::FromMass(1e300, {0., 0., 0.}), std::domain_error);
}

TEST(FourMomentum, RecordAccessReturnsValuesAndRequiresExplicitWrite) {
    using namespace siren::dataclasses;
    InteractionRecord record;
    record.primary_mass = 9.;
    record.primary_momentum = {1., 2., 3., 4.};
    record.secondary_momenta = {{-2., 0., 0., 0.}};
    auto value = PrimaryMomentum(record);
    value = value + value;
    EXPECT_DOUBLE_EQ(record.primary_momentum[0], 1.);
    SetPrimaryMomentum(record, value);
    EXPECT_DOUBLE_EQ(record.primary_momentum[0], 2.);
    EXPECT_DOUBLE_EQ(record.primary_mass, 9.);
    auto secondary = SecondaryMomentum(record, 0);
    SetSecondaryMomentum(record, 0, secondary + secondary);
    EXPECT_DOUBLE_EQ(record.secondary_momenta[0][0], -4.);
    EXPECT_THROW(SecondaryMomentum(record, 1), std::out_of_range);
    EXPECT_THROW(SetSecondaryMomentum(record, 1, value), std::out_of_range);
}
