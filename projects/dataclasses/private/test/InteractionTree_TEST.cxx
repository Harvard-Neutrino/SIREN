
#include <math.h>
#include <cmath>
#include <random>
#include <iostream>
#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/dataclasses/InteractionTree.h"

using namespace siren::dataclasses;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

int InteractionTreeDatum::depth() const {
    int depth = 0;
    if(parent==NULL) return depth;
    std::shared_ptr<InteractionTreeDatum> test = std::make_shared<InteractionTreeDatum>(*parent);
    while(true) {
        ++depth;
        if(test->parent==NULL) return depth;
        test = std::make_shared<InteractionTreeDatum>(*(test->parent));
    }
    return -1;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(std::shared_ptr<InteractionTreeDatum> datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    if (parent) {
        datum->parent = parent;
        parent->daughters.push_back(datum);
    }
    tree.push_back(datum);
    return datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionTreeDatum& datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> _datum = std::make_shared<InteractionTreeDatum>(datum);
    if (parent) {
        _datum->parent = parent;
        parent->daughters.push_back(_datum);
    }
    tree.push_back(_datum);
    return _datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionRecord& record,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    if (parent) {
        datum->parent = parent;
        parent->daughters.push_back(datum);
    }
    tree.push_back(datum);
    return datum;
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

TEST(TreeConstructor, Default)
{
    InteractionTree tree;
    EXPECT_EQ(tree.tree.size(), 0);
}

TEST(TreeAddEntry, Record)
{
    InteractionTree tree;
    InteractionRecord record;
    record.signature.primary_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.target_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);

    std::shared_ptr<InteractionTreeDatum> datum = tree.add_entry(record);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
}

TEST(TreeAddEntry, DatumReference)
{
    InteractionTree tree;
    InteractionRecord record;
    record.signature.primary_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.target_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);

    InteractionTreeDatum datum = InteractionTreeDatum(record);
    std::shared_ptr<InteractionTreeDatum> datum2 = tree.add_entry(datum);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum2->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
}

TEST(TreeAddEntry, DatumPointer)
{
    InteractionTree tree;
    InteractionRecord record;
    record.signature.primary_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.target_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);

    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    std::shared_ptr<InteractionTreeDatum> datum2 = tree.add_entry(datum);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum2->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
}

TEST(TreeAddEntry, Parent)
{
    InteractionTree tree;
    InteractionRecord record;
    record.signature.primary_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.target_type = siren::dataclasses::ParticleType::EPlus;
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    record.signature.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);

    auto datum = tree.add_entry(record);
    auto datum2 = tree.add_entry(record, datum);
    EXPECT_EQ(tree.tree.size(), 2);
    EXPECT_EQ(datum2->record.signature.primary_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.target_type, siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[0], siren::dataclasses::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.signature.secondary_types[1], siren::dataclasses::ParticleType::EMinus);
    EXPECT_EQ(datum2->parent, datum);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

