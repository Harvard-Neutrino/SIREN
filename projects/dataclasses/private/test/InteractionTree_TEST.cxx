
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

#include "LeptonInjector/dataclasses/InteractionTree.h"

using namespace LI::dataclasses;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

struct InteractionTreeDatum {
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
  std::set<std::shared_ptr<dataclasses::InteractionTreeDatum>> tree;
  std::shared_ptr<InteractionTreeDatum> add_entry(std::shared_prt<dataclasses::InteractionTreeDatum> datum,
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
    tree.insert(datum);
    return datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionTreeDatum& datum,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> _datum = std::make_shared<InteractionTreeDatum>(datum);
    if (parent) {
        _datum->parent = parent;
        parent->daughters.push_back(_datum);
    }
    tree.insert(_datum);
    return _datum;
}

std::shared_ptr<InteractionTreeDatum> InteractionTree::add_entry(InteractionRecord& record,
        std::shared_ptr<InteractionTreeDatum> parent) {
    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    if (parent) {
        datum->parent = parent;
        parent->daughters.push_back(datum);
    }
    tree.insert(datum);
    return datum;
}

TEST(DatumConstructor, InteractionRecord)
{
    InteractionRecord record;
    record.primary_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.target_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EPlus);
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EMinus);

    InteractionTreeDatum datum(record);
    EXPECT_EQ(datum.record.primary_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum.record.target_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum.record.secondary_types[0], LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum.record.secondary_types[1], LI::dataclasses::Particle::ParticleType::EMinus);
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
    record.primary_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.target_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EPlus);
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EMinus);

    std::shared_ptr<InteractionTreeDatum> datum = tree.add_entry(record);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum->record.primary_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum->record.target_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum->record.secondary_types[0], LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum->record.secondary_types[1], LI::dataclasses::Particle::ParticleType::EMinus);
}

TEST(TreeAddEntry, DatumReference)
{
    InteractionTree tree;
    InteractionRecord record;
    record.primary_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.target_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EPlus);
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EMinus);

    InteractionTreeDatum datum = InteractionTreeDatum(record);
    std::shared_ptr<InteractionTreeDatum> datum2 = tree.add_entry(datum);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum2->record.primary_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.target_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.secondary_types[0], LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.secondary_types[1], LI::dataclasses::Particle::ParticleType::EMinus);
}

TEST(TreeAddEntry, DatumPointer)
{
    InteractionTree tree;
    InteractionRecord record;
    record.primary_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.target_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EPlus);
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EMinus);

    std::shared_ptr<InteractionTreeDatum> datum = std::make_shared<InteractionTreeDatum>(record);
    std::shared_ptr<InteractionTreeDatum> datum2 = tree.add_entry(datum);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum2->record.primary_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.target_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.secondary_types[0], LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.secondary_types[1], LI::dataclasses::Particle::ParticleType::EMinus);
}

TEST(TreeAddEntry, Parent)
{
    InteractionTree tree;
    InteractionRecord record;
    record.primary_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.target_type = LI::dataclasses::Particle::ParticleType::EPlus;
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EPlus);
    record.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EMinus);

    auto datum = tree.add_entry(record);
    auto datum2 = tree.add_entry(record, datum);
    EXPECT_EQ(tree.tree.size(), 1);
    EXPECT_EQ(datum2->record.primary_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.target_type, LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.secondary_types[0], LI::dataclasses::Particle::ParticleType::EPlus);
    EXPECT_EQ(datum2->record.secondary_types[1], LI::dataclasses::Particle::ParticleType::EMinus);
    EXPECT_EQ(datum2->parent, datum);
}


TEST(Serialization, Save)
{
    InteractionTree A;
    A.primary_type = LI::dataclasses::Particle::ParticleType::EPlus;
    A.target_type = LI::dataclasses::Particle::ParticleType::EPlus;
    A.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EPlus);
    A.secondary_types.push_back(LI::dataclasses::Particle::ParticleType::EMinus);

    std::stringstream ss;
    {
        cereal::JSONOutputArchive oarchive(ss);
        oarchive(cereal::make_nvp("InteractionTree", A));
    }

    std::string expected = "{\n    \"InteractionTree\": {\n        \"cereal_class_version\": 0,\n        \"PrimaryType\": -11,\n        \"TargetType\": -11,\n        \"SecondaryTypes\": [\n            -11,\n            11\n        ]\n    }\n}";
    EXPECT_EQ(ss.str(), expected);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

