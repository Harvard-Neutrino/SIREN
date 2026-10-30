#include <gtest/gtest.h>
#include <sstream>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedAngularSectorChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"

using namespace siren::injection;
namespace {
using Channel = std::shared_ptr<PhaseSpaceChannel>;
template<class Archive> std::string Save(Channel const & channel) {
    std::stringstream stream;
    { Archive archive(stream); archive(channel); }
    return stream.str();
}
template<class Archive> Channel Load(std::string const & bytes) {
    std::stringstream stream(bytes);
    Channel channel;
    { Archive archive(stream); archive(channel); }
    return channel;
}
std::vector<Channel> Channels() {
    auto target = std::make_shared<siren::geometry::Sphere>(
        siren::geometry::Placement(siren::math::Vector3D(0,0,20)), 1., 0.);
    return {
        std::make_shared<DetectorDirected2BodyChannel>(target,1),
        std::make_shared<DetectorDirectedAngularSectorChannel>(target,.1,.9,-.5,.5),
        std::make_shared<DetectorDirected3BodyChannel>(target,0),
        std::make_shared<DetectorDirected3BodyChannel>(target,0,1,2,1),
        std::make_shared<DetectorDirectedScatteringChannel>(target,1)};
}
}
TEST(DirectedSerialization, EveryNativeFamilyRoundTripsJsonAndBinary) {
    for (auto const & channel : Channels()) {
        for (auto restored : {Load<cereal::JSONInputArchive>(Save<cereal::JSONOutputArchive>(channel)),
                              Load<cereal::BinaryInputArchive>(Save<cereal::BinaryOutputArchive>(channel))}) {
            EXPECT_EQ(restored->Name(),channel->Name());
            EXPECT_EQ(restored->Topology(),channel->Topology());
            EXPECT_TRUE(restored->Measure()==channel->Measure());
            EXPECT_EQ(Save<cereal::JSONOutputArchive>(restored),Save<cereal::JSONOutputArchive>(channel));
        }
    }
}
TEST(DirectedSerialization, ObsoleteDensityVersionsAreRefused) {
    for (auto const & channel : Channels()) {
        auto bytes=Save<cereal::JSONOutputArchive>(channel);
        int version = channel->Name()=="DetectorDirectedScattering" ? 1 : 2;
        std::string const key="\"cereal_class_version\": "+std::to_string(version);
        auto at=bytes.find(key);
        ASSERT_NE(at,std::string::npos);
        bytes.replace(at,key.size(),"\"cereal_class_version\": "+std::to_string(version-1));
        EXPECT_THROW(Load<cereal::JSONInputArchive>(bytes),std::runtime_error);
    }
}
TEST(DirectedSerialization, TamperedVolumeCannotChangeDensityOnLoad) {
    for (auto const & channel : Channels()) {
        if (channel->Name()=="DetectorDirectedAngularSector") continue;
        auto bytes=Save<cereal::JSONOutputArchive>(channel);
        std::string const key="\"TargetVolume\": ";
        auto at=bytes.find(key);
        ASSERT_NE(at,std::string::npos);
        at+=key.size();
        auto end=bytes.find_first_of(",\n}",at);
        bytes.replace(at,end-at,"1000000.0");
        EXPECT_THROW(Load<cereal::JSONInputArchive>(bytes),std::runtime_error);
    }
}

namespace {
std::string ReplaceNumber(std::string bytes, std::string const & field, std::string const & value) {
    std::string key="\""+field+"\": ";
    auto at=bytes.find(key);
    if(at==std::string::npos) throw std::runtime_error("Missing archive field");
    at+=key.size();
    auto end=bytes.find_first_of(",\n}",at);
    bytes.replace(at,end-at,value);
    return bytes;
}
}
TEST(DirectedSerialization, CorruptIndicesAndSectorBoundsAreRefused) {
    for (auto const & channel : Channels()) {
        auto bytes=Save<cereal::JSONOutputArchive>(channel);
        auto field=bytes.find("\"DaughterIndex\"")!=std::string::npos ? "DaughterIndex" : "DirectedIndex";
        EXPECT_ANY_THROW(Load<cereal::JSONInputArchive>(ReplaceNumber(bytes,field,"999")));
        if (channel->Name()=="DetectorDirectedAngularSector")
            EXPECT_ANY_THROW(Load<cereal::JSONInputArchive>(ReplaceNumber(bytes,"ULower","2")));
        if (channel->Name()=="DetectorDirected3Body")
            EXPECT_ANY_THROW(Load<cereal::JSONInputArchive>(ReplaceNumber(bytes,"PairFirstIndex","9")));
    }
}

TEST(DirectedSerialization, FailedLoadPreservesExistingObject) {
    auto channel=std::dynamic_pointer_cast<DetectorDirected2BodyChannel>(Channels()[0]);
    auto before=Save<cereal::JSONOutputArchive>(channel);
    std::stringstream encoded;
    { cereal::JSONOutputArchive archive(encoded); archive(*channel); }
    std::stringstream broken(ReplaceNumber(encoded.str(),"DaughterIndex","999"));
    EXPECT_ANY_THROW({ cereal::JSONInputArchive archive(broken); archive(*channel); });
    EXPECT_EQ(Save<cereal::JSONOutputArchive>(channel),before);
}
