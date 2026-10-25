#include <gtest/gtest.h>

#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

#include "SIREN/utilities/Random.h"

using siren::utilities::SIREN_random;

namespace {
std::string NativeState(uint64_t seed, unsigned int draws) {
    std::mt19937_64 engine(seed);
    engine.discard(draws);
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << engine;
    return stream.str();
}

std::string LegacyArchive(uint64_t seed, std::string const & state, unsigned int version) {
    std::ostringstream stream;
    {
        cereal::JSONOutputArchive archive(stream);
        archive(cereal::make_nvp("Seed", seed));
        if(version == 1) archive(cereal::make_nvp("EngineState", state));
    }
    return stream.str();
}

void Load(SIREN_random & random, std::string const & text, unsigned int version) {
    std::istringstream stream(text);
    cereal::JSONInputArchive archive(stream);
    random.load(archive, version);
}

void SameContinuation(SIREN_random left, SIREN_random right) {
    EXPECT_EQ(left.get_seed(), right.get_seed());
    for(int i = 0; i < 1000; ++i) EXPECT_EQ(left.Uniform(), right.Uniform());
}

class GroupedNumbers : public std::numpunct<char> {
    char do_thousands_sep() const override { return ','; }
    std::string do_grouping() const override { return "\3"; }
};

struct RestoreLocale {
    std::locale original = std::locale();
    ~RestoreLocale() { std::locale::global(original); }
};
}

TEST(RandomArchive, LegacySeedAndNativeStateContinueExactly) {
    for(uint64_t seed : {uint64_t(0), uint64_t(17), std::numeric_limits<uint64_t>::max()}) {
        for(unsigned int draws : {0, 1, 7, 311, 312, 313, 400, 10000}) {
            SIREN_random expected(seed), restored(999);
            for(unsigned int i = 0; i < draws; ++i) expected.Uniform();
            Load(restored, LegacyArchive(seed, NativeState(seed, draws), 1), 1);
            SameContinuation(expected, restored);
            Load(restored, LegacyArchive(seed, "", 0), 0);
            SameContinuation(SIREN_random(seed), restored);
        }
    }
}

TEST(RandomArchive, VersionTwoBinaryRoundtrip) {
    for(unsigned int draws : {0, 7, 312, 400, 10000}) {
        SIREN_random original(17), restored(999);
        for(unsigned int i = 0; i < draws; ++i) original.Uniform();
        std::stringstream stream;
        { cereal::BinaryOutputArchive archive(stream); archive(original); }
        { cereal::BinaryInputArchive archive(stream); archive(restored); }
        SameContinuation(original, restored);
    }
}

TEST(RandomArchive, ForeignAndCorruptLegacyStateLeaveDestinationUntouched) {
    auto state = NativeState(17, 7);
    auto end = state.find_last_of(' ');
    std::vector<std::string> invalid = {
        "", state + " 0", state + " garbage", state.substr(0, end),
        "-1" + state.substr(state.find(' ')),
        "+1" + state.substr(state.find(' ')),
        "18446744073709551616" + state.substr(state.find(' ')),
        "1x" + state.substr(state.find(' '))
    };
    std::string zero;
    for(int i = 0; i < 312; ++i) zero += "0 ";
#if defined(__GLIBCXX__)
    zero += "0";
    invalid.push_back(state.substr(0, end) + " 313");
#endif
    invalid.push_back(zero);
    for(auto const & bad : invalid) {
        SIREN_random target(51);
        for(int i = 0; i < 9; ++i) target.Uniform();
        auto before = target;
        EXPECT_THROW(Load(target, LegacyArchive(17, bad, 1), 1), std::runtime_error);
        SameContinuation(before, target);
    }
}

TEST(RandomArchive, ForeignVersionTwoFormatAndUnsupportedVersionAreAtomic) {
    SIREN_random original(17), target(51);
    std::stringstream stream;
    { cereal::JSONOutputArchive archive(stream); original.save(archive, 2); }
    auto text = stream.str();
    auto begin = text.find("mt19937_64/");
    ASSERT_NE(begin, std::string::npos);
    text.replace(begin, text.find('"', begin) - begin, "unknown-library");
    auto before = target;
    EXPECT_THROW(Load(target, text, 2), std::runtime_error);
    SameContinuation(before, target);
    EXPECT_THROW(Load(target, stream.str(), 3), std::runtime_error);
    SameContinuation(before, target);
}

TEST(RandomArchive, VersionTwoIgnoresGlobalNumberGrouping) {
    RestoreLocale restore;
    std::locale::global(std::locale(std::locale::classic(), new GroupedNumbers));
    SIREN_random original(17), restored(999);
    for(int i = 0; i < 400; ++i) original.Uniform();
    std::stringstream stream;
    { cereal::JSONOutputArchive archive(stream); archive(original); }
    { cereal::JSONInputArchive archive(stream); archive(restored); }
    SameContinuation(original, restored);
}
