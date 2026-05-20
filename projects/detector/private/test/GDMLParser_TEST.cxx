// Tests for GDML parser functionality:
//   - Expression evaluator (constants, arithmetic, precedence, parens)
//   - Density unit conversion
//   - Composite-of-composite material resolution
//   - Circular volume reference detection
//   - Forward-reference boolean solid parsing

#include <atomic>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <utility>

#include <unistd.h>

#include <gtest/gtest.h>

#include "SIREN/detector/GDMLParser.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/math/Vector3D.h"

using namespace siren::detector;

namespace {

// RAII temp file: writes content on construction, removes on destruction
// (even if the test throws inside EXPECT_THROW).
struct TempGDMLFile {
    std::string path;
    TempGDMLFile(std::string const & content) {
        static std::atomic<int> counter{0};
        std::string name = "siren_gdml_"
            + std::to_string(::getpid()) + "_"
            + std::to_string(counter.fetch_add(1)) + ".gdml";
        path = (std::filesystem::temp_directory_path() / name).string();
        std::ofstream f(path);
        f << content;
    }
    ~TempGDMLFile() { std::remove(path.c_str()); }
    TempGDMLFile(TempGDMLFile const &) = delete;
    TempGDMLFile & operator=(TempGDMLFile const &) = delete;
};

// Unique temp path for multi-file tests that manage their own lifecycle.
std::string TempPath(std::string const & tag) {
    static std::atomic<int> path_counter{0};
    std::string name = "siren_" + tag + "_"
        + std::to_string(::getpid()) + "_"
        + std::to_string(path_counter.fetch_add(1)) + ".gdml";
    return (std::filesystem::temp_directory_path() / name).string();
}

GDMLData ParseGDMLString(std::string const & gdml_content) {
    TempGDMLFile tmp(gdml_content);
    return ParseGDML(tmp.path);
}

std::string GDMLWithSolids(std::string const & solids, std::string const & world_solid_ref) {
    return std::string(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>)") + solids + R"(
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref=")" + world_solid_ref + R"("/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";
}

GDMLData ParseGDMLWithSolids(std::string const & solids, std::string const & world_solid_ref) {
    return ParseGDMLString(GDMLWithSolids(solids, world_solid_ref));
}

std::string GDMLWithSolidsAndDefine(std::string const & solids, std::string const & world_solid_ref,
                                    std::string const & define_content) {
    return std::string(R"(<?xml version="1.0"?>
<gdml>
  <define>)") + define_content + R"(
  </define>
  <materials/>
  <solids>)" + solids + R"(
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref=")" + world_solid_ref + R"("/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";
}

GDMLData ParseGDMLWithSolids(std::string const & solids, std::string const & world_solid_ref,
                             std::string const & define_content) {
    return ParseGDMLString(GDMLWithSolidsAndDefine(solids, world_solid_ref, define_content));
}

std::shared_ptr<siren::geometry::Geometry> RequireSolid(GDMLData const & data, std::string const & name) {
    auto it = data.solids.find(name);
    if(it == data.solids.end()) {
        ADD_FAILURE() << "Solid '" << name << "' should be parsed";
        return nullptr;
    }
    return it->second;
}

bool IsInside(std::shared_ptr<siren::geometry::Geometry> const & geo, double x, double y, double z) {
    static const siren::math::Vector3D dir(0, 0, 1);
    return geo->IsInside(siren::math::Vector3D(x, y, z), dir);
}

std::string GDMLWithDefine(std::string const & define_content) {
    return R"(<?xml version="1.0"?>
<gdml>
  <define>)" + define_content + R"(
  </define>
  <materials/>
  <solids><box name="w" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="w"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)";
}

GDMLData ParseGDMLWithDefine(std::string const & define_content) {
    return ParseGDMLString(GDMLWithDefine(define_content));
}

std::string GDMLFull(std::string const & define, std::string const & materials,
                     std::string const & solids, std::string const & structure,
                     std::string const & world_ref) {
    return R"(<?xml version="1.0"?>
<gdml>
  <define>)" + define + R"(</define>
  <materials>)" + materials + R"(</materials>
  <solids>)" + solids + R"(</solids>
  <structure>)" + structure + R"(</structure>
  <setup name="Default" version="1.0"><world ref=")" + world_ref + R"("/></setup>
</gdml>)";
}

std::unique_ptr<DetectorModel> LoadGDMLFromString(std::string const & gdml) {
    TempGDMLFile tmp(gdml);
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(tmp.path);
    return dm;
}

std::string MaterialAt(DetectorModel const & dm, double x, double y, double z) {
    auto sector = dm.GetContainingSector(
        DetectorPosition(siren::math::Vector3D(x, y, z)));
    return dm.GetMaterials().GetMaterialName(sector.material_id);
}

struct ProbePoint { double x, y, z; bool inside; };

void TestContainment(std::string const & solid_xml, std::string const & name,
                     std::vector<ProbePoint> const & probes) {
    GDMLData data = ParseGDMLWithSolids(solid_xml, name);
    auto geo = RequireSolid(data, name);
    ASSERT_NE(geo, nullptr);
    for(auto const & p : probes) {
        EXPECT_EQ(IsInside(geo, p.x, p.y, p.z), p.inside)
            << "at (" << p.x << ", " << p.y << ", " << p.z << ")";
    }
}

void ExpectNoWarningContaining(GDMLData const & data, std::string const & substr) {
    for(auto const & w : data.warnings) {
        EXPECT_EQ(w.find(substr), std::string::npos)
            << "Unexpected warning containing '" << substr << "': " << w;
    }
}

} // anonymous namespace


// Parser behavior groups. Each include is compiled in this translation unit so
// the shared ParseGDMLString helper above remains local to the test binary.
#include "gdml/GDMLParserExpressionUnits_TEST.inc"
#include "gdml/GDMLParserMaterials_TEST.inc"
#include "gdml/GDMLParserSolids_TEST.inc"
#include "gdml/GDMLParserBooleans_TEST.inc"
#include "gdml/GDMLParserStructure_TEST.inc"
#include "gdml/GDMLParserPreprocess_TEST.inc"
#include "gdml/GDMLParserMisc_TEST.inc"
