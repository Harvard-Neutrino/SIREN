// Tests for GDML parser functionality:
//   - Expression evaluator (constants, arithmetic, precedence, parens)
//   - Density unit conversion
//   - Composite-of-composite material resolution
//   - Circular volume reference detection
//   - Forward-reference boolean solid parsing

#include <cmath>
#include <set>
#include <string>
#include <fstream>
#include <cstdio>
#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "SIREN/detector/GDMLParser.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/math/Vector3D.h"

using namespace siren::detector;

namespace {

// Helper: write a GDML string to a temp file, parse it, clean up
GDMLData ParseGDMLString(std::string const & gdml_content) {
    std::string tmpfile = "/tmp/siren_gdml_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml_content;
    }
    GDMLData data = ParseGDML(tmpfile);
    std::remove(tmpfile.c_str());
    return data;
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
