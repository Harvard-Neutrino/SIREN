// Tests for GDML parser functionality:
//   - Expression evaluator (constants, arithmetic, precedence, parens)
//   - Density unit conversion
//   - Composite-of-composite material resolution
//   - Circular volume reference detection
//   - Forward-reference boolean solid parsing

#include <cmath>
#include <string>
#include <fstream>
#include <cstdio>

#include <gtest/gtest.h>

#include "SIREN/detector/GDMLParser.h"
#include "SIREN/detector/DetectorModel.h"

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

} // anonymous namespace


// =========================================================================
// Expression evaluator: constants, arithmetic, operator precedence
// =========================================================================
TEST(GDMLParser, ExpressionEvaluator) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="lit" value="3.14"/>
    <constant name="neg_lit" value="-2.5"/>
    <constant name="add_expr" value="1.0+2.0"/>
    <constant name="sub_expr" value="10.0-3.0"/>
    <constant name="mul_expr" value="3.0*4.0"/>
    <constant name="div_expr" value="10.0/4.0"/>
    <constant name="precedence" value="2.0+3.0*4.0"/>
    <constant name="parens" value="(2.0+3.0)*4.0"/>
    <constant name="ref_const" value="lit"/>
    <constant name="ref_arith" value="lit*2.0"/>
    <constant name="neg_ref" value="-lit"/>
    <constant name="chain" value="add_expr+1.0"/>
    <constant name="complex" value="(lit+neg_lit)*2.0"/>
  </define>
  <materials/>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    double tol = 1e-9;
    EXPECT_NEAR(data.constants["lit"], 3.14, tol);
    EXPECT_NEAR(data.constants["neg_lit"], -2.5, tol);
    EXPECT_NEAR(data.constants["add_expr"], 3.0, tol);
    EXPECT_NEAR(data.constants["sub_expr"], 7.0, tol);
    EXPECT_NEAR(data.constants["mul_expr"], 12.0, tol);
    EXPECT_NEAR(data.constants["div_expr"], 2.5, tol);
    // Precedence: 2 + 3*4 = 14, not (2+3)*4 = 20
    EXPECT_NEAR(data.constants["precedence"], 14.0, tol);
    EXPECT_NEAR(data.constants["parens"], 20.0, tol);
    // Constant references
    EXPECT_NEAR(data.constants["ref_const"], 3.14, tol);
    EXPECT_NEAR(data.constants["ref_arith"], 6.28, tol);
    EXPECT_NEAR(data.constants["neg_ref"], -3.14, tol);
    // Chained: add_expr=3, chain=3+1=4
    EXPECT_NEAR(data.constants["chain"], 4.0, tol);
    // Complex: (3.14 + (-2.5)) * 2 = 0.64 * 2 = 1.28
    EXPECT_NEAR(data.constants["complex"], 1.28, tol);
}


// =========================================================================
// Expression evaluator: solid dimensions using expressions
// =========================================================================
TEST(GDMLParser, ExpressionInSolidDimensions) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="half_size" value="50"/>
  </define>
  <materials/>
  <solids>
    <box name="test_box" x="half_size*2" y="half_size" z="half_size/2" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="test_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    // The box should exist as a solid
    ASSERT_TRUE(data.solids.count("test_box") > 0);
    // Verify dimensions via bounding box (in meters, since lunit=mm -> *0.001)
    // x = 50*2 = 100mm = 0.1m (full width = 0.2m since GDML x is half-width)
    // y = 50mm = 0.05m (full width = 0.1m)
    // z = 50/2 = 25mm = 0.025m (full width = 0.05m)
    auto bb = data.solids["test_box"]->GetBoundingBox();
    double tol = 1e-9;
    EXPECT_NEAR(bb.max_corner.GetX(), 0.1, tol);  // half of full width 0.2
    EXPECT_NEAR(bb.max_corner.GetY(), 0.05, tol); // half of full width 0.1
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.025, tol); // half of full width 0.05
}


// =========================================================================
// Density unit conversion
// =========================================================================
TEST(GDMLParser, DensityUnits) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="mat_gcm3">
      <D value="2.7" unit="g/cm3"/>
      <atom value="27" Z="13"/>
    </material>
    <material name="mat_kgm3">
      <D value="2700" unit="kg/m3"/>
      <atom value="27" Z="13"/>
    </material>
    <material name="mat_nounit">
      <D value="2.7"/>
      <atom value="27" Z="13"/>
    </material>
    <material name="mat_mgcm3">
      <D value="2700" unit="mg/cm3"/>
      <atom value="27" Z="13"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="mat_gcm3"/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    double tol = 1e-9;
    // All should resolve to 2.7 g/cm3
    EXPECT_NEAR(data.materials["mat_gcm3"].density, 2.7, tol);
    EXPECT_NEAR(data.materials["mat_kgm3"].density, 2.7, tol);
    EXPECT_NEAR(data.materials["mat_nounit"].density, 2.7, tol);
    EXPECT_NEAR(data.materials["mat_mgcm3"].density, 2.7, tol);
}


// =========================================================================
// Composite-of-composite material resolution via LoadGDML
// =========================================================================
TEST(GDMLParser, CompositeOfCompositeMaterial) {
    // Create a GDML file where Water is composed of elements H and O,
    // and a composite material "Mix" is composed of Water and another element.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <element name="Hydrogen" Z="1">
      <atom value="1.008"/>
    </element>
    <element name="Oxygen" Z="8">
      <atom value="15.999"/>
    </element>
    <element name="Carbon" Z="6">
      <atom value="12.011"/>
    </element>
    <material name="Water">
      <D value="1.0" unit="g/cm3"/>
      <fraction ref="Hydrogen" n="0.112"/>
      <fraction ref="Oxygen" n="0.888"/>
    </material>
    <material name="Mix">
      <D value="1.5" unit="g/cm3"/>
      <fraction ref="Water" n="0.8"/>
      <fraction ref="Carbon" n="0.2"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="inner_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="Inner">
      <materialref ref="Mix"/>
      <solidref ref="inner_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Water"/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="Inner"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_composite_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }

    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    // The "Mix" material should have been resolved recursively:
    // Water(0.8) -> H(0.8*0.112) + O(0.8*0.888)
    // Carbon(0.2)
    // So Mix = {H: 0.0896, O: 0.7104, C: 0.2}
    // All three should be present as PDG codes
    auto const & materials = dm.GetMaterials();
    ASSERT_TRUE(materials.HasMaterial("Mix"));
    int mix_id = materials.GetMaterialId("Mix");
    auto constituents = materials.GetMaterialConstituents(mix_id);
    // Should have 3 constituents (H, O, C)
    EXPECT_GE(constituents.size(), 3u)
        << "Mix should have at least 3 elemental constituents from recursive resolution";
}


// =========================================================================
// Forward-reference boolean solid (second operand defined before first)
// =========================================================================
TEST(GDMLParser, BooleanForwardReference) {
    // Define the boolean solid BEFORE one of its operands in the <solids> section.
    // With two-pass parsing, this should work.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <subtraction name="result">
      <first ref="big_sphere"/>
      <second ref="small_sphere"/>
    </subtraction>
    <sphere name="big_sphere" rmax="100" lunit="mm"/>
    <sphere name="small_sphere" rmax="50" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="result"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    // The boolean solid should have been parsed despite forward references
    ASSERT_TRUE(data.solids.count("result") > 0)
        << "Boolean solid with forward-referenced operands should be parsed";
    ASSERT_TRUE(data.solids.count("big_sphere") > 0);
    ASSERT_TRUE(data.solids.count("small_sphere") > 0);
}


// =========================================================================
// Circular volume reference detection
// =========================================================================
TEST(GDMLParser, CircularVolumeReference) {
    // Create a GDML where volume A contains physvol referencing A (self-loop)
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Air">
      <D value="0.00129" unit="g/cm3"/>
      <atom value="14" Z="7"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="World"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_cycle_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }

    // Should not crash or infinite-loop; cycle detection should break it
    DetectorModel dm;
    EXPECT_NO_THROW(dm.LoadGDML(tmpfile))
        << "LoadGDML should handle circular volume references without crashing";
    std::remove(tmpfile.c_str());
}
