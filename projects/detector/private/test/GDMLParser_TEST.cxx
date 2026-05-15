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


// =========================================================================
// Physvol inline position using defined constants
// =========================================================================
TEST(GDMLParser, PhysvolConstantsPosition) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="offset" value="500"/>
  </define>
  <materials>
    <element name="Hydrogen" Z="1">
      <atom value="1.008"/>
    </element>
    <element name="Oxygen" Z="8">
      <atom value="15.999"/>
    </element>
    <material name="Water">
      <D value="1.0" unit="g/cm3"/>
      <fraction ref="Hydrogen" n="0.112"/>
      <fraction ref="Oxygen" n="0.888"/>
    </material>
    <material name="Iron" Z="26">
      <D value="7.874" unit="g/cm3"/>
      <atom value="55.845"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="child_box" x="200" y="200" z="200" lunit="mm"/>
  </solids>
  <structure>
    <volume name="Child">
      <materialref ref="Iron"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Water"/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="Child"/>
        <position x="offset" unit="mm"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    // Verify via ParseGDMLString that the physvol position resolved the constant
    GDMLData data = ParseGDMLString(gdml);

    ASSERT_TRUE(data.volumes.count("World") > 0);
    auto const & world = data.volumes["World"];
    ASSERT_GE(world.children.size(), 1u);
    // offset=500, unit=mm -> 500mm = 0.5m
    double tol = 1e-9;
    EXPECT_NEAR(world.children[0].position.GetX(), 0.5, tol)
        << "Physvol position should resolve constant 'offset' to 500mm = 0.5m";
    EXPECT_NEAR(world.children[0].position.GetY(), 0.0, tol);
    EXPECT_NEAR(world.children[0].position.GetZ(), 0.0, tol);

    // Also verify via LoadGDML that the sector is placed correctly
    std::string tmpfile = "/tmp/siren_gdml_physvol_const_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    // A point at (0.5, 0, 0) meters should be inside the Iron child
    auto sector = dm.GetContainingSector(
        DetectorPosition(siren::math::Vector3D(0.5, 0.0, 0.0)));
    std::string mat_name = dm.GetMaterials().GetMaterialName(sector.material_id);
    EXPECT_EQ(mat_name, "Iron")
        << "Point at (0.5,0,0) should be in Iron child placed at offset=500mm";
}


// =========================================================================
// GDML default angle unit is degrees
// =========================================================================
TEST(GDMLParser, DefaultAngleUnitDegrees) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <rotation name="rot_nounit" z="90"/>
    <rotation name="rot_deg" z="90" unit="deg"/>
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

    ASSERT_TRUE(data.rotations.count("rot_nounit") > 0);
    ASSERT_TRUE(data.rotations.count("rot_deg") > 0);

    auto const & q_nounit = data.rotations["rot_nounit"];
    auto const & q_deg = data.rotations["rot_deg"];

    double tol = 1e-9;
    // Both should produce the same quaternion (default unit is degrees)
    EXPECT_NEAR(q_nounit.GetX(), q_deg.GetX(), tol);
    EXPECT_NEAR(q_nounit.GetY(), q_deg.GetY(), tol);
    EXPECT_NEAR(q_nounit.GetZ(), q_deg.GetZ(), tol);
    EXPECT_NEAR(q_nounit.GetW(), q_deg.GetW(), tol);

    // A z-rotation of 90 degrees = pi/2 radians.
    // Quaternion for z-rotation by angle a: (0, 0, sin(a/2), cos(a/2))
    // For a = pi/2: (0, 0, sin(pi/4), cos(pi/4)) = (0, 0, sqrt(2)/2, sqrt(2)/2)
    double s = std::sqrt(2.0) / 2.0;
    EXPECT_NEAR(q_nounit.GetX(), 0.0, tol);
    EXPECT_NEAR(q_nounit.GetY(), 0.0, tol);
    EXPECT_NEAR(q_nounit.GetZ(), s, tol)
        << "90-degree z-rotation should have z = sqrt(2)/2, not a radian-based value";
    EXPECT_NEAR(q_nounit.GetW(), s, tol);
}


// =========================================================================
// Length unit conversion: mm, cm, m all representing the same size
// =========================================================================
TEST(GDMLParser, LengthUnitConversion) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="box_mm" x="100" y="100" z="100" lunit="mm"/>
    <box name="box_cm" x="10" y="10" z="10" lunit="cm"/>
    <box name="box_m" x="0.1" y="0.1" z="0.1" lunit="m"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="box_mm"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    ASSERT_TRUE(data.solids.count("box_mm") > 0);
    ASSERT_TRUE(data.solids.count("box_cm") > 0);
    ASSERT_TRUE(data.solids.count("box_m") > 0);

    auto bb_mm = data.solids["box_mm"]->GetBoundingBox();
    auto bb_cm = data.solids["box_cm"]->GetBoundingBox();
    auto bb_m  = data.solids["box_m"]->GetBoundingBox();

    // GDML <box> x,y,z are HALF-widths. All three boxes specify the same half-width:
    // 100mm = 10cm = 0.1m. In SIREN meters, bounding box max_corner = 0.1m.
    double tol = 1e-9;
    double expected = 0.1; // half-width in meters

    // All three should produce the same bounding box
    EXPECT_NEAR(bb_mm.max_corner.GetX(), expected, tol);
    EXPECT_NEAR(bb_mm.max_corner.GetY(), expected, tol);
    EXPECT_NEAR(bb_mm.max_corner.GetZ(), expected, tol);

    EXPECT_NEAR(bb_cm.max_corner.GetX(), expected, tol);
    EXPECT_NEAR(bb_cm.max_corner.GetY(), expected, tol);
    EXPECT_NEAR(bb_cm.max_corner.GetZ(), expected, tol);

    EXPECT_NEAR(bb_m.max_corner.GetX(), expected, tol);
    EXPECT_NEAR(bb_m.max_corner.GetY(), expected, tol);
    EXPECT_NEAR(bb_m.max_corner.GetZ(), expected, tol);
}


// =========================================================================
// All solid types parse correctly
// =========================================================================
TEST(GDMLParser, AllSolidTypes) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="b1" x="100" y="100" z="100" lunit="mm"/>
    <sphere name="s1" rmax="50" lunit="mm"/>
    <tube name="t1" rmin="0" rmax="50" z="100" deltaphi="360" aunit="deg" lunit="mm"/>
    <cone name="c1" rmin1="0" rmax1="50" rmin2="0" rmax2="30" z="100" deltaphi="360" aunit="deg" lunit="mm"/>
    <polycone name="pc1" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-50"/>
      <zplane rmin="0" rmax="30" z="50"/>
    </polycone>
    <polyhedra name="ph1" numsides="6" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-40"/>
      <zplane rmin="0" rmax="50" z="40"/>
    </polyhedra>
    <trd name="trd1" x1="60" x2="40" y1="80" y2="50" z="100" lunit="mm"/>
    <xtru name="xt1" lunit="mm">
      <twoDimVertex x="-50" y="-50"/>
      <twoDimVertex x="50" y="-50"/>
      <twoDimVertex x="50" y="50"/>
      <twoDimVertex x="-50" y="50"/>
      <section zOrder="0" zPosition="-40" xOffset="0" yOffset="0" scalingFactor="1"/>
      <section zOrder="1" zPosition="40" xOffset="0" yOffset="0" scalingFactor="1"/>
    </xtru>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="b1"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    std::vector<std::string> names = {"b1", "s1", "t1", "c1", "pc1", "ph1", "trd1", "xt1"};
    for(auto const & name : names) {
        ASSERT_TRUE(data.solids.count(name) > 0)
            << "Solid '" << name << "' should be parsed";
        auto bb = data.solids[name]->GetBoundingBox();
        EXPECT_TRUE(bb.IsValid())
            << "Solid '" << name << "' should have a valid bounding box";
        EXPECT_TRUE(std::isfinite(bb.max_corner.GetX()))
            << "Solid '" << name << "' bounding box should have finite dimensions";
        EXPECT_TRUE(std::isfinite(bb.max_corner.GetY()))
            << "Solid '" << name << "' bounding box should have finite dimensions";
        EXPECT_TRUE(std::isfinite(bb.max_corner.GetZ()))
            << "Solid '" << name << "' bounding box should have finite dimensions";
    }
}


// =========================================================================
// Nested boolean solid resolution (forward references, multiple levels)
// =========================================================================
TEST(GDMLParser, NestedBooleanSolids) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <union name="final">
      <first ref="inner_sub"/>
      <second ref="some_box"/>
    </union>
    <subtraction name="inner_sub">
      <first ref="big_sphere"/>
      <second ref="small_sphere"/>
    </subtraction>
    <sphere name="big_sphere" rmax="100" lunit="mm"/>
    <sphere name="small_sphere" rmax="50" lunit="mm"/>
    <box name="some_box" x="40" y="40" z="40" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="final"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    // All composite/primitive solids should exist
    ASSERT_TRUE(data.solids.count("final") > 0)
        << "Outer union 'final' should be parsed despite nested forward references";
    ASSERT_TRUE(data.solids.count("inner_sub") > 0)
        << "Inner subtraction 'inner_sub' should be parsed";
    ASSERT_TRUE(data.solids.count("big_sphere") > 0);
    ASSERT_TRUE(data.solids.count("small_sphere") > 0);
    ASSERT_TRUE(data.solids.count("some_box") > 0);

    // Verify the boolean solids have valid bounding boxes
    auto bb_final = data.solids["final"]->GetBoundingBox();
    auto bb_inner = data.solids["inner_sub"]->GetBoundingBox();
    EXPECT_TRUE(bb_final.IsValid())
        << "Nested boolean 'final' should have a valid bounding box";
    EXPECT_TRUE(bb_inner.IsValid())
        << "Nested boolean 'inner_sub' should have a valid bounding box";
}


// =========================================================================
// Duplicate volume names: same logical volume placed twice
// =========================================================================
TEST(GDMLParser, DuplicateVolumeNames) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26">
      <D value="7.874" unit="g/cm3"/>
      <atom value="55.845"/>
    </material>
    <material name="Air" Z="7">
      <D value="0.00129" unit="g/cm3"/>
      <atom value="14.007"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="child_box" x="200" y="200" z="200" lunit="mm"/>
  </solids>
  <structure>
    <volume name="Child">
      <materialref ref="Iron"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="Child"/>
        <position x="100" unit="mm"/>
      </physvol>
      <physvol>
        <volumeref ref="Child"/>
        <position x="-100" unit="mm"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_dupvol_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }

    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    // World + 2 children + UNIVERSE = at least 4 sectors
    auto const & sectors = dm.GetSectors();
    EXPECT_GE(sectors.size(), 4u)
        << "Detector model should have at least 4 sectors (UNIVERSE + world + 2 duplicate children)";

    // Both child positions should resolve to Iron
    auto sector_pos = dm.GetContainingSector(
        DetectorPosition(siren::math::Vector3D(0.1, 0.0, 0.0)));
    std::string mat_pos = dm.GetMaterials().GetMaterialName(sector_pos.material_id);
    EXPECT_EQ(mat_pos, "Iron")
        << "Point at (0.1,0,0) should be in first Iron child";

    auto sector_neg = dm.GetContainingSector(
        DetectorPosition(siren::math::Vector3D(-0.1, 0.0, 0.0)));
    std::string mat_neg = dm.GetMaterials().GetMaterialName(sector_neg.material_id);
    EXPECT_EQ(mat_neg, "Iron")
        << "Point at (-0.1,0,0) should be in second Iron child";
}


// =========================================================================
// Rotation handling: child volume placed with rotation
// =========================================================================
TEST(GDMLParser, RotationHandling) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26">
      <D value="7.874" unit="g/cm3"/>
      <atom value="55.845"/>
    </material>
    <material name="Air" Z="7">
      <D value="0.00129" unit="g/cm3"/>
      <atom value="14.007"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="child_box" x="400" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="Child">
      <materialref ref="Iron"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="Child"/>
        <position x="200" unit="mm"/>
        <rotation z="90" unit="deg"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_rotation_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }

    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    // The child box is 400x100x100 mm (half-widths: 0.4x0.1x0.1 m)
    // Placed at x=200mm=0.2m with z-rotation of 90 degrees.
    // After 90-deg z-rotation, the long axis (originally x) becomes the y axis.
    // So the child occupies roughly: x in [0.1, 0.3], y in [-0.4, 0.4]

    // The child sector should exist
    auto const & sectors = dm.GetSectors();
    EXPECT_GE(sectors.size(), 3u)
        << "Should have at least 3 sectors (UNIVERSE + world + rotated child)";

    // A point at the child center (0.2, 0, 0) should be Iron
    auto sector_center = dm.GetContainingSector(
        DetectorPosition(siren::math::Vector3D(0.2, 0.0, 0.0)));
    std::string mat_center = dm.GetMaterials().GetMaterialName(sector_center.material_id);
    EXPECT_EQ(mat_center, "Iron")
        << "Center of rotated child at (0.2,0,0) should be Iron";

    // After rotation, the long axis is along y. A point at (0.2, 0.3, 0) should
    // still be inside (within the rotated half-width of 0.4 along y).
    auto sector_along_y = dm.GetContainingSector(
        DetectorPosition(siren::math::Vector3D(0.2, 0.3, 0.0)));
    std::string mat_along_y = dm.GetMaterials().GetMaterialName(sector_along_y.material_id);
    EXPECT_EQ(mat_along_y, "Iron")
        << "Point at (0.2, 0.3, 0) should be inside the rotated child (long axis now along y)";

    // A point far from the child should be in Air
    auto sector_far = dm.GetContainingSector(
        DetectorPosition(siren::math::Vector3D(1.0, 1.0, 0.0)));
    std::string mat_far = dm.GetMaterials().GetMaterialName(sector_far.material_id);
    EXPECT_EQ(mat_far, "Air")
        << "Point at (1.0,1.0,0) should be in Air (outside rotated child)";
}


// =========================================================================
// GDML error handling: unknown constant throws
// =========================================================================
TEST(GDMLParser, UnknownConstantThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="bad" value="nonexistent_name"/>
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

    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error);
}

// =========================================================================
// GDML error handling: unknown length unit throws
// =========================================================================
TEST(GDMLParser, UnknownLengthUnitThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="test_box" x="100" y="100" z="100" lunit="furlongs"/>
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

    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error);
}

// =========================================================================
// GDML error handling: unknown angle unit throws
// =========================================================================
TEST(GDMLParser, UnknownAngleUnitThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <rotation name="r1" x="45" y="0" z="0" unit="gradians"/>
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

    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error);
}

// =========================================================================
// GDML error handling: unsupported solid type in strict mode throws
// =========================================================================
TEST(GDMLParser, UnsupportedSolidStrictThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <hype name="h" rmin="1" rmax="2" inst="0.1" outst="0.2" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="h"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    // Strict mode: throws on unsupported solid
    std::string tmpfile = "/tmp/siren_gdml_test_strict.gdml";
    { std::ofstream f(tmpfile); f << gdml; }
    GDMLParseOptions strict_opts;
    strict_opts.strict = true;
    EXPECT_THROW(ParseGDML(tmpfile, strict_opts), std::runtime_error);
    std::remove(tmpfile.c_str());
}

// =========================================================================
// GDML: unsupported solid skipped with warning in default mode
// =========================================================================
TEST(GDMLParser, UnsupportedSolidSkipsWithWarning) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="good_box" x="100" y="100" z="100" lunit="mm"/>
    <hype name="h" rmin="1" rmax="2" inst="0.1" outst="0.2" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="good_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    // The box should have parsed, the hype should have been skipped
    EXPECT_EQ(data.solids.count("good_box"), 1u);
    EXPECT_EQ(data.solids.count("h"), 0u);
    // Warning should be recorded
    ASSERT_GE(data.warnings.size(), 1u);
    bool found_hype_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("hype") != std::string::npos) found_hype_warning = true;
    }
    EXPECT_TRUE(found_hype_warning) << "Expected warning about unsupported hype solid";
}

// =========================================================================
// GDML error handling: unparseable expression throws
// =========================================================================
TEST(GDMLParser, UnparseableExpressionThrows) {
    // Use custom delimiter to avoid )" in sin(3.14)" closing the raw string
    std::string gdml = R"GDML(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="bad" value="sin(3.14)"/>
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
</gdml>)GDML";

    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error);
}

// =========================================================================
// GDML: valid units should still work
// =========================================================================
TEST(GDMLParser, ValidUnitsStillWork) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <rotation name="r1" x="1.57" y="0" z="0" unit="rad"/>
    <position name="p1" x="10" y="20" z="30" unit="cm"/>
  </define>
  <materials/>
  <solids>
    <box name="box_mm" x="100" y="100" z="100" lunit="mm"/>
    <box name="box_cm" x="10" y="10" z="10" lunit="cm"/>
    <box name="box_m" x="0.1" y="0.1" z="0.1" lunit="m"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="box_mm"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data;
    EXPECT_NO_THROW(data = ParseGDMLString(gdml));
    // All three boxes should have the same size (0.1m half-width = 0.2m full)
    double tol = 1e-9;
    ASSERT_TRUE(data.solids.count("box_mm") > 0);
    ASSERT_TRUE(data.solids.count("box_cm") > 0);
    ASSERT_TRUE(data.solids.count("box_m") > 0);
    EXPECT_NEAR(data.solids["box_mm"]->GetBoundingBox().max_corner.GetX(),
                data.solids["box_cm"]->GetBoundingBox().max_corner.GetX(), tol);
    EXPECT_NEAR(data.solids["box_mm"]->GetBoundingBox().max_corner.GetX(),
                data.solids["box_m"]->GetBoundingBox().max_corner.GetX(), tol);
}


// =========================================================================
// Flexible root tag: <gdml_simple_extension> is accepted
// =========================================================================
TEST(GDMLParser, FlexibleRootTag) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml_simple_extension xmlns:gdml="http://www.example.org">
  <define/>
  <materials/>
  <solids>
    <box name="test_box" x="100" y="200" z="300" lunit="mm"/>
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
</gdml_simple_extension>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("test_box") > 0);
    double tol = 1e-9;
    // 100mm half-width = 0.1m, full width = 0.2m, AABB half = 0.1m
    EXPECT_NEAR(data.solids["test_box"]->GetBoundingBox().max_corner.GetX(), 0.1, tol);
    // Should have a warning about non-standard root tag
    ASSERT_GE(data.warnings.size(), 1u);
    bool found_root_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("non-standard root element") != std::string::npos) found_root_warning = true;
    }
    EXPECT_TRUE(found_root_warning) << "Expected warning about non-standard root element";
}


// =========================================================================
// Multiple <solids> sections: all are parsed
// =========================================================================
TEST(GDMLParser, MultipleSolidsSections) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="box_a" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <solids>
    <box name="box_b" x="200" y="200" z="200" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="box_a"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_EQ(data.solids.count("box_a"), 1u) << "box_a from first section should be parsed";
    EXPECT_EQ(data.solids.count("box_b"), 1u) << "box_b from second section should be parsed";
}


// =========================================================================
// Multiple sections: boolean references across sections
// =========================================================================
TEST(GDMLParser, MultipleSectionsBoolean) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="base_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <solids>
    <box name="cut_box" x="50" y="50" z="200" lunit="mm"/>
    <subtraction name="cut_result">
      <first ref="base_box"/>
      <second ref="cut_box"/>
    </subtraction>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="cut_result"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_EQ(data.solids.count("base_box"), 1u);
    EXPECT_EQ(data.solids.count("cut_box"), 1u);
    EXPECT_EQ(data.solids.count("cut_result"), 1u)
        << "Boolean solid referencing operand from different section should resolve";
}


// =========================================================================
// Multiple <define> sections: constants from all sections available
// =========================================================================
TEST(GDMLParser, MultipleDefineSections) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="half_x" value="50"/>
  </define>
  <define>
    <constant name="half_y" value="100"/>
  </define>
  <materials/>
  <solids>
    <box name="test_box" x="half_x" y="half_y" z="200" lunit="mm"/>
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
    EXPECT_NEAR(data.constants["half_x"], 50.0, 1e-9);
    EXPECT_NEAR(data.constants["half_y"], 100.0, 1e-9);
    ASSERT_TRUE(data.solids.count("test_box") > 0);
}


// =========================================================================
// Duplicate solid name: warns and overwrites
// =========================================================================
TEST(GDMLParser, DuplicateNameWarning) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <solids>
    <box name="mybox" x="200" y="200" z="200" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="mybox"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    // Last wins: mybox should be the 200mm version
    ASSERT_TRUE(data.solids.count("mybox") > 0);
    double tol = 1e-9;
    EXPECT_NEAR(data.solids["mybox"]->GetBoundingBox().max_corner.GetX(), 0.2, tol)
        << "Duplicate name should use last definition (200mm half-width = 0.2m)";
    // Should have a warning about duplicate
    bool found_dup_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("duplicate") != std::string::npos && w.find("mybox") != std::string::npos) {
            found_dup_warning = true;
        }
    }
    EXPECT_TRUE(found_dup_warning) << "Expected warning about duplicate solid name";
}


// =========================================================================
// Duplicate name in strict mode: throws
// =========================================================================
TEST(GDMLParser, DuplicateNameStrictThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <solids>
    <box name="mybox" x="200" y="200" z="200" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="mybox"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_test_dupstrict.gdml";
    { std::ofstream f(tmpfile); f << gdml; }
    GDMLParseOptions strict_opts;
    strict_opts.strict = true;
    EXPECT_THROW(ParseGDML(tmpfile, strict_opts), std::runtime_error);
    std::remove(tmpfile.c_str());
}


// =========================================================================
// Multiple <structure> sections: all volumes parsed
// =========================================================================
TEST(GDMLParser, MultipleStructureSections) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="box_a" x="100" y="100" z="100" lunit="mm"/>
    <box name="box_b" x="200" y="200" z="200" lunit="mm"/>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="VolA">
      <materialref ref=""/>
      <solidref ref="box_a"/>
    </volume>
  </structure>
  <structure>
    <volume name="VolB">
      <materialref ref=""/>
      <solidref ref="box_b"/>
    </volume>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="VolA"/>
      </physvol>
      <physvol>
        <volumeref ref="VolB"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_EQ(data.volumes.count("VolA"), 1u) << "VolA from first structure should be parsed";
    EXPECT_EQ(data.volumes.count("VolB"), 1u) << "VolB from second structure should be parsed";
    EXPECT_EQ(data.volumes.count("World"), 1u) << "World from second structure should be parsed";
}


// =========================================================================
// Torus integration: full-rotation torus parses correctly
// =========================================================================
TEST(GDMLParser, TorusIntegration) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <torus name="my_torus" rmin="10" rmax="20" rtor="100" startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="my_torus"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("my_torus") > 0);
    // Verify bounding box: rtor=100mm=0.1m, rmax=20mm=0.02m
    // xy extent = rtor + rmax = 0.12m, z extent = rmax = 0.02m
    double tol = 1e-6;
    auto bb = data.solids["my_torus"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 0.12, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.02, tol);
    // No warnings for a full-rotation torus
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("torus") == std::string::npos)
            << "Unexpected torus warning: " << w;
    }
}

// =========================================================================
// Torus: partial deltaphi emits warning
// =========================================================================
TEST(GDMLParser, PartialTorusWarning) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <torus name="elbow" rmin="10" rmax="20" rtor="100" startphi="0" deltaphi="90" lunit="mm" aunit="deg"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="elbow"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    // Should still parse (creates full torus) with a warning
    ASSERT_TRUE(data.solids.count("elbow") > 0);
    bool found_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("torus") != std::string::npos && w.find("deltaphi") != std::string::npos) {
            found_warning = true;
        }
    }
    EXPECT_TRUE(found_warning) << "Expected warning about partial torus angular extent";
}
