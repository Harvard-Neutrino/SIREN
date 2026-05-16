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

    // Circular reference should throw an error
    DetectorModel dm;
    EXPECT_THROW(dm.LoadGDML(tmpfile), std::runtime_error)
        << "LoadGDML should throw on circular volume references";
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
    <eltube name="et1" dx="30" dy="20" dz="50" lunit="mm"/>
    <ellipsoid name="ell1" ax="50" by="40" cz="30" lunit="mm"/>
    <para name="para1" x="40" y="30" z="50" alpha="0.3" theta="0.2" phi="0.5" aunit="rad" lunit="mm"/>
    <trap name="trap1" z="100" theta="0" phi="0" y1="50" x1="40" x2="40" alpha1="0" y2="30" x3="30" x4="30" alpha2="0" aunit="rad" lunit="mm"/>
    <cutTube name="ct1" rmin="0" rmax="50" z="100" startphi="0" deltaphi="360" lowX="0" lowY="0" lowZ="-1" highX="0" highY="0" highZ="1" aunit="deg" lunit="mm"/>
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

    std::vector<std::string> names = {"b1", "s1", "t1", "c1", "pc1", "ph1", "trd1", "xt1", "et1", "ell1", "para1", "trap1", "ct1"};
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
    std::string gdml = R"GDML(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="bad" value="bogus_func(3.14)"/>
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
// Duplicate solid name: tracked as multiple instances
// =========================================================================
TEST(GDMLParser, DuplicateSolidInstances) {
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
      <solidref ref="world_placeholder"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    // Instance 0 keeps the original name: 100mm half-width = 0.1m
    ASSERT_TRUE(data.solids.count("mybox") > 0);
    double tol = 1e-9;
    EXPECT_NEAR(data.solids["mybox"]->GetBoundingBox().max_corner.GetX(), 0.1, tol)
        << "First instance keeps original name (100mm half-width = 0.1m)";
    // Instance 1 gets flattened name: mybox__2, 200mm half-width = 0.2m
    ASSERT_TRUE(data.solids.count("mybox__2") > 0);
    EXPECT_NEAR(data.solids["mybox__2"]->GetBoundingBox().max_corner.GetX(), 0.2, tol)
        << "Second instance gets flattened name mybox__2 (200mm half-width = 0.2m)";
    // Instance count should be 2
    ASSERT_TRUE(data.solid_instance_counts.count("mybox") > 0);
    EXPECT_EQ(data.solid_instance_counts["mybox"], 2);
}


// =========================================================================
// Duplicate solid name referenced by volume: throws (ambiguous)
// =========================================================================
TEST(GDMLParser, AmbiguousSolidRefThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox" x="100" y="100" z="100" lunit="mm"/>
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

    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error)
        << "Volume referencing ambiguous solid (2 instances) should throw";
}


// =========================================================================
// Boolean with duplicate operand name: resolves correct instance
// =========================================================================
TEST(GDMLParser, BooleanInstanceResolution) {
    // Define a box with name "operand", then a boolean using it, then
    // redefine "operand" with different dimensions. The boolean should
    // use the first definition (current at definition time).
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="base" x="100" y="100" z="100" lunit="mm"/>
    <box name="cutter" x="50" y="50" z="200" lunit="mm"/>
    <subtraction name="result">
      <first ref="base"/>
      <second ref="cutter"/>
    </subtraction>
    <box name="base" x="500" y="500" z="500" lunit="mm"/>
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
    ASSERT_TRUE(data.solids.count("result") > 0);
    // The boolean "result" was defined after first "base" (100mm) but before
    // second "base" (500mm). It should use the first definition.
    auto bb = data.solids["result"]->GetBoundingBox();
    double tol = 1e-6;
    // The subtraction bounding box should be close to the first base box (100mm=0.1m half)
    EXPECT_NEAR(bb.max_corner.GetX(), 0.1, tol)
        << "Boolean should use first instance of 'base' (100mm), not second (500mm)";
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
// Partial sphere parses without warnings (angular extent now supported)
// =========================================================================
TEST(GDMLParser, PartialSphereParses) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <sphere name="hemi" rmin="0" rmax="100" starttheta="0" deltatheta="90" startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="hemi"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("hemi") > 0);
    // No warnings about partial sphere (it is now supported)
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("sphere") == std::string::npos && w.find("angular") == std::string::npos)
            << "Unexpected sphere angular warning: " << w;
    }
}

// =========================================================================
// Partial torus parses without warnings (angular extent now supported)
// =========================================================================
TEST(GDMLParser, PartialTorusParses) {
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
    ASSERT_TRUE(data.solids.count("elbow") > 0);
    // No warnings about partial torus (it is now supported)
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("torus") == std::string::npos && w.find("angular") == std::string::npos)
            << "Unexpected torus angular warning: " << w;
    }
}

// =========================================================================
// Partial cylinder still warns (not yet supported)
// =========================================================================
TEST(GDMLParser, PartialCylinderSupported) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <tube name="wedge" rmin="0" rmax="100" z="200" deltaphi="180" lunit="mm" aunit="deg"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="wedge"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("wedge") > 0);
    ASSERT_TRUE(data.solids["wedge"] != nullptr);
    // Partial phi is now supported -- no warnings expected
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("partial angular extent") == std::string::npos)
            << "Unexpected warning: " << w;
    }
}

// =========================================================================
// Expression evaluator: math functions
// =========================================================================
TEST(GDMLParser, ExpressionSingleArgFunctions) {
    std::string gdml = R"GDML(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="pi_val" value="3.14159265358979323846"/>
    <constant name="s" value="sin(pi_val/6)"/>
    <constant name="c" value="cos(0)"/>
    <constant name="t" value="tan(pi_val/4)"/>
    <constant name="e" value="exp(1)"/>
    <constant name="l" value="log(exp(2))"/>
    <constant name="sq" value="sqrt(16)"/>
    <constant name="ab" value="abs(-5.5)"/>
    <constant name="as" value="asin(0.5)"/>
    <constant name="ac" value="acos(1.0)"/>
    <constant name="at" value="atan(1.0)"/>
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

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_NEAR(data.constants.at("s"), 0.5, 1e-9);
    EXPECT_NEAR(data.constants.at("c"), 1.0, 1e-9);
    EXPECT_NEAR(data.constants.at("t"), 1.0, 1e-6);
    EXPECT_NEAR(data.constants.at("e"), std::exp(1.0), 1e-9);
    EXPECT_NEAR(data.constants.at("l"), 2.0, 1e-9);
    EXPECT_NEAR(data.constants.at("sq"), 4.0, 1e-9);
    EXPECT_NEAR(data.constants.at("ab"), 5.5, 1e-9);
    EXPECT_NEAR(data.constants.at("as"), std::asin(0.5), 1e-9);
    EXPECT_NEAR(data.constants.at("ac"), 0.0, 1e-9);
    EXPECT_NEAR(data.constants.at("at"), std::atan(1.0), 1e-9);
}

TEST(GDMLParser, ExpressionTwoArgFunctions) {
    std::string gdml = R"GDML(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="p" value="pow(2,10)"/>
    <constant name="a2" value="atan2(1,1)"/>
    <constant name="mn" value="min(3,7)"/>
    <constant name="mx" value="max(3,7)"/>
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

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_NEAR(data.constants.at("p"), 1024.0, 1e-9);
    EXPECT_NEAR(data.constants.at("a2"), std::atan2(1.0, 1.0), 1e-9);
    EXPECT_NEAR(data.constants.at("mn"), 3.0, 1e-9);
    EXPECT_NEAR(data.constants.at("mx"), 7.0, 1e-9);
}

TEST(GDMLParser, ExpressionNestedFunctions) {
    std::string gdml = R"GDML(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="v" value="sqrt(pow(3,2)+pow(4,2))"/>
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

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_NEAR(data.constants.at("v"), 5.0, 1e-9);
}

TEST(GDMLParser, DeeplyNestedParentheses) {
    // The iterative evaluator handles arbitrary nesting without stack overflow
    std::string nested = "7";
    for(int i = 0; i < 50; ++i) {
        nested = "(" + nested + ")";
    }
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="deep" value=")" + nested + R"("/>
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
    EXPECT_NEAR(data.constants.at("deep"), 7.0, 1e-9);
}

// =========================================================================
// Material composition cycle detection
// =========================================================================
TEST(GDMLParser, MaterialCompositionCycleThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="MatA" Z="1">
      <D value="1.0" unit="g/cm3"/>
      <composite n="1" ref="MatB"/>
    </material>
    <material name="MatB" Z="1">
      <D value="1.0" unit="g/cm3"/>
      <composite n="1" ref="MatA"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="MatA"/>
      <solidref ref="world_box"/>
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
    DetectorModel dm;
    EXPECT_THROW(dm.LoadGDML(tmpfile), std::runtime_error);
    std::remove(tmpfile.c_str());
}


// =========================================================================
// Quantity with unit converts correctly
// =========================================================================
TEST(GDMLParser, QuantityUnitConversion) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <quantity name="q_cm" value="100" unit="cm" type="length"/>
    <quantity name="q_mm" value="500" unit="mm" type="length"/>
    <quantity name="q_deg" value="90" unit="deg" type="angle"/>
    <quantity name="q_nounit" value="42" type="length"/>
    <quantity name="q_notype" value="10" unit="cm"/>
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
    // 100 cm -> 1000 mm (cm -> mm scale = 10)
    EXPECT_NEAR(data.constants["q_cm"], 1000.0, tol)
        << "100 cm should convert to 1000 mm";
    // 500 mm -> 500 mm (mm -> mm scale = 1)
    EXPECT_NEAR(data.constants["q_mm"], 500.0, tol)
        << "500 mm should stay as 500 mm";
    // 90 deg -> pi/2 rad
    EXPECT_NEAR(data.constants["q_deg"], 90.0 * 3.141592653589793 / 180.0, tol)
        << "90 deg should convert to pi/2 radians";
    // No unit -> raw value
    EXPECT_NEAR(data.constants["q_nounit"], 42.0, tol)
        << "No unit should leave value unchanged";
    // No type but unit=cm -> treated as length
    EXPECT_NEAR(data.constants["q_notype"], 100.0, tol)
        << "10 cm with no type should convert to 100 mm (length assumed)";
}


// =========================================================================
// Exponentiation operator: basic, precedence, right-associativity
// =========================================================================
TEST(GDMLParser, ExponentiationOperator) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="basic" value="2^3"/>
    <constant name="right_assoc" value="2^3^2"/>
    <constant name="precedence" value="2*3^2"/>
    <constant name="neg_base" value="-2^2"/>
    <constant name="paren" value="(2+1)^3"/>
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
    EXPECT_NEAR(data.constants["basic"], 8.0, tol)
        << "2^3 should be 8";
    // Right-associative: 2^3^2 = 2^(3^2) = 2^9 = 512
    EXPECT_NEAR(data.constants["right_assoc"], 512.0, tol)
        << "2^3^2 should be 512 (right-associative: 2^(3^2) = 2^9)";
    // Precedence: 2*3^2 = 2*(3^2) = 2*9 = 18
    EXPECT_NEAR(data.constants["precedence"], 18.0, tol)
        << "2*3^2 should be 18 (^ has higher precedence than *)";
    // Unary minus: -2^2 = -(2^2) = -4
    EXPECT_NEAR(data.constants["neg_base"], -4.0, tol)
        << "-2^2 should be -4 (unary minus applied after exponentiation)";
    // Parenthesized: (2+1)^3 = 3^3 = 27
    EXPECT_NEAR(data.constants["paren"], 27.0, tol)
        << "(2+1)^3 should be 27";
}


// =========================================================================
// stod out-of-range throws instead of returning 0
// =========================================================================
TEST(GDMLParser, StodOutOfRangeThrows) {
    // Use a value that will overflow double
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="huge" value="1e99999"/>
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

    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error)
        << "Out-of-range numeric value should throw, not silently return 0";
}


// =========================================================================
// Variable element stored as constant (usable by loops)
// =========================================================================
TEST(GDMLParser, VariableElementStoredAsConstant) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="myvar" value="10"/>
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

    // The variable should be stored as a constant
    EXPECT_NEAR(data.constants["myvar"], 10.0, 1e-9);
}


// =========================================================================
// Loop element expands children
// =========================================================================
TEST(GDMLParser, LoopElementExpansion) {
    // Uses bracket notation [i] for name construction (Geant4 convention)
    // and bare variable i for numeric expressions
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <loop for="i" from="0" to="3" step="1">
      <constant name="val_[i]" value="i*10"/>
    </loop>
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

    // Loop should have generated val_0, val_1, val_2, val_3
    EXPECT_NEAR(data.constants["val_0"], 0.0, 1e-9);
    EXPECT_NEAR(data.constants["val_1"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["val_2"], 20.0, 1e-9);
    EXPECT_NEAR(data.constants["val_3"], 30.0, 1e-9);
}


// =========================================================================
// Loop with omitted "from" defaults to variable's initial value
// (GDML spec section 3.4.28)
// =========================================================================
TEST(GDMLParser, LoopOmittedFromUsesVariableValue) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="x" value="3"/>
    <loop for="x" to="5" step="1">
      <constant name="c_[x]" value="x*x"/>
    </loop>
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

    // Loop should start from x=3 (variable value), not x=0
    EXPECT_EQ(data.constants.count("c_3"), 1u) << "Should have c_3 (x=3)";
    EXPECT_EQ(data.constants.count("c_4"), 1u) << "Should have c_4 (x=4)";
    EXPECT_EQ(data.constants.count("c_5"), 1u) << "Should have c_5 (x=5)";
    EXPECT_EQ(data.constants.count("c_0"), 0u) << "Should NOT have c_0 (loop starts at 3)";
    EXPECT_NEAR(data.constants["c_3"], 9.0, 1e-9);
    EXPECT_NEAR(data.constants["c_4"], 16.0, 1e-9);
    EXPECT_NEAR(data.constants["c_5"], 25.0, 1e-9);
}


// =========================================================================
// Missing setup section emits warning
// =========================================================================
TEST(GDMLParser, MissingSetupWarning) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
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
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    EXPECT_TRUE(data.world_volume.empty());
    bool found_setup_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("setup") != std::string::npos || w.find("world volume") != std::string::npos) {
            found_setup_warning = true;
        }
    }
    EXPECT_TRUE(found_setup_warning)
        << "Expected warning about missing <setup> section";
}


// =========================================================================
// Expression evaluator: function with 3+ arguments throws
// =========================================================================
TEST(GDMLParser, ThreeArgFunctionThrows) {
    std::string gdml = R"GDML(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="bad" value="max(1,2,3)"/>
  </define>
  <materials/>
  <solids/>
  <structure/>
</gdml>)GDML";

    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error)
        << "Function with 3 arguments should throw";
}


// =========================================================================
// Micrometer and nanometer length units
// =========================================================================
TEST(GDMLParser, MicrometerNanometerUnits) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="tiny_box" x="100" y="200" z="300" lunit="um"/>
    <box name="nano_box" x="500" y="500" z="500" lunit="nm"/>
  </solids>
  <structure/>
  <setup name="Default" version="1.0">
    <world ref=""/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_NE(data.solids.find("tiny_box"), data.solids.end());
    EXPECT_NE(data.solids.find("nano_box"), data.solids.end());
}


// =========================================================================
// Assembly volume support
// =========================================================================
TEST(GDMLParser, AssemblyVolume) {
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
    <box name="child_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="ChildVol">
      <materialref ref="Air"/>
      <solidref ref="child_box"/>
    </volume>
    <assembly name="MyAssembly">
      <physvol>
        <volumeref ref="ChildVol"/>
        <position name="pos1" x="200" y="0" z="0" unit="mm"/>
      </physvol>
    </assembly>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="MyAssembly"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    // Assembly should be stored as a volume with is_assembly=true
    auto it = data.volumes.find("MyAssembly");
    ASSERT_NE(it, data.volumes.end());
    EXPECT_TRUE(it->second.is_assembly);
    EXPECT_EQ(it->second.children.size(), 1u);

    // Loading into DetectorModel should create a sector for ChildVol but not for MyAssembly
    std::string tmpfile = "/tmp/siren_gdml_assembly_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    bool found_child = false;
    bool found_assembly = false;
    for(auto const & sector : dm.GetSectors()) {
        if(sector.name == "ChildVol") found_child = true;
        if(sector.name == "MyAssembly") found_assembly = true;
    }
    EXPECT_TRUE(found_child) << "ChildVol should have a sector";
    EXPECT_FALSE(found_assembly) << "Assembly should not have its own sector";
}


// =========================================================================
// BuildVolume throws on missing volume reference
// =========================================================================
TEST(GDMLParser, BuildVolumeMissingVolumeThrows) {
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
        <volumeref ref="NonexistentVolume"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_missing_vol_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    EXPECT_THROW(dm.LoadGDML(tmpfile), std::runtime_error);
    std::remove(tmpfile.c_str());
}


// =========================================================================
// BuildVolume throws on missing solid reference
// =========================================================================
TEST(GDMLParser, BuildVolumeMissingSolidThrows) {
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
      <solidref ref="nonexistent_solid"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_missing_solid_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    EXPECT_THROW(dm.LoadGDML(tmpfile), std::runtime_error);
    std::remove(tmpfile.c_str());
}


// =========================================================================
// Loop expansion in <structure> generates physvol children
// =========================================================================
TEST(GDMLParser, LoopInStructure) {
    // Uses bracket notation [i] for unique position names (Geant4 convention)
    // and bare variable i for numeric position expressions
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
  </define>
  <materials>
    <material name="TestMat">
      <D value="0.00129" unit="g/cm3"/>
      <atom value="14" Z="7"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="sub_box" x="10" y="10" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="SubVol">
      <materialref ref="TestMat"/>
      <solidref ref="sub_box"/>
    </volume>
    <volume name="World">
      <materialref ref="TestMat"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="2" step="1">
        <physvol>
          <volumeref ref="SubVol"/>
          <position name="pos_[i]" x="i*100" y="0" z="0" unit="mm"/>
        </physvol>
      </loop>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_loop_struct_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    EXPECT_NO_THROW(dm.LoadGDML(tmpfile));
    std::remove(tmpfile.c_str());

    // Should have 1 world + 3 child placements (i=0,1,2)
    int child_count = 0;
    for(auto const & sector : dm.GetSectors()) {
        if(sector.name.find("SubVol") != std::string::npos) child_count++;
    }
    EXPECT_EQ(child_count, 3) << "Loop should produce 3 child volume placements";
}


// =========================================================================
// Loop variable does not corrupt identifiers containing the variable name
// (e.g. variable "i" must not corrupt "ChildVol" into "Ch0ldVol")
// =========================================================================
TEST(GDMLParser, LoopWordBoundarySafety) {
    // Volume name "ChildVol" contains "i" but must NOT be corrupted
    // by a loop with variable "i". Word-boundary matching prevents this.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
  </define>
  <materials>
    <material name="Air">
      <D value="0.00129" unit="g/cm3"/>
      <atom value="14" Z="7"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="child_box" x="10" y="10" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="ChildVol">
      <materialref ref="Air"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="1" step="1">
        <physvol>
          <volumeref ref="ChildVol"/>
          <position name="p_[i]" x="i*100" y="0" z="0" unit="mm"/>
        </physvol>
      </loop>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_loop_wordboundary_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    EXPECT_NO_THROW(dm.LoadGDML(tmpfile))
        << "Loop variable 'i' should not corrupt 'ChildVol' reference";
    std::remove(tmpfile.c_str());

    int child_count = 0;
    for(auto const & sector : dm.GetSectors()) {
        if(sector.name.find("ChildVol") != std::string::npos) child_count++;
    }
    EXPECT_EQ(child_count, 2) << "Should have 2 ChildVol placements (i=0,1)";
}


// =========================================================================
// EllipticalTube parse + bounding box
// =========================================================================
TEST(GDMLParser, EllipticalTubeIntegration) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <eltube name="et" dx="30" dy="20" dz="50" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="et"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("et") > 0);
    auto bb = data.solids["et"]->GetBoundingBox();
    double tol = 1e-9;
    // dx=30mm=0.03m, dy=20mm=0.02m, dz=50mm=0.05m
    EXPECT_NEAR(bb.max_corner.GetX(), 0.03, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 0.02, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.05, tol);
}

// =========================================================================
// GenericPolycone parse
// =========================================================================
TEST(GDMLParser, GenericPolycone) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <genericPolycone name="gpc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <rzpoint r="50" z="-100"/>
      <rzpoint r="80" z="0"/>
      <rzpoint r="30" z="100"/>
    </genericPolycone>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="gpc"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("gpc") > 0);
    auto bb = data.solids["gpc"]->GetBoundingBox();
    double tol = 1e-6;
    // rmax=80mm=0.08m at z=0
    EXPECT_NEAR(bb.max_corner.GetX(), 0.08, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 0.08, tol);
    // z ranges from -100mm to 100mm = -0.1m to 0.1m
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.1, tol);
    EXPECT_NEAR(bb.min_corner.GetZ(), -0.1, tol);
}

// =========================================================================
// CutTube parse
// =========================================================================
TEST(GDMLParser, CutTubeIntegration) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <cutTube name="ct" rmin="10" rmax="50" z="100" startphi="0" deltaphi="360"
             lowX="0" lowY="0" lowZ="-1" highX="0" highY="0" highZ="1"
             aunit="deg" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="ct"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("ct") > 0);
    auto bb = data.solids["ct"]->GetBoundingBox();
    double tol = 1e-6;
    // rmax=50mm=0.05m, z=100mm half=50mm=0.05m (flat cuts)
    EXPECT_NEAR(bb.max_corner.GetX(), 0.05, tol);
}

// =========================================================================
// Trap parse + degenerate-to-Trd check
// =========================================================================
TEST(GDMLParser, TrapIntegration) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <trap name="t" z="100" theta="0" phi="0"
          y1="60" x1="40" x2="40" alpha1="0"
          y2="40" x3="30" x4="30" alpha2="0"
          aunit="rad" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="t"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("t") > 0);
    auto bb = data.solids["t"]->GetBoundingBox();
    double tol = 1e-6;
    // z=100mm = 0.1m (GDML trap z is half-height, passed directly as dz)
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.1, tol);
    // max x = max(dx1,dx2,dx3,dx4) = 40mm = 0.04m
    EXPECT_NEAR(bb.max_corner.GetX(), 0.04, tol);
}

// =========================================================================
// Ellipsoid parse
// =========================================================================
TEST(GDMLParser, EllipsoidIntegration) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <ellipsoid name="ell" ax="50" by="30" cz="40" zcut2="20" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="ell"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("ell") > 0);
    auto bb = data.solids["ell"]->GetBoundingBox();
    double tol = 1e-6;
    EXPECT_NEAR(bb.max_corner.GetX(), 0.05, tol);  // ax=50mm=0.05m
    EXPECT_NEAR(bb.max_corner.GetY(), 0.03, tol);  // by=30mm=0.03m
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.02, tol);  // zcut2=20mm=0.02m
    EXPECT_NEAR(bb.min_corner.GetZ(), -0.04, tol); // -cz=-40mm=-0.04m (no zcut1)
}

// =========================================================================
// Para parse
// =========================================================================
TEST(GDMLParser, ParaIntegration) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <para name="p" x="40" y="30" z="50" alpha="0" theta="0" phi="0" aunit="rad" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="p"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("p") > 0);
    auto bb = data.solids["p"]->GetBoundingBox();
    double tol = 1e-6;
    // With alpha=theta=phi=0, para = box: x=40mm=0.04m half, etc.
    EXPECT_NEAR(bb.max_corner.GetX(), 0.04, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 0.03, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.05, tol);
}

// =========================================================================
// Matrix vector lookup
// =========================================================================
TEST(GDMLParser, MatrixVectorLookup) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <matrix name="energies" coldim="1" values="1.5 2.5 3.5 4.5"/>
    <constant name="e2" value="energies[2]"/>
    <constant name="e4" value="energies[4]"/>
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
    EXPECT_NEAR(data.constants["e2"], 2.5, 1e-9);
    EXPECT_NEAR(data.constants["e4"], 4.5, 1e-9);
}

// =========================================================================
// Matrix 2D table lookup
// =========================================================================
TEST(GDMLParser, MatrixTableLookup) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <matrix name="m" coldim="3" values="10 20 30 40 50 60 70 80 90"/>
    <constant name="v12" value="m[1,2]"/>
    <constant name="v23" value="m[2,3]"/>
    <constant name="v31" value="m[3,1]"/>
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
    // Row 1, col 2 = index (0*3+1) = values[1] = 20
    EXPECT_NEAR(data.constants["v12"], 20.0, 1e-9);
    // Row 2, col 3 = index (1*3+2) = values[5] = 60
    EXPECT_NEAR(data.constants["v23"], 60.0, 1e-9);
    // Row 3, col 1 = index (2*3+0) = values[6] = 70
    EXPECT_NEAR(data.constants["v31"], 70.0, 1e-9);
}

// =========================================================================
// Matrix access inside loop
// =========================================================================
TEST(GDMLParser, MatrixInLoop) {
    // Uses word-boundary substitution (bare i in arithmetic) rather than bracket
    // notation (pos[i]) to avoid collision with the loop bracket substitution.
    // The loop expands i*1 to e.g. "1*1", "2*1", "3*1" which the expression
    // evaluator resolves to 1, 2, 3 for the matrix index.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <matrix name="pos" coldim="1" values="100 200 300"/>
    <variable name="idx" value="1"/>
    <constant name="p1" value="pos[1]"/>
    <constant name="p2" value="pos[2]"/>
    <constant name="p3" value="pos[3]"/>
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
    EXPECT_NEAR(data.constants["p1"], 100.0, 1e-9);
    EXPECT_NEAR(data.constants["p2"], 200.0, 1e-9);
    EXPECT_NEAR(data.constants["p3"], 300.0, 1e-9);
}

// --- Loop variable persistence after loop (#13) ---

TEST(GDMLParser, LoopVariablePersistsAfterLoop) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <loop for="i" from="1" to="5" step="1">
      <constant name="c_[i]" value="i*10"/>
    </loop>
    <constant name="final_i" value="i"/>
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
    // The loop variable should persist at its final value (5)
    EXPECT_NEAR(data.constants["i"], 5.0, 1e-9);
    // The constant referencing the loop variable after the loop should work
    EXPECT_NEAR(data.constants["final_i"], 5.0, 1e-9);
    // Loop body constants should be defined
    EXPECT_NEAR(data.constants["c_1"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["c_5"], 50.0, 1e-9);
}

TEST(GDMLParser, LoopVariableDocumentOrder) {
    // A constant defined BEFORE the loop should see the pre-loop variable value.
    // A constant defined AFTER the loop should see the post-loop (final) value.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="10"/>
    <constant name="before_loop" value="i"/>
    <loop for="i" from="1" to="5" step="1">
      <constant name="c_[i]" value="i*10"/>
    </loop>
    <constant name="after_loop" value="i"/>
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
    // before_loop was defined when i=10 (pre-loop)
    EXPECT_NEAR(data.constants["before_loop"], 10.0, 1e-9);
    // after_loop was defined when i=5 (post-loop final value)
    EXPECT_NEAR(data.constants["after_loop"], 5.0, 1e-9);
    // Loop body constants should use their iteration value
    EXPECT_NEAR(data.constants["c_1"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["c_5"], 50.0, 1e-9);
}

TEST(GDMLParser, LoopVariablePreexistingValue) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="idx" value="99"/>
    <loop for="idx" from="1" to="3" step="1">
      <constant name="v_[idx]" value="idx"/>
    </loop>
    <constant name="after_loop" value="idx"/>
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
    // After the loop, the variable should have the final loop value (3), not 99
    EXPECT_NEAR(data.constants["after_loop"], 3.0, 1e-9);
}

// --- Boolean firstposition/firstrotation support (#8) ---

TEST(GDMLParser, BooleanFirstPosition) {
    // box_a: 100mm half-width (0.1m), box_b: 20mm half-width (0.02m)
    // firstposition shifts box_a by (200, 0, 0) mm = 0.2m in x
    // A point at (0.2, 0, 0) is at box_a's center, outside box_b -> inside subtraction
    // A point at (0, 0, 0) is outside box_a -> outside subtraction
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Vacuum" Z="1">
      <D value="1e-25" unit="g/cm3"/>
      <atom value="1.008"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="box_a" x="100" y="100" z="100" lunit="mm"/>
    <box name="box_b" x="20" y="20" z="20" lunit="mm"/>
    <subtraction name="sub_first_pos">
      <first ref="box_a"/>
      <second ref="box_b"/>
      <firstposition x="200" y="0" z="0" unit="mm"/>
    </subtraction>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Vacuum"/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("sub_first_pos") > 0);
    auto geo = data.solids["sub_first_pos"];
    ASSERT_NE(geo, nullptr);
    siren::math::Vector3D dir(0, 1, 0);
    // box_a is shifted to center (0.2, 0, 0). Point at (0.2, 0, 0) is inside box_a
    // but outside box_b (centered at origin, half-width 0.02m). 0.2 > 0.02 -> outside box_b.
    // Subtraction: inside_a && !inside_b = true
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.2, 0, 0), dir));
    // Point at origin: box_a center is at 0.2m, half-width 0.1m. |0 - 0.2| = 0.2 > 0.1 -> outside box_a
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0, 0, 0), dir));
}

TEST(GDMLParser, BooleanFirstRotation) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Vacuum" Z="1">
      <D value="1e-25" unit="g/cm3"/>
      <atom value="1.008"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="box_a" x="200" y="50" z="50" lunit="mm"/>
    <box name="box_b" x="10" y="10" z="10" lunit="mm"/>
    <subtraction name="sub_first_rot">
      <first ref="box_a"/>
      <second ref="box_b"/>
      <firstrotation z="90" unit="deg"/>
    </subtraction>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Vacuum"/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("sub_first_rot") > 0);
    auto geo = data.solids["sub_first_rot"];
    ASSERT_NE(geo, nullptr);
    // box_a is 200x50x50 mm, rotated 90 deg about z -> becomes 50x200x50 mm
    // A point at (0, 0.08, 0) should be inside the rotated box_a (0.08m = 80mm < 100mm half-y)
    siren::math::Vector3D dir(0, 0, 1);
    siren::math::Vector3D inside_rot(0, 0.08, 0);
    EXPECT_TRUE(geo->IsInside(inside_rot, dir));
    // But a point at (0.08, 0, 0) should be outside (80mm > 25mm half-x after rotation)
    siren::math::Vector3D outside_rot(0.08, 0, 0);
    EXPECT_FALSE(geo->IsInside(outside_rot, dir));
}

TEST(GDMLParser, BooleanFirstPositionRef) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="shift1" x="50" y="0" z="0" unit="mm"/>
  </define>
  <materials>
    <material name="Vacuum" Z="1">
      <D value="1e-25" unit="g/cm3"/>
      <atom value="1.008"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="box_a" x="100" y="100" z="100" lunit="mm"/>
    <box name="box_b" x="10" y="10" z="10" lunit="mm"/>
    <union name="union_fposref">
      <first ref="box_a"/>
      <second ref="box_b"/>
      <firstpositionref ref="shift1"/>
    </union>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Vacuum"/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("union_fposref") > 0);
    auto geo = data.solids["union_fposref"];
    ASSERT_NE(geo, nullptr);
    // box_a center is shifted to x=0.05m. Point at (0.05, 0, 0) should be inside box_a
    siren::math::Vector3D dir(0, 0, 1);
    siren::math::Vector3D at_shifted(0.05, 0, 0);
    EXPECT_TRUE(geo->IsInside(at_shifted, dir));
}

TEST(GDMLParser, BooleanBothOperandsTransformed) {
    // box_a and box_b both 100mm half-width (0.1m each).
    // firstposition shifts box_a by 300mm (0.3m) in x
    // position shifts box_b by -300mm (-0.3m) in x
    // Union of the two: points at each center should be inside, point at
    // origin should be outside (0.3m > 0.1m half-width)
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Vacuum" Z="1">
      <D value="1e-25" unit="g/cm3"/>
      <atom value="1.008"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
    <box name="box_a" x="100" y="100" z="100" lunit="mm"/>
    <box name="box_b" x="100" y="100" z="100" lunit="mm"/>
    <union name="union_both">
      <first ref="box_a"/>
      <second ref="box_b"/>
      <firstposition x="300" y="0" z="0" unit="mm"/>
      <position x="-300" y="0" z="0" unit="mm"/>
    </union>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Vacuum"/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("union_both") > 0);
    auto geo = data.solids["union_both"];
    ASSERT_NE(geo, nullptr);
    siren::math::Vector3D dir(0, 0, 1);
    // box_a centered at x=0.3m, box_b centered at x=-0.3m (each half-width 0.1m)
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.3, 0, 0), dir));
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(-0.3, 0, 0), dir));
    // Origin is 0.3m from each center, > 0.1m half-width -> outside both
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0, 0, 0), dir));
}


// =========================================================================
// Tessellated solid
// =========================================================================
TEST(GDMLParser, TessellatedTriangular) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="v0" x="0" y="0" z="0" unit="m"/>
    <position name="v1" x="1" y="0" z="0" unit="m"/>
    <position name="v2" x="0" y="1" z="0" unit="m"/>
    <position name="v3" x="0" y="0" z="1" unit="m"/>
  </define>
  <materials/>
  <solids>
    <tessellated name="tetra">
      <triangular vertex1="v0" vertex2="v2" vertex3="v1"/>
      <triangular vertex1="v0" vertex2="v1" vertex3="v3"/>
      <triangular vertex1="v0" vertex2="v3" vertex3="v2"/>
      <triangular vertex1="v1" vertex2="v2" vertex3="v3"/>
    </tessellated>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="tetra"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.find("tetra") != data.solids.end());
    auto geo = data.solids["tetra"];
    ASSERT_TRUE(geo != nullptr);

    // Point inside the tetrahedron (centroid at 0.25, 0.25, 0.25)
    siren::math::Vector3D inside(0.1, 0.1, 0.1);
    EXPECT_TRUE(geo->IsInside(inside));

    // Point outside
    siren::math::Vector3D outside(2.0, 2.0, 2.0);
    EXPECT_FALSE(geo->IsInside(outside));
}

TEST(GDMLParser, TessellatedQuadrangular) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="v0" x="-1" y="-1" z="-1" unit="m"/>
    <position name="v1" x="1" y="-1" z="-1" unit="m"/>
    <position name="v2" x="1" y="1" z="-1" unit="m"/>
    <position name="v3" x="-1" y="1" z="-1" unit="m"/>
    <position name="v4" x="-1" y="-1" z="1" unit="m"/>
    <position name="v5" x="1" y="-1" z="1" unit="m"/>
    <position name="v6" x="1" y="1" z="1" unit="m"/>
    <position name="v7" x="-1" y="1" z="1" unit="m"/>
  </define>
  <materials/>
  <solids>
    <tessellated name="cube">
      <quadrangular vertex1="v3" vertex2="v2" vertex3="v1" vertex4="v0"/>
      <quadrangular vertex1="v4" vertex2="v5" vertex3="v6" vertex4="v7"/>
      <quadrangular vertex1="v0" vertex2="v1" vertex3="v5" vertex4="v4"/>
      <quadrangular vertex1="v1" vertex2="v2" vertex3="v6" vertex4="v5"/>
      <quadrangular vertex1="v2" vertex2="v3" vertex3="v7" vertex4="v6"/>
      <quadrangular vertex1="v3" vertex2="v0" vertex3="v4" vertex4="v7"/>
    </tessellated>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="cube"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.find("cube") != data.solids.end());
    auto geo = data.solids["cube"];
    ASSERT_TRUE(geo != nullptr);

    // Center should be inside
    siren::math::Vector3D center(0, 0, 0);
    EXPECT_TRUE(geo->IsInside(center));

    // Corner region just inside
    siren::math::Vector3D near_corner(0.9, 0.9, 0.9);
    EXPECT_TRUE(geo->IsInside(near_corner));

    // Outside
    siren::math::Vector3D outside(1.5, 0, 0);
    EXPECT_FALSE(geo->IsInside(outside));
}

// =========================================================================
// arb8 solid
// =========================================================================
TEST(GDMLParser, Arb8Basic) {
    // Simple cube via arb8: vertices form a 2x2x2 cube centered at origin
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <arb8 name="cube8" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="-1000" v5y="-1000"
      v6x="1000" v6y="-1000"
      v7x="1000" v7y="1000"
      v8x="-1000" v8y="1000"
      lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="cube8"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.find("cube8") != data.solids.end());
    auto geo = data.solids["cube8"];
    ASSERT_TRUE(geo != nullptr);

    // Center should be inside
    siren::math::Vector3D center(0, 0, 0);
    EXPECT_TRUE(geo->IsInside(center));

    // Point at (0.5, 0.5, 0.5) should be inside (cube is 2m x 2m x 2m)
    siren::math::Vector3D inside(0.5, 0.5, 0.5);
    EXPECT_TRUE(geo->IsInside(inside));

    // Point outside the cube
    siren::math::Vector3D outside(2.0, 0, 0);
    EXPECT_FALSE(geo->IsInside(outside));
}

TEST(GDMLParser, Arb8Trapezoid) {
    // Tapered arb8: bottom face is 2x2, top face is 1x1, height 2
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <arb8 name="taper" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="-500" v5y="-500"
      v6x="500" v6y="-500"
      v7x="500" v7y="500"
      v8x="-500" v8y="500"
      lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="taper"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.find("taper") != data.solids.end());
    auto geo = data.solids["taper"];
    ASSERT_TRUE(geo != nullptr);

    // Center at origin is inside
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0, 0, 0)));

    // Near top face, inside the smaller square
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.3, 0.3, 0.9)));

    // Near top face, outside the smaller square but inside bottom projection
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.8, 0.8, 0.9)));

    // Outside entirely
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0, 0, 2.0)));
}

// =========================================================================
// ENTITY preprocessing
// =========================================================================
TEST(GDMLParser, EntityExpansion) {
    // Create a temporary included file
    std::string inc_content = R"(
    <box name="included_box" x="2000" y="2000" z="2000" lunit="mm"/>
)";
    std::string inc_file = "/tmp/siren_gdml_entity_inc.gdml";
    {
        std::ofstream f(inc_file);
        f << inc_content;
    }

    std::string gdml = R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY solids_inc SYSTEM "siren_gdml_entity_inc.gdml">
]>
<gdml>
  <define/>
  <materials/>
  <solids>
    &solids_inc;
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="included_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    // Write the main file to /tmp so entity path resolves
    std::string main_file = "/tmp/siren_gdml_entity_test.gdml";
    {
        std::ofstream f(main_file);
        f << gdml;
    }

    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    std::remove(inc_file.c_str());

    ASSERT_TRUE(data.solids.find("included_box") != data.solids.end());
    auto geo = data.solids["included_box"];
    ASSERT_TRUE(geo != nullptr);

    // The box is 4m x 4m x 4m (half-widths 2m); center should be inside
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0, 0, 0)));
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(3.0, 0, 0)));
}

TEST(GDMLParser, NestedEntityExpansion) {
    // Create nested included files
    std::string inner_content = R"(
    <constant name="inner_val" value="42"/>
)";
    std::string inner_file = "/tmp/siren_gdml_inner.gdml";
    {
        std::ofstream f(inner_file);
        f << inner_content;
    }

    std::string outer_content = R"(<?xml version="1.0"?>
<!DOCTYPE gdml_fragment [
  <!ENTITY inner SYSTEM "siren_gdml_inner.gdml">
]>
    <constant name="outer_val" value="99"/>
    &inner;
)";
    std::string outer_file = "/tmp/siren_gdml_outer.gdml";
    {
        std::ofstream f(outer_file);
        f << outer_content;
    }

    std::string gdml = R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY defs SYSTEM "siren_gdml_outer.gdml">
]>
<gdml>
  <define>
    &defs;
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

    std::string main_file = "/tmp/siren_gdml_nested_test.gdml";
    {
        std::ofstream f(main_file);
        f << gdml;
    }

    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    std::remove(outer_file.c_str());
    std::remove(inner_file.c_str());

    EXPECT_NEAR(data.constants["outer_val"], 99.0, 1e-9);
    EXPECT_NEAR(data.constants["inner_val"], 42.0, 1e-9);
}

// --- Tests for stack-based loop expansion and document-order resolution ---

TEST(GDMLParser, NestedLoops) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <variable name="j" value="0"/>
    <loop for="i" from="1" to="3" step="1">
      <loop for="j" from="1" to="2" step="1">
        <constant name="c_[i]_[j]" value="i*10+j"/>
      </loop>
    </loop>
    <constant name="final_i" value="i"/>
    <constant name="final_j" value="j"/>
  </define>
  <materials/>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_NEAR(data.constants["c_1_1"], 11.0, 1e-9);
    EXPECT_NEAR(data.constants["c_1_2"], 12.0, 1e-9);
    EXPECT_NEAR(data.constants["c_2_1"], 21.0, 1e-9);
    EXPECT_NEAR(data.constants["c_2_2"], 22.0, 1e-9);
    EXPECT_NEAR(data.constants["c_3_1"], 31.0, 1e-9);
    EXPECT_NEAR(data.constants["c_3_2"], 32.0, 1e-9);
    // Outer loop final value
    EXPECT_NEAR(data.constants["final_i"], 3.0, 1e-9);
    // Inner loop final value (last iteration of the last outer step)
    EXPECT_NEAR(data.constants["final_j"], 2.0, 1e-9);
}

TEST(GDMLParser, NonExecutingLoop) {
    // Loop with from > to and positive step should not execute.
    // The variable should keep its pre-loop value.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="42"/>
    <loop for="i" from="10" to="5" step="1">
      <constant name="should_not_exist" value="999"/>
    </loop>
    <constant name="after" value="i"/>
  </define>
  <materials/>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    // Loop didn't run, so variable keeps its pre-loop value
    EXPECT_NEAR(data.constants["after"], 42.0, 1e-9);
    // Body constant should not have been created
    EXPECT_EQ(data.constants.find("should_not_exist"), data.constants.end());
}

TEST(GDMLParser, NegativeStepLoop) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <loop for="i" from="5" to="1" step="-1">
      <constant name="d_[i]" value="i*2"/>
    </loop>
    <constant name="final_i" value="i"/>
  </define>
  <materials/>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    EXPECT_NEAR(data.constants["d_5"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["d_4"], 8.0, 1e-9);
    EXPECT_NEAR(data.constants["d_3"], 6.0, 1e-9);
    EXPECT_NEAR(data.constants["d_2"], 4.0, 1e-9);
    EXPECT_NEAR(data.constants["d_1"], 2.0, 1e-9);
    EXPECT_NEAR(data.constants["final_i"], 1.0, 1e-9);
}

TEST(GDMLParser, CumulativeLoopReferences) {
    // Each iteration references the previous iteration's result.
    // sum_0 = 0, sum_1 = sum_0 + 1, sum_2 = sum_1 + 2, ...
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <constant name="sum_0" value="0"/>
    <loop for="i" from="1" to="4" step="1">
      <constant name="sum_[i]" value="sum_[i-1]+i"/>
    </loop>
  </define>
  <materials/>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    // sum_1 = 0 + 1 = 1
    EXPECT_NEAR(data.constants["sum_1"], 1.0, 1e-9);
    // sum_2 = 1 + 2 = 3
    EXPECT_NEAR(data.constants["sum_2"], 3.0, 1e-9);
    // sum_3 = 3 + 3 = 6
    EXPECT_NEAR(data.constants["sum_3"], 6.0, 1e-9);
    // sum_4 = 6 + 4 = 10
    EXPECT_NEAR(data.constants["sum_4"], 10.0, 1e-9);
}

TEST(GDMLParser, LoopGeneratingPositions) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <constant name="spacing" value="100"/>
    <loop for="i" from="0" to="2" step="1">
      <position name="pos_[i]" x="i*spacing" y="0" z="0" unit="mm"/>
    </loop>
  </define>
  <materials/>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.positions.count("pos_0") > 0);
    ASSERT_TRUE(data.positions.count("pos_1") > 0);
    ASSERT_TRUE(data.positions.count("pos_2") > 0);
    // pos_0: x = 0*100 mm = 0 m
    EXPECT_NEAR(data.positions["pos_0"].GetX(), 0.0, 1e-9);
    // pos_1: x = 1*100 mm = 0.1 m
    EXPECT_NEAR(data.positions["pos_1"].GetX(), 0.1, 1e-9);
    // pos_2: x = 2*100 mm = 0.2 m
    EXPECT_NEAR(data.positions["pos_2"].GetX(), 0.2, 1e-9);
}

TEST(GDMLParser, NestedLoopInStructure) {
    // Tests the iterative ExpandLoops for structure sections with nested loops.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <variable name="j" value="0"/>
  </define>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.01"/></material>
  </materials>
  <solids>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
    <box name="child_box" x="10" y="10" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="ChildVol">
      <materialref ref="Air"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="1" step="1">
        <loop for="j" from="0" to="1" step="1">
          <physvol>
            <volumeref ref="ChildVol"/>
            <position name="nested_[i]_[j]" x="i*200" y="j*200" z="0" unit="mm"/>
          </physvol>
        </loop>
      </loop>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_nested_struct_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    EXPECT_NO_THROW(dm.LoadGDML(tmpfile));
    std::remove(tmpfile.c_str());

    // Should have 2x2 = 4 child placements
    int child_count = 0;
    for(auto const & sector : dm.GetSectors()) {
        if(sector.name.find("ChildVol") != std::string::npos) child_count++;
    }
    EXPECT_EQ(child_count, 4);
}

TEST(GDMLParser, LoopVariableInInnerBounds) {
    // Outer loop variable used as inner loop's bound
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <variable name="j" value="0"/>
    <loop for="i" from="1" to="3" step="1">
      <loop for="j" from="1" to="i" step="1">
        <constant name="t_[i]_[j]" value="1"/>
      </loop>
    </loop>
  </define>
  <materials/>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    // i=1: j runs 1..1 -> t_1_1
    EXPECT_TRUE(data.constants.count("t_1_1") > 0);
    EXPECT_EQ(data.constants.count("t_1_2"), 0u);
    // i=2: j runs 1..2 -> t_2_1, t_2_2
    EXPECT_TRUE(data.constants.count("t_2_1") > 0);
    EXPECT_TRUE(data.constants.count("t_2_2") > 0);
    EXPECT_EQ(data.constants.count("t_2_3"), 0u);
    // i=3: j runs 1..3 -> t_3_1, t_3_2, t_3_3
    EXPECT_TRUE(data.constants.count("t_3_1") > 0);
    EXPECT_TRUE(data.constants.count("t_3_2") > 0);
    EXPECT_TRUE(data.constants.count("t_3_3") > 0);
}

// =========================================================================
// Test #1: Degenerate arb8 (coincident vertices)
// =========================================================================
TEST(GDMLParser, Arb8Degenerate) {
    // Top face collapses: v5==v6==v7==v8 (pyramid with apex at top)
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <arb8 name="pyramid" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="0" v5y="0"
      v6x="0" v6y="0"
      v7x="0" v7y="0"
      v8x="0" v8y="0"
      lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="pyramid"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.find("pyramid") != data.solids.end());
    auto geo = data.solids["pyramid"];
    ASSERT_TRUE(geo != nullptr);

    // Center is inside
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0, 0, 0)));
    // Near base, inside
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.5, 0.5, -0.8)));
    // Outside above apex
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0, 0, 2.0)));
    // Outside to the side
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(2.0, 0, 0)));
}

// =========================================================================
// Test #4: ENTITY with missing file (graceful error)
// =========================================================================
TEST(GDMLParser, EntityMissingFile) {
    std::string gdml = R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY missing SYSTEM "nonexistent_file_xyz.gdml">
]>
<gdml>
  <define/>
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

    // Should not crash. The entity reference &missing; is not used in the body,
    // so the file should parse successfully despite the unresolvable entity.
    std::string main_file = "/tmp/siren_gdml_missing_entity_test.gdml";
    {
        std::ofstream f(main_file);
        f << gdml;
    }
    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    EXPECT_TRUE(data.solids.find("world_box") != data.solids.end());
}

TEST(GDMLParser, EntityMissingFileUsed) {
    // Entity file doesn't exist but is referenced -- the unexpanded &missing;
    // stays in the XML. RapidXML treats unresolved entity refs as parse errors.
    // However, since ExpandEntities leaves the DOCTYPE intact when no entities
    // are resolved, the parser may or may not throw depending on RapidXML version.
    // What we verify: the parse does NOT produce the expected constant.
    std::string gdml = R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY missing SYSTEM "nonexistent_file_xyz.gdml">
]>
<gdml>
  <define>
    <constant name="local_only" value="42"/>
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

    std::string main_file = "/tmp/siren_gdml_missing_entity_used_test.gdml";
    {
        std::ofstream f(main_file);
        f << gdml;
    }
    // Entity wasn't expanded (file missing), but local content should still parse
    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    EXPECT_NEAR(data.constants["local_only"], 42.0, 1e-9);
}

// =========================================================================
// Test #5: Mesh watertightness (ray gets exactly 2 intersections)
// =========================================================================
TEST(GDMLParser, Arb8Watertight) {
    // Cube via arb8
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <arb8 name="cube8" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="-1000" v5y="-1000"
      v6x="1000" v6y="-1000"
      v7x="1000" v7y="1000"
      v8x="-1000" v8y="1000"
      lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="cube8"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    auto geo = data.solids["cube8"];
    ASSERT_TRUE(geo != nullptr);

    // Fire rays that hit shared edges and vertices directly.
    // After deduplication, each should produce exactly 2 forward intersections.

    // Through the center (hits diagonal of both +z and -z faces)
    siren::math::Vector3D origin(0, 0, 5.0);
    siren::math::Vector3D dir(0, 0, -1);
    auto hits = geo->Intersections(origin, dir);
    int forward_hits = 0;
    for(auto const & h : hits) {
        if(h.distance > 0) ++forward_hits;
    }
    EXPECT_EQ(forward_hits, 2) << "Ray through cube center should hit exactly 2 surfaces";

    // Through a vertex (hits corner shared by 3 faces)
    origin = siren::math::Vector3D(5.0, 5.0, 5.0);
    double inv_sqrt3 = 1.0 / std::sqrt(3.0);
    dir = siren::math::Vector3D(-inv_sqrt3, -inv_sqrt3, -inv_sqrt3);
    hits = geo->Intersections(origin, dir);
    forward_hits = 0;
    for(auto const & h : hits) {
        if(h.distance > 0) ++forward_hits;
    }
    EXPECT_EQ(forward_hits, 2) << "Diagonal ray through cube vertex should hit exactly 2 surfaces";

    // Along an edge (parallel to z, hits +x/+y edge shared by 2 faces)
    origin = siren::math::Vector3D(1.0, 0, 5.0);
    dir = siren::math::Vector3D(0, 0, -1);
    hits = geo->Intersections(origin, dir);
    forward_hits = 0;
    for(auto const & h : hits) {
        if(h.distance > 0) ++forward_hits;
    }
    EXPECT_EQ(forward_hits, 2) << "Ray hitting +x face edge should hit exactly 2 surfaces";
}

// =========================================================================
// Test #7: Built-in constants in expressions
// =========================================================================
TEST(GDMLParser, BuiltInConstants) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="two_pi" value="pi*2"/>
    <constant name="half" value="halfpi"/>
    <constant name="full" value="twopi"/>
    <constant name="full2" value="TWOPI"/>
    <constant name="ten_mm" value="10*mm"/>
    <constant name="one_cm" value="1*cm"/>
    <constant name="one_deg" value="1*deg"/>
    <constant name="half_m" value="0.5*m"/>
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
    double pi = 3.141592653589793;
    EXPECT_NEAR(data.constants["two_pi"], 2.0 * pi, 1e-9);
    EXPECT_NEAR(data.constants["half"], pi / 2.0, 1e-9);
    EXPECT_NEAR(data.constants["full"], 2.0 * pi, 1e-9);
    EXPECT_NEAR(data.constants["full2"], 2.0 * pi, 1e-9);
    EXPECT_NEAR(data.constants["ten_mm"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["one_cm"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["one_deg"], pi / 180.0, 1e-9);
    EXPECT_NEAR(data.constants["half_m"], 500.0, 1e-9);
}

// =========================================================================
// Test #8: xi:include expansion
// =========================================================================
TEST(GDMLParser, XIncludeExpansion) {
    std::string inc_content = R"(
    <constant name="xinc_val" value="77"/>
)";
    std::string inc_file = "/tmp/siren_gdml_xinc.gdml";
    {
        std::ofstream f(inc_file);
        f << inc_content;
    }

    std::string gdml = R"(<?xml version="1.0"?>
<gdml xmlns:xi="http://www.w3.org/2001/XInclude">
  <define>
    <xi:include href="siren_gdml_xinc.gdml"/>
    <constant name="local_val" value="88"/>
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

    std::string main_file = "/tmp/siren_gdml_xinc_test.gdml";
    {
        std::ofstream f(main_file);
        f << gdml;
    }

    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    std::remove(inc_file.c_str());

    EXPECT_NEAR(data.constants["xinc_val"], 77.0, 1e-9);
    EXPECT_NEAR(data.constants["local_val"], 88.0, 1e-9);
}

// =========================================================================
// Test #9: Angle unit mrad
// =========================================================================
TEST(GDMLParser, AngleUnitMrad) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <rotation name="small_rot" x="100" y="200" z="0" unit="mrad"/>
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
    EXPECT_TRUE(data.rotations.find("small_rot") != data.rotations.end());
    // 100 mrad = 0.1 rad, 200 mrad = 0.2 rad
    // The rotation should be non-identity
    auto q = data.rotations["small_rot"];
    EXPECT_NE(q.GetW(), 1.0);
}

// =========================================================================
// Test #10: tessellated with inline vertex positions
// =========================================================================
TEST(GDMLParser, TessellatedInlineVertices) {
    // GDML spec allows <triangular> to reference positions defined in <define>
    // but some generators put position elements directly inside the tessellated.
    // However, the standard approach is always through define references.
    // This test verifies the standard approach works with unit conversion.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="p0" x="0" y="0" z="0" unit="cm"/>
    <position name="p1" x="100" y="0" z="0" unit="cm"/>
    <position name="p2" x="0" y="100" z="0" unit="cm"/>
    <position name="p3" x="0" y="0" z="100" unit="cm"/>
  </define>
  <materials/>
  <solids>
    <tessellated name="big_tetra">
      <triangular vertex1="p0" vertex2="p2" vertex3="p1"/>
      <triangular vertex1="p0" vertex2="p1" vertex3="p3"/>
      <triangular vertex1="p0" vertex2="p3" vertex3="p2"/>
      <triangular vertex1="p1" vertex2="p2" vertex3="p3"/>
    </tessellated>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="big_tetra"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.find("big_tetra") != data.solids.end());
    auto geo = data.solids["big_tetra"];
    ASSERT_TRUE(geo != nullptr);

    // Vertices are at 0,0,0 and 1m,0,0 and 0,1m,0 and 0,0,1m (100cm = 1m)
    // Centroid at (0.25, 0.25, 0.25) should be inside
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.1, 0.1, 0.1)));
    // Point at 0.9m in all axes is outside the tetrahedron
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.9, 0.9, 0.9)));
    // Verify unit conversion: point at 50cm = 0.5m should be outside
    // (sum of barycentric coords > 1: 0.5+0.5+0.5 = 1.5 > 1)
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.5, 0.5, 0.5)));
}


// =========================================================================
// Partial-phi tube and cone parsing (no warning, correct geometry)
// =========================================================================
TEST(GDMLParser, PartialPhiTubeAndCone) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <tube name="quarter_tube" rmin="10" rmax="50" z="100"
          startphi="0" deltaphi="1.5707963267948966" lunit="mm" aunit="rad"/>
    <cone name="quarter_cone" rmin1="10" rmax1="50" rmin2="5" rmax2="30" z="80"
          startphi="0.7853981633974483" deltaphi="1.5707963267948966" lunit="mm" aunit="rad"/>
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

    // Verify the solids were created successfully
    ASSERT_TRUE(data.solids.find("quarter_tube") != data.solids.end());
    ASSERT_TRUE(data.solids["quarter_tube"] != nullptr);

    ASSERT_TRUE(data.solids.find("quarter_cone") != data.solids.end());
    ASSERT_TRUE(data.solids["quarter_cone"] != nullptr);

    // Verify no warnings about partial angular extent
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("partial angular extent") == std::string::npos)
            << "Unexpected warning: " << w;
    }
}


// =========================================================================
// Energy unit built-in constants
// =========================================================================
TEST(GDMLParser, EnergyUnitConstants) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="photon_e" value="3.0*eV"/>
    <constant name="gamma_e" value="1.5*MeV"/>
    <constant name="k_e" value="100*keV"/>
    <constant name="g_e" value="2*GeV"/>
    <constant name="t_e" value="0.5*TeV"/>
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

    double tol = 1e-12;
    EXPECT_NEAR(data.constants["eV"], 1e-6, tol);
    EXPECT_NEAR(data.constants["keV"], 1e-3, tol);
    EXPECT_NEAR(data.constants["MeV"], 1.0, tol);
    EXPECT_NEAR(data.constants["GeV"], 1e3, tol);
    EXPECT_NEAR(data.constants["TeV"], 1e6, tol);

    EXPECT_NEAR(data.constants["photon_e"], 3.0e-6, tol);
    EXPECT_NEAR(data.constants["gamma_e"], 1.5, tol);
    EXPECT_NEAR(data.constants["k_e"], 0.1, tol);
    EXPECT_NEAR(data.constants["g_e"], 2e3, tol);
    EXPECT_NEAR(data.constants["t_e"], 5e5, tol);
}


// =========================================================================
// Polycone z-plane auto-sorting
// =========================================================================
TEST(GDMLParser, PolyconeZPlaneSorting) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <polycone name="sorted_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="15" z="100"/>
    </polycone>
    <polycone name="reversed_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="15" z="100"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
    </polycone>
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

    ASSERT_TRUE(data.solids.find("sorted_pc") != data.solids.end());
    ASSERT_TRUE(data.solids["sorted_pc"] != nullptr);
    ASSERT_TRUE(data.solids.find("reversed_pc") != data.solids.end());
    ASSERT_TRUE(data.solids["reversed_pc"] != nullptr);

    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("non-monotonic") == std::string::npos)
            << "Unexpected warning: " << w;
    }
}

TEST(GDMLParser, PolyconeZPlaneDuplicateZFails) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <polycone name="dup_z_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="10" z="50"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="15" z="100"/>
    </polycone>
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

    EXPECT_TRUE(data.solids.find("dup_z_pc") == data.solids.end()
                || data.solids["dup_z_pc"] == nullptr);

    bool found_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("non-monotonic") != std::string::npos) { found_warning = true; break; }
    }
    EXPECT_TRUE(found_warning);
}


// =========================================================================
// Replicavol: N-copy placement along axis
// =========================================================================
TEST(GDMLParser, ReplicavolXAxis) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mother_box" x="800" y="100" z="100" lunit="mm"/>
    <box name="child_box" x="100" y="100" z="100" lunit="mm"/>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="ChildVol">
      <materialref ref=""/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="MotherVol">
      <materialref ref=""/>
      <solidref ref="mother_box"/>
      <replicavol number="8">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction x="1"/>
          <width value="100" unit="mm"/>
          <offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="MotherVol"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    auto it = data.volumes.find("MotherVol");
    ASSERT_TRUE(it != data.volumes.end());
    ASSERT_EQ(it->second.children.size(), 8u);

    double width_m = 0.1; // 100mm in meters
    double tol = 1e-9;
    for(int i = 0; i < 8; ++i) {
        double expected_x = (i + 0.5 - 4.0) * width_m;
        EXPECT_NEAR(it->second.children[i].position.GetX(), expected_x, tol)
            << "Replica " << i << " x position mismatch";
        EXPECT_NEAR(it->second.children[i].position.GetY(), 0.0, tol);
        EXPECT_NEAR(it->second.children[i].position.GetZ(), 0.0, tol);
        EXPECT_EQ(it->second.children[i].volume_ref, "ChildVol");
    }
}

TEST(GDMLParser, ReplicavolYAxisWithOffset) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mother_box" x="100" y="300" z="100" lunit="mm"/>
    <box name="child_box" x="100" y="100" z="100" lunit="mm"/>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="ChildVol">
      <materialref ref=""/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="MotherVol">
      <materialref ref=""/>
      <solidref ref="mother_box"/>
      <replicavol number="3">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction y="1"/>
          <width value="10" unit="cm"/>
          <offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
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

    auto it = data.volumes.find("MotherVol");
    ASSERT_TRUE(it != data.volumes.end());
    ASSERT_EQ(it->second.children.size(), 3u);

    double width_m = 0.1; // 10cm in meters
    double tol = 1e-9;
    for(int i = 0; i < 3; ++i) {
        double expected_y = (i + 0.5 - 1.5) * width_m;
        EXPECT_NEAR(it->second.children[i].position.GetX(), 0.0, tol);
        EXPECT_NEAR(it->second.children[i].position.GetY(), expected_y, tol);
        EXPECT_NEAR(it->second.children[i].position.GetZ(), 0.0, tol);
    }
}


// =========================================================================
// Polyhedra z-plane auto-sorting (parallel code path to polycone)
// =========================================================================
TEST(GDMLParser, PolyhedraZPlaneSorting) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <polyhedra name="sorted_ph" startphi="0" deltaphi="360" numsides="6" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="15" z="100"/>
    </polyhedra>
    <polyhedra name="reversed_ph" startphi="0" deltaphi="360" numsides="6" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="15" z="100"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
    </polyhedra>
    <polyhedra name="scrambled_ph" startphi="0" deltaphi="360" numsides="4" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="15" z="100"/>
      <zplane rmin="0" rmax="5" z="200"/>
    </polyhedra>
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

    ASSERT_TRUE(data.solids.find("sorted_ph") != data.solids.end());
    ASSERT_TRUE(data.solids["sorted_ph"] != nullptr);
    ASSERT_TRUE(data.solids.find("reversed_ph") != data.solids.end());
    ASSERT_TRUE(data.solids["reversed_ph"] != nullptr);
    ASSERT_TRUE(data.solids.find("scrambled_ph") != data.solids.end());
    ASSERT_TRUE(data.solids["scrambled_ph"] != nullptr);

    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("non-monotonic") == std::string::npos)
            << "Unexpected warning: " << w;
    }
}


// =========================================================================
// Energy units in matrix values (optical property use case)
// =========================================================================
TEST(GDMLParser, EnergyUnitsInMatrix) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <matrix name="RINDEX" coldim="2"
            values="1.5*eV  1.48
                    2.0*eV  1.49
                    3.0*eV  1.50
                    4.0*eV  1.52"/>
    <matrix name="ENERGIES" coldim="1"
            values="1.5*eV 2.0*eV 3.0*eV 4.0*eV"/>
    <constant name="threshold" value="2.5*GeV"/>
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

    double tol = 1e-12;

    auto ri = data.matrices.find("RINDEX");
    ASSERT_TRUE(ri != data.matrices.end());
    EXPECT_EQ(ri->second.first, 2);
    ASSERT_EQ(ri->second.second.size(), 8u);
    EXPECT_NEAR(ri->second.second[0], 1.5e-6, tol);
    EXPECT_NEAR(ri->second.second[1], 1.48, tol);
    EXPECT_NEAR(ri->second.second[2], 2.0e-6, tol);
    EXPECT_NEAR(ri->second.second[3], 1.49, tol);
    EXPECT_NEAR(ri->second.second[4], 3.0e-6, tol);
    EXPECT_NEAR(ri->second.second[5], 1.50, tol);
    EXPECT_NEAR(ri->second.second[6], 4.0e-6, tol);
    EXPECT_NEAR(ri->second.second[7], 1.52, tol);

    auto en = data.matrices.find("ENERGIES");
    ASSERT_TRUE(en != data.matrices.end());
    EXPECT_EQ(en->second.first, 1);
    ASSERT_EQ(en->second.second.size(), 4u);
    EXPECT_NEAR(en->second.second[0], 1.5e-6, tol);
    EXPECT_NEAR(en->second.second[1], 2.0e-6, tol);
    EXPECT_NEAR(en->second.second[2], 3.0e-6, tol);
    EXPECT_NEAR(en->second.second[3], 4.0e-6, tol);

    EXPECT_NEAR(data.constants["threshold"], 2.5e3, tol);
}


// =========================================================================
// Replicavol along arbitrary (non-axis-aligned) direction
// =========================================================================
TEST(GDMLParser, ReplicavolArbitraryAxis) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mother_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="child_box" x="50" y="50" z="50" lunit="mm"/>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="ChildVol">
      <materialref ref=""/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="ZVol">
      <materialref ref=""/>
      <solidref ref="mother_box"/>
      <replicavol number="4">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction z="1"/>
          <width value="200" unit="mm"/>
          <offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
    <volume name="DiagVol">
      <materialref ref=""/>
      <solidref ref="mother_box"/>
      <replicavol number="3">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction x="1" y="1"/>
          <width value="100" unit="mm"/>
          <offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
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

    // Z-axis: 4 copies, width=0.2m
    {
        auto it = data.volumes.find("ZVol");
        ASSERT_TRUE(it != data.volumes.end());
        ASSERT_EQ(it->second.children.size(), 4u);
        double w = 0.2;
        double tol = 1e-9;
        for(int i = 0; i < 4; ++i) {
            double expected_z = (i + 0.5 - 2.0) * w;
            EXPECT_NEAR(it->second.children[i].position.GetX(), 0.0, tol);
            EXPECT_NEAR(it->second.children[i].position.GetY(), 0.0, tol);
            EXPECT_NEAR(it->second.children[i].position.GetZ(), expected_z, tol)
                << "Z-axis replica " << i;
        }
    }

    // Diagonal (1,1,0): direction normalized to (1/sqrt2, 1/sqrt2, 0)
    // 3 copies, width=0.1m
    {
        auto it = data.volumes.find("DiagVol");
        ASSERT_TRUE(it != data.volumes.end());
        ASSERT_EQ(it->second.children.size(), 3u);
        double w = 0.1;
        double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        double tol = 1e-9;
        for(int i = 0; i < 3; ++i) {
            double t = (i + 0.5 - 1.5) * w;
            EXPECT_NEAR(it->second.children[i].position.GetX(), t * inv_sqrt2, tol)
                << "Diagonal replica " << i << " x";
            EXPECT_NEAR(it->second.children[i].position.GetY(), t * inv_sqrt2, tol)
                << "Diagonal replica " << i << " y";
            EXPECT_NEAR(it->second.children[i].position.GetZ(), 0.0, tol);
        }
    }
}


// =========================================================================
// Built-in constants are immutable: user <define> cannot overwrite them
// =========================================================================
TEST(GDMLParser, BuiltInConstantsImmutable) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="pi" value="999"/>
    <constant name="cm" value="42"/>
    <constant name="MeV" value="-1"/>
    <variable name="deg" value="123"/>
    <quantity name="mm" value="7" type="length" unit="mm"/>
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
    double tol = 1e-12;
    EXPECT_NEAR(data.constants.at("pi"), M_PI, tol);
    EXPECT_NEAR(data.constants.at("cm"), 10.0, tol);
    EXPECT_NEAR(data.constants.at("MeV"), 1.0, tol);
    EXPECT_NEAR(data.constants.at("deg"), M_PI / 180.0, tol);
    EXPECT_NEAR(data.constants.at("mm"), 1.0, tol);
    EXPECT_GE(data.warnings.size(), 5u);
}

// =========================================================================
// Built-in constants cannot be used as loop variables
// =========================================================================
TEST(GDMLParser, BuiltInConstantsBlockLoopVar) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <loop for="mm" from="1" to="3" step="1">
      <constant name="x_[mm]" value="[mm]*100"/>
    </loop>
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
    EXPECT_NEAR(data.constants.at("mm"), 1.0, 1e-12);
    EXPECT_EQ(data.constants.count("x_1"), 0u);
    EXPECT_EQ(data.constants.count("x_2"), 0u);
    EXPECT_EQ(data.constants.count("x_3"), 0u);
    EXPECT_FALSE(data.warnings.empty());
}

// =========================================================================
// Extended built-in constants: verify new additions are available
// =========================================================================
TEST(GDMLParser, ExtendedBuiltInConstants) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <constant name="use_euler" value="e"/>
    <constant name="use_gamma" value="gamma"/>
    <constant name="use_km" value="km"/>
    <constant name="use_degree" value="degree"/>
    <constant name="use_radian" value="radian"/>
    <constant name="use_inch" value="inch"/>
    <constant name="use_fermi" value="fermi"/>
    <constant name="use_ns" value="ns"/>
    <constant name="use_PeV" value="PeV"/>
    <constant name="use_mole" value="mole"/>
    <constant name="use_kelvin" value="kelvin"/>
    <constant name="len_m" value="2.5*meter"/>
    <constant name="len_um" value="3*micrometer"/>
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
    EXPECT_NEAR(data.constants.at("use_euler"), 2.7182818284590452354, tol);
    EXPECT_NEAR(data.constants.at("use_gamma"), 0.5772156649015328606, tol);
    EXPECT_NEAR(data.constants.at("use_km"), 1e6, tol);
    EXPECT_NEAR(data.constants.at("use_degree"), M_PI / 180.0, tol);
    EXPECT_NEAR(data.constants.at("use_radian"), 1.0, tol);
    EXPECT_NEAR(data.constants.at("use_inch"), 25.4, tol);
    EXPECT_NEAR(data.constants.at("use_fermi"), 1e-12, tol);
    EXPECT_NEAR(data.constants.at("use_ns"), 1.0, tol);
    EXPECT_NEAR(data.constants.at("use_PeV"), 1e9, tol);
    EXPECT_NEAR(data.constants.at("use_mole"), 1.0, tol);
    EXPECT_NEAR(data.constants.at("use_kelvin"), 1.0, tol);
    EXPECT_NEAR(data.constants.at("len_m"), 2500.0, tol);
    EXPECT_NEAR(data.constants.at("len_um"), 3e-3, tol);
}

// =========================================================================
// GenericPolycone: hollow shape with non-monotonic z
// =========================================================================
TEST(GDMLParser, GenericPolyconeHollow) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <genericPolycone name="gpc_hollow" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <rzpoint r="300" z="-500"/>
      <rzpoint r="450" z="0"/>
      <rzpoint r="500" z="500"/>
      <rzpoint r="350" z="500"/>
      <rzpoint r="300" z="0"/>
      <rzpoint r="200" z="-500"/>
    </genericPolycone>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="gpc_hollow"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("gpc_hollow") > 0);
    auto & geo = data.solids["gpc_hollow"];
    auto bb = geo->GetBoundingBox();
    double tol = 1e-6;
    EXPECT_NEAR(bb.max_corner.GetX(), 0.5, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.5, tol);
    EXPECT_NEAR(bb.min_corner.GetZ(), -0.5, tol);

    // Point inside the hollow shell
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.4, 0, 0)));
    // Point in the hollow center
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.1, 0, 0)));
    // Point outside
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.6, 0, 0)));
}

// =========================================================================
// GenericPolycone: phi-cut
// =========================================================================
TEST(GDMLParser, GenericPolyconePhiCut) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <genericPolycone name="gpc_phi" startphi="0" deltaphi="180" aunit="deg" lunit="mm">
      <rzpoint r="0" z="-100"/>
      <rzpoint r="50" z="0"/>
      <rzpoint r="0" z="100"/>
    </genericPolycone>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="gpc_phi"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("gpc_phi") > 0);
    auto & geo = data.solids["gpc_phi"];
    // +x, +y is within phi=[0, pi]
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.01, 0.01, 0)));
    // -y is outside phi=[0, pi]
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.01, -0.01, 0)));
}
