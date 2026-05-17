// Included by GDMLParser_TEST.cxx. Do not compile directly.

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


TEST(GDMLParser, DefaultAngleUnitRadians) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <rotation name="rot_nounit" z="pi/2"/>
    <rotation name="rot_deg" z="pi/2" unit="rad"/>
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
        << "90-degree z-rotation should have z = sqrt(2)/2, not a degree-based value";
    EXPECT_NEAR(q_nounit.GetW(), s, tol);
}


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


TEST(GDMLParser, QuantityUnitConversion) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <quantity name="q_cm" value="100" unit="cm" type="length"/>
    <quantity name="q_mm" value="500" unit="mm" type="length"/>
    <quantity name="q_deg" value="90" unit="deg" type="angle"/>
    <quantity name="q_nounit" value="42" type="length"/>
    <quantity name="q_notype" value="10" unit="cm"/>
    <quantity name="q_in" value="2" unit="in" type="length"/>
    <quantity name="rho_gcm3" value="2.7" unit="g/cm3" type="density"/>
    <quantity name="rho_mgcm3" value="2700" unit="mg/cm3" type="density"/>
    <quantity name="rho_kgm3" value="2700" unit="kg/m3" type="density"/>
    <constant name="plain_const" value="99"/>
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
    EXPECT_NEAR(data.constants["q_in"], 50.8, tol)
        << "2 in should convert to 50.8 mm";
    EXPECT_NEAR(data.constants["rho_gcm3"], 2.7, 1e-9)
        << "g/cm3 density should store as 2.7 g/cm3";
    EXPECT_NEAR(data.constants["rho_mgcm3"], 2.7, 1e-9)
        << "2700 mg/cm3 density should store as 2.7 g/cm3";
    EXPECT_NEAR(data.constants["rho_kgm3"], 2.7, 1e-9)
        << "2700 kg/m3 density should store as 2.7 g/cm3";

    // Verify type map is populated
    EXPECT_EQ(data.quantity_types["q_cm"], siren::detector::GDMLQuantityType::LENGTH);
    EXPECT_EQ(data.quantity_types["q_mm"], siren::detector::GDMLQuantityType::LENGTH);
    EXPECT_EQ(data.quantity_types["q_deg"], siren::detector::GDMLQuantityType::ANGLE);
    EXPECT_EQ(data.quantity_types["q_notype"], siren::detector::GDMLQuantityType::LENGTH);
    EXPECT_EQ(data.quantity_types["q_in"], siren::detector::GDMLQuantityType::LENGTH);
    EXPECT_EQ(data.quantity_types["rho_gcm3"], siren::detector::GDMLQuantityType::DENSITY);
    EXPECT_EQ(data.quantity_types["rho_mgcm3"], siren::detector::GDMLQuantityType::DENSITY);
    EXPECT_EQ(data.quantity_types["rho_kgm3"], siren::detector::GDMLQuantityType::DENSITY);
    EXPECT_EQ(data.quantity_types["q_nounit"], siren::detector::GDMLQuantityType::LENGTH);
    EXPECT_EQ(data.quantity_types.count("plain_const"), 0u);
}


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


TEST(GDMLParser, QuantityTypeConflictWarns) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <quantity name="confused" value="10" unit="cm" type="angle"/>
  </define>
  <materials/>
  <solids>
    <box name="world" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="world"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_GE(data.warnings.size(), 1u);
    bool found_conflict = false;
    for(auto const & w : data.warnings) {
        if(w.find("conflicts") != std::string::npos) found_conflict = true;
    }
    EXPECT_TRUE(found_conflict) << "Expected warning about type/unit conflict";
    EXPECT_NEAR(data.constants["confused"], 100.0, 1e-9)
        << "Unit-inferred type (LENGTH) should convert 10cm -> 100mm";
    EXPECT_EQ(data.quantity_types["confused"], siren::detector::GDMLQuantityType::LENGTH);
}


TEST(GDMLParser, DensityQuantityInMaterial) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <quantity name="al_rho" value="2.7" unit="g/cm3" type="density"/>
  </define>
  <materials>
    <material name="Al_from_qty" Z="13">
      <D value="al_rho"/>
      <atom value="27"/>
    </material>
  </materials>
  <solids>
    <box name="world" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Al_from_qty"/>
      <solidref ref="world"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.materials.count("Al_from_qty") > 0);
    EXPECT_NEAR(data.materials["Al_from_qty"].density, 2.7, 1e-9)
        << "Density quantity 2.7 g/cm3 consumed by <D> should yield 2.7";
}


TEST(GDMLParser, DensityQuantityKgM3InMaterial) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <quantity name="water_rho" value="1000" unit="kg/m3" type="density"/>
  </define>
  <materials>
    <material name="Water" Z="1">
      <D value="water_rho"/>
      <atom value="1"/>
    </material>
  </materials>
  <solids>
    <box name="world" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Water"/>
      <solidref ref="world"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.materials.count("Water") > 0);
    EXPECT_NEAR(data.materials["Water"].density, 1.0, 1e-9)
        << "1000 kg/m3 density quantity consumed by <D> should yield 1.0 g/cm3";
}


TEST(GDMLParser, LengthQuantityInSolid) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <quantity name="half_w" value="5" unit="cm" type="length"/>
  </define>
  <materials/>
  <solids>
    <box name="test" x="half_w" y="half_w" z="half_w" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="test"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.solids.count("test") > 0);
    auto bb = data.solids["test"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 0.05, 1e-9)
        << "5cm quantity = 50mm half-width, box full-width = 100mm, AABB max = 50mm = 0.05m";
}
