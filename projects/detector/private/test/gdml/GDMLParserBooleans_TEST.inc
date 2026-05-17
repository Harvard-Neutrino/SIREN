// Included by GDMLParser_TEST.cxx. Do not compile directly.

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
    // Subtraction probe: origin is inside rotated box_a AND inside box_b
    // (box_b half-width = 10mm = 0.01m), so subtraction should exclude it
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0, 0, 0), dir));
    // Just outside box_b: z=0.012 > 0.01m half-width, but inside rotated box_a
    // (z half-width = 25mm = 0.025m), so subtraction should include it
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0, 0, 0.012), dir));
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


TEST(GDMLParser, BooleanSubtractionHollow) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="outer" x="100" y="100" z="100" lunit="mm"/>
    <box name="inner" x="60" y="60" z="60" lunit="mm"/>
    <subtraction name="hollow_box">
      <first ref="outer"/>
      <second ref="inner"/>
    </subtraction>
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
    ASSERT_TRUE(data.solids.count("hollow_box") > 0);
    auto geo = data.solids["hollow_box"];
    ASSERT_NE(geo, nullptr);
    siren::math::Vector3D dir(0, 0, 1);

    // Center is inside both boxes, so subtraction should exclude it
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0, 0, 0), dir));
    // In the wall: x=0.08m is outside inner (half-width 0.06m) but inside outer (0.1m)
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.08, 0, 0), dir));
    // In the wall off-axis
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0, 0.08, 0), dir));
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0, 0, 0.08), dir));
    // Outside the outer box entirely (half-width 0.1m)
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.12, 0, 0), dir));
    // Inside the inner box: subtracted away
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.03, 0.03, 0.03), dir));
}


TEST(GDMLParser, NestedBooleanSolidsContainment) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <sphere name="big_sphere" rmax="100" lunit="mm"/>
    <sphere name="small_sphere" rmax="60" lunit="mm"/>
    <box name="some_box" x="40" y="40" z="40" lunit="mm"/>
    <subtraction name="inner_sub">
      <first ref="big_sphere"/>
      <second ref="small_sphere"/>
    </subtraction>
    <union name="final">
      <first ref="inner_sub"/>
      <second ref="some_box"/>
    </union>
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
    ASSERT_TRUE(data.solids.count("final") > 0);
    auto geo = data.solids["final"];
    ASSERT_NE(geo, nullptr);
    siren::math::Vector3D dir(0, 0, 1);

    // In the shell: r=0.08m is between rmin=0.06 and rmax=0.1
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0.08, 0, 0), dir));
    // In the hollow but inside some_box: (0,0,0) is in the box (half-width=0.04m)
    EXPECT_TRUE(geo->IsInside(siren::math::Vector3D(0, 0, 0), dir));
    // In the hollow but outside some_box: r=0.05 < rmin=0.06, outside box (0.05 > 0.04)
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.05, 0, 0), dir));
    // Outside everything: r=0.15 > rmax=0.1
    EXPECT_FALSE(geo->IsInside(siren::math::Vector3D(0.15, 0, 0), dir));
}


TEST(GDMLParser, UnresolvedBooleanOperandWarns) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="existing_box" x="100" y="100" z="100" lunit="mm"/>
    <subtraction name="bad_sub">
      <first ref="existing_box"/>
      <second ref="nonexistent_solid"/>
    </subtraction>
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

    // The boolean solid should not be created
    EXPECT_EQ(data.solids.count("bad_sub"), 0u);
    // A warning should have been emitted
    bool found_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("bad_sub") != std::string::npos && w.find("unresolved") != std::string::npos) {
            found_warning = true;
            break;
        }
    }
    EXPECT_TRUE(found_warning)
        << "Should warn about unresolved operand in boolean solid 'bad_sub'";
}
