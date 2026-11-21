// Included by GDMLParser_TEST.cxx. Do not compile directly.

// ---------------------------------------------------------------------------
// Integration smoke test: a single GDML exercising define (with expressions),
// materials, solids, and a two-level structure with a named-position physvol.
// ---------------------------------------------------------------------------

const char* kSmokeGDML = R"GDML(<?xml version="1.0" encoding="UTF-8"?>
<gdml>
  <define>
    <constant name="hw" value="5.0"/>
    <constant name="neg" value="-2+3"/>
    <position name="ctr" x="1" y="2" z="3" unit="m"/>
  </define>
  <materials>
    <material name="Water" Z="7.42">
      <D value="1.0" unit="g/cm3"/>
      <atom value="18.0"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="hw" y="hw" z="hw" lunit="m"/>
    <orb name="ball" r="2" lunit="m"/>
  </solids>
  <structure>
    <volume name="ballVol">
      <materialref ref="Water"/>
      <solidref ref="ball"/>
    </volume>
    <volume name="worldVol">
      <materialref ref="Water"/>
      <solidref ref="world_box"/>
      <physvol>
        <volumeref ref="ballVol"/>
        <positionref ref="ctr"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="worldVol"/>
  </setup>
</gdml>
)GDML";

TEST(GDMLParser, ParsesDefineSection) {
    TempGDMLFile f(kSmokeGDML);
    GDMLData data = ParseGDML(f.path);

    ASSERT_EQ(data.constants.count("hw"), 1u);
    EXPECT_NEAR(data.constants.at("hw"), 5.0, 1e-12);

    ASSERT_EQ(data.constants.count("neg"), 1u);
    EXPECT_NEAR(data.constants.at("neg"), 1.0, 1e-12);

    ASSERT_EQ(data.positions.count("ctr"), 1u);
    siren::math::Vector3D const & ctr = data.positions.at("ctr");
    EXPECT_NEAR(ctr.GetX(), 1.0, 1e-9);
    EXPECT_NEAR(ctr.GetY(), 2.0, 1e-9);
    EXPECT_NEAR(ctr.GetZ(), 3.0, 1e-9);
}

TEST(GDMLParser, ParsesMaterials) {
    TempGDMLFile f(kSmokeGDML);
    GDMLData data = ParseGDML(f.path);

    ASSERT_EQ(data.materials.count("Water"), 1u);
    GDMLMaterial const & w = data.materials.at("Water");
    EXPECT_NEAR(w.density, 1.0, 1e-9);
    EXPECT_NEAR(w.Z, 7.42, 1e-9);
}

TEST(GDMLParser, ParsesSolids) {
    TempGDMLFile f(kSmokeGDML);
    GDMLData data = ParseGDML(f.path);

    EXPECT_EQ(data.solids.count("world_box"), 1u);
    EXPECT_EQ(data.solids.count("ball"), 1u);
    ASSERT_NE(data.solids.at("world_box"), nullptr);
    ASSERT_NE(data.solids.at("ball"), nullptr);
}

TEST(GDMLParser, ParsesStructureAndPlacement) {
    TempGDMLFile f(kSmokeGDML);
    GDMLData data = ParseGDML(f.path);

    EXPECT_EQ(data.world_volume, "worldVol");

    ASSERT_EQ(data.volumes.count("ballVol"), 1u);
    ASSERT_EQ(data.volumes.count("worldVol"), 1u);

    GDMLVolume const & world = data.volumes.at("worldVol");
    EXPECT_EQ(world.solid_ref, "world_box");
    EXPECT_EQ(world.material_ref, "Water");
    ASSERT_EQ(world.children.size(), 1u);
    EXPECT_EQ(world.children[0].volume_ref, "ballVol");
    EXPECT_NEAR(world.children[0].position.GetX(), 1.0, 1e-9);
    EXPECT_NEAR(world.children[0].position.GetY(), 2.0, 1e-9);
    EXPECT_NEAR(world.children[0].position.GetZ(), 3.0, 1e-9);
}

TEST(GDMLParser, MissingFileThrows) {
    EXPECT_THROW(ParseGDML("this_file_does_not_exist_12345.gdml"), std::runtime_error);
}

// ---------------------------------------------------------------------------

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
