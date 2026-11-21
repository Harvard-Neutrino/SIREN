// Included by GDMLParser_TEST.cxx. Do not compile directly.

static const std::string kAirIron = R"(
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>)";

TEST(GDMLParser, CircularVolumeReference) {
    // Create a GDML where volume A contains physvol referencing A (self-loop)
    std::string gdml = GDMLFull("", R"(
    <material name="Air"><D value="0.00129" unit="g/cm3"/><atom value="14" Z="7"/></material>)",
        R"(<box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>)",
        R"(<volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol><volumeref ref="World"/></physvol>
    </volume>)", "World");

    // Circular reference should throw an error
    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error)
        << "LoadGDML should throw on circular volume references";
}


TEST(GDMLParser, PhysvolConstantsPosition) {
    std::string gdml = GDMLFull(
        R"(<constant name="offset" value="500"/>)",
        R"(<element name="Hydrogen" Z="1"><atom value="1.008"/></element>
    <element name="Oxygen" Z="8"><atom value="15.999"/></element>
    <material name="Water"><D value="1.0" unit="g/cm3"/>
      <fraction ref="Hydrogen" n="0.112"/><fraction ref="Oxygen" n="0.888"/>
    </material>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>)",
        R"(<box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="child_box" x="200" y="200" z="200" lunit="mm"/>)",
        R"(<volume name="Child">
      <materialref ref="Iron"/><solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Water"/><solidref ref="world_box"/>
      <physvol><volumeref ref="Child"/><position x="offset" unit="mm"/></physvol>
    </volume>)", "World");

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
    auto dm = LoadGDMLFromString(gdml);

    // A point at (0.5, 0, 0) meters should be inside the Iron child
    EXPECT_EQ(MaterialAt(*dm, 0.5, 0.0, 0.0), "Iron")
        << "Point at (0.5,0,0) should be in Iron child placed at offset=500mm";
}


TEST(GDMLParser, DuplicateVolumeNames) {
    std::string gdml = GDMLFull("", kAirIron,
        R"(<box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="child_box" x="200" y="200" z="200" lunit="mm"/>)",
        R"(<volume name="Child">
      <materialref ref="Iron"/><solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Air"/><solidref ref="world_box"/>
      <physvol><volumeref ref="Child"/><position x="100" unit="mm"/></physvol>
      <physvol><volumeref ref="Child"/><position x="-100" unit="mm"/></physvol>
    </volume>)", "World");

    auto dm = LoadGDMLFromString(gdml);

    // World + 2 children + UNIVERSE = at least 4 sectors
    auto const & sectors = dm->GetSectors();
    EXPECT_GE(sectors.size(), 4u)
        << "Detector model should have at least 4 sectors (UNIVERSE + world + 2 duplicate children)";

    // Both child positions should resolve to Iron
    EXPECT_EQ(MaterialAt(*dm, 0.1, 0.0, 0.0), "Iron")
        << "Point at (0.1,0,0) should be in first Iron child";
    EXPECT_EQ(MaterialAt(*dm, -0.1, 0.0, 0.0), "Iron")
        << "Point at (-0.1,0,0) should be in second Iron child";
}


TEST(GDMLParser, RotationHandling) {
    std::string gdml = GDMLFull("", kAirIron,
        R"(<box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="child_box" x="400" y="100" z="100" lunit="mm"/>)",
        R"(<volume name="Child">
      <materialref ref="Iron"/><solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Air"/><solidref ref="world_box"/>
      <physvol><volumeref ref="Child"/>
        <position x="200" unit="mm"/><rotation z="90" unit="deg"/>
      </physvol>
    </volume>)", "World");

    auto dm = LoadGDMLFromString(gdml);

    // The child box is 400x100x100 mm full-widths (half-extents 0.2x0.05x0.05 m).
    // Placed at x=200mm=0.2m with z-rotation of 90 degrees.
    // After 90-deg z-rotation, the long axis (originally x) becomes the y axis.
    // So the child occupies roughly: x in [0.15, 0.25], y in [-0.2, 0.2]

    // The child sector should exist
    auto const & sectors = dm->GetSectors();
    EXPECT_GE(sectors.size(), 3u)
        << "Should have at least 3 sectors (UNIVERSE + world + rotated child)";

    // A point at the child center (0.2, 0, 0) should be Iron
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.0, 0.0), "Iron")
        << "Center of rotated child at (0.2,0,0) should be Iron";

    // After rotation, the long axis is along y. A point at (0.2, 0.15, 0) should
    // still be inside (within the rotated half-extent of 0.2 along y).
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.15, 0.0), "Iron")
        << "Point at (0.2, 0.15, 0) should be inside the rotated child (long axis now along y)";

    // A point far from the child should be in Air
    EXPECT_EQ(MaterialAt(*dm, 1.0, 1.0, 0.0), "Air")
        << "Point at (1.0,1.0,0) should be in Air (outside rotated child)";
}


TEST(GDMLParser, CompoundRotationGroundTruth) {
    // child_box GDML full-widths 1.0x0.4x0.2 m -> half-extents 0.5x0.2x0.1 m,
    // matching the local probe margins (0.1/0.2/0.1 at local (0.4,0,0)) used below.
    std::string gdml = GDMLFull("", kAirIron,
        R"(<box name="world_box" x="5" y="5" z="5" lunit="m"/>
    <box name="child_box" x="1.0" y="0.4" z="0.2" lunit="m"/>)",
        R"(<volume name="Child">
      <materialref ref="Iron"/><solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Air"/><solidref ref="world_box"/>
      <physvol><volumeref ref="Child"/>
        <rotation x="30" z="60" unit="deg"/>
      </physvol>
    </volume>)", "World");

    auto dm = LoadGDMLFromString(gdml);

    // Box center maps to local origin under any convention: sanity that the
    // rotated child sector was built and is reachable.
    EXPECT_EQ(MaterialAt(*dm, 0.0, 0.0, 0.0), "Iron")
        << "Rotated child center must be Iron";

    // P1 = 0.4 * (first row of M). Geant4 local = (0.4, 0, 0): comfortably
    // inside (margins 0.1/0.2/0.1). FAILS as Air under the QFromXYZs sign
    // bug and under the wrong (active) handedness.
    EXPECT_EQ(MaterialAt(*dm, 0.2, -0.3, 0.17320508075688773), "Iron")
        << "P1 is inside the box under the Geant4-true passive mapping; "
           "Air here indicates the QFromXYZs sign bug or wrong handedness";

    // P3 = 0.4 * (first column of M). Geant4 local = (-0.1598, 0.3232,
    // 0.1732): outside (|y|,|z| exceed half-extents). The active-handedness
    // bug instead maps P3 to (0.4,0,0) -> Iron.
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.3464101615137755, 0.0), "Air")
        << "P3 is outside the box under the Geant4-true passive mapping; "
           "Iron here indicates the rotation is applied with wrong handedness";

    // Far point: outside the child under every convention.
    EXPECT_EQ(MaterialAt(*dm, 1.0, 1.0, 1.0), "Air")
        << "Point far from the rotated child must be Air";
}


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


TEST(GDMLParser, AssemblyVolume) {
    std::string gdml = GDMLFull("", R"(
    <material name="Air"><D value="0.00129" unit="g/cm3"/><atom value="14" Z="7"/></material>)",
        R"(<box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="child_box" x="100" y="100" z="100" lunit="mm"/>)",
        R"(<volume name="ChildVol">
      <materialref ref="Air"/><solidref ref="child_box"/>
    </volume>
    <assembly name="MyAssembly">
      <physvol><volumeref ref="ChildVol"/>
        <position name="pos1" x="200" y="0" z="0" unit="mm"/>
      </physvol>
    </assembly>
    <volume name="World">
      <materialref ref="Air"/><solidref ref="world_box"/>
      <physvol><volumeref ref="MyAssembly"/></physvol>
    </volume>)", "World");

    GDMLData data = ParseGDMLString(gdml);

    // Assembly should be stored as a volume with is_assembly=true
    auto it = data.volumes.find("MyAssembly");
    ASSERT_NE(it, data.volumes.end());
    EXPECT_TRUE(it->second.is_assembly);
    EXPECT_EQ(it->second.children.size(), 1u);

    // Loading into DetectorModel should create a sector for ChildVol but not for MyAssembly
    auto dm = LoadGDMLFromString(gdml);

    bool found_child = false;
    bool found_assembly = false;
    for(auto const & sector : dm->GetSectors()) {
        if(sector.name == "ChildVol") found_child = true;
        if(sector.name == "MyAssembly") found_assembly = true;
    }
    EXPECT_TRUE(found_child) << "ChildVol should have a sector";
    EXPECT_FALSE(found_assembly) << "Assembly should not have its own sector";
}


TEST(GDMLParser, BuildVolumeMissingVolumeThrows) {
    std::string gdml = GDMLFull("", R"(
    <material name="Air"><D value="0.00129" unit="g/cm3"/><atom value="14" Z="7"/></material>)",
        R"(<box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>)",
        R"(<volume name="World">
      <materialref ref="Air"/><solidref ref="world_box"/>
      <physvol><volumeref ref="NonexistentVolume"/></physvol>
    </volume>)", "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error);
}


TEST(GDMLParser, BuildVolumeMissingSolidThrows) {
    std::string gdml = GDMLFull("", R"(
    <material name="Air"><D value="0.00129" unit="g/cm3"/><atom value="14" Z="7"/></material>)",
        R"(<box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>)",
        R"(<volume name="World">
      <materialref ref="Air"/><solidref ref="nonexistent_solid"/>
    </volume>)", "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error);
}


TEST(GDMLParser, ReplicavolXAxis) {
    std::string gdml = GDMLFull("", "",
        R"(<box name="mother_box" x="800" y="100" z="100" lunit="mm"/>
    <box name="child_box" x="100" y="100" z="100" lunit="mm"/>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>)",
        R"(<volume name="ChildVol">
      <materialref ref=""/><solidref ref="child_box"/>
    </volume>
    <volume name="MotherVol">
      <materialref ref=""/><solidref ref="mother_box"/>
      <replicavol number="8">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction x="1"/><width value="100" unit="mm"/><offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
    <volume name="World">
      <materialref ref=""/><solidref ref="world_box"/>
      <physvol><volumeref ref="MotherVol"/></physvol>
    </volume>)", "World");

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
    std::string gdml = GDMLFull("", "",
        R"(<box name="mother_box" x="100" y="300" z="100" lunit="mm"/>
    <box name="child_box" x="100" y="100" z="100" lunit="mm"/>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>)",
        R"(<volume name="ChildVol">
      <materialref ref=""/><solidref ref="child_box"/>
    </volume>
    <volume name="MotherVol">
      <materialref ref=""/><solidref ref="mother_box"/>
      <replicavol number="3">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction y="1"/><width value="10" unit="cm"/><offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
    <volume name="World">
      <materialref ref=""/><solidref ref="world_box"/>
    </volume>)", "World");

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


TEST(GDMLParser, ReplicavolArbitraryAxis) {
    std::string gdml = GDMLFull("", "",
        R"(<box name="mother_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="child_box" x="50" y="50" z="50" lunit="mm"/>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>)",
        R"(<volume name="ChildVol">
      <materialref ref=""/><solidref ref="child_box"/>
    </volume>
    <volume name="ZVol">
      <materialref ref=""/><solidref ref="mother_box"/>
      <replicavol number="4">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction z="1"/><width value="200" unit="mm"/><offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
    <volume name="DiagVol">
      <materialref ref=""/><solidref ref="mother_box"/>
      <replicavol number="3">
        <volumeref ref="ChildVol"/>
        <replicate_along_axis>
          <direction x="1" y="1"/><width value="100" unit="mm"/><offset value="0" unit="mm"/>
        </replicate_along_axis>
      </replicavol>
    </volume>
    <volume name="World">
      <materialref ref=""/><solidref ref="world_box"/>
    </volume>)", "World");

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


TEST(GDMLParser, AssemblyPlacementGeometry) {
    // Assembly placed at (200,0,0) mm in World, containing a child box
    // at local position (100,0,0) mm. The child's global position should
    // be (300,0,0) mm = (0.3,0,0) m.
    // Child box half-width = 50mm = 0.05m, so its extent centered at (0.3,0,0)
    // is [0.25, 0.35] x [-0.05, 0.05] x [-0.05, 0.05].
    std::string gdml = GDMLFull("", kAirIron,
        R"(<box name="world_box" x="4000" y="4000" z="4000" lunit="mm"/>
    <box name="child_box" x="50" y="50" z="50" lunit="mm"/>)",
        R"(<volume name="IronBlock">
      <materialref ref="Iron"/><solidref ref="child_box"/>
    </volume>
    <assembly name="Assy">
      <physvol><volumeref ref="IronBlock"/>
        <position name="assy_child_pos" x="100" y="0" z="0" unit="mm"/>
      </physvol>
    </assembly>
    <volume name="World">
      <materialref ref="Air"/><solidref ref="world_box"/>
      <physvol><volumeref ref="Assy"/>
        <position name="assy_pos" x="200" y="0" z="0" unit="mm"/>
      </physvol>
    </volume>)", "World");

    auto dm = LoadGDMLFromString(gdml);

    // Child center at global (0.3, 0, 0) m should be Iron
    EXPECT_EQ(MaterialAt(*dm, 0.3, 0.0, 0.0), "Iron")
        << "Assembly child center at (200+100, 0, 0) mm = (0.3, 0, 0) m should be Iron";

    // Assembly position (0.2, 0, 0) m has no child solid there -- should be Air
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.0, 0.0), "Air")
        << "Assembly origin at (0.2, 0, 0) m has no child geometry, should be Air";

    // World center should be Air
    EXPECT_EQ(MaterialAt(*dm, 0.0, 0.0, 0.0), "Air")
        << "World center should be Air";
}


TEST(GDMLParser, AssemblyPlacementMultiAxis) {
    std::string gdml = GDMLFull("", kAirIron,
        R"(<box name="world_box" x="4000" y="4000" z="4000" lunit="mm"/>
    <box name="child_box" x="50" y="50" z="50" lunit="mm"/>)",
        R"(<volume name="IronBlock">
      <materialref ref="Iron"/><solidref ref="child_box"/>
    </volume>
    <assembly name="Assy">
      <physvol><volumeref ref="IronBlock"/>
        <position name="assy_child_pos" x="0" y="100" z="0" unit="mm"/>
      </physvol>
    </assembly>
    <volume name="World">
      <materialref ref="Air"/><solidref ref="world_box"/>
      <physvol><volumeref ref="Assy"/>
        <position name="assy_pos" x="200" y="0" z="0" unit="mm"/>
      </physvol>
    </volume>)", "World");

    auto dm = LoadGDMLFromString(gdml);

    // Child global position: assy at (0.2,0,0) + child at (0,0.1,0) = (0.2,0.1,0) m
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.1, 0.0), "Iron")
        << "Child at (200, 100, 0) mm should be Iron";

    // Assembly origin has no child: just Air
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.0, 0.0), "Air")
        << "Assembly origin should be Air (no child there)";

    // The original x-only position from old test should be Air
    EXPECT_EQ(MaterialAt(*dm, 0.3, 0.0, 0.0), "Air")
        << "(300, 0, 0) mm should be Air -- child is at y=100mm";
}


TEST(GDMLParser, AssemblyAsWorldVolume) {
    // Use an assembly as the world volume -- children placed, no bounding box
    std::string gdml = GDMLFull("",
        kAirIron,
        R"(<box name="block" x="100" y="100" z="100" lunit="mm"/>)",
        R"(<volume name="IronBlock">
             <materialref ref="Iron"/>
             <solidref ref="block"/>
           </volume>
           <assembly name="World">
             <physvol>
               <volumeref ref="IronBlock"/>
               <position name="p1" x="200" y="0" z="0" unit="mm"/>
             </physvol>
             <physvol>
               <volumeref ref="IronBlock"/>
               <position name="p2" x="-200" y="0" z="0" unit="mm"/>
             </physvol>
           </assembly>)",
        "World");

    auto dm = LoadGDMLFromString(gdml);

    // Children placed at +/-200mm = +/-0.2m should be Iron
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0, 0), "Iron");
    EXPECT_EQ(MaterialAt(*dm, -0.2, 0, 0), "Iron");

    // Points between the blocks have no geometry (assembly has no solid)
    EXPECT_NE(MaterialAt(*dm, 0.2, 0, 0), MaterialAt(*dm, 0, 0, 0))
        << "Assembly world should not fill the space between children";
}


// Missing positionref/rotationref should use
// identity defaults (origin position, identity rotation), not throw.
// A physvol with an undefined ref should still parse and place the
// child at the parent origin.
TEST(GDMLParser, MissingPositionRefThrows) {
    std::string gdml = GDMLFull("",
        kAirIron,
        R"(
        <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="inner_box" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="Inner">
          <materialref ref="Iron"/>
          <solidref ref="inner_box"/>
        </volume>
        <volume name="World">
          <materialref ref="Air"/>
          <solidref ref="world_box"/>
          <physvol>
            <volumeref ref="Inner"/>
            <positionref ref="undefined_pos"/>
          </physvol>
        </volume>
        )", "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error)
        << "Unresolved positionref should throw, not silently default to origin";
}

TEST(GDMLParser, MissingRotationRefThrows) {
    std::string gdml = GDMLFull("",
        kAirIron,
        R"(
        <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="inner_box" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="Inner">
          <materialref ref="Iron"/>
          <solidref ref="inner_box"/>
        </volume>
        <volume name="World">
          <materialref ref="Air"/>
          <solidref ref="world_box"/>
          <physvol>
            <volumeref ref="Inner"/>
            <rotationref ref="undefined_rot"/>
          </physvol>
        </volume>
        )", "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error)
        << "Unresolved rotationref should throw, not silently default to identity";
}


// A misspelled child element in <physvol> (e.g. <positoinref> instead of
// <positionref>) should throw, not silently produce identity placement.
TEST(GDMLParser, UnrecognizedPhysvolChildThrows) {
    std::string gdml = GDMLFull("",
        kAirIron,
        R"(
        <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="inner_box" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="Inner">
          <materialref ref="Iron"/>
          <solidref ref="inner_box"/>
        </volume>
        <volume name="World">
          <materialref ref="Air"/>
          <solidref ref="world_box"/>
          <physvol>
            <volumeref ref="Inner"/>
            <positoinref ref="some_pos"/>
          </physvol>
        </volume>
        )", "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error)
        << "Misspelled <positoinref> should throw, not silently use identity";
}

TEST(GDMLParser, UnrecognizedPhysvolRotationChildThrows) {
    std::string gdml = GDMLFull("",
        kAirIron,
        R"(
        <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="inner_box" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="Inner">
          <materialref ref="Iron"/>
          <solidref ref="inner_box"/>
        </volume>
        <volume name="World">
          <materialref ref="Air"/>
          <solidref ref="world_box"/>
          <physvol>
            <volumeref ref="Inner"/>
            <rotatoinref ref="some_rot"/>
          </physvol>
        </volume>
        )", "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error)
        << "Misspelled <rotatoinref> should throw, not silently use identity";
}

// Geant4 GDML exports may include <auxiliary> metadata inside <physvol>.
// These are valid extension elements that we skip, not typos.
TEST(GDMLParser, PhysvolAuxiliaryElementAccepted) {
    std::string gdml = GDMLFull("",
        kAirIron,
        R"(
        <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="inner_box" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="Inner">
          <materialref ref="Iron"/>
          <solidref ref="inner_box"/>
        </volume>
        <volume name="World">
          <materialref ref="Air"/>
          <solidref ref="world_box"/>
          <physvol>
            <volumeref ref="Inner"/>
            <auxiliary auxtype="SensDet" auxvalue="Tracker"/>
          </physvol>
        </volume>
        )", "World");

    EXPECT_NO_THROW(LoadGDMLFromString(gdml))
        << "<auxiliary> inside physvol should be accepted (Geant4 extension)";
}

// Physvol <file> references: the child volume is defined in an external
// GDML file. The parser should load and merge that file's definitions.
TEST(GDMLParser, PhysvolFileReference) {
    // Write a sub-GDML file that defines a box volume
    std::string sub_path = TempPath("siren_gdml_sub.gdml");
    {
        std::ofstream f(sub_path);
        f << R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids>
    <box name="sub_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="SubVol">
      <materialref ref="Iron"/>
      <solidref ref="sub_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubVol"/></setup>
</gdml>)";
    }

    // Main GDML references the sub-file via <file>
    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol>
        <file name=")" + sub_path + R"(" volname="SubVol"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_gdml_fileref.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(main_path);
    std::remove(main_path.c_str());
    std::remove(sub_path.c_str());

    // The sub-volume's iron box should be at the origin
    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "Iron")
        << "Physvol <file> should load the external volume";
}

// When volname is omitted, the sub-file's world volume (from <setup>)
// should be used as the child. This is the common case: 446 of 447
// <file> tags in real GDML corpora omit volname.
TEST(GDMLParser, PhysvolFileReferenceNoVolname) {
    std::string sub_path = TempPath("siren_gdml_sub_novolname.gdml");
    {
        std::ofstream f(sub_path);
        f << R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids>
    <box name="sub_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="SubWorld">
      <materialref ref="Iron"/>
      <solidref ref="sub_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubWorld"/></setup>
</gdml>)";
    }

    // No volname attribute -- should use the sub-file's world volume
    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol>
        <file name=")" + sub_path + R"("/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_gdml_fileref_novolname.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(main_path);
    std::remove(main_path.c_str());
    std::remove(sub_path.c_str());

    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "Iron")
        << "Physvol <file> without volname should use the sub-file world volume";
}

// When both <element name="Si"> and <material name="Si"> exist, volumes
// referencing <materialref ref="Si"/> must resolve to the material, not
// the element.
TEST(GDMLParser, MaterialElementNameCollision) {
    std::string gdml = GDMLFull("",
        R"(
        <element name="Si" Z="14"><atom value="28.085"/></element>
        <material name="Si" Z="14"><D value="2.33" unit="g/cm3"/><atom value="28.085"/></material>
        <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
        )",
        R"(
        <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="si_box" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="SiVol">
          <materialref ref="Si"/>
          <solidref ref="si_box"/>
        </volume>
        <volume name="World">
          <materialref ref="Air"/>
          <solidref ref="world_box"/>
          <physvol><volumeref ref="SiVol"/></physvol>
        </volume>
        )", "World");

    auto dm = LoadGDMLFromString(gdml);
    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "Si")
        << "materialref 'Si' should resolve to the material, not the element";
}

// A GDML file without a <setup> section (materials-only include file)
// should load without error. The world volume is empty so no sectors
// are built, but parsing succeeds.
TEST(GDMLParser, NoSetupSectionLoads) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids/>
  <structure/>
</gdml>)";

    EXPECT_NO_THROW(LoadGDMLFromString(gdml))
        << "GDML without <setup> should load (materials-only include file)";
}


// ======================================================================
// <file> element extensions: as_assembly and auto-collision handling
// ======================================================================

// Helper: write a sub-GDML file to a temp path and return the path.
static std::string WriteSubGDML(std::string const & name, std::string const & content) {
    std::string path = TempPath(name);
    { std::ofstream f(path); f << content; }
    return path;
}

// Two sub-files define a material with the same name but DIFFERENT definitions.
// The auto-collision handler should rename the second one so both are usable.
TEST(GDMLParser, FileAutoRenameOnCollision) {
    // Sub A: "Mat" is iron (Z=26)
    std::string sub_a = WriteSubGDML("siren_gdml_collision_a.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Mat" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids><box name="box_a" x="100" y="100" z="100" lunit="mm"/></solids>
  <structure>
    <volume name="VolA"><materialref ref="Mat"/><solidref ref="box_a"/></volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="VolA"/></setup>
</gdml>)");

    // Sub B: "Mat" is aluminum (Z=13) -- different definition, same name
    std::string sub_b = WriteSubGDML("siren_gdml_collision_b.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Mat" Z="13"><D value="2.70" unit="g/cm3"/><atom value="26.982"/></material>
  </materials>
  <solids><box name="box_b" x="100" y="100" z="100" lunit="mm"/></solids>
  <structure>
    <volume name="VolB"><materialref ref="Mat"/><solidref ref="box_b"/></volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="VolB"/></setup>
</gdml>)");

    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="pos_a" x="-200" y="0" z="0" unit="mm"/>
    <position name="pos_b" x="200" y="0" z="0" unit="mm"/>
  </define>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
  </materials>
  <solids><box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/></solids>
  <structure>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol><file name=")" + sub_a + R"("/><positionref ref="pos_a"/></physvol>
      <physvol><file name=")" + sub_b + R"("/><positionref ref="pos_b"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_gdml_collision_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(main_path);
    std::remove(main_path.c_str());
    std::remove(sub_a.c_str());
    std::remove(sub_b.c_str());

    // Sub A uses the original "Mat" name (first wins for identical keys)
    // Sub B's "Mat" was auto-renamed because it differs
    std::string mat_a = MaterialAt(*dm, -0.2, 0, 0);
    std::string mat_b = MaterialAt(*dm, 0.2, 0, 0);
    EXPECT_EQ(mat_a, "Mat") << "First sub-file keeps the original name";
    EXPECT_NE(mat_b, "Mat") << "Second sub-file's conflicting name should be renamed";
    EXPECT_NE(mat_a, mat_b) << "The two materials should be distinct";
}

// Two sub-files define a material with the same name and IDENTICAL definitions.
// The auto-collision handler should silently reuse the first definition.
TEST(GDMLParser, FileIdenticalDefinitionReused) {
    // Both sub-files define "Iron" identically
    std::string iron_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids><box name="the_box" x="100" y="100" z="100" lunit="mm"/></solids>
  <structure>
    <volume name="TheVol"><materialref ref="Iron"/><solidref ref="the_box"/></volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="TheVol"/></setup>
</gdml>)";

    std::string sub_a = WriteSubGDML("siren_gdml_ident_a.gdml", iron_gdml);
    std::string sub_b = WriteSubGDML("siren_gdml_ident_b.gdml", iron_gdml);

    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="pos_a" x="-200" y="0" z="0" unit="mm"/>
    <position name="pos_b" x="200" y="0" z="0" unit="mm"/>
  </define>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
  </materials>
  <solids><box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/></solids>
  <structure>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol><file name=")" + sub_a + R"("/><positionref ref="pos_a"/></physvol>
      <physvol><file name=")" + sub_b + R"("/><positionref ref="pos_b"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_gdml_ident_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(main_path);
    std::remove(main_path.c_str());
    std::remove(sub_a.c_str());
    std::remove(sub_b.c_str());

    // Both should resolve to the same "Iron" material (no rename needed)
    EXPECT_EQ(MaterialAt(*dm, -0.2, 0, 0), "Iron");
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0, 0), "Iron");
}

// as_assembly="true" on <file> should treat the referenced volume as an
// assembly: no sector created for it, but its children are placed.
TEST(GDMLParser, FileAsAssembly) {
    std::string sub = WriteSubGDML("siren_gdml_asasm_sub.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="SubAir" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids>
    <box name="sub_world_box" x="500" y="500" z="500" lunit="mm"/>
    <box name="iron_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="IronVol">
      <materialref ref="Iron"/>
      <solidref ref="iron_box"/>
    </volume>
    <volume name="SubWorld">
      <materialref ref="SubAir"/>
      <solidref ref="sub_world_box"/>
      <physvol><volumeref ref="IronVol"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubWorld"/></setup>
</gdml>)");

    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Rock" Z="14"><D value="2.65" unit="g/cm3"/><atom value="28.085"/></material>
  </materials>
  <solids><box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/></solids>
  <structure>
    <volume name="World">
      <materialref ref="Rock"/>
      <solidref ref="world_box"/>
      <physvol><file name=")" + sub + R"(" as_assembly="true"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_gdml_asasm_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(main_path);
    std::remove(main_path.c_str());
    std::remove(sub.c_str());

    // Iron child should be placed
    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "Iron")
        << "Iron child should be placed through the assembly";
    // Outside iron box but inside sub_world_box: Rock (no SubAir sector)
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0, 0), "Rock")
        << "as_assembly should prevent SubAir sector from being created";
}

// Without as_assembly, the sub-file's world volume creates a sector
// (the air bubble). Verify this so the as_assembly test is meaningful.
TEST(GDMLParser, FileWithoutAsAssemblyCreatesWorldSector) {
    std::string sub = WriteSubGDML("siren_gdml_noasm_sub.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="SubAir" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids>
    <box name="sub_world_box" x="500" y="500" z="500" lunit="mm"/>
    <box name="iron_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="IronVol">
      <materialref ref="Iron"/>
      <solidref ref="iron_box"/>
    </volume>
    <volume name="SubWorld">
      <materialref ref="SubAir"/>
      <solidref ref="sub_world_box"/>
      <physvol><volumeref ref="IronVol"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubWorld"/></setup>
</gdml>)");

    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Rock" Z="14"><D value="2.65" unit="g/cm3"/><atom value="28.085"/></material>
  </materials>
  <solids><box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/></solids>
  <structure>
    <volume name="World">
      <materialref ref="Rock"/>
      <solidref ref="world_box"/>
      <physvol><file name=")" + sub + R"("/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_gdml_noasm_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(main_path);
    std::remove(main_path.c_str());
    std::remove(sub.c_str());

    EXPECT_EQ(MaterialAt(*dm, 0.2, 0, 0), "SubAir")
        << "Without as_assembly, the sub-world air sector should exist";
}

// Combine as_assembly with a collision: sub-file defines "Air" differently
// from the parent. The sub's Air should be auto-renamed AND the sub-world
// should be treated as assembly.
TEST(GDMLParser, FileAsAssemblyWithCollision) {
    // Sub-file: "Air" has different density than parent
    std::string sub = WriteSubGDML("siren_gdml_asmcoll_sub.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Air" Z="7"><D value="0.00100" unit="g/cm3"/><atom value="14.007"/></material>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids>
    <box name="sub_world" x="500" y="500" z="500" lunit="mm"/>
    <box name="child_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="Child">
      <materialref ref="Iron"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="SubWorld">
      <materialref ref="Air"/>
      <solidref ref="sub_world"/>
      <physvol><volumeref ref="Child"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubWorld"/></setup>
</gdml>)");

    // Parent defines "Air" with different density
    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
    <material name="Rock" Z="14"><D value="2.65" unit="g/cm3"/><atom value="28.085"/></material>
  </materials>
  <solids><box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/></solids>
  <structure>
    <volume name="World">
      <materialref ref="Rock"/>
      <solidref ref="world_box"/>
      <physvol><file name=")" + sub + R"(" as_assembly="true"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_gdml_asmcoll_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    dm->LoadGDML(main_path);
    std::remove(main_path.c_str());
    std::remove(sub.c_str());

    // Iron child placed
    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "Iron");
    // Outside child: Rock (assembly suppresses the sub-world sector)
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0, 0), "Rock");
}

// ======================================================================
// <file> scope isolation: parent must not reference sub-file internals
// ======================================================================

// Parent directly references a sub-file internal volume via <volumeref>.
// This should be rejected: only the sub-file's world volume is accessible.
TEST(GDMLParser, FileScopeIsolationVolumeref) {
    std::string sub = WriteSubGDML("siren_scope_sub.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids>
    <box name="sub_world_box" x="500" y="500" z="500" lunit="mm"/>
    <box name="internal_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="InternalVol">
      <materialref ref="Iron"/>
      <solidref ref="internal_box"/>
    </volume>
    <volume name="SubWorld">
      <materialref ref="Iron"/>
      <solidref ref="sub_world_box"/>
      <physvol><volumeref ref="InternalVol"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubWorld"/></setup>
</gdml>)");

    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
  </materials>
  <solids><box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/></solids>
  <structure>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol><file name=")" + sub + R"("/></physvol>
      <physvol><volumeref ref="InternalVol"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_scope_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    EXPECT_THROW(dm->LoadGDML(main_path), std::runtime_error)
        << "Parent should not reference sub-file internal volume";
    std::remove(main_path.c_str());
    std::remove(sub.c_str());
}

// Same violation but with the <volumeref> BEFORE the <file> tag in XML
// order. The post-parse validation must still catch it.
TEST(GDMLParser, FileScopeIsolationReverseOrder) {
    std::string sub = WriteSubGDML("siren_scope_rev_sub.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids>
    <box name="sub_world_box" x="500" y="500" z="500" lunit="mm"/>
    <box name="internal_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="InternalVol">
      <materialref ref="Iron"/>
      <solidref ref="internal_box"/>
    </volume>
    <volume name="SubWorld">
      <materialref ref="Iron"/>
      <solidref ref="sub_world_box"/>
      <physvol><volumeref ref="InternalVol"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubWorld"/></setup>
</gdml>)");

    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
  </materials>
  <solids><box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/></solids>
  <structure>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <physvol><volumeref ref="InternalVol"/></physvol>
      <physvol><file name=")" + sub + R"("/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_scope_rev_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    EXPECT_THROW(dm->LoadGDML(main_path), std::runtime_error)
        << "Scope isolation should catch references regardless of XML order";
    std::remove(main_path.c_str());
    std::remove(sub.c_str());
}

// Parent defines a volume with the same name as a sub-file internal
// volume. The parent's definition should shadow the file-internal one:
// no scope error, and the reference resolves to the parent's volume.
TEST(GDMLParser, FileScopeParentShadowsInternal) {
    // Sub-file defines "SharedVol" internally (not placed anywhere)
    std::string sub = WriteSubGDML("siren_scope_shadow_sub.gdml", R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
    <material name="Copper" Z="29"><D value="8.96" unit="g/cm3"/><atom value="63.546"/></material>
  </materials>
  <solids>
    <box name="sub_world_box" x="500" y="500" z="500" lunit="mm"/>
    <box name="shared_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="SharedVol">
      <materialref ref="Copper"/>
      <solidref ref="shared_box"/>
    </volume>
    <volume name="SubWorld">
      <materialref ref="Iron"/>
      <solidref ref="sub_world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="SubWorld"/></setup>
</gdml>)");

    // Parent also defines "SharedVol" (Lead) and references it
    std::string main_gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="pos_sv" x="300" y="0" z="0" unit="mm"/>
  </define>
  <materials>
    <material name="Rock" Z="14"><D value="2.65" unit="g/cm3"/><atom value="28.085"/></material>
    <material name="Lead" Z="82"><D value="11.35" unit="g/cm3"/><atom value="207.2"/></material>
  </materials>
  <solids>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
    <box name="parent_box" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="SharedVol">
      <materialref ref="Lead"/>
      <solidref ref="parent_box"/>
    </volume>
    <volume name="World">
      <materialref ref="Rock"/>
      <solidref ref="world_box"/>
      <physvol><file name=")" + sub + R"("/></physvol>
      <physvol><volumeref ref="SharedVol"/><positionref ref="pos_sv"/></physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)";

    std::string main_path = TempPath("siren_scope_shadow_main.gdml");
    { std::ofstream f(main_path); f << main_gdml; }
    auto dm = std::make_unique<DetectorModel>();
    EXPECT_NO_THROW(dm->LoadGDML(main_path))
        << "Parent-defined volume should shadow file-internal name";
    std::remove(main_path.c_str());
    std::remove(sub.c_str());

    EXPECT_EQ(MaterialAt(*dm, 0.3, 0, 0), "Lead")
        << "SharedVol should resolve to the parent's definition";
}
