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

    // The child box is 400x100x100 mm (half-widths: 0.4x0.1x0.1 m)
    // Placed at x=200mm=0.2m with z-rotation of 90 degrees.
    // After 90-deg z-rotation, the long axis (originally x) becomes the y axis.
    // So the child occupies roughly: x in [0.1, 0.3], y in [-0.4, 0.4]

    // The child sector should exist
    auto const & sectors = dm->GetSectors();
    EXPECT_GE(sectors.size(), 3u)
        << "Should have at least 3 sectors (UNIVERSE + world + rotated child)";

    // A point at the child center (0.2, 0, 0) should be Iron
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.0, 0.0), "Iron")
        << "Center of rotated child at (0.2,0,0) should be Iron";

    // After rotation, the long axis is along y. A point at (0.2, 0.3, 0) should
    // still be inside (within the rotated half-width of 0.4 along y).
    EXPECT_EQ(MaterialAt(*dm, 0.2, 0.3, 0.0), "Iron")
        << "Point at (0.2, 0.3, 0) should be inside the rotated child (long axis now along y)";

    // A point far from the child should be in Air
    EXPECT_EQ(MaterialAt(*dm, 1.0, 1.0, 0.0), "Air")
        << "Point at (1.0,1.0,0) should be in Air (outside rotated child)";
}


TEST(GDMLParser, CompoundRotationGroundTruth) {
    std::string gdml = GDMLFull("", kAirIron,
        R"(<box name="world_box" x="5" y="5" z="5" lunit="m"/>
    <box name="child_box" x="0.5" y="0.2" z="0.1" lunit="m"/>)",
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
    // Use an assembly as the world volume. This is the expected integration
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
