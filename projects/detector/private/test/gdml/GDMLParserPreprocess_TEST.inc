// Included by GDMLParser_TEST.cxx. Do not compile directly.

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
    EXPECT_NEAR(data.solids["test_box"]->GetBoundingBox().max_corner.GetX(), 0.1, tol);
    ASSERT_GE(data.warnings.size(), 1u);
    bool found_root_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("non-standard root element") != std::string::npos) found_root_warning = true;
    }
    EXPECT_TRUE(found_root_warning) << "Expected warning about non-standard root element";
}

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

TEST(GDMLParser, VariableElementStoredAsConstant) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="myvar" value="10"/>)");
    EXPECT_NEAR(data.constants["myvar"], 10.0, 1e-9);
}

TEST(GDMLParser, LoopElementExpansion) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="0"/>
    <loop for="i" from="0" to="3" step="1">
      <constant name="val_[i]" value="i*10"/>
    </loop>)");

    EXPECT_NEAR(data.constants["val_0"], 0.0, 1e-9);
    EXPECT_NEAR(data.constants["val_1"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["val_2"], 20.0, 1e-9);
    EXPECT_NEAR(data.constants["val_3"], 30.0, 1e-9);
}

TEST(GDMLParser, LoopOmittedFromUsesVariableValue) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="x" value="3"/>
    <loop for="x" to="5" step="1">
      <constant name="c_[x]" value="x*x"/>
    </loop>)");

    EXPECT_EQ(data.constants.count("c_3"), 1u) << "Should have c_3 (x=3)";
    EXPECT_EQ(data.constants.count("c_4"), 1u) << "Should have c_4 (x=4)";
    EXPECT_EQ(data.constants.count("c_5"), 1u) << "Should have c_5 (x=5)";
    EXPECT_EQ(data.constants.count("c_0"), 0u) << "Should NOT have c_0 (loop starts at 3)";
    EXPECT_NEAR(data.constants["c_3"], 9.0, 1e-9);
    EXPECT_NEAR(data.constants["c_4"], 16.0, 1e-9);
    EXPECT_NEAR(data.constants["c_5"], 25.0, 1e-9);
}

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

TEST(GDMLParser, LoopInStructure) {
    // Merged: checks both DetectorModel material lookup AND parsed positions.
    std::string gdml = GDMLFull(
        R"(
    <variable name="i" value="0"/>)",
        R"(
    <material name="TestMat" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.01"/></material>)",
        R"(
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
    <box name="child_box" x="10" y="10" z="10" lunit="mm"/>)",
        R"(
    <volume name="Child"><materialref ref="TestMat"/><solidref ref="child_box"/></volume>
    <volume name="World">
      <materialref ref="TestMat"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="3" step="1">
        <physvol>
          <volumeref ref="Child"/>
          <position name="lpos_[i]" x="i*100" y="0" z="0" unit="mm"/>
        </physvol>
      </loop>
    </volume>)",
        "World");

    auto dm = LoadGDMLFromString(gdml);
    int child_count = 0;
    for(auto const & sector : dm->GetSectors()) {
        if(sector.name.find("Child") != std::string::npos) child_count++;
    }
    EXPECT_EQ(child_count, 4) << "Loop should produce 4 child volume placements";

    GDMLData data = ParseGDMLString(gdml);
    auto it = data.volumes.find("World");
    ASSERT_NE(it, data.volumes.end());
    ASSERT_EQ(it->second.children.size(), 4u)
        << "Loop i=0..3 should produce exactly 4 physvol placements";
    for(size_t idx = 0; idx < 4; ++idx) {
        double expected_x = idx * 0.1; // i*100 mm -> i*0.1 m
        auto const & child = it->second.children[idx];
        EXPECT_NEAR(child.position.GetX(), expected_x, 1e-9)
            << "Child " << idx << " x position should be " << expected_x << " m";
        EXPECT_NEAR(child.position.GetY(), 0.0, 1e-9);
        EXPECT_NEAR(child.position.GetZ(), 0.0, 1e-9);
    }
}

TEST(GDMLParser, LoopWordBoundarySafety) {
    // Volume name "ChildVol" contains "i" but must NOT be corrupted
    // by a loop with variable "i". Word-boundary matching prevents this.
    std::string gdml = GDMLFull(
        R"(
    <variable name="i" value="0"/>)",
        R"(
    <material name="Air"><D value="0.00129" unit="g/cm3"/><atom value="14" Z="7"/></material>)",
        R"(
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <box name="child_box" x="10" y="10" z="10" lunit="mm"/>)",
        R"(
    <volume name="ChildVol"><materialref ref="Air"/><solidref ref="child_box"/></volume>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="1" step="1">
        <physvol>
          <volumeref ref="ChildVol"/>
          <position name="p_[i]" x="i*100" y="0" z="0" unit="mm"/>
        </physvol>
      </loop>
    </volume>)",
        "World");

    auto dm = LoadGDMLFromString(gdml);
    int child_count = 0;
    for(auto const & sector : dm->GetSectors()) {
        if(sector.name.find("ChildVol") != std::string::npos) child_count++;
    }
    EXPECT_EQ(child_count, 2) << "Should have 2 ChildVol placements (i=0,1)";
}

TEST(GDMLParser, LoopVariablePersistsAfterLoop) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="0"/>
    <loop for="i" from="1" to="5" step="1">
      <constant name="c_[i]" value="i*10"/>
    </loop>
    <constant name="final_i" value="i"/>)");

    EXPECT_NEAR(data.constants["i"], 5.0, 1e-9);
    EXPECT_NEAR(data.constants["final_i"], 5.0, 1e-9);
    EXPECT_NEAR(data.constants["c_1"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["c_5"], 50.0, 1e-9);
}

TEST(GDMLParser, LoopVariableDocumentOrder) {
    // Before-loop sees pre-loop value; after-loop sees post-loop (final) value.
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="10"/>
    <constant name="before_loop" value="i"/>
    <loop for="i" from="1" to="5" step="1">
      <constant name="c_[i]" value="i*10"/>
    </loop>
    <constant name="after_loop" value="i"/>)");

    EXPECT_NEAR(data.constants["before_loop"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["after_loop"], 5.0, 1e-9);
    EXPECT_NEAR(data.constants["c_1"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["c_5"], 50.0, 1e-9);
}

TEST(GDMLParser, LoopVariablePreexistingValue) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="idx" value="99"/>
    <loop for="idx" from="1" to="3" step="1">
      <constant name="v_[idx]" value="idx"/>
    </loop>
    <constant name="after_loop" value="idx"/>)");

    EXPECT_NEAR(data.constants["after_loop"], 3.0, 1e-9);
}

TEST(GDMLParser, EntityExpansion) {
    std::string inc_file = "/tmp/siren_gdml_entity_inc.gdml";
    { std::ofstream f(inc_file); f << R"(
    <box name="included_box" x="2000" y="2000" z="2000" lunit="mm"/>
)"; }

    std::string main_file = "/tmp/siren_gdml_entity_test.gdml";
    { std::ofstream f(main_file); f << R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY solids_inc SYSTEM "siren_gdml_entity_inc.gdml">
]>
<gdml>
  <define/><materials/>
  <solids>&solids_inc;</solids>
  <structure><volume name="World"><materialref ref=""/><solidref ref="included_box"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)"; }

    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    std::remove(inc_file.c_str());

    auto geo = RequireSolid(data, "included_box");
    ASSERT_NE(geo, nullptr);
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));
    EXPECT_FALSE(IsInside(geo, 3.0, 0, 0));
}

TEST(GDMLParser, NestedEntityExpansion) {
    std::string inner_file = "/tmp/siren_gdml_inner.gdml";
    { std::ofstream f(inner_file); f << R"(
    <constant name="inner_val" value="42"/>
)"; }

    std::string outer_file = "/tmp/siren_gdml_outer.gdml";
    { std::ofstream f(outer_file); f << R"(<?xml version="1.0"?>
<!DOCTYPE gdml_fragment [
  <!ENTITY inner SYSTEM "siren_gdml_inner.gdml">
]>
    <constant name="outer_val" value="99"/>
    &inner;
)"; }

    std::string main_file = "/tmp/siren_gdml_nested_test.gdml";
    { std::ofstream f(main_file); f << R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY defs SYSTEM "siren_gdml_outer.gdml">
]>
<gdml>
  <define>&defs;</define>
  <materials/>
  <solids><box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/></solids>
  <structure><volume name="World"><materialref ref=""/><solidref ref="world_box"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)"; }

    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    std::remove(outer_file.c_str());
    std::remove(inner_file.c_str());

    EXPECT_NEAR(data.constants["outer_val"], 99.0, 1e-9);
    EXPECT_NEAR(data.constants["inner_val"], 42.0, 1e-9);
}

// --- Tests for stack-based loop expansion and document-order resolution ---

TEST(GDMLParser, NestedLoops) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="0"/>
    <variable name="j" value="0"/>
    <loop for="i" from="1" to="3" step="1">
      <loop for="j" from="1" to="2" step="1">
        <constant name="c_[i]_[j]" value="i*10+j"/>
      </loop>
    </loop>
    <constant name="final_i" value="i"/>
    <constant name="final_j" value="j"/>)");

    EXPECT_NEAR(data.constants["c_1_1"], 11.0, 1e-9);
    EXPECT_NEAR(data.constants["c_1_2"], 12.0, 1e-9);
    EXPECT_NEAR(data.constants["c_2_1"], 21.0, 1e-9);
    EXPECT_NEAR(data.constants["c_2_2"], 22.0, 1e-9);
    EXPECT_NEAR(data.constants["c_3_1"], 31.0, 1e-9);
    EXPECT_NEAR(data.constants["c_3_2"], 32.0, 1e-9);
    EXPECT_NEAR(data.constants["final_i"], 3.0, 1e-9);
    EXPECT_NEAR(data.constants["final_j"], 2.0, 1e-9);
}

TEST(GDMLParser, NonExecutingLoop) {
    // Loop with from > to and positive step should not execute.
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="42"/>
    <loop for="i" from="10" to="5" step="1">
      <constant name="should_not_exist" value="999"/>
    </loop>
    <constant name="after" value="i"/>)");

    EXPECT_NEAR(data.constants["after"], 42.0, 1e-9);
    EXPECT_EQ(data.constants.find("should_not_exist"), data.constants.end());
}

TEST(GDMLParser, NegativeStepLoop) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="0"/>
    <loop for="i" from="5" to="1" step="-1">
      <constant name="d_[i]" value="i*2"/>
    </loop>
    <constant name="final_i" value="i"/>)");

    EXPECT_NEAR(data.constants["d_5"], 10.0, 1e-9);
    EXPECT_NEAR(data.constants["d_4"], 8.0, 1e-9);
    EXPECT_NEAR(data.constants["d_3"], 6.0, 1e-9);
    EXPECT_NEAR(data.constants["d_2"], 4.0, 1e-9);
    EXPECT_NEAR(data.constants["d_1"], 2.0, 1e-9);
    EXPECT_NEAR(data.constants["final_i"], 1.0, 1e-9);
}

TEST(GDMLParser, CumulativeLoopReferences) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="0"/>
    <constant name="sum_0" value="0"/>
    <loop for="i" from="1" to="4" step="1">
      <constant name="sum_[i]" value="sum_[i-1]+i"/>
    </loop>)");

    EXPECT_NEAR(data.constants["sum_1"], 1.0, 1e-9);
    EXPECT_NEAR(data.constants["sum_2"], 3.0, 1e-9);
    EXPECT_NEAR(data.constants["sum_3"], 6.0, 1e-9);
    EXPECT_NEAR(data.constants["sum_4"], 10.0, 1e-9);
}

TEST(GDMLParser, LoopGeneratingPositions) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="0"/>
    <constant name="spacing" value="100"/>
    <loop for="i" from="0" to="2" step="1">
      <position name="pos_[i]" x="i*spacing" y="0" z="0" unit="mm"/>
    </loop>)");

    ASSERT_TRUE(data.positions.count("pos_0") > 0);
    ASSERT_TRUE(data.positions.count("pos_1") > 0);
    ASSERT_TRUE(data.positions.count("pos_2") > 0);
    EXPECT_NEAR(data.positions["pos_0"].GetX(), 0.0, 1e-9);
    EXPECT_NEAR(data.positions["pos_1"].GetX(), 0.1, 1e-9);
    EXPECT_NEAR(data.positions["pos_2"].GetX(), 0.2, 1e-9);
}

TEST(GDMLParser, NestedLoopInStructure) {
    // Merged: checks both DetectorModel child count AND parsed positions.
    std::string gdml = GDMLFull(
        R"(
    <variable name="i" value="0"/>
    <variable name="j" value="0"/>)",
        R"(
    <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.01"/></material>)",
        R"(
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
    <box name="child_box" x="10" y="10" z="10" lunit="mm"/>)",
        R"(
    <volume name="Child"><materialref ref="Air"/><solidref ref="child_box"/></volume>
    <volume name="World">
      <materialref ref="Air"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="1" step="1">
        <loop for="j" from="0" to="1" step="1">
          <physvol>
            <volumeref ref="Child"/>
            <position name="nlpos_[i]_[j]" x="i*200" y="j*200" z="0" unit="mm"/>
          </physvol>
        </loop>
      </loop>
    </volume>)",
        "World");

    auto dm = LoadGDMLFromString(gdml);
    int child_count = 0;
    for(auto const & sector : dm->GetSectors()) {
        if(sector.name.find("Child") != std::string::npos) child_count++;
    }
    EXPECT_EQ(child_count, 4);

    GDMLData data = ParseGDMLString(gdml);
    auto it = data.volumes.find("World");
    ASSERT_NE(it, data.volumes.end());
    ASSERT_EQ(it->second.children.size(), 4u)
        << "Nested loops i=0..1, j=0..1 should produce 4 physvol placements";

    std::set<std::pair<double, double>> actual_positions;
    for(auto const & child : it->second.children) {
        double x = std::round(child.position.GetX() * 1e6) / 1e6;
        double y = std::round(child.position.GetY() * 1e6) / 1e6;
        actual_positions.insert({x, y});
        EXPECT_NEAR(child.position.GetZ(), 0.0, 1e-9);
    }

    std::set<std::pair<double, double>> expected_positions = {
        {0.0, 0.0}, {0.0, 0.2}, {0.2, 0.0}, {0.2, 0.2}
    };
    EXPECT_EQ(actual_positions, expected_positions)
        << "Nested loop should produce positions at (0,0), (0,0.2), (0.2,0), (0.2,0.2) meters";
}

TEST(GDMLParser, LoopVariableInInnerBounds) {
    GDMLData data = ParseGDMLWithDefine(R"(
    <variable name="i" value="0"/>
    <variable name="j" value="0"/>
    <loop for="i" from="1" to="3" step="1">
      <loop for="j" from="1" to="i" step="1">
        <constant name="t_[i]_[j]" value="1"/>
      </loop>
    </loop>)");

    EXPECT_TRUE(data.constants.count("t_1_1") > 0);
    EXPECT_EQ(data.constants.count("t_1_2"), 0u);
    EXPECT_TRUE(data.constants.count("t_2_1") > 0);
    EXPECT_TRUE(data.constants.count("t_2_2") > 0);
    EXPECT_EQ(data.constants.count("t_2_3"), 0u);
    EXPECT_TRUE(data.constants.count("t_3_1") > 0);
    EXPECT_TRUE(data.constants.count("t_3_2") > 0);
    EXPECT_TRUE(data.constants.count("t_3_3") > 0);
}

TEST(GDMLParser, EntityMissingFile) {
    // Merged: tests both unused and used missing entity references.

    // Part 1: Entity not referenced in body -- should parse successfully.
    std::string main_file = "/tmp/siren_gdml_missing_entity_test.gdml";
    { std::ofstream f(main_file); f << R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY missing SYSTEM "nonexistent_file_xyz.gdml">
]>
<gdml>
  <define/><materials/>
  <solids><box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/></solids>
  <structure><volume name="World"><materialref ref=""/><solidref ref="world_box"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)"; }

    GDMLData data_unused = ParseGDML(main_file);
    std::remove(main_file.c_str());
    EXPECT_TRUE(data_unused.solids.find("world_box") != data_unused.solids.end());

    // Part 2: Entity file missing but local content should still parse.
    std::string main_file2 = "/tmp/siren_gdml_missing_entity_used_test.gdml";
    { std::ofstream f(main_file2); f << R"(<?xml version="1.0"?>
<!DOCTYPE gdml [
  <!ENTITY missing SYSTEM "nonexistent_file_xyz.gdml">
]>
<gdml>
  <define><constant name="local_only" value="42"/></define>
  <materials/>
  <solids><box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/></solids>
  <structure><volume name="World"><materialref ref=""/><solidref ref="world_box"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)"; }

    GDMLData data_used = ParseGDML(main_file2);
    std::remove(main_file2.c_str());
    EXPECT_NEAR(data_used.constants["local_only"], 42.0, 1e-9);
}

TEST(GDMLParser, XIncludeExpansion) {
    std::string inc_file = "/tmp/siren_gdml_xinc.gdml";
    { std::ofstream f(inc_file); f << R"(
    <constant name="xinc_val" value="77"/>
)"; }

    std::string main_file = "/tmp/siren_gdml_xinc_test.gdml";
    { std::ofstream f(main_file); f << R"(<?xml version="1.0"?>
<gdml xmlns:xi="http://www.w3.org/2001/XInclude">
  <define>
    <xi:include href="siren_gdml_xinc.gdml"/>
    <constant name="local_val" value="88"/>
  </define>
  <materials/>
  <solids><box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/></solids>
  <structure><volume name="World"><materialref ref=""/><solidref ref="world_box"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)"; }

    GDMLData data = ParseGDML(main_file);
    std::remove(main_file.c_str());
    std::remove(inc_file.c_str());

    EXPECT_NEAR(data.constants["xinc_val"], 77.0, 1e-9);
    EXPECT_NEAR(data.constants["local_val"], 88.0, 1e-9);
}

TEST(GDMLParser, RealWorldGDMLSmoke) {
    std::string test_dir = std::string(__FILE__);
    auto pos = test_dir.rfind("projects/");
    ASSERT_NE(pos, std::string::npos)
        << "Could not locate repo root from __FILE__";
    test_dir = test_dir.substr(0, pos) + "test_gdml/dark_matter/cats/";

    struct TestCase { std::string filename; };
    std::vector<TestCase> cases = {
        {"simpleLArTPC.gdml"},
        {"DRSingleCrystal.gdml"},
        {"zen.gdml"},
    };

    for(auto const & tc : cases) {
        std::string path = test_dir + tc.filename;
        GDMLData data;
        EXPECT_NO_THROW(data = ParseGDML(path))
            << "ParseGDML should not throw for " << tc.filename;
        EXPECT_GE(data.solids.size(), 1u)
            << tc.filename << " should contain at least one solid";
        bool has_world = false;
        for(auto const & kv : data.volumes) {
            if(!kv.second.children.empty() || !kv.second.solid_ref.empty()) {
                has_world = true;
                break;
            }
        }
        EXPECT_TRUE(has_world)
            << tc.filename << " should have at least one non-empty volume";
    }
}
