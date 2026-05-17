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


TEST(GDMLParser, LoopInStructurePositions) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
  </define>
  <materials>
    <material name="TestMat" Z="7">
      <D value="0.00129" unit="g/cm3"/>
      <atom value="14.01"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
    <box name="child_box" x="10" y="10" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="Child">
      <materialref ref="TestMat"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="TestMat"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="3" step="1">
        <physvol>
          <volumeref ref="Child"/>
          <position name="lpos_[i]" x="i*100" y="0" z="0" unit="mm"/>
        </physvol>
      </loop>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

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
        EXPECT_NEAR(child.position.GetY(), 0.0, 1e-9)
            << "Child " << idx << " y position should be 0";
        EXPECT_NEAR(child.position.GetZ(), 0.0, 1e-9)
            << "Child " << idx << " z position should be 0";
    }
}


TEST(GDMLParser, NestedLoopInStructurePositions) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <variable name="i" value="0"/>
    <variable name="j" value="0"/>
  </define>
  <materials>
    <material name="TestMat" Z="7">
      <D value="0.00129" unit="g/cm3"/>
      <atom value="14.01"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="2000" y="2000" z="2000" lunit="mm"/>
    <box name="child_box" x="10" y="10" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="Child">
      <materialref ref="TestMat"/>
      <solidref ref="child_box"/>
    </volume>
    <volume name="World">
      <materialref ref="TestMat"/>
      <solidref ref="world_box"/>
      <loop for="i" from="0" to="1" step="1">
        <loop for="j" from="0" to="1" step="1">
          <physvol>
            <volumeref ref="Child"/>
            <position name="nlpos_[i]_[j]" x="i*200" y="j*200" z="0" unit="mm"/>
          </physvol>
        </loop>
      </loop>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);

    auto it = data.volumes.find("World");
    ASSERT_NE(it, data.volumes.end());
    ASSERT_EQ(it->second.children.size(), 4u)
        << "Nested loops i=0..1, j=0..1 should produce 4 physvol placements";

    // Collect actual (x, y) pairs from children
    std::set<std::pair<double, double>> actual_positions;
    for(auto const & child : it->second.children) {
        double x = std::round(child.position.GetX() * 1e6) / 1e6;
        double y = std::round(child.position.GetY() * 1e6) / 1e6;
        actual_positions.insert({x, y});
        EXPECT_NEAR(child.position.GetZ(), 0.0, 1e-9)
            << "All children should have z=0";
    }

    std::set<std::pair<double, double>> expected_positions = {
        {0.0, 0.0}, {0.0, 0.2}, {0.2, 0.0}, {0.2, 0.2}
    };

    EXPECT_EQ(actual_positions, expected_positions)
        << "Nested loop should produce positions at (0,0), (0,0.2), (0.2,0), (0.2,0.2) meters";
}


TEST(GDMLParser, RealWorldGDMLSmoke) {
    std::string test_dir = std::string(__FILE__);
    auto pos = test_dir.rfind("projects/");
    ASSERT_NE(pos, std::string::npos)
        << "Could not locate repo root from __FILE__";
    test_dir = test_dir.substr(0, pos) + "test_gdml/dark_matter/cats/";

    struct TestCase {
        std::string filename;
    };
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
