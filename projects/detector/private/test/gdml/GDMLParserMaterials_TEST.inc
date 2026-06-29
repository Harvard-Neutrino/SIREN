// Included by GDMLParser_TEST.cxx. Do not compile directly.

// Shared material fragments used by multiple tests
static const std::string kHOC_Elements = R"(
    <element name="Hydrogen" Z="1"><atom value="1.008"/></element>
    <element name="Oxygen" Z="8"><atom value="15.999"/></element>
    <element name="Carbon" Z="6"><atom value="12.011"/></element>)";

static const std::string kWaterFrac =
    R"(<material name="Water"><D value="1.0" unit="g/cm3"/>
      <fraction ref="Hydrogen" n="0.112"/><fraction ref="Oxygen" n="0.888"/></material>)";

// Solids/structure for a two-volume world (inner box inside world box)
static const std::string kTwoBoxSolids =
    R"(<box name="world_box" x="5000" y="5000" z="5000" lunit="mm"/>
    <box name="inner_box" x="1000" y="1000" z="1000" lunit="mm"/>)";

static std::string TwoVolumeStructure(std::string const & inner_mat,
                                      std::string const & outer_mat) {
    return R"(<volume name="Inner"><materialref ref=")" + inner_mat +
           R"("/><solidref ref="inner_box"/></volume>
    <volume name="World"><materialref ref=")" + outer_mat +
           R"("/><solidref ref="world_box"/>
      <physvol><volumeref ref="Inner"/></physvol></volume>)";
}

// Solids/structure for a single-volume world
static const std::string kOneBoxSolid =
    R"(<box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>)";

static std::string OneVolumeStructure(std::string const & mat) {
    return R"(<volume name="World"><materialref ref=")" + mat +
           R"("/><solidref ref="world_box"/></volume>)";
}

// Helper: extract the set of nuclear Z values from material constituents
static std::set<int> NuclearZSet(MaterialModel const & materials, int mat_id) {
    auto constituents = materials.GetMaterialConstituents(mat_id);
    std::set<int> nuclear_Z;
    for(auto pt : constituents) {
        long long pdg = static_cast<long long>(pt);
        if(pdg < 0) pdg = -pdg;
        if(pdg >= 1000000000LL) {
            int Z = static_cast<int>((pdg / 10000) % 1000);
            if(Z >= 1) nuclear_Z.insert(Z);
        }
    }
    return nuclear_Z;
}


TEST(GDMLParser, CompositeOfCompositeMaterial) {
    std::string gdml = GDMLFull("",
        kHOC_Elements + "\n    " + kWaterFrac + R"(
    <material name="Mix"><D value="1.5" unit="g/cm3"/>
      <fraction ref="Water" n="0.8"/><fraction ref="Carbon" n="0.2"/></material>)",
        kTwoBoxSolids, TwoVolumeStructure("Mix", "Water"), "World");

    auto dm = LoadGDMLFromString(gdml);

    auto const & materials = dm->GetMaterials();
    ASSERT_TRUE(materials.HasMaterial("Mix"));
    int mix_id = materials.GetMaterialId("Mix");
    auto constituents = materials.GetMaterialConstituents(mix_id);
    EXPECT_GE(constituents.size(), 3u) << "Mix should resolve to >= H, C, O";

    // Check that the nuclear constituents include H, C, O after recursive resolution
    std::set<int> expected_Z = {1, 6, 8};        // H, C, O
    EXPECT_EQ(NuclearZSet(materials, mix_id), expected_Z)
        << "recursive composite resolution must yield exactly H, C, O nuclei";
}


TEST(GDMLParser, MaterialCompositionCycleThrows) {
    std::string gdml = GDMLFull("", R"(
    <material name="MatA" Z="1"><D value="1.0" unit="g/cm3"/>
      <composite n="1" ref="MatB"/></material>
    <material name="MatB" Z="1"><D value="1.0" unit="g/cm3"/>
      <composite n="1" ref="MatA"/></material>)",
        kOneBoxSolid, OneVolumeStructure("MatA"), "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error);
}


TEST(GDMLParser, CompositeMaterialMassFractions) {
    // Water = H(0.112) + O(0.888), Mix = Water(0.8) + Carbon(0.2)
    // Mix after recursive resolution:
    //   H: 0.8 * 0.112 = 0.0896
    //   O: 0.8 * 0.888 = 0.7104
    //   C: 0.2
    std::string gdml = GDMLFull("",
        kHOC_Elements + "\n    " + kWaterFrac + R"(
    <material name="Mix"><D value="1.5" unit="g/cm3"/>
      <fraction ref="Water" n="0.8"/><fraction ref="Carbon" n="0.2"/></material>)",
        kTwoBoxSolids, TwoVolumeStructure("Mix", "Water"), "World");

    auto dm = LoadGDMLFromString(gdml);

    auto const & materials = dm->GetMaterials();
    ASSERT_TRUE(materials.HasMaterial("Mix"));
    int mix_id = materials.GetMaterialId("Mix");

    // Verify nuclear constituents are present
    std::set<int> expected_Z = {1, 6, 8};
    ASSERT_EQ(NuclearZSet(materials, mix_id), expected_Z)
        << "Mix should resolve to H (Z=1), C (Z=6), O (Z=8)";

    // Verify mass fractions via GetTargetMassFraction
    using PT = siren::dataclasses::ParticleType;
    double frac_H = materials.GetTargetMassFraction(mix_id, PT::HNucleus);
    double frac_C = materials.GetTargetMassFraction(mix_id, PT::C12Nucleus);
    double frac_O = materials.GetTargetMassFraction(mix_id, PT::O16Nucleus);

    EXPECT_NEAR(frac_H, 0.0896, 1e-3)
        << "H mass fraction: 0.8 * 0.112 = 0.0896";
    EXPECT_NEAR(frac_O, 0.7104, 1e-3)
        << "O mass fraction: 0.8 * 0.888 = 0.7104";
    EXPECT_NEAR(frac_C, 0.2, 1e-3)
        << "C mass fraction: 0.2";

    // Nuclear mass fractions should sum to approximately 1.0
    double nuclear_sum = frac_H + frac_C + frac_O;
    EXPECT_NEAR(nuclear_sum, 1.0, 0.05)
        << "Nuclear mass fractions should sum close to 1.0";

    // Also verify at GDMLData level that the composition map is correct
    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.materials.count("Mix") > 0);
    auto const & mix_comp = data.materials["Mix"].composition;
    ASSERT_EQ(mix_comp.size(), 2u) << "Mix has 2 direct constituents: Water and Carbon";
    ASSERT_TRUE(mix_comp.count("Water") > 0);
    ASSERT_TRUE(mix_comp.count("Carbon") > 0);
    EXPECT_NEAR(mix_comp.at("Water"), 0.8, 1e-9);
    EXPECT_NEAR(mix_comp.at("Carbon"), 0.2, 1e-9);
}


TEST(GDMLParser, CompositeAtomCountMaterial) {
    // Water defined via <composite> atom counts: H=2, O=1
    std::string gdml = GDMLFull("", R"(
    <element name="Hydrogen" Z="1"><atom value="1.008"/></element>
    <element name="Oxygen" Z="8"><atom value="15.999"/></element>
    <material name="Water" formula="H2O"><D value="1.0" unit="g/cm3"/>
      <composite n="2" ref="Hydrogen"/><composite n="1" ref="Oxygen"/></material>)",
        kOneBoxSolid, OneVolumeStructure("Water"), "World");

    GDMLData data = ParseGDMLString(gdml);
    ASSERT_TRUE(data.materials.count("Water") > 0);
    auto const & comp = data.materials["Water"].composition;

    double total_mass = 2.0 * 1.008 + 1.0 * 15.999;
    double expected_H = 2.0 * 1.008 / total_mass;
    double expected_O = 1.0 * 15.999 / total_mass;

    ASSERT_EQ(comp.size(), 2u);
    EXPECT_NEAR(comp.at("Hydrogen"), expected_H, 1e-4)
        << "H mass fraction from atom count conversion";
    EXPECT_NEAR(comp.at("Oxygen"), expected_O, 1e-4)
        << "O mass fraction from atom count conversion";

    // Fractions must sum to 1
    double sum = 0;
    for(auto const & kv : comp) sum += kv.second;
    EXPECT_NEAR(sum, 1.0, 1e-9);
}


TEST(GDMLParser, CompositeAtomCountThroughDetectorModel) {
    std::string gdml = GDMLFull("", R"(
    <element name="Carbon" Z="6"><atom value="12.011"/></element>
    <element name="Hydrogen" Z="1"><atom value="1.008"/></element>
    <material name="Polyethylene" formula="C2H4"><D value="0.96" unit="g/cm3"/>
      <composite n="2" ref="Carbon"/><composite n="4" ref="Hydrogen"/></material>)",
        kOneBoxSolid, OneVolumeStructure("Polyethylene"), "World");

    auto dm = LoadGDMLFromString(gdml);

    auto const & materials = dm->GetMaterials();
    ASSERT_TRUE(materials.HasMaterial("Polyethylene"));
    int mat_id = materials.GetMaterialId("Polyethylene");

    using PT = siren::dataclasses::ParticleType;
    double frac_C = materials.GetTargetMassFraction(mat_id, PT::C12Nucleus);
    double frac_H = materials.GetTargetMassFraction(mat_id, PT::HNucleus);

    double total_mass = 2.0 * 12.011 + 4.0 * 1.008;
    double expected_C = 2.0 * 12.011 / total_mass;
    double expected_H = 4.0 * 1.008 / total_mass;

    EXPECT_NEAR(frac_C, expected_C, 1e-3)
        << "C mass fraction in polyethylene";
    EXPECT_NEAR(frac_H, expected_H, 1e-3)
        << "H mass fraction in polyethylene";
    EXPECT_NEAR(frac_C + frac_H, 1.0, 0.05)
        << "Nuclear mass fractions should sum close to 1.0";
}


TEST(GDMLParser, RecursiveCompositeAtomCount) {
    // Water via atom counts, then a Mix using mass fractions
    // that references the Water material.
    std::string gdml = GDMLFull("",
        kHOC_Elements + R"(
    <material name="Water" formula="H2O"><D value="1.0" unit="g/cm3"/>
      <composite n="2" ref="Hydrogen"/><composite n="1" ref="Oxygen"/></material>
    <material name="Mix"><D value="1.2" unit="g/cm3"/>
      <fraction ref="Water" n="0.7"/><fraction ref="Carbon" n="0.3"/></material>)",
        kTwoBoxSolids, TwoVolumeStructure("Mix", "Water"), "World");

    auto dm = LoadGDMLFromString(gdml);

    auto const & materials = dm->GetMaterials();
    ASSERT_TRUE(materials.HasMaterial("Mix"));
    int mix_id = materials.GetMaterialId("Mix");

    std::set<int> expected_Z = {1, 6, 8};
    EXPECT_EQ(NuclearZSet(materials, mix_id), expected_Z)
        << "Mix should resolve to H, C, O";

    double total_water = 2.0 * 1.008 + 15.999;
    double water_H = 2.0 * 1.008 / total_water;
    double water_O = 15.999 / total_water;

    using PT = siren::dataclasses::ParticleType;
    double frac_H = materials.GetTargetMassFraction(mix_id, PT::HNucleus);
    double frac_C = materials.GetTargetMassFraction(mix_id, PT::C12Nucleus);
    double frac_O = materials.GetTargetMassFraction(mix_id, PT::O16Nucleus);

    EXPECT_NEAR(frac_H, 0.7 * water_H, 1e-3)
        << "H: 0.7 * water_H_fraction";
    EXPECT_NEAR(frac_O, 0.7 * water_O, 1e-3)
        << "O: 0.7 * water_O_fraction";
    EXPECT_NEAR(frac_C, 0.3, 1e-3)
        << "C: 0.3 direct mass fraction";
}

// Geant4 exports append pointer addresses to names (e.g. "Si0x7f3a4c").
// Identical definitions with different addresses should dedup to one entry,
// and references using the pointer-suffixed names should resolve correctly.
TEST(GDMLParser, PointerSuffixStrippingAndDedup) {
    GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <element name="Hydrogen0x1111111" Z="1"><atom value="1.008"/></element>
    <element name="Hydrogen0x2222222" Z="1"><atom value="1.008"/></element>
    <element name="Oxygen0xAAAAAAA" Z="8"><atom value="15.999"/></element>
    <material name="Water0xBBBBBBB">
      <D value="1.0" unit="g/cm3"/>
      <fraction ref="Hydrogen0x2222222" n="0.112"/>
      <fraction ref="Oxygen0xAAAAAAA" n="0.888"/>
    </material>
  </materials>
  <solids>
    <box name="world0xCCCCCCC" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World0xDDDDDDD">
      <materialref ref="Water0xBBBBBBB"/>
      <solidref ref="world0xCCCCCCC"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World0xDDDDDDD"/></setup>
</gdml>)");

    // Identical Hydrogen definitions should dedup to one entry
    EXPECT_EQ(data.elements.count("Hydrogen"), 1u)
        << "Two identical Hydrogen0x... should dedup to one 'Hydrogen'";
    EXPECT_EQ(data.elements.count("Hydrogen__2"), 0u)
        << "Should not have a second instance";

    // Oxygen and Water should have clean names
    EXPECT_EQ(data.elements.count("Oxygen"), 1u);
    EXPECT_EQ(data.materials.count("Water"), 1u);

    // Solid should have clean name
    EXPECT_EQ(data.solids.count("world"), 1u);

    // Volume should have clean name
    EXPECT_EQ(data.volumes.count("World"), 1u);

    // Composition references should resolve through stripped names
    auto const & water = data.materials.at("Water");
    EXPECT_TRUE(water.composition.count("Hydrogen") > 0)
        << "Water composition should reference 'Hydrogen', not 'Hydrogen__1' or 'Hydrogen0x...'";
    EXPECT_TRUE(water.composition.count("Oxygen") > 0)
        << "Water composition should reference 'Oxygen', not 'Oxygen__1' or 'Oxygen0x...'";
    EXPECT_EQ(water.composition.count("Hydrogen0x1111111"), 0u);
    EXPECT_EQ(water.composition.count("Hydrogen0x2222222"), 0u);
    EXPECT_EQ(water.composition.count("Hydrogen__1"), 0u);
}

// Identical materials with different pointer suffixes should dedup,
// and the composition resolution through DetectorModel should work.
TEST(GDMLParser, PointerSuffixMaterialDedupThroughModel) {
    std::string gdml = GDMLFull("",
        R"(
        <element name="Hydrogen0x1111111" Z="1"><atom value="1.008"/></element>
        <element name="Oxygen0x2222222" Z="8"><atom value="15.999"/></element>
        <material name="Water0xAAAAAAA">
          <D value="1.0" unit="g/cm3"/>
          <fraction ref="Hydrogen0x1111111" n="0.112"/>
          <fraction ref="Oxygen0x2222222" n="0.888"/>
        </material>
        <material name="Air0xBBBBBBB" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
        )",
        R"(
        <box name="world0xCCCCCCC" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="inner0xDDDDDDD" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="Inner0xEEEEEEE">
          <materialref ref="Water0xAAAAAAA"/>
          <solidref ref="inner0xDDDDDDD"/>
        </volume>
        <volume name="World0xFFFFFFF">
          <materialref ref="Air0xBBBBBBB"/>
          <solidref ref="world0xCCCCCCC"/>
          <physvol><volumeref ref="Inner0xEEEEEEE"/></physvol>
        </volume>
        )", "World0xFFFFFFF");

    auto dm = LoadGDMLFromString(gdml);
    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "Water")
        << "Volume with pointer-suffixed materialref should resolve to stripped name";
}

// Non-identical solids with the same base name but different pointer
// suffixes should keep both instances and resolve correctly through aliases.
TEST(GDMLParser, NonIdenticalSolidsWithPointerSuffixes) {
    GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox0x1111111" x="100" y="100" z="100" lunit="mm"/>
    <box name="mybox0x2222222" x="200" y="200" z="200" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="mybox0x1111111"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)");

    // Two non-identical boxes: both kept as separate instances
    EXPECT_EQ(data.solid_instance_counts["mybox"], 2);
    EXPECT_EQ(data.solids.count("mybox__1"), 1u);
    EXPECT_EQ(data.solids.count("mybox__2"), 1u);

    // solidref="mybox0x1111111" should resolve to the first (smaller) instance
    auto const & vol = data.volumes.at("World");
    EXPECT_EQ(vol.solid_ref, "mybox__1")
        << "solidref with pointer suffix should resolve to specific instance";
    double tol = 1e-9;
    EXPECT_NEAR(data.solids["mybox__1"]->GetBoundingBox().max_corner.GetX(), 0.05, tol)
        << "First instance should be x=100mm full-width -> half-extent 0.05m";
    EXPECT_NEAR(data.solids["mybox__2"]->GetBoundingBox().max_corner.GetX(), 0.1, tol)
        << "Second instance should be x=200mm full-width -> half-extent 0.1m";
}

// Identical materials with different pointer suffixes should collapse
// to a single entry, not create __1 and __2.
TEST(GDMLParser, IdenticalMaterialsCollapse) {
    GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <material name="Iron0x1111111" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
    <material name="Iron0x2222222" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
    <material name="Iron0x3333333" Z="26"><D value="7.874" unit="g/cm3"/><atom value="55.845"/></material>
  </materials>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref=""/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)");

    // All three are identical, should collapse to one "Iron"
    EXPECT_EQ(data.materials.count("Iron"), 1u)
        << "Three identical Iron0x... should collapse to one 'Iron'";
    EXPECT_EQ(data.materials.count("Iron__2"), 0u);
    EXPECT_EQ(data.materials.count("Iron__3"), 0u);
    EXPECT_EQ(data.material_instance_counts["Iron"], 1);
}


// Mixed bare and pointer-suffixed names that overlap after stripping.
// A hand-written <box name="mybox"/> followed by a Geant4-exported
// <box name="mybox0x1234ABC"/> should coexist if they differ, or dedup
// if identical.
TEST(GDMLParser, MixedBareAndPointerSuffixedSolids) {
    // Case 1: bare and suffixed are DIFFERENT sizes -> two instances
    {
        GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox" x="100" y="100" z="100" lunit="mm"/>
    <box name="mybox0xABCDABC" x="200" y="200" z="200" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="mybox0xABCDABC"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)");

        EXPECT_EQ(data.solid_instance_counts["mybox"], 2)
            << "Bare 'mybox' and stripped 'mybox0xABCDABC' are different, so 2 instances";
        EXPECT_EQ(data.solids.count("mybox__1"), 1u);
        EXPECT_EQ(data.solids.count("mybox__2"), 1u);

        // The pointer-suffixed ref should resolve to the second instance
        auto const & vol = data.volumes.at("World");
        EXPECT_EQ(vol.solid_ref, "mybox__2")
            << "solidref='mybox0xABCDABC' should alias to the second instance";
    }

    // Case 2: bare and suffixed are IDENTICAL -> dedup to one
    {
        GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox" x="100" y="100" z="100" lunit="mm"/>
    <box name="mybox0xABCDABC" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="mybox0xABCDABC"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)");

        EXPECT_EQ(data.solid_instance_counts["mybox"], 1)
            << "Bare and stripped-identical should dedup to 1 instance";
        EXPECT_EQ(data.solids.count("mybox"), 1u)
            << "Single instance should have bare name (no __1 suffix)";

        // Both the bare name and the pointer-suffixed name should work
        auto const & vol = data.volumes.at("World");
        EXPECT_EQ(vol.solid_ref, "mybox")
            << "solidref='mybox0xABCDABC' should alias to 'mybox' after dedup";
    }
}

TEST(GDMLParser, MixedBareAndPointerSuffixedMaterials) {
    // Bare element "Hydrogen" and pointer-suffixed identical "Hydrogen0x1234ABC"
    // should dedup. A material referencing the suffixed name should resolve.
    GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <element name="Hydrogen" Z="1"><atom value="1.008"/></element>
    <element name="Hydrogen0xABCDABC" Z="1"><atom value="1.008"/></element>
    <element name="Oxygen" Z="8"><atom value="15.999"/></element>
    <material name="Water">
      <D value="1.0" unit="g/cm3"/>
      <fraction ref="Hydrogen0xABCDABC" n="0.112"/>
      <fraction ref="Oxygen" n="0.888"/>
    </material>
  </materials>
  <solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
  <structure><volume name="W"><materialref ref="Water"/><solidref ref="b"/></volume></structure>
  <setup name="Default" version="1.0"><world ref="W"/></setup>
</gdml>)");

    EXPECT_EQ(data.elements.count("Hydrogen"), 1u)
        << "Bare + identical pointer-suffixed should dedup to one 'Hydrogen'";
    EXPECT_EQ(data.elements.count("Hydrogen__2"), 0u);

    // Water's composition should reference "Hydrogen" (the stripped name)
    auto const & water = data.materials.at("Water");
    EXPECT_TRUE(water.composition.count("Hydrogen") > 0)
        << "Water should reference 'Hydrogen', not 'Hydrogen0xABCDABC' or 'Hydrogen__1'";
}


// Same as MixedBareAndPointerSuffixedSolids but with the pointer-suffixed
// definition FIRST and the bare one second. Order should not matter.
TEST(GDMLParser, MixedPointerThenBareSolids) {
    // Different sizes: pointer-suffixed first, bare second -> two instances
    {
        GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox0xABCDABC" x="200" y="200" z="200" lunit="mm"/>
    <box name="mybox" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="mybox0xABCDABC"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)");

        EXPECT_EQ(data.solid_instance_counts["mybox"], 2);
        // Pointer-suffixed ref should resolve to the first (larger) instance
        auto const & vol = data.volumes.at("World");
        EXPECT_EQ(vol.solid_ref, "mybox__1");
        double tol = 1e-9;
        EXPECT_NEAR(data.solids["mybox__1"]->GetBoundingBox().max_corner.GetX(), 0.1, tol)
            << "First instance (from pointer-suffixed) should be the larger box (200mm->0.1m)";
        EXPECT_NEAR(data.solids["mybox__2"]->GetBoundingBox().max_corner.GetX(), 0.05, tol)
            << "Second instance (from bare) should be the smaller box (100mm->0.05m)";
    }

    // Identical: pointer-suffixed first, bare second -> dedup to one
    {
        GDMLData data = ParseGDMLString(R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <box name="mybox0xABCDABC" x="100" y="100" z="100" lunit="mm"/>
    <box name="mybox" x="100" y="100" z="100" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="mybox0xABCDABC"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0"><world ref="World"/></setup>
</gdml>)");

        EXPECT_EQ(data.solid_instance_counts["mybox"], 1)
            << "Pointer-first then identical bare should dedup to 1";
        EXPECT_EQ(data.solids.count("mybox"), 1u);
        auto const & vol = data.volumes.at("World");
        EXPECT_EQ(vol.solid_ref, "mybox")
            << "Alias should resolve to bare name after dedup cleanup";
    }
}


// Zero-density materials (common for GDML vacuum world volumes) should
// parse without error, converting to a tiny vacuum density with a warning.
TEST(GDMLParser, ZeroDensityVacuumMaterial) {
    std::string gdml = GDMLFull("",
        R"(
        <material name="Vacuum" Z="1">
          <D value="0" unit="g/cm3"/>
          <atom value="1.008"/>
        </material>
        <material name="Iron" Z="26">
          <D value="7.874" unit="g/cm3"/>
          <atom value="55.845"/>
        </material>
        )",
        kTwoBoxSolids,
        TwoVolumeStructure("Iron", "Vacuum"),
        "World");

    EXPECT_NO_THROW(LoadGDMLFromString(gdml))
        << "Zero-density vacuum material should parse without error";

    auto dm = LoadGDMLFromString(gdml);
    // The vacuum world volume should still be present
    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "Iron");
}

TEST(GDMLParser, NegativeDensityThrows) {
    std::string gdml = GDMLFull("",
        R"(
        <material name="BadMat" Z="1">
          <D value="-1.0" unit="g/cm3"/>
          <atom value="1.008"/>
        </material>
        )",
        R"(<box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>)",
        R"(
        <volume name="World">
          <materialref ref="BadMat"/>
          <solidref ref="world_box"/>
        </volume>
        )", "World");

    EXPECT_THROW(LoadGDMLFromString(gdml), std::runtime_error)
        << "Negative density should still throw";
}

// When an element and material share the same name (common in Geant4
// exports), composition references should resolve to the element, not
// the material. This prevents circular resolution.
TEST(GDMLParser, ElementMaterialSameNameComposition) {
    std::string gdml = GDMLFull("",
        R"(
        <element name="C" Z="6"><atom value="12.011"/></element>
        <material name="C" state="solid">
          <D value="2.0" unit="g/cm3"/>
          <fraction ref="C" n="1"/>
        </material>
        <material name="Air" Z="7"><D value="0.00129" unit="g/cm3"/><atom value="14.007"/></material>
        )",
        R"(
        <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
        <box name="c_box" x="100" y="100" z="100" lunit="mm"/>
        )",
        R"(
        <volume name="CarbonVol">
          <materialref ref="C"/>
          <solidref ref="c_box"/>
        </volume>
        <volume name="World">
          <materialref ref="Air"/>
          <solidref ref="world_box"/>
          <physvol><volumeref ref="CarbonVol"/></physvol>
        </volume>
        )", "World");

    // Should not throw with circular composition error
    EXPECT_NO_THROW(LoadGDMLFromString(gdml))
        << "Element/material 'C' name collision should not cause circular resolution";
    auto dm = LoadGDMLFromString(gdml);
    EXPECT_EQ(MaterialAt(*dm, 0, 0, 0), "C");
}

// ---------------------------------------------------------------------------
// <D> value-expression guard
// ---------------------------------------------------------------------------
//
// <D value="<expr>"/> with no unit= attribute is evaluated against the
// constants table, which is seeded with CLHEP unit identifiers (gram,
// cm, eV, ...) at their CLHEP-internal scale. A value like "1*gram"
// would silently store density = 6.24e21 because the unit-attribute
// conversion path below the SafeParseDouble call has no way to know
// the value already absorbed a unit scale. The materials parser
// rejects any value-expression that references a dimensional built-in
// (math-only built-ins like pi/e/gamma remain allowed).

namespace {
std::string DensityGDML(std::string const & d_attr,
                        std::string const & defines = "") {
    return R"(<?xml version="1.0"?><gdml>
<define>)" + defines + R"(</define>
<materials>
  <material name="TestMat" Z="1">
    <D )" + d_attr + R"(/>
    <atom value="1.0"/>
  </material>
</materials>
<solids><box name="b" x="1" y="1" z="1" lunit="mm"/></solids>
<structure>
  <volume name="World"><materialref ref="TestMat"/><solidref ref="b"/></volume>
</structure>
<setup name="default" version="1.0"><world ref="World"/></setup>
</gdml>)";
}
}

TEST(GDMLParser, DensityStandardFormWorks) {
    // Sanity: <D value="N" unit="g/cm3"/> is the canonical form and
    // must continue to work after the inline-unit guard lands.
    auto data = ParseGDMLString(DensityGDML(R"(value="2.7" unit="g/cm3")"));
    ASSERT_EQ(data.materials.size(), 1u);
    EXPECT_NEAR(data.materials.begin()->second.density, 2.7, 1e-9);
}

TEST(GDMLParser, DensityBareNumberWorks) {
    auto data = ParseGDMLString(DensityGDML(R"(value="2.7")"));
    ASSERT_EQ(data.materials.size(), 1u);
    EXPECT_NEAR(data.materials.begin()->second.density, 2.7, 1e-9);
}

TEST(GDMLParser, DensityPureMathExpressionWorks) {
    // pi/e/gamma are dimensionless math constants and must NOT trip
    // the guard. Use pi/3.14 just to confirm a meaningful expression
    // round-trips.
    auto data = ParseGDMLString(DensityGDML(R"(value="pi*0.5")"));
    ASSERT_EQ(data.materials.size(), 1u);
    EXPECT_NEAR(data.materials.begin()->second.density, M_PI * 0.5, 1e-9);
}

TEST(GDMLParser, DensityUserConstantWorks) {
    auto data = ParseGDMLString(DensityGDML(
        R"(value="my_density")",
        R"(<constant name="my_density" value="2.7"/>)"));
    ASSERT_EQ(data.materials.size(), 1u);
    EXPECT_NEAR(data.materials.begin()->second.density, 2.7, 1e-9);
}

TEST(GDMLParser, DensityTypedQuantityReferenceWorks) {
    // <quantity type="density"> routes through QuantityUnitScale, so
    // referencing one by name from <D> must remain allowed even though
    // the quantity name lives in the same constants table.
    auto data = ParseGDMLString(DensityGDML(
        R"(value="my_qty")",
        R"(<quantity name="my_qty" type="density" value="1000" unit="kg/m3"/>)"));
    ASSERT_EQ(data.materials.size(), 1u);
    EXPECT_NEAR(data.materials.begin()->second.density, 1.0, 1e-9)
        << "1000 kg/m3 = 1 g/cm3 via QuantityUnitScale";
}

TEST(GDMLParser, DensityInlineUnitGramThrows) {
    // The bug class: '1*gram' would silently store ~6.24e21 g/cm3.
    EXPECT_THROW(ParseGDMLString(DensityGDML(R"(value="1*gram")")),
                 std::runtime_error);
}

TEST(GDMLParser, DensityInlineUnitCmThrows) {
    // Length identifier inside a density expression - also nonsense.
    EXPECT_THROW(ParseGDMLString(DensityGDML(R"(value="3*cm")")),
                 std::runtime_error);
}

TEST(GDMLParser, DensityInlineUnitEnergyThrows) {
    // Energy identifier - eV evaluates to 1e-6 in CLHEP-internal
    // units, silently storing 1e-6 g/cm3 without the guard.
    EXPECT_THROW(ParseGDMLString(DensityGDML(R"(value="1*eV")")),
                 std::runtime_error);
}

TEST(GDMLParser, DensityInlineUnitErrorMessageHelpful) {
    // Diagnostic must identify the offending unit and point the user
    // at the standard form (unit attribute or typed quantity).
    try {
        ParseGDMLString(DensityGDML(R"(value="1*gram")"));
        FAIL() << "Expected runtime_error";
    } catch(std::runtime_error const & e) {
        std::string msg = e.what();
        EXPECT_NE(msg.find("TestMat"), std::string::npos)
            << "Error should name the material: " << msg;
        EXPECT_NE(msg.find("gram"), std::string::npos)
            << "Error should name the offending unit: " << msg;
        EXPECT_NE(msg.find("unit=\"g/cm3\""), std::string::npos)
            << "Error should suggest the standard form: " << msg;
    }
}
