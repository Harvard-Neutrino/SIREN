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
