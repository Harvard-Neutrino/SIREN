// Included by GDMLParser_TEST.cxx. Do not compile directly.

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
    // The original test only counted constituents (>=3). That passes even
    // if recursive resolution silently drops Carbon or fails to expand
    // Water into H+O. SIREN decomposes a material into scattering targets:
    // the element nuclei plus free p/n/e. Value-check that the resolved
    // NUCLEAR set is exactly {H (Z=1), C (Z=6), O (Z=8)} -- nuclear PDG
    // codes have the form 10LZZZAAAI (>= 1000000000).
    EXPECT_GE(constituents.size(), 3u) << "Mix should resolve to >= H, C, O";
    std::set<int> nuclear_Z;
    for(auto pt : constituents) {
        long long pdg = static_cast<long long>(pt);
        if(pdg < 0) pdg = -pdg;
        if(pdg >= 1000000000LL) {                 // nuclear PDG 10LZZZAAAI
            int Z = static_cast<int>((pdg / 10000) % 1000);
            if(Z >= 1) nuclear_Z.insert(Z);       // skip SIREN pseudo-targets
        }
    }
    std::set<int> expected_Z = {1, 6, 8};        // H, C, O
    EXPECT_EQ(nuclear_Z, expected_Z)
        << "recursive composite resolution must yield exactly H, C, O nuclei";
}


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


TEST(GDMLParser, CompositeMaterialMassFractions) {
    // Water = H(0.112) + O(0.888), Mix = Water(0.8) + Carbon(0.2)
    // After recursive resolution:
    //   H: 0.8 * 0.112 = 0.0896
    //   O: 0.8 * 0.888 = 0.7104
    //   C: 0.2
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

    std::string tmpfile = "/tmp/siren_gdml_massfrac_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }

    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    auto const & materials = dm.GetMaterials();
    ASSERT_TRUE(materials.HasMaterial("Mix"));
    int mix_id = materials.GetMaterialId("Mix");

    // Verify nuclear constituents are present
    auto constituents = materials.GetMaterialConstituents(mix_id);
    std::set<int> nuclear_Z;
    for(auto pt : constituents) {
        long long pdg = static_cast<long long>(pt);
        if(pdg < 0) pdg = -pdg;
        if(pdg >= 1000000000LL) {
            int Z = static_cast<int>((pdg / 10000) % 1000);
            if(Z >= 1) nuclear_Z.insert(Z);
        }
    }
    std::set<int> expected_Z = {1, 6, 8};
    ASSERT_EQ(nuclear_Z, expected_Z)
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
    // Expected mass fractions: H = 2*1.008 / (2*1.008 + 15.999) = 0.11191
    //                          O = 15.999 / (2*1.008 + 15.999) = 0.88809
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
    <material name="Water" formula="H2O">
      <D value="1.0" unit="g/cm3"/>
      <composite n="2" ref="Hydrogen"/>
      <composite n="1" ref="Oxygen"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Water"/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

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
    // Polyethylene: C2H4
    //   C: 2 * 12.011 = 24.022
    //   H: 4 * 1.008  =  4.032
    //   total = 28.054
    //   C frac = 24.022/28.054 = 0.8563
    //   H frac =  4.032/28.054 = 0.1437
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <element name="Carbon" Z="6">
      <atom value="12.011"/>
    </element>
    <element name="Hydrogen" Z="1">
      <atom value="1.008"/>
    </element>
    <material name="Polyethylene" formula="C2H4">
      <D value="0.96" unit="g/cm3"/>
      <composite n="2" ref="Carbon"/>
      <composite n="4" ref="Hydrogen"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref="Polyethylene"/>
      <solidref ref="world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    std::string tmpfile = "/tmp/siren_gdml_composite_atom_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    auto const & materials = dm.GetMaterials();
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
    // Water: H2O -> H frac = 2*1.008/(2*1.008+15.999), O frac = 15.999/(...)
    // Mix: Water(0.7) + Carbon(0.3)
    // After resolution: H = 0.7 * H_frac, O = 0.7 * O_frac, C = 0.3
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
    <material name="Water" formula="H2O">
      <D value="1.0" unit="g/cm3"/>
      <composite n="2" ref="Hydrogen"/>
      <composite n="1" ref="Oxygen"/>
    </material>
    <material name="Mix">
      <D value="1.2" unit="g/cm3"/>
      <fraction ref="Water" n="0.7"/>
      <fraction ref="Carbon" n="0.3"/>
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

    std::string tmpfile = "/tmp/siren_gdml_recursive_composite_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml;
    }
    DetectorModel dm;
    dm.LoadGDML(tmpfile);
    std::remove(tmpfile.c_str());

    auto const & materials = dm.GetMaterials();
    ASSERT_TRUE(materials.HasMaterial("Mix"));
    int mix_id = materials.GetMaterialId("Mix");

    auto constituents = materials.GetMaterialConstituents(mix_id);
    std::set<int> nuclear_Z;
    for(auto pt : constituents) {
        long long pdg = static_cast<long long>(pt);
        if(pdg < 0) pdg = -pdg;
        if(pdg >= 1000000000LL) {
            int Z = static_cast<int>((pdg / 10000) % 1000);
            if(Z >= 1) nuclear_Z.insert(Z);
        }
    }
    std::set<int> expected_Z = {1, 6, 8};
    EXPECT_EQ(nuclear_Z, expected_Z)
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
