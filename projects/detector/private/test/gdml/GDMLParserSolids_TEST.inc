// Included by GDMLParser_TEST.cxx. Do not compile directly.

TEST(GDMLParser, AllSolidTypes) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <box name="b1" x="100" y="100" z="100" lunit="mm"/>
    <sphere name="s1" rmax="50" lunit="mm"/>
    <tube name="t1" rmin="0" rmax="50" z="100" deltaphi="360" aunit="deg" lunit="mm"/>
    <cone name="c1" rmin1="0" rmax1="50" rmin2="0" rmax2="30" z="100" deltaphi="360" aunit="deg" lunit="mm"/>
    <polycone name="pc1" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-50"/>
      <zplane rmin="0" rmax="30" z="50"/>
    </polycone>
    <polyhedra name="ph1" numsides="6" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-40"/>
      <zplane rmin="0" rmax="50" z="40"/>
    </polyhedra>
    <trd name="trd1" x1="60" x2="40" y1="80" y2="50" z="100" lunit="mm"/>
    <xtru name="xt1" lunit="mm">
      <twoDimVertex x="-50" y="-50"/>
      <twoDimVertex x="50" y="-50"/>
      <twoDimVertex x="50" y="50"/>
      <twoDimVertex x="-50" y="50"/>
      <section zOrder="0" zPosition="-40" xOffset="0" yOffset="0" scalingFactor="1"/>
      <section zOrder="1" zPosition="40" xOffset="0" yOffset="0" scalingFactor="1"/>
    </xtru>
    <eltube name="et1" dx="30" dy="20" dz="50" lunit="mm"/>
    <ellipsoid name="ell1" ax="50" by="40" cz="30" lunit="mm"/>
    <para name="para1" x="40" y="30" z="50" alpha="0.3" theta="0.2" phi="0.5" aunit="rad" lunit="mm"/>
    <trap name="trap1" z="100" theta="0" phi="0" y1="50" x1="40" x2="40" alpha1="0" y2="30" x3="30" x4="30" alpha2="0" aunit="rad" lunit="mm"/>
    <cutTube name="ct1" rmin="0" rmax="50" z="100" startphi="0" deltaphi="360" lowX="0" lowY="0" lowZ="-1" highX="0" highY="0" highZ="1" aunit="deg" lunit="mm"/>
)", "b1");

    std::vector<std::string> names = {"b1", "s1", "t1", "c1", "pc1", "ph1", "trd1", "xt1", "et1", "ell1", "para1", "trap1", "ct1"};
    for(auto const & name : names) {
        ASSERT_TRUE(data.solids.count(name) > 0)
            << "Solid '" << name << "' should be parsed";
        auto bb = data.solids[name]->GetBoundingBox();
        EXPECT_TRUE(bb.IsValid())
            << "Solid '" << name << "' should have a valid bounding box";
        EXPECT_TRUE(std::isfinite(bb.max_corner.GetX()))
            << "Solid '" << name << "' bounding box should have finite dimensions";
        EXPECT_TRUE(std::isfinite(bb.max_corner.GetY()))
            << "Solid '" << name << "' bounding box should have finite dimensions";
        EXPECT_TRUE(std::isfinite(bb.max_corner.GetZ()))
            << "Solid '" << name << "' bounding box should have finite dimensions";
    }

    double tol = 1e-6;
    // GDML box x/y/z are full-widths: 100mm full-width -> 0.05m half-extent
    auto bb_b1 = data.solids["b1"]->GetBoundingBox();
    EXPECT_NEAR(bb_b1.max_corner.GetX(), 0.05, tol);
    EXPECT_NEAR(bb_b1.max_corner.GetZ(), 0.05, tol);

    // sphere: rmax=50mm = 0.05m
    auto bb_s1 = data.solids["s1"]->GetBoundingBox();
    EXPECT_NEAR(bb_s1.max_corner.GetX(), 0.05, tol);

    // tube: rmax=50mm=0.05m, GDML z=100mm full-height -> 0.05m half-extent
    auto bb_t1 = data.solids["t1"]->GetBoundingBox();
    EXPECT_NEAR(bb_t1.max_corner.GetX(), 0.05, tol);
    EXPECT_NEAR(bb_t1.max_corner.GetZ(), 0.05, tol);

    // trd: GDML full-widths halved -> max(x1=60,x2=40)/2=30mm=0.03m, z=100mm/2=0.05m
    auto bb_trd = data.solids["trd1"]->GetBoundingBox();
    EXPECT_NEAR(bb_trd.max_corner.GetX(), 0.03, tol);
    EXPECT_NEAR(bb_trd.max_corner.GetZ(), 0.05, tol);

    // elliptical tube: dx=30mm=0.03m, dy=20mm=0.02m, dz=50mm half-z=0.05m
    auto bb_et = data.solids["et1"]->GetBoundingBox();
    EXPECT_NEAR(bb_et.max_corner.GetX(), 0.03, tol);
    EXPECT_NEAR(bb_et.max_corner.GetY(), 0.02, tol);
    EXPECT_NEAR(bb_et.max_corner.GetZ(), 0.05, tol);

    // ellipsoid: ax=50mm=0.05m, by=40mm=0.04m, cz=30mm=0.03m
    auto bb_ell = data.solids["ell1"]->GetBoundingBox();
    EXPECT_NEAR(bb_ell.max_corner.GetX(), 0.05, tol);
    EXPECT_NEAR(bb_ell.max_corner.GetY(), 0.04, tol);
    EXPECT_NEAR(bb_ell.max_corner.GetZ(), 0.03, tol);
}


TEST(GDMLParser, UnsupportedSolidStrictThrows) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define/>
  <materials/>
  <solids>
    <hype name="h" rmin="1" rmax="2" inst="0.1" outst="0.2" z="10" lunit="mm"/>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="h"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    // Strict mode: throws on unsupported solid
    std::string tmpfile = TempPath("siren_gdml_test_strict.gdml");
    { std::ofstream f(tmpfile); f << gdml; }
    GDMLParseOptions strict_opts;
    strict_opts.strict = true;
    EXPECT_THROW(ParseGDML(tmpfile, strict_opts), std::runtime_error);
    std::remove(tmpfile.c_str());
}


TEST(GDMLParser, UnsupportedSolidSkipsWithWarning) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <box name="good_box" x="100" y="100" z="100" lunit="mm"/>
    <hype name="h" rmin="1" rmax="2" inst="0.1" outst="0.2" z="10" lunit="mm"/>
)", "good_box");
    // The box should have parsed, the hype should have been skipped
    EXPECT_EQ(data.solids.count("good_box"), 1u);
    EXPECT_EQ(data.solids.count("h"), 0u);
    // Warning should be recorded
    ASSERT_GE(data.warnings.size(), 1u);
    bool found_hype_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("hype") != std::string::npos) found_hype_warning = true;
    }
    EXPECT_TRUE(found_hype_warning) << "Expected warning about unsupported hype solid";
}


TEST(GDMLParser, TorusIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <torus name="my_torus" rmin="10" rmax="20" rtor="100" startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>
)", "my_torus");
    ASSERT_TRUE(data.solids.count("my_torus") > 0);
    // Verify bounding box: rtor=100mm=0.1m, rmax=20mm=0.02m
    // xy extent = rtor + rmax = 0.12m, z extent = rmax = 0.02m
    double tol = 1e-6;
    auto bb = data.solids["my_torus"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 0.12, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.02, tol);
    // No warnings for a full-rotation torus
    ExpectNoWarningContaining(data, "torus");
}


TEST(GDMLParser, PartialSphereParses) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <sphere name="hemi" rmin="0" rmax="100" starttheta="0" deltatheta="90" startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>
)", "hemi");
    ASSERT_TRUE(data.solids.count("hemi") > 0);
    // No warnings about partial sphere
    ExpectNoWarningContaining(data, "sphere");
    ExpectNoWarningContaining(data, "angular");
}


TEST(GDMLParser, PartialTorusParses) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <torus name="elbow" rmin="10" rmax="20" rtor="100" startphi="0" deltaphi="90" lunit="mm" aunit="deg"/>
)", "elbow");
    ASSERT_TRUE(data.solids.count("elbow") > 0);
    // No warnings about partial torus
    ExpectNoWarningContaining(data, "torus");
    ExpectNoWarningContaining(data, "angular");
}


TEST(GDMLParser, PartialCylinderSupported) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <tube name="wedge" rmin="0" rmax="100" z="200" deltaphi="180" lunit="mm" aunit="deg"/>
)", "wedge");
    ASSERT_TRUE(data.solids.count("wedge") > 0);
    ASSERT_TRUE(data.solids["wedge"] != nullptr);
    // Partial phi should produce no warnings
    ExpectNoWarningContaining(data, "partial angular extent");
}


TEST(GDMLParser, GenericPolycone) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <genericPolycone name="gpc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <rzpoint r="50" z="-100"/>
      <rzpoint r="80" z="0"/>
      <rzpoint r="30" z="100"/>
    </genericPolycone>
)", "gpc");
    auto geo = RequireSolid(data, "gpc");
    auto bb = geo->GetBoundingBox();
    double tol = 1e-6;
    // rmax=80mm=0.08m at z=0
    EXPECT_NEAR(bb.max_corner.GetX(), 0.08, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 0.08, tol);
    // z ranges from -100mm to 100mm = -0.1m to 0.1m
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.1, tol);
    EXPECT_NEAR(bb.min_corner.GetZ(), -0.1, tol);
}


TEST(GDMLParser, CutTubeIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <cutTube name="ct" rmin="10" rmax="50" z="100" startphi="0" deltaphi="360"
             lowX="0" lowY="0" lowZ="-1" highX="0" highY="0" highZ="1"
             aunit="deg" lunit="mm"/>
)", "ct");
    auto geo = RequireSolid(data, "ct");
    auto bb = geo->GetBoundingBox();
    double tol = 1e-6;
    // rmax=50mm=0.05m; GDML z=100mm is full-height -> dz=0.05m half (flat cuts)
    EXPECT_NEAR(bb.max_corner.GetX(), 0.05, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.05, tol);
}


TEST(GDMLParser, TrapIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <trap name="t" z="100" theta="0" phi="0"
          y1="60" x1="40" x2="40" alpha1="0"
          y2="40" x3="30" x4="30" alpha2="0"
          aunit="rad" lunit="mm"/>
)", "t");
    auto geo = RequireSolid(data, "t");
    auto bb = geo->GetBoundingBox();
    double tol = 1e-6;
    // GDML trap z=100mm is a full-length -> dz=0.05m half-height
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.05, tol);
    // max x = max(dx1..dx4) = 40mm full / 2 = 0.02m
    EXPECT_NEAR(bb.max_corner.GetX(), 0.02, tol);
}


TEST(GDMLParser, EllipsoidIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <ellipsoid name="ell" ax="50" by="30" cz="40" zcut2="20" lunit="mm"/>
)", "ell");
    auto geo = RequireSolid(data, "ell");
    auto bb = geo->GetBoundingBox();
    double tol = 1e-6;
    EXPECT_NEAR(bb.max_corner.GetX(), 0.05, tol);  // ax=50mm=0.05m
    EXPECT_NEAR(bb.max_corner.GetY(), 0.03, tol);  // by=30mm=0.03m
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.02, tol);  // zcut2=20mm=0.02m
    EXPECT_NEAR(bb.min_corner.GetZ(), -0.04, tol); // -cz=-40mm=-0.04m (no zcut1)
}


TEST(GDMLParser, ParaIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <para name="p" x="40" y="30" z="50" alpha="0" theta="0" phi="0" aunit="rad" lunit="mm"/>
)", "p");
    auto geo = RequireSolid(data, "p");
    auto bb = geo->GetBoundingBox();
    double tol = 1e-6;
    // With alpha=theta=phi=0, para = box. GDML x,y,z are full-lengths:
    // x=40mm->0.02m, y=30mm->0.015m, z=50mm->0.025m half-extents.
    EXPECT_NEAR(bb.max_corner.GetX(), 0.02, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 0.015, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.025, tol);
}


TEST(GDMLParser, TessellatedTriangular) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="v0" x="0" y="0" z="0" unit="m"/>
    <position name="v1" x="1" y="0" z="0" unit="m"/>
    <position name="v2" x="0" y="1" z="0" unit="m"/>
    <position name="v3" x="0" y="0" z="1" unit="m"/>
  </define>
  <materials/>
  <solids>
    <tessellated name="tetra">
      <triangular vertex1="v0" vertex2="v2" vertex3="v1"/>
      <triangular vertex1="v0" vertex2="v1" vertex3="v3"/>
      <triangular vertex1="v0" vertex2="v3" vertex3="v2"/>
      <triangular vertex1="v1" vertex2="v2" vertex3="v3"/>
    </tessellated>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="tetra"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    auto geo = RequireSolid(data, "tetra");
    ASSERT_NE(geo, nullptr);

    // Point inside the tetrahedron (centroid at 0.25, 0.25, 0.25)
    siren::math::Vector3D inside(0.1, 0.1, 0.1);
    EXPECT_TRUE(geo->IsInside(inside));

    // Point outside
    siren::math::Vector3D outside(2.0, 2.0, 2.0);
    EXPECT_FALSE(geo->IsInside(outside));
}


TEST(GDMLParser, TessellatedQuadrangular) {
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="v0" x="-1" y="-1" z="-1" unit="m"/>
    <position name="v1" x="1" y="-1" z="-1" unit="m"/>
    <position name="v2" x="1" y="1" z="-1" unit="m"/>
    <position name="v3" x="-1" y="1" z="-1" unit="m"/>
    <position name="v4" x="-1" y="-1" z="1" unit="m"/>
    <position name="v5" x="1" y="-1" z="1" unit="m"/>
    <position name="v6" x="1" y="1" z="1" unit="m"/>
    <position name="v7" x="-1" y="1" z="1" unit="m"/>
  </define>
  <materials/>
  <solids>
    <tessellated name="cube">
      <quadrangular vertex1="v3" vertex2="v2" vertex3="v1" vertex4="v0"/>
      <quadrangular vertex1="v4" vertex2="v5" vertex3="v6" vertex4="v7"/>
      <quadrangular vertex1="v0" vertex2="v1" vertex3="v5" vertex4="v4"/>
      <quadrangular vertex1="v1" vertex2="v2" vertex3="v6" vertex4="v5"/>
      <quadrangular vertex1="v2" vertex2="v3" vertex3="v7" vertex4="v6"/>
      <quadrangular vertex1="v3" vertex2="v0" vertex3="v4" vertex4="v7"/>
    </tessellated>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="cube"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    auto geo = RequireSolid(data, "cube");
    ASSERT_NE(geo, nullptr);

    // Center should be inside
    siren::math::Vector3D center(0, 0, 0);
    EXPECT_TRUE(geo->IsInside(center));

    // Corner region just inside
    siren::math::Vector3D near_corner(0.9, 0.9, 0.9);
    EXPECT_TRUE(geo->IsInside(near_corner));

    // Outside
    siren::math::Vector3D outside(1.5, 0, 0);
    EXPECT_FALSE(geo->IsInside(outside));
}


TEST(GDMLParser, Arb8Basic) {
    // Simple cube via arb8: vertices form a 2x2x2 cube centered at origin
    GDMLData data = ParseGDMLWithSolids(R"(
    <arb8 name="cube8" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="-1000" v5y="-1000"
      v6x="1000" v6y="-1000"
      v7x="1000" v7y="1000"
      v8x="-1000" v8y="1000"
      lunit="mm"/>
)", "cube8");
    auto geo = RequireSolid(data, "cube8");
    ASSERT_NE(geo, nullptr);

    // Center should be inside
    siren::math::Vector3D center(0, 0, 0);
    EXPECT_TRUE(geo->IsInside(center));

    // Point at (0.5, 0.5, 0.5) should be inside (cube is 2m x 2m x 2m)
    siren::math::Vector3D inside(0.5, 0.5, 0.5);
    EXPECT_TRUE(geo->IsInside(inside));

    // Point outside the cube
    siren::math::Vector3D outside(2.0, 0, 0);
    EXPECT_FALSE(geo->IsInside(outside));
}


TEST(GDMLParser, Arb8Trapezoid) {
    // Tapered arb8: bottom face is 2x2, top face is 1x1, height 2
    GDMLData data = ParseGDMLWithSolids(R"(
    <arb8 name="taper" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="-500" v5y="-500"
      v6x="500" v6y="-500"
      v7x="500" v7y="500"
      v8x="-500" v8y="500"
      lunit="mm"/>
)", "taper");
    auto geo = RequireSolid(data, "taper");
    ASSERT_NE(geo, nullptr);

    // Center at origin is inside
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));

    // Near top face, inside the smaller square
    EXPECT_TRUE(IsInside(geo, 0.3, 0.3, 0.9));

    // Near top face, outside the smaller square but inside bottom projection
    EXPECT_FALSE(IsInside(geo, 0.8, 0.8, 0.9));

    // Outside entirely
    EXPECT_FALSE(IsInside(geo, 0, 0, 2.0));
}


TEST(GDMLParser, Arb8Degenerate) {
    // Top face collapses: v5==v6==v7==v8 (pyramid with apex at top)
    GDMLData data = ParseGDMLWithSolids(R"(
    <arb8 name="pyramid" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="0" v5y="0"
      v6x="0" v6y="0"
      v7x="0" v7y="0"
      v8x="0" v8y="0"
      lunit="mm"/>
)", "pyramid");
    auto geo = RequireSolid(data, "pyramid");
    ASSERT_NE(geo, nullptr);

    // Center is inside
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));
    // Near base, inside
    EXPECT_TRUE(IsInside(geo, 0.5, 0.5, -0.8));
    // Outside above apex
    EXPECT_FALSE(IsInside(geo, 0, 0, 2.0));
    // Outside to the side
    EXPECT_FALSE(IsInside(geo, 2.0, 0, 0));
}


TEST(GDMLParser, Arb8Watertight) {
    // Cube via arb8
    GDMLData data = ParseGDMLWithSolids(R"(
    <arb8 name="cube8" dz="1000"
      v1x="-1000" v1y="-1000"
      v2x="1000" v2y="-1000"
      v3x="1000" v3y="1000"
      v4x="-1000" v4y="1000"
      v5x="-1000" v5y="-1000"
      v6x="1000" v6y="-1000"
      v7x="1000" v7y="1000"
      v8x="-1000" v8y="1000"
      lunit="mm"/>
)", "cube8");
    auto geo = data.solids["cube8"];
    ASSERT_TRUE(geo != nullptr);

    // Fire rays that hit shared edges and vertices directly.
    // After deduplication, each should produce exactly 2 forward intersections.

    // Through the center (hits diagonal of both +z and -z faces)
    siren::math::Vector3D origin(0, 0, 5.0);
    siren::math::Vector3D dir(0, 0, -1);
    auto hits = geo->Intersections(origin, dir);
    int forward_hits = 0;
    for(auto const & h : hits) {
        if(h.distance > 0) ++forward_hits;
    }
    EXPECT_EQ(forward_hits, 2) << "Ray through cube center should hit exactly 2 surfaces";

    // Through a vertex (hits corner shared by 3 faces)
    origin = siren::math::Vector3D(5.0, 5.0, 5.0);
    double inv_sqrt3 = 1.0 / std::sqrt(3.0);
    dir = siren::math::Vector3D(-inv_sqrt3, -inv_sqrt3, -inv_sqrt3);
    hits = geo->Intersections(origin, dir);
    forward_hits = 0;
    for(auto const & h : hits) {
        if(h.distance > 0) ++forward_hits;
    }
    EXPECT_EQ(forward_hits, 2) << "Diagonal ray through cube vertex should hit exactly 2 surfaces";

    // Along an edge (parallel to z, hits +x/+y edge shared by 2 faces)
    origin = siren::math::Vector3D(1.0, 0, 5.0);
    dir = siren::math::Vector3D(0, 0, -1);
    hits = geo->Intersections(origin, dir);
    forward_hits = 0;
    for(auto const & h : hits) {
        if(h.distance > 0) ++forward_hits;
    }
    EXPECT_EQ(forward_hits, 2) << "Ray hitting +x face edge should hit exactly 2 surfaces";
}


TEST(GDMLParser, TessellatedCmUnitConversion) {
    // Tessellated solid referencing positions defined in <define> with
    // non-default units (cm). Verifies unit conversion to SIREN meters.
    std::string gdml = R"(<?xml version="1.0"?>
<gdml>
  <define>
    <position name="p0" x="0" y="0" z="0" unit="cm"/>
    <position name="p1" x="100" y="0" z="0" unit="cm"/>
    <position name="p2" x="0" y="100" z="0" unit="cm"/>
    <position name="p3" x="0" y="0" z="100" unit="cm"/>
  </define>
  <materials/>
  <solids>
    <tessellated name="big_tetra">
      <triangular vertex1="p0" vertex2="p2" vertex3="p1"/>
      <triangular vertex1="p0" vertex2="p1" vertex3="p3"/>
      <triangular vertex1="p0" vertex2="p3" vertex3="p2"/>
      <triangular vertex1="p1" vertex2="p2" vertex3="p3"/>
    </tessellated>
  </solids>
  <structure>
    <volume name="World">
      <materialref ref=""/>
      <solidref ref="big_tetra"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="World"/>
  </setup>
</gdml>)";

    GDMLData data = ParseGDMLString(gdml);
    auto geo = RequireSolid(data, "big_tetra");
    ASSERT_NE(geo, nullptr);

    // Vertices are at 0,0,0 and 1m,0,0 and 0,1m,0 and 0,0,1m (100cm = 1m)
    // Centroid at (0.25, 0.25, 0.25) should be inside
    EXPECT_TRUE(IsInside(geo, 0.1, 0.1, 0.1));
    // Point at 0.9m in all axes is outside the tetrahedron
    EXPECT_FALSE(IsInside(geo, 0.9, 0.9, 0.9));
    // Verify unit conversion: point at 50cm = 0.5m should be outside
    // (sum of barycentric coords > 1: 0.5+0.5+0.5 = 1.5 > 1)
    EXPECT_FALSE(IsInside(geo, 0.5, 0.5, 0.5));
}


TEST(GDMLParser, PartialPhiTubeAndCone) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <tube name="quarter_tube" rmin="10" rmax="50" z="100"
          startphi="0" deltaphi="1.5707963267948966" lunit="mm" aunit="rad"/>
    <cone name="quarter_cone" rmin1="10" rmax1="50" rmin2="5" rmax2="30" z="80"
          startphi="0.7853981633974483" deltaphi="1.5707963267948966" lunit="mm" aunit="rad"/>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
)", "world_box");

    // Verify the solids were created successfully
    ASSERT_TRUE(data.solids.find("quarter_tube") != data.solids.end());
    ASSERT_TRUE(data.solids["quarter_tube"] != nullptr);

    ASSERT_TRUE(data.solids.find("quarter_cone") != data.solids.end());
    ASSERT_TRUE(data.solids["quarter_cone"] != nullptr);

    // Verify no warnings about partial angular extent
    ExpectNoWarningContaining(data, "partial angular extent");
}


TEST(GDMLParser, PolyconeZPlaneReversed) {
    // Descending z-planes should be accepted and reversed by the constructor
    GDMLData data = ParseGDMLWithSolids(R"(
    <polycone name="sorted_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="15" z="100"/>
    </polycone>
    <polycone name="reversed_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="15" z="100"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
    </polycone>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
)", "world_box");

    auto sorted = RequireSolid(data, "sorted_pc");
    auto reversed = RequireSolid(data, "reversed_pc");
    ASSERT_NE(sorted, nullptr);
    ASSERT_NE(reversed, nullptr);
    // Both should produce the same geometry
    EXPECT_TRUE(IsInside(sorted, 0, 0, 0.025));
    EXPECT_TRUE(IsInside(reversed, 0, 0, 0.025));
}


TEST(GDMLParser, PolyconeZPlaneDuplicateZ) {
    // Duplicate z-planes create step profiles and should be accepted
    GDMLData data = ParseGDMLWithSolids(R"(
    <polycone name="step_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="30" z="-50"/>
      <zplane rmin="0" rmax="30" z="-20"/>
      <zplane rmin="0" rmax="50" z="-20"/>
      <zplane rmin="0" rmax="50" z="20"/>
      <zplane rmin="0" rmax="30" z="20"/>
      <zplane rmin="0" rmax="30" z="50"/>
    </polycone>
)", "step_pc");

    auto geo = RequireSolid(data, "step_pc");
    ASSERT_NE(geo, nullptr);
    // In the flange: r=0.04 at z=0 should be inside (rmax=0.05 in flange)
    EXPECT_TRUE(IsInside(geo, 0.04, 0, 0));
    // In the narrow section: r=0.04 at z=-0.04 should be outside (rmax=0.03)
    EXPECT_FALSE(IsInside(geo, 0.04, 0, -0.04));
}


TEST(GDMLParser, PolyconeNonMonotonicThrows) {
    // Non-monotonic z-planes are malformed and should throw
    std::string gdml = GDMLWithSolids(R"(
    <polycone name="bad_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="15" z="100"/>
    </polycone>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    )", "world_box");
    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error);
}


TEST(GDMLParser, PolyhedraZPlaneReversed) {
    // Descending z-planes should be accepted for polyhedra too
    GDMLData data = ParseGDMLWithSolids(R"(
    <polyhedra name="sorted_ph" startphi="0" deltaphi="360" numsides="6" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="15" z="100"/>
    </polyhedra>
    <polyhedra name="reversed_ph" startphi="0" deltaphi="360" numsides="6" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="15" z="100"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
    </polyhedra>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
)", "world_box");

    auto sorted = RequireSolid(data, "sorted_ph");
    auto reversed = RequireSolid(data, "reversed_ph");
    ASSERT_NE(sorted, nullptr);
    ASSERT_NE(reversed, nullptr);
    EXPECT_TRUE(IsInside(sorted, 0, 0, 0.025));
    EXPECT_TRUE(IsInside(reversed, 0, 0, 0.025));
}


TEST(GDMLParser, PolyhedraNonMonotonicThrows) {
    // Non-monotonic z-planes should throw
    std::string gdml = GDMLWithSolids(R"(
    <polyhedra name="bad_ph" startphi="0" deltaphi="360" numsides="4" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="15" z="100"/>
      <zplane rmin="0" rmax="5" z="200"/>
    </polyhedra>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    )", "world_box");
    EXPECT_THROW(ParseGDMLString(gdml), std::runtime_error);
}


TEST(GDMLParser, GenericPolyconeHollow) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <genericPolycone name="gpc_hollow" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <rzpoint r="300" z="-500"/>
      <rzpoint r="450" z="0"/>
      <rzpoint r="500" z="500"/>
      <rzpoint r="350" z="500"/>
      <rzpoint r="300" z="0"/>
      <rzpoint r="200" z="-500"/>
    </genericPolycone>
)", "gpc_hollow");
    auto geo = RequireSolid(data, "gpc_hollow");
    auto bb = geo->GetBoundingBox();
    double tol = 1e-6;
    EXPECT_NEAR(bb.max_corner.GetX(), 0.5, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.5, tol);
    EXPECT_NEAR(bb.min_corner.GetZ(), -0.5, tol);

    // Point inside the hollow shell
    EXPECT_TRUE(IsInside(geo, 0.4, 0, 0));
    // Point in the hollow center
    EXPECT_FALSE(IsInside(geo, 0.1, 0, 0));
    // Point outside
    EXPECT_FALSE(IsInside(geo, 0.6, 0, 0));
}


TEST(GDMLParser, GenericPolyconePhiCut) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <genericPolycone name="gpc_phi" startphi="0" deltaphi="180" aunit="deg" lunit="mm">
      <rzpoint r="0" z="-100"/>
      <rzpoint r="50" z="0"/>
      <rzpoint r="0" z="100"/>
    </genericPolycone>
)", "gpc_phi");
    auto geo = RequireSolid(data, "gpc_phi");
    // +x, +y is within phi=[0, pi]
    EXPECT_TRUE(IsInside(geo, 0.01, 0.01, 0));
    // -y is outside phi=[0, pi]
    EXPECT_FALSE(IsInside(geo, 0.01, -0.01, 0));
}


TEST(GDMLParser, BoxContainmentThroughParser) {
    // box x,y,z are GDML full-widths: 50mm, 30mm, 20mm
    // Half-extents in meters: 0.025, 0.015, 0.01
    TestContainment(
        R"(<box name="b" x="50" y="30" z="20" lunit="mm"/>)",
        "b", {
            {0, 0, 0, true},             // Interior: center
            {0.024, 0.014, 0.009, true}, // Interior: just inside corners
            {0.026, 0, 0, false},        // Exterior: just past x half-width
            {0, 0.016, 0, false},        // Exterior: just past y half-width
            {0, 0, 0.011, false},        // Exterior: just past z half-width
        });
}


TEST(GDMLParser, SphereContainmentThroughParser) {
    // Hollow sphere: rmin=20mm=0.02m, rmax=50mm=0.05m
    TestContainment(
        R"(<sphere name="s" rmin="20" rmax="50" lunit="mm"/>)",
        "s", {
            {0.035, 0, 0, true},        // Interior: between inner and outer radius
            {0.015, 0, 0, false},       // Exterior: inside hollow core (r<rmin)
            {0.051, 0, 0, false},       // Exterior: past outer radius
            {0.02, 0.02, 0.02, true},   // Off-axis interior: r=0.0346 in [0.02, 0.05]
            {0.04, 0.03, 0.01, false},  // Off-axis exterior: r~0.05 ~ rmax boundary
        });
}


TEST(GDMLParser, TubeContainmentThroughParser) {
    // Hollow tube: rmin=10mm=0.01m, rmax=50mm=0.05m
    // GDML z=100mm is full-height -> half-height=0.05m
    TestContainment(
        R"(<tube name="t" rmin="10" rmax="50" z="100" deltaphi="360" aunit="deg" lunit="mm"/>)",
        "t", {
            {0.03, 0, 0, true},         // Interior: in shell, at center height
            {0, 0.03, 0.04, true},      // Interior: in shell, near top (|z|<0.05)
            {0.005, 0, 0, false},       // Exterior: inside hollow (r<rmin)
            {0.051, 0, 0, false},       // Exterior: past outer radius
            {0.03, 0, 0.051, false},    // Exterior: past z half-height
        });
}


TEST(GDMLParser, ConeContainmentThroughParser) {
    // Solid cone: rmax1=50mm at bottom, rmax2=20mm at top
    // GDML z=100mm is full-height -> spans z=-0.05 to z=0.05
    // rmax(z) = 0.035 - 0.3*z. At z=0: 0.035m; z=-0.04: 0.047m; z=+0.04: 0.023m
    TestContainment(
        R"(<cone name="c" rmin1="0" rmax1="50" rmin2="0" rmax2="20" z="100" deltaphi="360" aunit="deg" lunit="mm"/>)",
        "c", {
            {0, 0, 0, true},            // Interior: center
            {0.04, 0, -0.04, true},     // Interior: near bottom (r=0.04<rmax~0.047)
            {0.04, 0, 0.04, false},     // Exterior: near top (r=0.04>rmax~0.023)
            {0.02, 0.02, 0, true},      // Off-axis interior: r=0.0283 < 0.035
            {0.03, 0.03, 0, false},     // Off-axis exterior: r=0.0424 > 0.035
        });
}


TEST(GDMLParser, TrdContainmentThroughParser) {
    // Trd: GDML full-widths x1=60,x2=30,y1=40,y2=20,z=50 mm -> half-lengths
    // dx1=0.03,dx2=0.015,dy1=0.02,dy2=0.01,dz=0.025 (spans z=-0.025..0.025)
    // At z=0 (center): x half-width=0.0225m, y half-width=0.015m
    TestContainment(
        R"(<trd name="tr" x1="60" x2="30" y1="40" y2="20" z="50" lunit="mm"/>)",
        "tr", {
            {0.02, 0.012, 0, true},     // Interior: center region
            {0.025, 0, 0, false},       // Exterior: x>interpolated half-width at z=0
            {0.025, 0, -0.02, true},    // Interior: near bottom, x<~0.0285m
            {0.025, 0, 0.02, false},    // Exterior: near top, x>~0.0165m
        });
}


TEST(GDMLParser, EllipticalTubeContainmentThroughParser) {
    // eltube dx=40mm=0.04m, dy=20mm=0.02m, dz=50mm=0.05m
    // Cross-section: (x/dx)^2 + (y/dy)^2 < 1, |z| < dz
    TestContainment(
        R"(<eltube name="et" dx="40" dy="20" dz="50" lunit="mm"/>)",
        "et", {
            {0.03, 0, 0, true},         // Interior: (0.03/0.04)^2 = 0.5625 < 1
            {0.03, 0.015, 0, false},    // Exterior: 0.5625 + 0.5625 = 1.125 > 1
            {0, 0, 0.051, false},       // Exterior: past z extent
        });
}


TEST(GDMLParser, TorusContainmentThroughParser) {
    // Solid torus: rmin=0, rmax=20mm=0.02m, rtor=100mm=0.1m
    TestContainment(
        R"(<torus name="tor" rmin="0" rmax="20" rtor="100" startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>)",
        "tor", {
            {0.1, 0, 0, true},          // Interior: on major radius circle
            {0.11, 0, 0, true},         // Interior: 0.01m from ring center < rmax
            {0, 0, 0, false},           // Exterior: origin is 0.1m from ring > rmax
            {0.1, 0, 0.025, false},     // Exterior: z=0.025 from ring center > rmax
            // Off-axis interior: ring at phi=45, center at (0.1/sqrt2, 0.1/sqrt2, 0)
            {0.1/1.41421356, 0.1/1.41421356, 0, true},
        });
}


TEST(GDMLParser, TrapContainmentThroughParser) {
    // Trap: GDML full-lengths z=50,y1=30,x1=x2=20,y2=15,x3=x4=10 mm -> half-lengths
    // dz=0.025m, theta=phi=0. At z=-dz: dy1=0.015m, dx1=dx2=0.01m
    // At z=+dz: dy2=0.0075m, dx3=dx4=0.005m
    TestContainment(
        R"(<trap name="trp" z="50" theta="0" phi="0" y1="30" x1="20" x2="20" alpha1="0" y2="15" x3="10" x4="10" alpha2="0" aunit="rad" lunit="mm"/>)",
        "trp", {
            {0, 0, 0, true},            // Interior: center
            // Interior: near bottom face (z=-0.02), x_half~0.0095, y_half~0.0143
            {0.008, 0.012, -0.02, true},
            // Exterior: near top face (z=+0.02), x_half~0.0055, x=0.008>0.0055
            {0.008, 0, 0.02, false},
        });
}


TEST(GDMLParser, PolyconeContainmentThroughParser) {
    // Polycone with 3 z-planes:
    //   z=-50mm rmin=0 rmax=50mm
    //   z=0     rmin=0 rmax=80mm
    //   z=+50mm rmin=0 rmax=30mm
    TestContainment(
        R"(<polycone name="pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-50"/>
      <zplane rmin="0" rmax="80" z="0"/>
      <zplane rmin="0" rmax="30" z="50"/>
    </polycone>)",
        "pc", {
            {0.06, 0, 0, true},         // Interior: at z=0, r=0.06 < rmax=0.08
            {0.04, 0, -0.04, true},     // Interior: at z=-0.04m, r<rmax~0.056
            {0.06, 0, 0.04, false},     // Exterior: at z=0.04m, r>rmax~0.04m
            {0.03, 0.04, 0, true},      // Off-axis interior: r=0.05 < rmax=0.08 at z=0
        });
}


TEST(GDMLParser, EllipsoidContainmentThroughParser) {
    // Ellipsoid: ax=50mm=0.05m, by=30mm=0.03m, cz=40mm=0.04m
    // Equation: (x/0.05)^2 + (y/0.03)^2 + (z/0.04)^2 < 1
    TestContainment(
        R"(<ellipsoid name="ell" ax="50" by="30" cz="40" lunit="mm"/>)",
        "ell", {
            {0, 0, 0, true},            // Interior: origin
            // Interior: 0.36 + 0.1111 + 0.0625 = 0.5336 < 1
            {0.03, 0.01, 0.01, true},
            // Exterior: 0.64 + 0.4444 + 0.5625 = 1.647 > 1
            {0.04, 0.02, 0.03, false},
        });
}


TEST(GDMLParser, PartialSphereContainmentThroughParser) {
    // Upper hemisphere: starttheta=0, deltatheta=90 deg
    // rmax = 100mm = 0.1m
    TestContainment(
        R"(<sphere name="hemi" rmax="100" starttheta="0" deltatheta="90" startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>)",
        "hemi", {
            {0, 0, 0.05, true},         // Interior: upper hemisphere (z>0, r<0.1)
            {0, 0, -0.05, false},       // Exterior: lower hemisphere (theta > 90 deg)
        });
}


TEST(GDMLParser, PartialTorusContainmentThroughParser) {
    // Quarter torus: deltaphi=90 deg (phi from 0 to 90 deg)
    // rmin=0, rmax=20mm=0.02m, rtor=100mm=0.1m
    TestContainment(
        R"(<torus name="qt" rmin="0" rmax="20" rtor="100" startphi="0" deltaphi="90" lunit="mm" aunit="deg"/>)",
        "qt", {
            {0.1, 0, 0, true},          // Interior: at phi=0, on the ring
            {0.0707, 0.0707, 0, true},  // Interior: at phi=45 deg
            {-0.1, 0, 0, false},        // Exterior: at phi=180 deg
            {0, -0.1, 0, false},        // Exterior: at phi=270 deg
        });
}


TEST(GDMLParser, PartialCylinderContainmentThroughParser) {
    // Half cylinder: deltaphi=180 deg (phi from 0 to 180 deg, y >= 0 half)
    // rmax=50mm=0.05m, z=100mm half-height=0.05m
    TestContainment(
        R"(<tube name="hc" rmin="0" rmax="50" z="100" startphi="0" deltaphi="180" lunit="mm" aunit="deg"/>)",
        "hc", {
            {0.03, 0.03, 0, true},      // Interior: phi~45 deg, within 180 deg cut
            {0.03, -0.03, 0, false},    // Exterior: phi~315 deg, outside 180 deg cut
        });
}


TEST(GDMLParser, TrapShearedContainment) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <trap name="sheared_trap" z="100" theta="20" phi="0"
          y1="30" x1="20" x2="20" alpha1="10"
          y2="30" x3="20" x4="20" alpha2="10"
          aunit="deg" lunit="mm"/>
)", "world_box");
    auto geo = RequireSolid(data, "sheared_trap");
    ASSERT_NE(geo, nullptr);

    // GDML full-lengths -> dz=0.05m, dx=0.01m, dy=0.015m half-lengths.
    // Center of the solid: should be inside regardless of theta
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));

    // Near top face (z close to +dz=0.05m): theta=20deg shifts the face
    // center to x ~ +0.0178m, x half-width 0.01m -> x range [0.0078, 0.0278]
    EXPECT_TRUE(IsInside(geo, 0.018, 0, 0.049));
    // Same x at bottom face (z=-0.049m): face center ~ x=-0.0178m,
    // range [-0.0278, -0.0078], so x=0.018 is outside
    EXPECT_FALSE(IsInside(geo, 0.018, 0, -0.049));
}


TEST(GDMLParser, PartialPolyhedraContainment) {
    // 6-sided polyhedra with deltaphi=180deg (half-shell)
    // rmax=50mm = 0.05m
    TestContainment(
        R"(<polyhedra name="half_hex" numsides="6" startphi="0" deltaphi="180" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-40"/>
      <zplane rmin="0" rmax="50" z="40"/>
    </polyhedra>)",
        "half_hex", {
            {0, 0.03, 0, true},         // Interior: positive-y half (phi~90, within [0,180])
            {0.03, 0, 0, true},         // Interior: along +x axis (phi=0, boundary)
            {0, -0.03, 0, false},       // Exterior: negative-y half (phi~270, outside [0,180])
            {0, 0.06, 0, false},        // Exterior: beyond the solid radius
        });
}


TEST(GDMLParser, SolidHalfLengthConventions) {
    // Regression guard for the GDML->SIREN dimension conventions (Geant4 ground
    // truth). GDML box/tube/cone/para/trd/trap/cutTube linear attributes are
    // FULL lengths, so each SIREN bounding-box half-extent must equal half the
    // GDML value. eltube/ellipsoid attributes are semi-axes (no halving).
    // Polyhedra rmin/rmax are apothem (tangent) distances, which SIREN stores and
    // uses directly (pass-through); it only scales by 1/cos(pi/n) when it needs a
    // circumradius (e.g. for vertices/AABB). All dimensions in meters for legibility.
    GDMLData data = ParseGDMLWithSolids(R"(
    <box name="cv_box" x="2" y="4" z="6" lunit="m"/>
    <tube name="cv_tube" rmin="0" rmax="1" z="4" deltaphi="360" aunit="deg" lunit="m"/>
    <cone name="cv_cone" rmin1="0" rmax1="1" rmin2="0" rmax2="1" z="4" deltaphi="360" aunit="deg" lunit="m"/>
    <para name="cv_para" x="2" y="2" z="2" alpha="0" theta="0" phi="0" aunit="rad" lunit="m"/>
    <trd name="cv_trd" x1="2" x2="2" y1="2" y2="2" z="2" lunit="m"/>
    <trap name="cv_trap" z="2" theta="0" phi="0" y1="2" x1="2" x2="2" alpha1="0" y2="2" x3="2" x4="2" alpha2="0" aunit="rad" lunit="m"/>
    <cutTube name="cv_cuttube" rmin="0" rmax="1" z="4" startphi="0" deltaphi="360" lowX="0" lowY="0" lowZ="-1" highX="0" highY="0" highZ="1" aunit="deg" lunit="m"/>
    <eltube name="cv_eltube" dx="1" dy="2" dz="3" lunit="m"/>
    <ellipsoid name="cv_ellipsoid" ax="1" by="2" cz="3" lunit="m"/>
    <polyhedra name="cv_polyhedra" numsides="4" startphi="45" deltaphi="360" aunit="deg" lunit="m">
      <zplane rmin="0" rmax="1" z="-1"/>
      <zplane rmin="0" rmax="1" z="1"/>
    </polyhedra>
)", "cv_box");

    double tol = 1e-6;

    // box: GDML full-widths 2,4,6 m -> half-extents 1,2,3 m
    auto bb = data.solids["cv_box"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 2.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 3.0, tol);

    // tube: GDML z=4 m full-height -> half-height 2 m; rmax=1 m actual radius
    bb = data.solids["cv_tube"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 2.0, tol);

    // cone: GDML z=4 m full-height -> half-height 2 m; rmax=1 m actual radius
    bb = data.solids["cv_cone"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 2.0, tol);

    // para: GDML full-lengths 2,2,2 m (axis-aligned) -> half-extents 1,1,1 m
    bb = data.solids["cv_para"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 1.0, tol);

    // trd: GDML full-widths 2,2,2,2,2 m -> half-extents 1,1,1 m
    bb = data.solids["cv_trd"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 1.0, tol);

    // trap: GDML full-lengths (orthogonal) -> half-extents 1,1,1 m
    bb = data.solids["cv_trap"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 1.0, tol);

    // cutTube: GDML z=4 m full-height -> half-height 2 m (flat cuts); rmax=1 m
    bb = data.solids["cv_cuttube"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 2.0, tol);

    // eltube: GDML dx,dy,dz are semi-axes (half) -> bbox 1,2,3 m (no halving)
    bb = data.solids["cv_eltube"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 2.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 3.0, tol);

    // ellipsoid: GDML ax,by,cz are semi-axes (half) -> bbox 1,2,3 m (no halving)
    bb = data.solids["cv_ellipsoid"]->GetBoundingBox();
    EXPECT_NEAR(bb.max_corner.GetX(), 1.0, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 2.0, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 3.0, tol);

    // polyhedra: 4-sided, startphi=45deg so the flat faces are axis-aligned at
    // x,y = +/- the apothem. GDML rmax=1 m is the apothem (tangent distance to
    // the flat face) and SIREN's Polyhedra uses the same apothem convention, so
    // the right face sits at x=1 m: a point at x=0.9 is inside, x=1.1 is outside.
    // This pins the rmin/rmax pass-through (no spurious apothem<->corner
    // rescaling). Containment is checked rather than the bounding box, which is
    // a conservative circumscribed-square bound.
    auto poly = data.solids["cv_polyhedra"];
    EXPECT_TRUE(IsInside(poly, 0.9, 0, 0));
    EXPECT_FALSE(IsInside(poly, 1.1, 0, 0));
}
