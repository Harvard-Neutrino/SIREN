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
    // GDML box x/y/z are half-widths: 100mm half-width = 0.1m
    auto bb_b1 = data.solids["b1"]->GetBoundingBox();
    EXPECT_NEAR(bb_b1.max_corner.GetX(), 0.1, tol);
    EXPECT_NEAR(bb_b1.max_corner.GetZ(), 0.1, tol);

    // sphere: rmax=50mm = 0.05m
    auto bb_s1 = data.solids["s1"]->GetBoundingBox();
    EXPECT_NEAR(bb_s1.max_corner.GetX(), 0.05, tol);

    // tube: rmax=50mm=0.05m, GDML z=100mm is half-length=0.1m
    auto bb_t1 = data.solids["t1"]->GetBoundingBox();
    EXPECT_NEAR(bb_t1.max_corner.GetX(), 0.05, tol);
    EXPECT_NEAR(bb_t1.max_corner.GetZ(), 0.1, tol);

    // trd: GDML half-widths: max(x1=60,x2=40)=60mm=0.06m, z=100mm half-height=0.1m
    auto bb_trd = data.solids["trd1"]->GetBoundingBox();
    EXPECT_NEAR(bb_trd.max_corner.GetX(), 0.06, tol);
    EXPECT_NEAR(bb_trd.max_corner.GetZ(), 0.1, tol);

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
    std::string tmpfile = "/tmp/siren_gdml_test_strict.gdml";
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
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("torus") == std::string::npos)
            << "Unexpected torus warning: " << w;
    }
}


TEST(GDMLParser, PartialSphereParses) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <sphere name="hemi" rmin="0" rmax="100" starttheta="0" deltatheta="90" startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>
)", "hemi");
    ASSERT_TRUE(data.solids.count("hemi") > 0);
    // No warnings about partial sphere (it is now supported)
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("sphere") == std::string::npos && w.find("angular") == std::string::npos)
            << "Unexpected sphere angular warning: " << w;
    }
}


TEST(GDMLParser, PartialTorusParses) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <torus name="elbow" rmin="10" rmax="20" rtor="100" startphi="0" deltaphi="90" lunit="mm" aunit="deg"/>
)", "elbow");
    ASSERT_TRUE(data.solids.count("elbow") > 0);
    // No warnings about partial torus (it is now supported)
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("torus") == std::string::npos && w.find("angular") == std::string::npos)
            << "Unexpected torus angular warning: " << w;
    }
}


TEST(GDMLParser, PartialCylinderSupported) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <tube name="wedge" rmin="0" rmax="100" z="200" deltaphi="180" lunit="mm" aunit="deg"/>
)", "wedge");
    ASSERT_TRUE(data.solids.count("wedge") > 0);
    ASSERT_TRUE(data.solids["wedge"] != nullptr);
    // Partial phi is now supported -- no warnings expected
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("partial angular extent") == std::string::npos)
            << "Unexpected warning: " << w;
    }
}


TEST(GDMLParser, EllipticalTubeIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <eltube name="et" dx="30" dy="20" dz="50" lunit="mm"/>
)", "et");
    ASSERT_TRUE(data.solids.count("et") > 0);
    auto bb = data.solids["et"]->GetBoundingBox();
    double tol = 1e-9;
    // dx=30mm=0.03m, dy=20mm=0.02m, dz=50mm=0.05m
    EXPECT_NEAR(bb.max_corner.GetX(), 0.03, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 0.02, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.05, tol);
}


TEST(GDMLParser, GenericPolycone) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <genericPolycone name="gpc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <rzpoint r="50" z="-100"/>
      <rzpoint r="80" z="0"/>
      <rzpoint r="30" z="100"/>
    </genericPolycone>
)", "gpc");
    ASSERT_TRUE(data.solids.count("gpc") > 0);
    auto bb = data.solids["gpc"]->GetBoundingBox();
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
    ASSERT_TRUE(data.solids.count("ct") > 0);
    auto bb = data.solids["ct"]->GetBoundingBox();
    double tol = 1e-6;
    // rmax=50mm=0.05m, z=100mm half=50mm=0.05m (flat cuts)
    EXPECT_NEAR(bb.max_corner.GetX(), 0.05, tol);
}


TEST(GDMLParser, TrapIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <trap name="t" z="100" theta="0" phi="0"
          y1="60" x1="40" x2="40" alpha1="0"
          y2="40" x3="30" x4="30" alpha2="0"
          aunit="rad" lunit="mm"/>
)", "t");
    ASSERT_TRUE(data.solids.count("t") > 0);
    auto bb = data.solids["t"]->GetBoundingBox();
    double tol = 1e-6;
    // z=100mm = 0.1m (GDML trap z is half-height, passed directly as dz)
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.1, tol);
    // max x = max(dx1,dx2,dx3,dx4) = 40mm = 0.04m
    EXPECT_NEAR(bb.max_corner.GetX(), 0.04, tol);
}


TEST(GDMLParser, EllipsoidIntegration) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <ellipsoid name="ell" ax="50" by="30" cz="40" zcut2="20" lunit="mm"/>
)", "ell");
    ASSERT_TRUE(data.solids.count("ell") > 0);
    auto bb = data.solids["ell"]->GetBoundingBox();
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
    ASSERT_TRUE(data.solids.count("p") > 0);
    auto bb = data.solids["p"]->GetBoundingBox();
    double tol = 1e-6;
    // With alpha=theta=phi=0, para = box: x=40mm=0.04m half, etc.
    EXPECT_NEAR(bb.max_corner.GetX(), 0.04, tol);
    EXPECT_NEAR(bb.max_corner.GetY(), 0.03, tol);
    EXPECT_NEAR(bb.max_corner.GetZ(), 0.05, tol);
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


TEST(GDMLParser, TessellatedInlineVertices) {
    // GDML spec allows <triangular> to reference positions defined in <define>
    // but some generators put position elements directly inside the tessellated.
    // However, the standard approach is always through define references.
    // This test verifies the standard approach works with unit conversion.
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
    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("partial angular extent") == std::string::npos)
            << "Unexpected warning: " << w;
    }
}


TEST(GDMLParser, PolyconeZPlaneSorting) {
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

    ASSERT_TRUE(data.solids.find("sorted_pc") != data.solids.end());
    ASSERT_TRUE(data.solids["sorted_pc"] != nullptr);
    ASSERT_TRUE(data.solids.find("reversed_pc") != data.solids.end());
    ASSERT_TRUE(data.solids["reversed_pc"] != nullptr);

    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("non-monotonic") == std::string::npos)
            << "Unexpected warning: " << w;
    }
}


TEST(GDMLParser, PolyconeZPlaneDuplicateZFails) {
    GDMLData data = ParseGDMLWithSolids(R"(
    <polycone name="dup_z_pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="10" z="50"/>
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="15" z="100"/>
    </polycone>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
)", "world_box");

    EXPECT_TRUE(data.solids.find("dup_z_pc") == data.solids.end()
                || data.solids["dup_z_pc"] == nullptr);

    bool found_warning = false;
    for(auto const & w : data.warnings) {
        if(w.find("non-monotonic") != std::string::npos) { found_warning = true; break; }
    }
    EXPECT_TRUE(found_warning);
}


TEST(GDMLParser, PolyhedraZPlaneSorting) {
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
    <polyhedra name="scrambled_ph" startphi="0" deltaphi="360" numsides="4" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="20" z="50"/>
      <zplane rmin="0" rmax="10" z="0"/>
      <zplane rmin="0" rmax="15" z="100"/>
      <zplane rmin="0" rmax="5" z="200"/>
    </polyhedra>
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
)", "world_box");

    ASSERT_TRUE(data.solids.find("sorted_ph") != data.solids.end());
    ASSERT_TRUE(data.solids["sorted_ph"] != nullptr);
    ASSERT_TRUE(data.solids.find("reversed_ph") != data.solids.end());
    ASSERT_TRUE(data.solids["reversed_ph"] != nullptr);
    ASSERT_TRUE(data.solids.find("scrambled_ph") != data.solids.end());
    ASSERT_TRUE(data.solids["scrambled_ph"] != nullptr);

    for(auto const & w : data.warnings) {
        EXPECT_TRUE(w.find("non-monotonic") == std::string::npos)
            << "Unexpected warning: " << w;
    }
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
    ASSERT_TRUE(data.solids.count("gpc_hollow") > 0);
    auto & geo = data.solids["gpc_hollow"];
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
    ASSERT_TRUE(data.solids.count("gpc_phi") > 0);
    auto & geo = data.solids["gpc_phi"];
    // +x, +y is within phi=[0, pi]
    EXPECT_TRUE(IsInside(geo, 0.01, 0.01, 0));
    // -y is outside phi=[0, pi]
    EXPECT_FALSE(IsInside(geo, 0.01, -0.01, 0));
}


TEST(GDMLParser, BoxContainmentThroughParser) {
    // box x,y,z are GDML half-widths: 50mm, 30mm, 20mm
    // Parser doubles them; SIREN stores full-widths internally.
    // Half-extents in meters: 0.05, 0.03, 0.02
    GDMLData data = ParseGDMLWithSolids(R"(
    <box name="b" x="50" y="30" z="20" lunit="mm"/>
)", "b");
    auto geo = RequireSolid(data, "b");
    ASSERT_NE(geo, nullptr);

    // Interior: center
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));
    // Interior: just inside corners
    EXPECT_TRUE(IsInside(geo, 0.049, 0.029, 0.019));
    // Exterior: just past x half-width
    EXPECT_FALSE(IsInside(geo, 0.051, 0, 0));
    // Exterior: just past y half-width
    EXPECT_FALSE(IsInside(geo, 0, 0.031, 0));
    // Exterior: just past z half-width
    EXPECT_FALSE(IsInside(geo, 0, 0, 0.021));
}


TEST(GDMLParser, SphereContainmentThroughParser) {
    // Hollow sphere: rmin=20mm=0.02m, rmax=50mm=0.05m
    GDMLData data = ParseGDMLWithSolids(R"(
    <sphere name="s" rmin="20" rmax="50" lunit="mm"/>
)", "s");
    auto geo = RequireSolid(data, "s");
    ASSERT_NE(geo, nullptr);

    // Interior: between inner and outer radius
    EXPECT_TRUE(IsInside(geo, 0.035, 0, 0));
    // Exterior: inside hollow core (r=0.015 < rmin=0.02)
    EXPECT_FALSE(IsInside(geo, 0.015, 0, 0));
    // Exterior: past outer radius (r=0.051 > rmax=0.05)
    EXPECT_FALSE(IsInside(geo, 0.051, 0, 0));
    // Off-axis interior: r=sqrt(0.02^2+0.02^2+0.02^2)=0.0346 in [0.02, 0.05]
    EXPECT_TRUE(IsInside(geo, 0.02, 0.02, 0.02));
    // Off-axis exterior: r=sqrt(0.04^2+0.03^2)=0.05 ~ rmax boundary
    EXPECT_FALSE(IsInside(geo, 0.04, 0.03, 0.01));
}


TEST(GDMLParser, TubeContainmentThroughParser) {
    // Hollow tube: rmin=10mm=0.01m, rmax=50mm=0.05m
    // z=100mm is GDML half-height; parser doubles it -> full height 200mm, half=0.1m
    GDMLData data = ParseGDMLWithSolids(R"(
    <tube name="t" rmin="10" rmax="50" z="100" deltaphi="360" aunit="deg" lunit="mm"/>
)", "t");
    auto geo = RequireSolid(data, "t");
    ASSERT_NE(geo, nullptr);

    // Interior: in shell, at center height
    EXPECT_TRUE(IsInside(geo, 0.03, 0, 0));
    // Interior: in shell, near top (z=0.09 < half-height 0.1)
    EXPECT_TRUE(IsInside(geo, 0, 0.03, 0.09));
    // Exterior: inside hollow (r=0.005 < rmin=0.01)
    EXPECT_FALSE(IsInside(geo, 0.005, 0, 0));
    // Exterior: past outer radius
    EXPECT_FALSE(IsInside(geo, 0.051, 0, 0));
    // Exterior: past z half-height (0.101 > 0.1)
    EXPECT_FALSE(IsInside(geo, 0.03, 0, 0.101));
}


TEST(GDMLParser, ConeContainmentThroughParser) {
    // Solid cone: rmax1=50mm at z=-0.05m, rmax2=20mm at z=+0.05m
    // z=100mm is GDML half-height, so full height = 0.2m, spanning z=-0.1 to z=0.1
    // Wait -- parser doubles z for cone, so z goes from -0.1 to +0.1
    // rmax linearly interpolates: rmax(z) = 50 + (20-50)*(z+0.1)/0.2 mm
    //   = 50 - 150*(z+0.1) mm  (z in meters)
    // At z=0: rmax = 50 - 150*0.1 = 35mm = 0.035m
    // At z=-0.09: rmax = 50 - 150*0.01 = 48.5mm = 0.0485m
    // At z=+0.09: rmax = 50 - 150*0.19 = 21.5mm = 0.0215m
    GDMLData data = ParseGDMLWithSolids(R"(
    <cone name="c" rmin1="0" rmax1="50" rmin2="0" rmax2="20" z="100"
          deltaphi="360" aunit="deg" lunit="mm"/>
)", "c");
    auto geo = RequireSolid(data, "c");
    ASSERT_NE(geo, nullptr);

    // Interior: center (r=0 < rmax=0.035m at z=0)
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));
    // Interior: near bottom (r=0.04 < rmax~0.0485m at z=-0.09)
    EXPECT_TRUE(IsInside(geo, 0.04, 0, -0.09));
    // Exterior: near top (r=0.04 > rmax~0.0215m at z=+0.09)
    EXPECT_FALSE(IsInside(geo, 0.04, 0, 0.09));
    // Off-axis interior: r=sqrt(0.02^2+0.02^2)=0.0283 < rmax~0.035 at z=0
    EXPECT_TRUE(IsInside(geo, 0.02, 0.02, 0));
    // Off-axis exterior: r=sqrt(0.03^2+0.03^2)=0.0424 > rmax~0.035 at z=0
    EXPECT_FALSE(IsInside(geo, 0.03, 0.03, 0));
}


TEST(GDMLParser, TrdContainmentThroughParser) {
    // Trd: x1=60, x2=30, y1=40, y2=20, z=50 (all GDML half-widths/half-height, mm)
    // At z=-0.05m (bottom face): x half-width=0.06m, y half-width=0.04m
    // At z=+0.05m (top face): x half-width=0.03m, y half-width=0.02m
    // At z=0 (center): x half-width=0.045m, y half-width=0.03m
    GDMLData data = ParseGDMLWithSolids(R"(
    <trd name="tr" x1="60" x2="30" y1="40" y2="20" z="50" lunit="mm"/>
)", "tr");
    auto geo = RequireSolid(data, "tr");
    ASSERT_NE(geo, nullptr);

    // Interior: center region (|x|=0.04 < 0.045, |y|=0.025 < 0.03)
    EXPECT_TRUE(IsInside(geo, 0.04, 0.025, 0));
    // Exterior: at center height, x=0.046 > 0.045 (interpolated half-width)
    EXPECT_FALSE(IsInside(geo, 0.046, 0, 0));
    // Interior: near bottom z=-0.045, x half-width ~ 0.057m, 0.055 < 0.057
    EXPECT_TRUE(IsInside(geo, 0.055, 0, -0.045));
    // Exterior: near top z=+0.045, x half-width ~ 0.033m, 0.055 > 0.033
    EXPECT_FALSE(IsInside(geo, 0.055, 0, 0.045));
}


TEST(GDMLParser, EllipticalTubeContainmentThroughParser) {
    // eltube dx=40mm=0.04m, dy=20mm=0.02m, dz=50mm=0.05m
    // Cross-section: (x/dx)^2 + (y/dy)^2 < 1, |z| < dz
    GDMLData data = ParseGDMLWithSolids(R"(
    <eltube name="et" dx="40" dy="20" dz="50" lunit="mm"/>
)", "et");
    auto geo = RequireSolid(data, "et");
    ASSERT_NE(geo, nullptr);

    // Interior: (0.03/0.04)^2 = 0.5625 < 1
    EXPECT_TRUE(IsInside(geo, 0.03, 0, 0));
    // Exterior: (0.03/0.04)^2 + (0.015/0.02)^2 = 0.5625 + 0.5625 = 1.125 > 1
    EXPECT_FALSE(IsInside(geo, 0.03, 0.015, 0));
    // Exterior: past z extent
    EXPECT_FALSE(IsInside(geo, 0, 0, 0.051));
}


TEST(GDMLParser, TorusContainmentThroughParser) {
    // Solid torus: rmin=0, rmax=20mm=0.02m, rtor=100mm=0.1m
    // Points on the major radius circle (distance 0.1m from z-axis) are
    // within the tube cross-section.
    GDMLData data = ParseGDMLWithSolids(R"(
    <torus name="tor" rmin="0" rmax="20" rtor="100" startphi="0" deltaphi="360"
           lunit="mm" aunit="deg"/>
)", "tor");
    auto geo = RequireSolid(data, "tor");
    ASSERT_NE(geo, nullptr);

    // Interior: on major radius circle (distance from ring center = 0)
    EXPECT_TRUE(IsInside(geo, 0.1, 0, 0));
    // Interior: 0.01m from ring center (0.11 - 0.1 = 0.01 < rmax=0.02)
    EXPECT_TRUE(IsInside(geo, 0.11, 0, 0));
    // Exterior: origin is 0.1m from ring, 0.1 > rmax=0.02
    EXPECT_FALSE(IsInside(geo, 0, 0, 0));
    // Exterior: z=0.025m from ring center at (0.1,0,0): 0.025 > rmax=0.02
    EXPECT_FALSE(IsInside(geo, 0.1, 0, 0.025));
    // Off-axis interior: ring at phi=45, center at (0.1/sqrt2, 0.1/sqrt2, 0)
    double s2 = 1.0 / std::sqrt(2.0);
    EXPECT_TRUE(IsInside(geo, 0.1*s2, 0.1*s2, 0));
}


TEST(GDMLParser, TrapContainmentThroughParser) {
    // Trap: z=50mm dz=0.05m. theta=phi=0.
    // At z=-dz: dy1=0.03m, dx1=dx2=0.02m
    // At z=+dz: dy2=0.015m, dx3=dx4=0.01m
    // At z=0: interpolated dy~0.0225m, dx~0.015m
    GDMLData data = ParseGDMLWithSolids(R"(
    <trap name="trp" z="50" theta="0" phi="0"
          y1="30" x1="20" x2="20" alpha1="0"
          y2="15" x3="10" x4="10" alpha2="0"
          aunit="rad" lunit="mm"/>
)", "trp");
    auto geo = RequireSolid(data, "trp");
    ASSERT_NE(geo, nullptr);

    // Interior: center
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));
    // Interior: near bottom face (z=-0.04). Interpolated at z=-0.04:
    //   t = (-0.04+0.05)/0.1 = 0.1, dx = 0.02*(1-0.1)+0.01*0.1 = 0.019
    //   dy = 0.03*(1-0.1)+0.015*0.1 = 0.0285
    //   x=0.015 < 0.019, y=0.02 < 0.0285 -> inside
    EXPECT_TRUE(IsInside(geo, 0.015, 0.02, -0.04));
    // Exterior: near top face (z=+0.04). Interpolated:
    //   t = (0.04+0.05)/0.1 = 0.9, dx = 0.02*0.1+0.01*0.9 = 0.011
    //   x=0.015 > 0.011 -> outside
    EXPECT_FALSE(IsInside(geo, 0.015, 0, 0.04));
}


TEST(GDMLParser, PolyconeContainmentThroughParser) {
    // Polycone with 3 z-planes:
    //   z=-50mm rmin=0 rmax=50mm
    //   z=0     rmin=0 rmax=80mm
    //   z=+50mm rmin=0 rmax=30mm
    // At z=0: rmax=0.08m
    // At z=40mm (between 0 and 50mm): rmax = 80 + (30-80)*(40/50) = 80-40 = 40mm = 0.04m
    // At z=-40mm (between -50 and 0mm): rmax = 50 + (80-50)*(10/50) = 56mm = 0.056m
    GDMLData data = ParseGDMLWithSolids(R"(
    <polycone name="pc" startphi="0" deltaphi="360" aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-50"/>
      <zplane rmin="0" rmax="80" z="0"/>
      <zplane rmin="0" rmax="30" z="50"/>
    </polycone>
)", "pc");
    auto geo = RequireSolid(data, "pc");
    ASSERT_NE(geo, nullptr);

    // Interior: at z=0, r=0.06 < rmax=0.08
    EXPECT_TRUE(IsInside(geo, 0.06, 0, 0));
    // Interior: at z=-0.04m, r=0.04 < rmax~0.056
    EXPECT_TRUE(IsInside(geo, 0.04, 0, -0.04));
    // Exterior: at z=0.04m, r=0.06 > rmax~0.04m (interpolated between 80mm and 30mm)
    EXPECT_FALSE(IsInside(geo, 0.06, 0, 0.04));
    // Off-axis interior: r=sqrt(0.03^2+0.04^2)=0.05 < rmax=0.08 at z=0
    EXPECT_TRUE(IsInside(geo, 0.03, 0.04, 0));
}


TEST(GDMLParser, EllipsoidContainmentThroughParser) {
    // Ellipsoid: ax=50mm=0.05m, by=30mm=0.03m, cz=40mm=0.04m
    // Equation: (x/0.05)^2 + (y/0.03)^2 + (z/0.04)^2 < 1
    GDMLData data = ParseGDMLWithSolids(R"(
    <ellipsoid name="ell" ax="50" by="30" cz="40" lunit="mm"/>
)", "ell");
    auto geo = RequireSolid(data, "ell");
    ASSERT_NE(geo, nullptr);

    // Interior: origin
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));
    // Interior: (0.03/0.05)^2 + (0.01/0.03)^2 + (0.01/0.04)^2
    //         = 0.36 + 0.1111 + 0.0625 = 0.5336 < 1
    EXPECT_TRUE(IsInside(geo, 0.03, 0.01, 0.01));
    // Exterior: (0.04/0.05)^2 + (0.02/0.03)^2 + (0.03/0.04)^2
    //         = 0.64 + 0.4444 + 0.5625 = 1.647 > 1
    EXPECT_FALSE(IsInside(geo, 0.04, 0.02, 0.03));
}


TEST(GDMLParser, PartialSphereContainmentThroughParser) {
    // Upper hemisphere: starttheta=0, deltatheta=90 deg
    // Only theta in [0, 90deg] is included (z >= 0 hemisphere)
    // rmax = 100mm = 0.1m
    GDMLData data = ParseGDMLWithSolids(R"(
    <sphere name="hemi" rmax="100" starttheta="0" deltatheta="90"
            startphi="0" deltaphi="360" lunit="mm" aunit="deg"/>
)", "hemi");
    auto geo = RequireSolid(data, "hemi");
    ASSERT_NE(geo, nullptr);

    // Interior: upper hemisphere (z>0, r=0.05 < 0.1)
    EXPECT_TRUE(IsInside(geo, 0, 0, 0.05));
    // Exterior: lower hemisphere (z<0, theta > 90 deg)
    EXPECT_FALSE(IsInside(geo, 0, 0, -0.05));
}


TEST(GDMLParser, PartialTorusContainmentThroughParser) {
    // Quarter torus: deltaphi=90 deg (phi from 0 to 90 deg)
    // rmin=0, rmax=20mm=0.02m, rtor=100mm=0.1m
    GDMLData data = ParseGDMLWithSolids(R"(
    <torus name="qt" rmin="0" rmax="20" rtor="100" startphi="0" deltaphi="90"
           lunit="mm" aunit="deg"/>
)", "qt");
    auto geo = RequireSolid(data, "qt");
    ASSERT_NE(geo, nullptr);

    // Interior: at phi=0, on the ring
    EXPECT_TRUE(IsInside(geo, 0.1, 0, 0));
    // Interior: at phi=45 deg (x=y=0.1/sqrt(2)~0.0707)
    EXPECT_TRUE(IsInside(geo, 0.0707, 0.0707, 0));
    // Exterior: at phi=180 deg (negative x axis, outside 90 deg cut)
    EXPECT_FALSE(IsInside(geo, -0.1, 0, 0));
    // Exterior: at phi=270 deg (negative y axis)
    EXPECT_FALSE(IsInside(geo, 0, -0.1, 0));
}


TEST(GDMLParser, PartialCylinderContainmentThroughParser) {
    // Half cylinder: deltaphi=180 deg (phi from 0 to 180 deg, y >= 0 half)
    // rmax=50mm=0.05m, z=100mm half-height=0.05m
    GDMLData data = ParseGDMLWithSolids(R"(
    <tube name="hc" rmin="0" rmax="50" z="100" startphi="0" deltaphi="180"
          lunit="mm" aunit="deg"/>
)", "hc");
    auto geo = RequireSolid(data, "hc");
    ASSERT_NE(geo, nullptr);

    // Interior: phi~45 deg (first quadrant, within 180 deg cut)
    EXPECT_TRUE(IsInside(geo, 0.03, 0.03, 0));
    // Exterior: phi~315 deg (fourth quadrant, outside 180 deg cut)
    EXPECT_FALSE(IsInside(geo, 0.03, -0.03, 0));
}


TEST(GDMLParser, TrapShearedContainment) {
    // theta=20deg tilts the trapezoid so center faces shift in x-y.
    // phi=0 means the tilt is purely in the x-direction.
    // alpha1=alpha2=10deg shears x within each face.
    // dz=50mm, dy1=dy2=30mm, dx1=dx2=dx3=dx4=20mm
    // Face center offset at z=+dz: x_off = +dz*tan(theta)*cos(phi)
    //   = 0.05*tan(20deg)*1.0 ~ 0.05*0.364 = 0.0182m
    // So at z=+0.05, the face center is near x=+0.0182, and dx=0.02m.
    // A point at x=0 should be outside the top face since 0 < 0.0182-0.02=-0.0018? No,
    // let's be more careful. With alpha, the x-extent is sheared by y*tan(alpha).
    // At y=0: x range is [xc-dx, xc+dx] = [0.0182-0.02, 0.0182+0.02] = [-0.0018, 0.0382]
    // So x=0 at z=+0.05, y=0 is barely inside (0 > -0.0018).
    // But x=-0.01 at z=+0.05, y=0 would be outside (-0.01 < -0.0018).
    //
    // At z=-0.05: x_off = -0.05*tan(20)*cos(0) ~ -0.0182m
    // x range = [-0.0182-0.02, -0.0182+0.02] = [-0.0382, -0.0018]
    // So x=0 at z=-0.05, y=0 is outside (0 > -0.0018, i.e., outside).
    // Actually wait: 0 > -0.0018 means 0 is NOT in [-0.0382, -0.0018]. Correct, outside.
    //
    // Let's use a clear test: point (0.03, 0, 0.05) should be inside the top face
    // (0.03 in [-0.0018, 0.0382]), and point (0.03, 0, -0.05) should be outside the
    // bottom face (0.03 not in [-0.0382, -0.0018]).
    // This asymmetry only exists because theta != 0 -- a zero-theta trap would
    // have the same x-range at both faces.
    GDMLData data = ParseGDMLWithSolids(R"(
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <trap name="sheared_trap" z="100" theta="20" phi="0"
          y1="30" x1="20" x2="20" alpha1="10"
          y2="30" x3="20" x4="20" alpha2="10"
          aunit="deg" lunit="mm"/>
)", "world_box");
    auto geo = RequireSolid(data, "sheared_trap");
    ASSERT_NE(geo, nullptr);

    // Center of the solid: should be inside regardless of theta
    EXPECT_TRUE(IsInside(geo, 0, 0, 0));

    // Near top face (z close to +dz=0.05m): point shifted in +x direction
    // should be inside, while unshifted point may be outside
    // At z=+0.049m (near top), face center ~ x=+0.0178m
    EXPECT_TRUE(IsInside(geo, 0.03, 0, 0.049));
    // Same x at bottom face (z=-0.049m): face center ~ x=-0.0178m
    // x=0.03 is far from [-0.0378, 0.0022] range, so outside
    EXPECT_FALSE(IsInside(geo, 0.03, 0, -0.049));
}


TEST(GDMLParser, PartialPolyhedraContainment) {
    // 6-sided polyhedra with deltaphi=180deg (half-shell)
    // rmax=50mm = 0.05m
    GDMLData data = ParseGDMLWithSolids(R"(
    <box name="world_box" x="1000" y="1000" z="1000" lunit="mm"/>
    <polyhedra name="half_hex" numsides="6" startphi="0" deltaphi="180"
               aunit="deg" lunit="mm">
      <zplane rmin="0" rmax="50" z="-40"/>
      <zplane rmin="0" rmax="50" z="40"/>
    </polyhedra>
)", "world_box");
    auto geo = RequireSolid(data, "half_hex");
    ASSERT_NE(geo, nullptr);

    // Interior: point in the positive-y half (phi ~ 90deg, within [0,180])
    EXPECT_TRUE(IsInside(geo, 0, 0.03, 0));
    // Interior: point along +x axis (phi=0, boundary of sector)
    EXPECT_TRUE(IsInside(geo, 0.03, 0, 0));
    // Exterior: point in the negative-y half (phi ~ 270deg, outside [0,180])
    EXPECT_FALSE(IsInside(geo, 0, -0.03, 0));
    // Exterior: beyond the solid radius
    EXPECT_FALSE(IsInside(geo, 0, 0.06, 0));
    // Bounding box should be tighter than full rotation
    auto bb = geo->GetBoundingBox();
    EXPECT_LT(bb.min_corner.GetY(), 0.001);
    EXPECT_GT(bb.max_corner.GetY(), 0.04);
}
