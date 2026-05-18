// Included by GDMLParser_TEST.cxx. Do not compile directly.

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
