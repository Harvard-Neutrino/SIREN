// Tests for GDML parser functionality:
//   - Expression evaluator (constants, arithmetic, precedence, parens)
//   - Density unit conversion
//   - Composite-of-composite material resolution
//   - Circular volume reference detection
//   - Forward-reference boolean solid parsing

#include <cmath>
#include <set>
#include <string>
#include <fstream>
#include <cstdio>
#include <utility>

#include <gtest/gtest.h>

#include "SIREN/detector/GDMLParser.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/dataclasses/ParticleType.h"

using namespace siren::detector;

namespace {

// Helper: write a GDML string to a temp file, parse it, clean up
GDMLData ParseGDMLString(std::string const & gdml_content) {
    std::string tmpfile = "/tmp/siren_gdml_test.gdml";
    {
        std::ofstream f(tmpfile);
        f << gdml_content;
    }
    GDMLData data = ParseGDML(tmpfile);
    std::remove(tmpfile.c_str());
    return data;
}

} // anonymous namespace


// Parser behavior groups. Each include is compiled in this translation unit so
// the shared ParseGDMLString helper above remains local to the test binary.
#include "gdml/GDMLParserExpressionUnits_TEST.inc"
#include "gdml/GDMLParserMaterials_TEST.inc"
#include "gdml/GDMLParserSolids_TEST.inc"
#include "gdml/GDMLParserBooleans_TEST.inc"
#include "gdml/GDMLParserStructure_TEST.inc"
#include "gdml/GDMLParserPreprocess_TEST.inc"
#include "gdml/GDMLParserMisc_TEST.inc"
