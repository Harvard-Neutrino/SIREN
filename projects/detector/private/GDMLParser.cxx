#include "SIREN/detector/GDMLParser.h"

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <map>
#include <cctype>
#include <functional>
#include <unordered_set>

#include <cereal/external/rapidxml/rapidxml.hpp>

namespace rapidxml = cereal::rapidxml;

#include "SIREN/utilities/Constants.h"

#include "SIREN/math/Vector3D.h"
#include "SIREN/math/Quaternion.h"
#include "SIREN/math/EulerQuaternionConversions.h"

#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/geometry/ExtrPoly.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Polyhedra.h"
#include "SIREN/geometry/Torus.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/EllipticalTube.h"
#include "SIREN/geometry/CutTube.h"
#include "SIREN/geometry/Trap.h"
#include "SIREN/geometry/Ellipsoid.h"
#include "SIREN/geometry/Para.h"
#include "SIREN/geometry/GenericPolycone.h"
#include "SIREN/geometry/GeometryMesh.h"

using namespace siren::math;
using namespace siren::geometry;

namespace siren {
namespace detector {

namespace {

// Safe attribute value accessor: returns empty string if attribute is missing
const char* SafeAttrVal(rapidxml::xml_node<>* node, const char* attr_name) {
    if(!node) return "";
    rapidxml::xml_attribute<>* attr = node->first_attribute(attr_name);
    if(!attr) return "";
    return attr->value();
}

// Trim leading and trailing whitespace from a string
std::string TrimWhitespace(std::string const & s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    if(start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

// ---- Built-in constants (CLHEP/Geant4 unit system) ----
// Single source of truth for all constants seeded before user-defined <define>.
// Values are in CLHEP base units: mm=1 (length), MeV=1 (energy), radian=1
// (angle), nanosecond=1 (time), kelvin=1 (temp), mole=1 (amount).
// ParseLength/ParseAngle handle conversion from CLHEP units to SIREN units.
//
// Unit values derived from CLHEP SystemOfUnits.h and PhysicalConstants.h:
//   https://gitlab.cern.ch/CLHEP/CLHEP/-/blob/master/Units/Units/SystemOfUnits.h
// Math constants from CLHEP Evaluator (setStdMath):
//   https://gitlab.cern.ch/CLHEP/CLHEP/-/blob/master/Evaluator/src/Evaluator.cc
//
// GDML schema and element semantics per GDML User's Guide v2.9:
//   https://gdml.web.cern.ch/doc/GDMLmanual.pdf

struct BuiltInEntry {
    const char* name;
    double value;
};

static const BuiltInEntry kBuiltInConstants[] = {
    // Math (from CLHEP Evaluator::setStdMath)
    {"pi", M_PI},
    {"twopi", 2.0 * M_PI},
    {"TWOPI", 2.0 * M_PI},
    {"halfpi", M_PI / 2.0},
    {"e", 2.7182818284590452354},
    {"gamma", 0.5772156649015328606},

    // Angle (radian = 1)
    {"rad", 1.0},
    {"radian", 1.0},
    {"mrad", 1e-3},
    {"milliradian", 1e-3},
    {"deg", M_PI / 180.0},
    {"degree", M_PI / 180.0},

    // Length (mm = 1)
    {"mm", 1.0},
    {"millimeter", 1.0},
    {"cm", 10.0},
    {"centimeter", 10.0},
    {"m", 1000.0},
    {"meter", 1000.0},
    {"km", 1e6},
    {"kilometer", 1e6},
    {"micrometer", 1e-3},
    {"nanometer", 1e-6},
    {"angstrom", 1e-7},
    {"fermi", 1e-12},
    {"inch", 25.4},
    {"in", 25.4},

    // Energy (MeV = 1)
    {"eV", 1e-6},
    {"keV", 1e-3},
    {"MeV", 1.0},
    {"GeV", 1e3},
    {"TeV", 1e6},
    {"PeV", 1e9},
    {"joule", 6.241509074460763e12},

    // Time (nanosecond = 1)
    {"ns", 1.0},
    {"nanosecond", 1.0},
    {"ps", 1e-3},
    {"picosecond", 1e-3},
    {"us", 1e3},
    {"microsecond", 1e3},
    {"ms", 1e6},
    {"millisecond", 1e6},
    {"second", 1e9},

    // Mass (derived: kg = joule * s^2 / m^2)
    {"milligram", 6.241509074460763e18},
    {"gram", 6.241509074460763e21},
    {"kg", 6.241509074460763e24},
    {"kilogram", 6.241509074460763e24},

    // Volumic mass (CLHEP mass / mm^3)
    {"g/cm3", 6.241509074460763e18},
    {"mg/cm3", 6.241509074460763e15},
    {"kg/m3", 6.241509074460763e15},

    // Amount (mole = 1)
    {"mole", 1.0},
    {"mol", 1.0},

    // Temperature (kelvin = 1)
    {"kelvin", 1.0},
};

static const std::map<std::string, double>& GetBuiltInMap() {
    static const std::map<std::string, double> m = []() {
        std::map<std::string, double> m;
        for(auto const & e : kBuiltInConstants) m[e.name] = e.value;
        return m;
    }();
    return m;
}

static const std::unordered_set<std::string>& GetBuiltInNames() {
    static const std::unordered_set<std::string> names = []() {
        std::unordered_set<std::string> s;
        for(auto const & e : kBuiltInConstants) s.insert(e.name);
        return s;
    }();
    return names;
}

static bool IsBuiltInConstant(std::string const & name) {
    return GetBuiltInNames().count(name) != 0;
}

static void SeedBuiltInConstants(GDMLData & data) {
    for(auto const & e : kBuiltInConstants) {
        data.constants[e.name] = e.value;
    }
}

// Iterative expression evaluator using the shunting-yard algorithm.
// Tokenizes the input, converts infix to postfix (RPN), and evaluates
// the RPN with an explicit value stack. No recursion, O(1) call-stack depth.

enum class TokenKind { NUMBER, OP, LPAREN, RPAREN, COMMA, FUNC, UNARY_NEG };

struct Token {
    TokenKind kind;
    double value;
    char op;
    std::string name;
    int argc;
};

static int OpPrecedence(char op) {
    if(op == '+' || op == '-') return 1;
    if(op == '*' || op == '/') return 2;
    if(op == '^') return 3;
    return 0;
}

static bool IsRightAssociative(char op) {
    return op == '^';
}

static bool IsIdentStart(char c) {
    return std::isalpha(static_cast<unsigned char>(c)) || c == '_';
}

static bool IsIdentChar(char c) {
    return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}

static double ApplyFunc(std::string const & name, double const * args, int argc, std::string const & expr) {
    if(argc == 1) {
        double a = args[0];
        if(name == "sin") return std::sin(a);
        if(name == "cos") return std::cos(a);
        if(name == "tan") return std::tan(a);
        if(name == "asin") return std::asin(a);
        if(name == "acos") return std::acos(a);
        if(name == "atan") return std::atan(a);
        if(name == "exp") return std::exp(a);
        if(name == "log") return std::log(a);
        if(name == "sqrt") return std::sqrt(a);
        if(name == "abs") return std::fabs(a);
    } else if(argc == 2) {
        double a = args[0], b = args[1];
        if(name == "pow") return std::pow(a, b);
        if(name == "atan2") return std::atan2(a, b);
        if(name == "min") return std::min(a, b);
        if(name == "max") return std::max(a, b);
    }
    throw std::runtime_error("GDML expression error: unknown function '" + name + "' in '" + expr + "'");
}

static double EvalExpression(std::string const & expr, std::map<std::string, double> const & constants,
                             std::map<std::string, std::pair<int, std::vector<double>>> const & matrices = {});

static std::vector<Token> Tokenize(std::string const & s, std::map<std::string, double> const & constants, std::string const & expr,
                                   std::map<std::string, std::pair<int, std::vector<double>>> const & matrices = {}) {
    std::vector<Token> tokens;
    size_t i = 0;
    while(i < s.size()) {
        char c = s[i];
        if(c == ' ' || c == '\t' || c == '\n' || c == '\r') { ++i; continue; }

        if(c == '(') { tokens.push_back({TokenKind::LPAREN}); ++i; continue; }
        if(c == ')') { tokens.push_back({TokenKind::RPAREN}); ++i; continue; }
        if(c == ',') { tokens.push_back({TokenKind::COMMA}); ++i; continue; }

        if(c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
            bool is_unary = (c == '-' || c == '+') &&
                (tokens.empty() || tokens.back().kind == TokenKind::LPAREN ||
                 tokens.back().kind == TokenKind::OP || tokens.back().kind == TokenKind::COMMA ||
                 tokens.back().kind == TokenKind::UNARY_NEG);
            if(is_unary && c == '-') {
                tokens.push_back({TokenKind::UNARY_NEG});
            } else if(is_unary && c == '+') {
                // unary plus: skip
            } else {
                Token t; t.kind = TokenKind::OP; t.op = c;
                tokens.push_back(t);
            }
            ++i; continue;
        }

        if(std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            size_t pos = 0;
            double val = std::stod(s.substr(i), &pos);
            Token t; t.kind = TokenKind::NUMBER; t.value = val;
            tokens.push_back(t);
            i += pos; continue;
        }

        if(IsIdentStart(c)) {
            size_t start = i;
            while(i < s.size() && IsIdentChar(s[i])) ++i;
            std::string name = s.substr(start, i - start);
            size_t j = i;
            while(j < s.size() && (s[j] == ' ' || s[j] == '\t')) ++j;
            if(j < s.size() && s[j] == '(') {
                Token t; t.kind = TokenKind::FUNC; t.name = name; t.argc = 1;
                tokens.push_back(t);
            } else if(j < s.size() && s[j] == '[') {
                // Matrix/vector element access: name[idx] or name[row,col]
                auto mit = matrices.find(name);
                if(mit == matrices.end()) {
                    throw std::runtime_error("GDML expression error: unknown matrix '" + name + "' in '" + expr + "'");
                }
                i = j + 1; // skip '['
                // Parse first index expression (up to ',' or ']')
                int depth = 1;
                size_t idx_start = i;
                size_t comma_pos = std::string::npos;
                while(i < s.size() && depth > 0) {
                    if(s[i] == '[') ++depth;
                    else if(s[i] == ']') --depth;
                    else if(s[i] == ',' && depth == 1) comma_pos = i;
                    if(depth > 0) ++i;
                }
                if(depth != 0) {
                    throw std::runtime_error("GDML expression error: unmatched '[' in '" + expr + "'");
                }
                size_t bracket_end = i;
                ++i; // skip ']'

                int coldim = mit->second.first;
                auto const & vals = mit->second.second;
                double result;

                if(comma_pos != std::string::npos) {
                    // Two indices: name[row, col] (1-based)
                    std::string row_expr = s.substr(idx_start, comma_pos - idx_start);
                    std::string col_expr = s.substr(comma_pos + 1, bracket_end - comma_pos - 1);
                    // Recursively evaluate the index expressions
                    double row_val = EvalExpression(row_expr, constants);
                    double col_val = EvalExpression(col_expr, constants);
                    int row = static_cast<int>(row_val + 0.5);
                    int col = static_cast<int>(col_val + 0.5);
                    int flat_idx = (row - 1) * coldim + (col - 1);
                    if(flat_idx < 0 || flat_idx >= static_cast<int>(vals.size())) {
                        throw std::runtime_error("GDML expression error: matrix index out of range for '" + name + "[" + std::to_string(row) + "," + std::to_string(col) + "]' in '" + expr + "'");
                    }
                    result = vals[flat_idx];
                } else {
                    // Single index: name[idx] (1-based)
                    std::string idx_expr = s.substr(idx_start, bracket_end - idx_start);
                    double idx_val = EvalExpression(idx_expr, constants);
                    int idx = static_cast<int>(idx_val + 0.5);
                    int flat_idx;
                    if(coldim == 1) {
                        flat_idx = idx - 1;
                    } else {
                        flat_idx = idx - 1;
                    }
                    if(flat_idx < 0 || flat_idx >= static_cast<int>(vals.size())) {
                        throw std::runtime_error("GDML expression error: matrix index out of range for '" + name + "[" + std::to_string(idx) + "]' in '" + expr + "'");
                    }
                    result = vals[flat_idx];
                }
                Token t; t.kind = TokenKind::NUMBER; t.value = result;
                tokens.push_back(t);
            } else {
                auto it = constants.find(name);
                if(it == constants.end()) {
                    throw std::runtime_error("GDML expression error: unknown constant '" + name + "' in '" + expr + "'");
                }
                Token t; t.kind = TokenKind::NUMBER; t.value = it->second;
                tokens.push_back(t);
            }
            continue;
        }

        throw std::runtime_error("GDML expression error: unexpected character '" + std::string(1, c) + "' in '" + expr + "'");
    }
    return tokens;
}

double EvalExpression(std::string const & expr, std::map<std::string, double> const & constants,
                      std::map<std::string, std::pair<int, std::vector<double>>> const & matrices) {
    std::string s = TrimWhitespace(expr);
    if(s.empty()) return 0.0;

    // Fast path: bare number
    {
        size_t pos = 0;
        try {
            double val = std::stod(s, &pos);
            if(pos == s.size()) return val;
        } catch(std::invalid_argument &) {
        } catch(std::out_of_range &) {
            throw std::runtime_error("GDML expression error: numeric value out of range in '" + expr + "'");
        }
    }

    // Fast path: bare constant name (only if no bracket follows)
    {
        bool bare = true;
        for(size_t i = 0; i < s.size(); ++i) {
            if(!IsIdentChar(s[i])) { bare = false; break; }
        }
        if(bare && IsIdentStart(s[0])) {
            auto it = constants.find(s);
            if(it != constants.end()) return it->second;
            if(matrices.find(s) == matrices.end()) {
                throw std::runtime_error("GDML expression error: unknown constant '" + s + "'");
            }
        }
    }

    std::vector<Token> tokens = Tokenize(s, constants, expr, matrices);

    // Shunting-yard: convert to RPN
    std::vector<Token> output;
    std::vector<Token> op_stack;
    output.reserve(tokens.size());
    op_stack.reserve(tokens.size());

    for(size_t i = 0; i < tokens.size(); ++i) {
        Token const & tok = tokens[i];
        switch(tok.kind) {
        case TokenKind::NUMBER:
            output.push_back(tok);
            break;
        case TokenKind::FUNC:
            op_stack.push_back(tok);
            break;
        case TokenKind::COMMA:
            while(!op_stack.empty() && op_stack.back().kind != TokenKind::LPAREN) {
                output.push_back(op_stack.back());
                op_stack.pop_back();
            }
            if(!op_stack.empty() && op_stack.size() >= 2) {
                Token & func_below = op_stack[op_stack.size() - 2];
                if(func_below.kind == TokenKind::FUNC) func_below.argc++;
            }
            break;
        case TokenKind::OP:
            while(!op_stack.empty() && op_stack.back().kind == TokenKind::OP) {
                int stack_prec = OpPrecedence(op_stack.back().op);
                int cur_prec = OpPrecedence(tok.op);
                // Right-associative: pop only when stack precedence is strictly greater
                // Left-associative: pop when stack precedence is greater or equal
                if(IsRightAssociative(tok.op) ? (stack_prec > cur_prec) : (stack_prec >= cur_prec)) {
                    output.push_back(op_stack.back());
                    op_stack.pop_back();
                } else {
                    break;
                }
            }
            op_stack.push_back(tok);
            break;
        case TokenKind::UNARY_NEG:
            op_stack.push_back(tok);
            break;
        case TokenKind::LPAREN:
            op_stack.push_back(tok);
            break;
        case TokenKind::RPAREN:
            while(!op_stack.empty() && op_stack.back().kind != TokenKind::LPAREN) {
                output.push_back(op_stack.back());
                op_stack.pop_back();
            }
            if(op_stack.empty()) {
                throw std::runtime_error("GDML expression error: mismatched parentheses in '" + expr + "'");
            }
            op_stack.pop_back(); // pop LPAREN
            if(!op_stack.empty() && op_stack.back().kind == TokenKind::FUNC) {
                output.push_back(op_stack.back());
                op_stack.pop_back();
            }
            break;
        }
    }
    while(!op_stack.empty()) {
        if(op_stack.back().kind == TokenKind::LPAREN) {
            throw std::runtime_error("GDML expression error: mismatched parentheses in '" + expr + "'");
        }
        output.push_back(op_stack.back());
        op_stack.pop_back();
    }

    // Evaluate RPN
    std::vector<double> val_stack;
    val_stack.reserve(output.size());
    for(auto const & tok : output) {
        switch(tok.kind) {
        case TokenKind::NUMBER:
            val_stack.push_back(tok.value);
            break;
        case TokenKind::OP: {
            if(val_stack.size() < 2) {
                throw std::runtime_error("GDML expression error: malformed expression '" + expr + "'");
            }
            double b = val_stack.back(); val_stack.pop_back();
            double a = val_stack.back(); val_stack.pop_back();
            switch(tok.op) {
            case '+': val_stack.push_back(a + b); break;
            case '-': val_stack.push_back(a - b); break;
            case '*': val_stack.push_back(a * b); break;
            case '/':
                if(b == 0.0) throw std::runtime_error("GDML expression error: division by zero in '" + expr + "'");
                val_stack.push_back(a / b);
                break;
            case '^':
                val_stack.push_back(std::pow(a, b));
                break;
            }
            break;
        }
        case TokenKind::UNARY_NEG: {
            if(val_stack.empty()) {
                throw std::runtime_error("GDML expression error: malformed expression '" + expr + "'");
            }
            val_stack.back() = -val_stack.back();
            break;
        }
        case TokenKind::FUNC: {
            int argc = tok.argc;
            if(argc > 2) {
                throw std::runtime_error("GDML expression error: function '" + tok.name + "' called with " + std::to_string(argc) + " arguments (max 2) in '" + expr + "'");
            }
            if(static_cast<int>(val_stack.size()) < argc) {
                throw std::runtime_error("GDML expression error: not enough arguments for '" + tok.name + "' in '" + expr + "'");
            }
            double args[2];
            for(int j = argc - 1; j >= 0; --j) {
                args[j] = val_stack.back();
                val_stack.pop_back();
            }
            val_stack.push_back(ApplyFunc(tok.name, args, argc, expr));
            break;
        }
        default:
            throw std::runtime_error("GDML expression error: unexpected token in RPN for '" + expr + "'");
        }
    }
    if(val_stack.size() != 1) {
        throw std::runtime_error("GDML expression error: malformed expression '" + expr + "'");
    }
    return val_stack[0];
}

// Parse a double from a string, returning 0 if empty.
// Optionally evaluates expressions using the provided constants map.
double SafeParseDouble(const char* val, std::map<std::string, double> const & constants = {},
                       std::map<std::string, std::pair<int, std::vector<double>>> const & matrices = {}) {
    if(!val || val[0] == '\0') return 0.0;
    // Fast path: try simple numeric literal
    try {
        size_t pos = 0;
        double result = std::stod(std::string(val), &pos);
        std::string s(val);
        if(pos == s.size()) return result;
    } catch(std::invalid_argument &) {
        // Not a simple number
    } catch(std::out_of_range &) {
        throw std::runtime_error("GDML expression error: numeric value out of range in '" + std::string(val) + "'");
    }
    // Fall back to expression evaluation
    return EvalExpression(std::string(val), constants, matrices);
}

// Get the length scale factor for a given unit string
double LengthScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        return siren::utilities::Constants::mm; // GDML default is mm
    }
    std::string u(unit);
    if(u == "mm") return siren::utilities::Constants::mm;
    if(u == "cm") return siren::utilities::Constants::cm;
    if(u == "m")  return siren::utilities::Constants::m;
    if(u == "um") return siren::utilities::Constants::um;
    if(u == "nm") return siren::utilities::Constants::nm;
    if(u == "km") return siren::utilities::Constants::km;
    throw std::runtime_error("GDML error: unrecognized length unit '" + u + "'");
}

// Convert a GDML length value+unit to meters (SIREN base unit)
// GDML default length unit is mm
double ParseLength(const char* value, const char* unit,
                   std::map<std::string, double> const & constants = {}) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);

    return val * LengthScale(unit);
}

// Get the angle scale factor for a given unit string
double AngleScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        // GDML default angle unit is radians
        return siren::utilities::Constants::radian;
    }
    std::string u(unit);
    if(u == "deg" || u == "degree") return siren::utilities::Constants::degrees;
    if(u == "rad" || u == "radian") return siren::utilities::Constants::radian;
    if(u == "mrad") return siren::utilities::Constants::mrad;
    throw std::runtime_error("GDML error: unrecognized angle unit '" + u + "'");
}

// Convert a GDML angle value+unit to radians
// GDML default angle unit is degrees
double ParseAngle(const char* value, const char* unit,
                  std::map<std::string, double> const & constants = {}) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);
    return val * AngleScale(unit);
}

// Build a quaternion from GDML rotation convention:
// Extrinsic rotations: rotate around X by rx, then Y by ry, then Z by rz
// Frame is static, so XYZs is the appropriate convention
Quaternion QuatFromGDMLRotation(double rx, double ry, double rz) {
    return siren::math::QFromXYZs(rx, ry, rz);
}

// Convert a value from the given unit to the GDML default for that unit type.
// For length: convert to mm (GDML default length)
// For angle: convert to rad (GDML default angle)
// For energy: convert to MeV (common default)
// Returns the scale factor to multiply the raw value by.
// Forward declaration (defined after ExpandEntities)
static void EmitWarning(GDMLData & data, GDMLParseOptions const & options, std::string const & msg);

// Infer quantity type from a unit string.
static GDMLQuantityType TypeFromUnit(std::string const & u) {
    if(u.empty()) return GDMLQuantityType::NONE;
    if(u == "mm" || u == "cm" || u == "m" || u == "km" ||
       u == "um" || u == "nm" || u == "in" ||
       u == "millimeter" || u == "centimeter" || u == "meter" ||
       u == "kilometer" || u == "micrometer" || u == "nanometer" ||
       u == "angstrom" || u == "fermi" || u == "inch")
        return GDMLQuantityType::LENGTH;
    if(u == "rad" || u == "deg" || u == "mrad" ||
       u == "radian" || u == "degree" || u == "milliradian")
        return GDMLQuantityType::ANGLE;
    if(u == "eV" || u == "keV" || u == "MeV" || u == "GeV" ||
       u == "TeV" || u == "PeV" || u == "joule")
        return GDMLQuantityType::ENERGY;
    if(u == "ns" || u == "ps" || u == "us" || u == "ms" ||
       u == "nanosecond" || u == "picosecond" || u == "microsecond" ||
       u == "millisecond" || u == "second")
        return GDMLQuantityType::TIME;
    if(u == "g/cm3" || u == "mg/cm3" || u == "kg/m3")
        return GDMLQuantityType::DENSITY;
    if(u == "kelvin")
        return GDMLQuantityType::TEMPERATURE;
    if(u == "mole" || u == "mol")
        return GDMLQuantityType::AMOUNT;
    return GDMLQuantityType::NONE;
}

// Infer quantity type from the type attribute string.
static GDMLQuantityType TypeFromTypeAttr(std::string const & type_str) {
    if(type_str.empty()) return GDMLQuantityType::NONE;
    std::string t = type_str;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    if(t == "length" || t == "lenght" || t == "coordinate")
        return GDMLQuantityType::LENGTH;
    if(t == "angle") return GDMLQuantityType::ANGLE;
    if(t == "energy" || t == "threshold") return GDMLQuantityType::ENERGY;
    if(t == "time") return GDMLQuantityType::TIME;
    if(t == "density") return GDMLQuantityType::DENSITY;
    if(t == "temperature") return GDMLQuantityType::TEMPERATURE;
    return GDMLQuantityType::NONE;
}

// Resolve the quantity type from both the type attribute and the unit string.
// If both are present and disagree, warn and prefer the unit-inferred type.
static GDMLQuantityType ResolveQuantityType(
        const char* type_attr, const char* unit,
        std::string const & name,
        GDMLData & data, GDMLParseOptions const & options) {
    GDMLQuantityType from_attr = TypeFromTypeAttr(type_attr ? type_attr : "");
    GDMLQuantityType from_unit = TypeFromUnit(unit ? unit : "");

    if(from_attr != GDMLQuantityType::NONE && from_unit != GDMLQuantityType::NONE) {
        if(from_attr != from_unit) {
            EmitWarning(data, options, "quantity '" + name + "': type '"
                + std::string(type_attr) + "' conflicts with unit '"
                + std::string(unit) + "'; using unit-inferred type");
        }
        return from_unit;
    }
    if(from_unit != GDMLQuantityType::NONE) return from_unit;
    return from_attr;
}

// Convert a quantity value from its declared unit to the canonical storage
// representation for its type.
//
// LENGTH/ANGLE/ENERGY/TIME: stored in CLHEP base units (mm, rad, MeV, ns).
//   These equal the GDML default units, so downstream consumption via
//   lscale/ascale works correctly.
// DENSITY: stored in g/cm3. The <D> parser expects g/cm3, so we must NOT
//   use CLHEP density units here.
static double QuantityUnitScale(const char* unit, GDMLQuantityType resolved_type) {
    if(!unit || unit[0] == '\0') return 1.0;

    std::string u(unit);

    if(resolved_type == GDMLQuantityType::DENSITY) {
        if(u == "g/cm3") return 1.0;
        if(u == "mg/cm3") return 1.0e-3;
        if(u == "kg/m3") return 1.0e-3;
        throw std::runtime_error("GDML error: unrecognized density unit '" + u + "'");
    }

    auto it = GetBuiltInMap().find(u);
    if(it != GetBuiltInMap().end()) return it->second;
    throw std::runtime_error("GDML error: unrecognized unit '" + u + "'");
}

// Parse SYSTEM entity declarations from a DOCTYPE internal subset.
// Returns a list of (entity_name, resolved_file_path) pairs.
static std::vector<std::pair<std::string, std::string>> ParseEntityDeclarations(
        std::string const & content, std::string const & base_dir) {
    std::vector<std::pair<std::string, std::string>> decls;

    size_t doctype_start = content.find("<!DOCTYPE");
    if(doctype_start == std::string::npos)
        doctype_start = content.find("<!doctype");
    if(doctype_start == std::string::npos) return decls;

    size_t bracket_open = content.find('[', doctype_start);
    if(bracket_open == std::string::npos) return decls;
    size_t bracket_close = content.find(']', bracket_open);
    if(bracket_close == std::string::npos) return decls;

    std::string dtd = content.substr(bracket_open + 1, bracket_close - bracket_open - 1);
    size_t pos = 0;
    while((pos = dtd.find("<!ENTITY", pos)) != std::string::npos) {
        pos += 8;
        while(pos < dtd.size() && std::isspace(static_cast<unsigned char>(dtd[pos]))) ++pos;

        size_t name_start = pos;
        while(pos < dtd.size() && !std::isspace(static_cast<unsigned char>(dtd[pos]))) ++pos;
        std::string ent_name = dtd.substr(name_start, pos - name_start);

        while(pos < dtd.size() && std::isspace(static_cast<unsigned char>(dtd[pos]))) ++pos;

        if(dtd.substr(pos, 6) != "SYSTEM") continue;
        pos += 6;
        while(pos < dtd.size() && std::isspace(static_cast<unsigned char>(dtd[pos]))) ++pos;

        if(pos >= dtd.size()) break;
        char quote = dtd[pos];
        if(quote != '"' && quote != '\'') continue;
        ++pos;
        size_t path_start = pos;
        while(pos < dtd.size() && dtd[pos] != quote) ++pos;
        std::string filepath = dtd.substr(path_start, pos - path_start);

        std::string full_path;
        if(!filepath.empty() && filepath[0] == '/') {
            full_path = filepath;
        } else {
            full_path = base_dir + "/" + filepath;
        }
        decls.push_back({ent_name, full_path});
    }
    return decls;
}

// Strip the DOCTYPE declaration from content (RapidXML cannot parse it).
static std::string StripDoctype(std::string const & content) {
    size_t doctype_start = content.find("<!DOCTYPE");
    if(doctype_start == std::string::npos)
        doctype_start = content.find("<!doctype");
    if(doctype_start == std::string::npos) return content;

    size_t bracket_close = content.find(']', doctype_start);
    if(bracket_close == std::string::npos) return content;
    size_t doctype_end = content.find('>', bracket_close);
    if(doctype_end == std::string::npos) return content;

    std::string result = content;
    result.erase(doctype_start, doctype_end - doctype_start + 1);
    return result;
}

// Substitute entity references (&name;) in content using the given entity map.
static std::string SubstituteEntities(std::string const & content,
        std::map<std::string, std::string> const & entities) {
    std::string result = content;
    for(auto const & ent : entities) {
        std::string ref = "&" + ent.first + ";";
        size_t pos = 0;
        while((pos = result.find(ref, pos)) != std::string::npos) {
            result.replace(pos, ref.size(), ent.second);
            pos += ent.second.size();
        }
    }
    return result;
}

// Preprocess XML content to expand ENTITY references.
// Reads the DOCTYPE for ENTITY SYSTEM declarations, loads the referenced files,
// and substitutes all &name; references with file contents.
// Uses an iterative stack to handle nested entity files without recursion.
// Each stack frame owns its own scope: entities declared in a file's DOCTYPE
// are only visible for substitution within that file's content.
static std::string ExpandEntities(std::string const & content, std::string const & base_dir) {
    struct Frame {
        std::string content;
        std::string base_dir;
        std::vector<std::pair<std::string, std::string>> decls; // name, file_path
        std::map<std::string, std::string> resolved;            // name -> expanded content
        int next_decl = 0;
    };

    std::vector<Frame> stack;
    stack.push_back({content, base_dir, ParseEntityDeclarations(content, base_dir), {}, 0});

    static constexpr int MAX_DEPTH = 16;

    while(!stack.empty()) {
        Frame & top = stack.back();

        if(top.next_decl < (int)top.decls.size() && (int)stack.size() <= MAX_DEPTH) {
            auto & [ent_name, ent_path] = top.decls[top.next_decl];
            top.next_decl++;

            std::ifstream ent_file(ent_path.c_str(), std::ios::binary);
            if(!ent_file.is_open()) continue;

            std::string file_content((std::istreambuf_iterator<char>(ent_file)),
                                      std::istreambuf_iterator<char>());
            ent_file.close();

            size_t last_slash = ent_path.rfind('/');
            std::string ent_dir = (last_slash != std::string::npos) ?
                ent_path.substr(0, last_slash) : top.base_dir;

            auto child_decls = ParseEntityDeclarations(file_content, ent_dir);
            if(child_decls.empty()) {
                top.resolved[ent_name] = file_content;
            } else {
                stack.push_back({file_content, ent_dir, std::move(child_decls), {}, 0});
                // Tag the child frame so we know which entity it resolves
                // Store the entity name in the parent's next slot (already incremented)
                // We'll retrieve it when the child pops
            }
        } else {
            // All declarations resolved (or depth limit reached) — substitute and pop
            std::string expanded = StripDoctype(top.content);
            if(!top.resolved.empty()) {
                expanded = SubstituteEntities(expanded, top.resolved);
            }

            if(stack.size() == 1) {
                return expanded;
            }

            // Pop this frame and store result in parent as the resolved entity
            std::string result = std::move(expanded);
            stack.pop_back();

            Frame & parent = stack.back();
            // The entity name is from the declaration we just finished processing
            // (parent.next_decl - 1, since we incremented before pushing)
            std::string const & ent_name = parent.decls[parent.next_decl - 1].first;
            parent.resolved[ent_name] = std::move(result);
        }
    }

    return content;
}

// Record a warning. In strict mode, throws instead of continuing.
void EmitWarning(GDMLData & data, GDMLParseOptions const & options, std::string const & msg) {
    if(options.strict) {
        throw std::runtime_error("GDML error (strict mode): " + msg);
    }
    data.warnings.push_back(msg);
    std::cerr << "GDML warning: " << msg << std::endl;
}

// Substitute loop variable references in an attribute value string.
// Two substitution modes following Geant4 convention:
//   1. Bracket notation: [expr] is evaluated and replaced with the integer
//      result. Used in name attributes to construct unique identifiers
//      (e.g. "layer_[i]" -> "layer_0"). The expression inside brackets
//      can be any evaluable expression using the loop variable.
//   2. Word-boundary replacement: bare occurrences of the variable name
//      that form a complete identifier token are replaced with the numeric
//      value. Used in numeric expressions (e.g. "i*10" -> "0*10").
//      Does NOT replace the variable inside longer identifiers
//      (e.g. "ChildVol" is untouched when the variable is "i").
static std::string SubstLoopVar(
    std::string const & val,
    std::string const & var_name,
    std::string const & var_value_str,
    std::map<std::string, double> const & constants) {
    std::string result = val;

    // Pass 1: bracket substitution [expr] -> evaluated integer
    size_t pos = 0;
    while((pos = result.find('[', pos)) != std::string::npos) {
        size_t end = result.find(']', pos + 1);
        if(end == std::string::npos) break;
        std::string expr = result.substr(pos + 1, end - pos - 1);
        try {
            double v = EvalExpression(expr, constants);
            long long iv = static_cast<long long>(v);
            result.replace(pos, end - pos + 1, std::to_string(iv));
        } catch(...) {
            pos = end + 1;
        }
    }

    // Pass 2: word-boundary replacement of the bare variable name
    pos = 0;
    while((pos = result.find(var_name, pos)) != std::string::npos) {
        bool left_ok = (pos == 0 || !IsIdentChar(result[pos - 1]));
        bool right_ok = (pos + var_name.size() >= result.size() || !IsIdentChar(result[pos + var_name.size()]));
        if(left_ok && right_ok) {
            result.replace(pos, var_name.size(), var_value_str);
            pos += var_value_str.size();
        } else {
            pos += var_name.size();
        }
    }

    return result;
}

// Deep-clone an XML node and all its children/attributes into a document.
// Performs loop variable substitution in all attribute values.
static rapidxml::xml_node<>* CloneNodeWithSubst(
    rapidxml::xml_document<> & doc,
    rapidxml::xml_node<>* src,
    std::string const & var_name,
    std::string const & var_value_str,
    std::map<std::string, double> const & constants) {

    auto* node = doc.allocate_node(src->type());
    node->name(src->name(), src->name_size());

    for(auto* attr = src->first_attribute(); attr; attr = attr->next_attribute()) {
        std::string val(attr->value(), attr->value_size());
        val = SubstLoopVar(val, var_name, var_value_str, constants);
        char* aname = doc.allocate_string(attr->name(), attr->name_size() + 1);
        char* aval = doc.allocate_string(val.c_str(), val.size() + 1);
        node->append_attribute(doc.allocate_attribute(aname, aval));
    }

    for(auto* child = src->first_node(); child; child = child->next_sibling()) {
        node->append_node(CloneNodeWithSubst(doc, child, var_name, var_value_str, constants));
    }
    return node;
}

// Expand all <loop> elements in a DOM subtree (iterative, stack-based).
// Clones loop body nodes for each iteration with text substitution, inserts
// them into the parent, and removes the loop node. Used for non-define
// sections (structure, solids) where the DOM must be modified for downstream
// parsers that iterate over it.
static void ExpandLoops(rapidxml::xml_document<> & doc,
                        rapidxml::xml_node<>* root,
                        std::map<std::string, double> & constants) {
    if(!root) return;

    struct Frame {
        rapidxml::xml_node<>* parent;
        rapidxml::xml_node<>* current; // next child to process
    };

    std::vector<Frame> stack;
    stack.push_back({root, root->first_node()});

    static constexpr int MAX_DEPTH = 64;

    while(!stack.empty()) {
        if((int)stack.size() > MAX_DEPTH) {
            stack.pop_back();
            continue;
        }

        Frame & top = stack.back();
        if(!top.current) {
            stack.pop_back();
            continue;
        }

        auto* child = top.current;
        top.current = child->next_sibling();

        if(std::string(child->name()) != "loop") {
            // Descend into non-loop nodes to find nested loops
            if(child->first_node()) {
                stack.push_back({child, child->first_node()});
            }
            continue;
        }

        // Process the loop: expand inline
        std::string var_name = SafeAttrVal(child, "for");
        if(IsBuiltInConstant(var_name)) {
            top.parent->remove_node(child);
            continue;
        }
        const char* from_attr = SafeAttrVal(child, "from");
        double from_val;
        if(from_attr[0] != '\0') {
            from_val = SafeParseDouble(from_attr, constants);
        } else {
            auto it = constants.find(var_name);
            from_val = (it != constants.end()) ? it->second : 0.0;
        }
        double to_val = SafeParseDouble(SafeAttrVal(child, "to"), constants);
        double step_val = SafeParseDouble(SafeAttrVal(child, "step"), constants);
        if(step_val == 0) step_val = 1.0;

        // Track the first cloned node so we can resume processing from there
        // (cloned nodes may contain nested loops)
        rapidxml::xml_node<>* first_cloned = nullptr;
        int count = 0;
        for(double v = from_val;
            (step_val > 0) ? (v <= to_val + 1e-9) : (v >= to_val - 1e-9);
            v += step_val) {
            if(++count > 10000) break;
            constants[var_name] = v;

            std::string v_str;
            if(v == std::floor(v) && std::fabs(v) < 1e15) {
                v_str = std::to_string(static_cast<long long>(v));
            } else {
                std::ostringstream oss;
                oss << v;
                v_str = oss.str();
            }

            for(auto* body = child->first_node(); body; body = body->next_sibling()) {
                auto* cloned = CloneNodeWithSubst(doc, body, var_name, v_str, constants);
                top.parent->insert_node(child, cloned);
                if(!first_cloned) first_cloned = cloned;
            }
        }

        top.parent->remove_node(child);

        // Resume processing from the first cloned node (to catch nested loops)
        if(first_cloned) {
            top.current = first_cloned;
        }
    }
}

// Resolve all <define> sections in document order with inline loop handling.
// Uses an explicit stack to iterate loop bodies without recursion or DOM
// modification. Constants/variables are evaluated immediately when encountered,
// giving document-order semantics (a constant before a loop sees the pre-loop
// variable value; a constant after sees the post-loop value).
static void ResolveDefineInOrder(rapidxml::xml_node<>* gdml_node,
                                 GDMLData & data,
                                 GDMLParseOptions const & options) {
    SeedBuiltInConstants(data);

    // Evaluate bracket expressions [expr] in an attribute string using current
    // constants. Transforms "sum_[i+1]" into "sum_4" (if i=3 in constants).
    // Used for both name and value attributes within loop bodies.
    auto resolveBrackets = [&](const char* raw) -> std::string {
        std::string s(raw);
        size_t pos = 0;
        while((pos = s.find('[', pos)) != std::string::npos) {
            size_t end = s.find(']', pos + 1);
            if(end == std::string::npos) break;
            std::string expr = s.substr(pos + 1, end - pos - 1);
            try {
                double v = EvalExpression(expr, data.constants, data.matrices);
                long long iv = static_cast<long long>(v);
                s.replace(pos, end - pos + 1, std::to_string(iv));
            } catch(...) {
                pos = end + 1;
            }
        }
        return s;
    };

    // Stack frame: either a sequence of sibling nodes, or a loop iteration context.
    struct Frame {
        rapidxml::xml_node<>* current; // next node to process
        bool is_loop = false;
        std::string var_name;
        double to_val = 0;
        double step_val = 0;
        rapidxml::xml_node<>* body_first = nullptr; // first child of loop (for restart)
    };

    std::vector<Frame> stack;
    // Push all define sections (in order) as top-level frames
    for(auto* def = gdml_node->first_node("define"); def; def = def->next_sibling("define")) {
        stack.push_back({def->first_node(), false, {}, 0, 0, nullptr});
    }

    static constexpr int MAX_DEPTH = 64;

    while(!stack.empty()) {
        Frame & top = stack.back();

        if(!top.current) {
            // End of this scope
            if(top.is_loop) {
                // Advance loop iteration
                double next_val = data.constants[top.var_name] + top.step_val;
                bool in_range = (top.step_val > 0)
                    ? (next_val <= top.to_val + 1e-9)
                    : (next_val >= top.to_val - 1e-9);
                if(in_range) {
                    data.constants[top.var_name] = next_val;
                    top.current = top.body_first;
                    continue;
                }
            }
            stack.pop_back();
            continue;
        }

        auto* node = top.current;
        top.current = node->next_sibling();
        std::string tag(node->name());

        if(tag == "loop") {
            if((int)stack.size() >= MAX_DEPTH) continue;
            std::string var_name = SafeAttrVal(node, "for");
            if(IsBuiltInConstant(var_name)) {
                EmitWarning(data, options, "loop variable '" + var_name + "' shadows a built-in constant, skipping loop");
                continue;
            }
            const char* from_attr = SafeAttrVal(node, "from");
            double from_val;
            if(from_attr[0] != '\0') {
                from_val = SafeParseDouble(from_attr, data.constants, data.matrices);
            } else {
                auto it = data.constants.find(var_name);
                from_val = (it != data.constants.end()) ? it->second : 0.0;
            }
            double to_val = SafeParseDouble(SafeAttrVal(node, "to"), data.constants, data.matrices);
            double step_val = SafeParseDouble(SafeAttrVal(node, "step"), data.constants, data.matrices);
            if(step_val == 0) step_val = 1.0;

            // Check if the loop will execute at all
            bool will_run = (step_val > 0) ? (from_val <= to_val + 1e-9) : (from_val >= to_val - 1e-9);
            if(!will_run || !node->first_node()) continue;

            data.constants[var_name] = from_val;
            stack.push_back({node->first_node(), true, var_name, to_val, step_val, node->first_node()});
        }
        else if(tag == "constant" || tag == "variable") {
            std::string name = resolveBrackets(SafeAttrVal(node, "name"));
            if(name.empty()) continue;
            if(IsBuiltInConstant(name)) {
                EmitWarning(data, options, "cannot redefine built-in constant '" + name + "', ignoring");
                continue;
            }
            bool in_loop = std::any_of(stack.begin(), stack.end(),
                [](Frame const & f) { return f.is_loop; });
            double value;
            if(in_loop) {
                std::string val_str = resolveBrackets(SafeAttrVal(node, "value"));
                value = SafeParseDouble(val_str.c_str(), data.constants, data.matrices);
            } else {
                value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants, data.matrices);
            }
            data.constants[name] = value;
        }
        else if(tag == "quantity") {
            std::string name = resolveBrackets(SafeAttrVal(node, "name"));
            if(name.empty()) continue;
            if(IsBuiltInConstant(name)) {
                EmitWarning(data, options, "cannot redefine built-in constant '" + name + "', ignoring");
                continue;
            }
            bool in_loop = std::any_of(stack.begin(), stack.end(),
                [](Frame const & f) { return f.is_loop; });
            double value;
            if(in_loop) {
                std::string val_str = resolveBrackets(SafeAttrVal(node, "value"));
                value = SafeParseDouble(val_str.c_str(), data.constants, data.matrices);
            } else {
                value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants, data.matrices);
            }
            const char* unit = SafeAttrVal(node, "unit");
            const char* type_attr = SafeAttrVal(node, "type");
            GDMLQuantityType resolved_type = ResolveQuantityType(type_attr, unit, name, data, options);
            value *= QuantityUnitScale(unit, resolved_type);
            data.constants[name] = value;
            if(resolved_type != GDMLQuantityType::NONE) {
                data.quantity_types[name] = resolved_type;
            }
        }
        else if(tag == "position") {
            std::string name = resolveBrackets(SafeAttrVal(node, "name"));
            if(name.empty()) continue;
            bool in_loop = std::any_of(stack.begin(), stack.end(),
                [](Frame const & f) { return f.is_loop; });
            const char* unit = SafeAttrVal(node, "unit");
            double x, y, z;
            if(in_loop) {
                x = ParseLength(resolveBrackets(SafeAttrVal(node, "x")).c_str(), unit, data.constants);
                y = ParseLength(resolveBrackets(SafeAttrVal(node, "y")).c_str(), unit, data.constants);
                z = ParseLength(resolveBrackets(SafeAttrVal(node, "z")).c_str(), unit, data.constants);
            } else {
                x = ParseLength(SafeAttrVal(node, "x"), unit, data.constants);
                y = ParseLength(SafeAttrVal(node, "y"), unit, data.constants);
                z = ParseLength(SafeAttrVal(node, "z"), unit, data.constants);
            }
            data.positions[name] = Vector3D(x, y, z);
        }
        else if(tag == "rotation") {
            std::string name = resolveBrackets(SafeAttrVal(node, "name"));
            if(name.empty()) continue;
            bool in_loop = std::any_of(stack.begin(), stack.end(),
                [](Frame const & f) { return f.is_loop; });
            const char* unit = SafeAttrVal(node, "unit");
            double rx, ry, rz;
            if(in_loop) {
                rx = ParseAngle(resolveBrackets(SafeAttrVal(node, "x")).c_str(), unit, data.constants);
                ry = ParseAngle(resolveBrackets(SafeAttrVal(node, "y")).c_str(), unit, data.constants);
                rz = ParseAngle(resolveBrackets(SafeAttrVal(node, "z")).c_str(), unit, data.constants);
            } else {
                rx = ParseAngle(SafeAttrVal(node, "x"), unit, data.constants);
                ry = ParseAngle(SafeAttrVal(node, "y"), unit, data.constants);
                rz = ParseAngle(SafeAttrVal(node, "z"), unit, data.constants);
            }
            data.rotations[name] = QuatFromGDMLRotation(rx, ry, rz);
        }
        else if(tag == "matrix") {
            std::string name = resolveBrackets(SafeAttrVal(node, "name"));
            if(name.empty()) continue;
            int coldim = 1;
            const char* cd = SafeAttrVal(node, "coldim");
            if(cd[0] != '\0') { try { coldim = std::stoi(std::string(cd)); } catch(...) { coldim = 1; } }
            if(coldim < 1) coldim = 1;
            std::string values_str = SafeAttrVal(node, "values");
            std::vector<double> values;
            std::istringstream iss(values_str);
            std::string tok;
            while(iss >> tok) {
                values.push_back(EvalExpression(tok, data.constants, data.matrices));
            }
            if(!values.empty()) data.matrices[name] = {coldim, values};
        }
        else if(tag == "scale") {
            std::string sname = SafeAttrVal(node, "name");
            EmitWarning(data, options, "GDML <scale> element '" + sname + "' is not supported and will be ignored");
        }
    }
}

} // anonymous namespace


// ---- Section parsers ----

// Parse the <define> section: constants, positions, rotations
static void ParseDefine(rapidxml::xml_node<>* define_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!define_node) return;

    for(auto* node = define_node->first_node(); node; node = node->next_sibling()) {
        std::string tag(node->name());

        if(tag == "constant") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants, data.matrices);
            if(!name.empty()) {
                if(IsBuiltInConstant(name)) {
                    EmitWarning(data, options, "cannot redefine built-in constant '" + name + "', ignoring");
                } else {
                    if(data.constants.find(name) != data.constants.end()) {
                        EmitWarning(data, options, "duplicate constant name '" + name + "', overwriting");
                    }
                    data.constants[name] = value;
                }
            }
        }
        else if(tag == "quantity") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants, data.matrices);
            const char* unit = SafeAttrVal(node, "unit");
            const char* type_attr = SafeAttrVal(node, "type");
            GDMLQuantityType resolved_type = ResolveQuantityType(type_attr, unit, name, data, options);
            value *= QuantityUnitScale(unit, resolved_type);
            if(!name.empty()) {
                if(IsBuiltInConstant(name)) {
                    EmitWarning(data, options, "cannot redefine built-in constant '" + name + "', ignoring");
                } else {
                    if(data.constants.find(name) != data.constants.end()) {
                        EmitWarning(data, options, "duplicate constant name '" + name + "', overwriting");
                    }
                    data.constants[name] = value;
                    if(resolved_type != GDMLQuantityType::NONE) {
                        data.quantity_types[name] = resolved_type;
                    }
                }
            }
        }
        else if(tag == "position") {
            std::string name = SafeAttrVal(node, "name");
            const char* unit = SafeAttrVal(node, "unit");
            double x = ParseLength(SafeAttrVal(node, "x"), unit, data.constants);
            double y = ParseLength(SafeAttrVal(node, "y"), unit, data.constants);
            double z = ParseLength(SafeAttrVal(node, "z"), unit, data.constants);
            if(!name.empty()) {
                if(data.positions.find(name) != data.positions.end()) {
                    EmitWarning(data, options, "duplicate position name '" + name + "', overwriting");
                }
                data.positions[name] = Vector3D(x, y, z);
            }
        }
        else if(tag == "rotation") {
            std::string name = SafeAttrVal(node, "name");
            const char* unit = SafeAttrVal(node, "unit");
            double rx = ParseAngle(SafeAttrVal(node, "x"), unit, data.constants);
            double ry = ParseAngle(SafeAttrVal(node, "y"), unit, data.constants);
            double rz = ParseAngle(SafeAttrVal(node, "z"), unit, data.constants);
            if(!name.empty()) {
                if(data.rotations.find(name) != data.rotations.end()) {
                    EmitWarning(data, options, "duplicate rotation name '" + name + "', overwriting");
                }
                data.rotations[name] = QuatFromGDMLRotation(rx, ry, rz);
            }
        }
        else if(tag == "scale") {
            std::string sname = SafeAttrVal(node, "name");
            EmitWarning(data, options, "GDML <scale> element '" + sname + "' is not supported and will be ignored");
        }
        else if(tag == "variable") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants, data.matrices);
            if(!name.empty()) {
                if(IsBuiltInConstant(name)) {
                    EmitWarning(data, options, "cannot redefine built-in constant '" + name + "', ignoring");
                } else {
                    if(data.constants.find(name) != data.constants.end()) {
                        EmitWarning(data, options, "duplicate constant name '" + name + "', overwriting");
                    }
                    data.constants[name] = value;
                }
            }
        }
        else if(tag == "matrix") {
            std::string name = SafeAttrVal(node, "name");
            int coldim = 1;
            const char* cd = SafeAttrVal(node, "coldim");
            if(cd[0] != '\0') {
                try { coldim = std::stoi(std::string(cd)); }
                catch(...) { coldim = 1; }
            }
            if(coldim < 1) coldim = 1;
            std::string values_str = SafeAttrVal(node, "values");
            std::vector<double> values;
            std::istringstream iss(values_str);
            std::string tok;
            while(iss >> tok) {
                values.push_back(EvalExpression(tok, data.constants, data.matrices));
            }
            if(!name.empty() && !values.empty()) {
                data.matrices[name] = {coldim, values};
            }
        }
    }
}


// Parse all <materials> sections with instance-scoped name resolution.
// Each name can have multiple definitions (instances). Composition references
// resolve to the current instance of the referenced name at definition time.
// After parsing, multi-instance names are flattened to unique names.
static void ParseAllMaterials(rapidxml::xml_node<>* root_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!root_node) return;

    // Collect all <materials> section nodes
    std::vector<rapidxml::xml_node<>*> mat_sections;
    for(auto* node = root_node->first_node("materials"); node; node = node->next_sibling("materials")) {
        mat_sections.push_back(node);
    }
    if(mat_sections.empty()) return;

    // Instance-aware storage, separated by GDML type.
    // In GDML, isotopes, elements, and materials are separate namespaces.
    // An <element name="Iron"> and <material name="Iron"> are different things.
    // We use separate instance stacks so they don't collide.
    struct MatTypeStacks {
        std::map<std::string, std::vector<GDMLMaterial>> instances;
        std::map<std::string, int> count;
    };
    MatTypeStacks isotope_stacks, element_stacks, material_stacks;

    auto flattenedName = [](std::string const & name, int instance_idx) -> std::string {
        if(instance_idx == 0) return name;
        return name + "__" + std::to_string(instance_idx + 1);
    };

    // Store into the appropriate type stack. Dedup against most recent instance.
    auto storeInStack = [&](MatTypeStacks & stacks, GDMLMaterial const & mat) {
        auto & vec = stacks.instances[mat.name];
        if(!vec.empty() && vec.back() == mat) {
            return;
        }
        vec.push_back(mat);
        stacks.count[mat.name] = (int)vec.size();
    };

    // Resolve a composition reference: look up in isotope, then element, then
    // material stacks (in that order) to find the current instance.
    auto resolveCompRef = [&](std::string const & ref_name) -> std::string {
        // Try isotopes first (most specific)
        auto it = isotope_stacks.count.find(ref_name);
        if(it != isotope_stacks.count.end()) {
            return flattenedName(ref_name, it->second - 1);
        }
        // Then elements
        it = element_stacks.count.find(ref_name);
        if(it != element_stacks.count.end()) {
            return flattenedName(ref_name, it->second - 1);
        }
        // Then materials
        it = material_stacks.count.find(ref_name);
        if(it != material_stacks.count.end()) {
            return flattenedName(ref_name, it->second - 1);
        }
        // Forward reference: assume instance 0
        return ref_name;
    };

    // Process all sections in order
    for(auto* mat_node : mat_sections) {
        for(auto* node = mat_node->first_node(); node; node = node->next_sibling()) {
            std::string tag(node->name());

            if(tag == "isotope") {
                GDMLMaterial mat;
                mat.name = SafeAttrVal(node, "name");
                mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"), data.constants);
                mat.A = SafeParseDouble(SafeAttrVal(node, "N"), data.constants);
                mat.density = 0.0;

                auto* atom_node = node->first_node("atom");
                if(atom_node) {
                    mat.A = SafeParseDouble(SafeAttrVal(atom_node, "value"), data.constants);
                }

                if(!mat.name.empty()) {
                    storeInStack(isotope_stacks, mat);
                }
            }
            else if(tag == "element") {
                GDMLMaterial mat;
                mat.name = SafeAttrVal(node, "name");
                mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"), data.constants);
                mat.density = 0.0;
                mat.A = 0.0;

                auto* atom_node = node->first_node("atom");
                if(atom_node) {
                    const char* a_val = SafeAttrVal(atom_node, "value");
                    mat.A = SafeParseDouble(a_val, data.constants);
                }

                // Fraction children (element defined as isotope mix)
                // Resolve each ref to the current instance
                for(auto* frac = node->first_node("fraction"); frac; frac = frac->next_sibling("fraction")) {
                    std::string ref = SafeAttrVal(frac, "ref");
                    double n = SafeParseDouble(SafeAttrVal(frac, "n"), data.constants);
                    if(!ref.empty()) {
                        mat.composition[resolveCompRef(ref)] = n;
                    }
                }

                if(!mat.name.empty()) {
                    storeInStack(element_stacks, mat);
                }
            }
            else if(tag == "material") {
                GDMLMaterial mat;
                mat.name = SafeAttrVal(node, "name");
                mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"), data.constants);
                mat.density = 0.0;
                mat.A = 0.0;

                auto* d_node = node->first_node("D");
                if(d_node) {
                    const char* d_val_str = SafeAttrVal(d_node, "value");
                    mat.density = SafeParseDouble(d_val_str, data.constants);

                    std::string d_val_trimmed = TrimWhitespace(d_val_str);
                    auto qty_it = data.quantity_types.find(d_val_trimmed);
                    bool is_density_qty = (qty_it != data.quantity_types.end()
                                           && qty_it->second == GDMLQuantityType::DENSITY);

                    if(!is_density_qty) {
                        const char* d_unit = SafeAttrVal(d_node, "unit");
                        if(d_unit && d_unit[0] != '\0') {
                            std::string du(d_unit);
                            if(du == "kg/m3") {
                                mat.density /= 1000.0;
                            } else if(du == "mg/cm3") {
                                mat.density /= 1000.0;
                            } else if(du == "g/cm3") {
                                // No conversion needed
                            } else {
                                EmitWarning(data, options, "unrecognized density unit '" + du + "' for material '" + mat.name + "', assuming g/cm3");
                            }
                        }
                    }
                }

                auto* atom_node = node->first_node("atom");
                if(atom_node) {
                    mat.A = SafeParseDouble(SafeAttrVal(atom_node, "value"), data.constants);
                }

                // Fraction children: resolve each ref to the current instance
                for(auto* frac = node->first_node("fraction"); frac; frac = frac->next_sibling("fraction")) {
                    std::string ref = SafeAttrVal(frac, "ref");
                    double n = SafeParseDouble(SafeAttrVal(frac, "n"), data.constants);
                    if(!ref.empty()) {
                        mat.composition[resolveCompRef(ref)] = n;
                    }
                }

                // Composite children: resolve each ref to the current instance
                for(auto* comp = node->first_node("composite"); comp; comp = comp->next_sibling("composite")) {
                    std::string ref = SafeAttrVal(comp, "ref");
                    double n = SafeParseDouble(SafeAttrVal(comp, "n"), data.constants);
                    if(!ref.empty()) {
                        mat.composition[resolveCompRef(ref)] = n;
                    }
                }

                if(!mat.name.empty()) {
                    storeInStack(material_stacks, mat);
                }
            }
        }
    }

    // Flatten all three type stacks into data.materials with unique names.
    // Each type stack is independent, so <element name="Iron"> and
    // <material name="Iron"> get separate entries.
    // Within each stack, multi-instance names get __2, __3 suffixes.
    // The material_instance_counts tracks the material stack only (since
    // that is what volumes reference via materialref).
    auto flattenStack = [&](MatTypeStacks const & stacks) {
        for(auto const & pair : stacks.instances) {
            std::string const & name = pair.first;
            std::vector<GDMLMaterial> const & instances = pair.second;
            for(int i = 0; i < (int)instances.size(); ++i) {
                GDMLMaterial mat = instances[i];
                std::string flat_name = flattenedName(name, i);
                // If this flat_name collides with one already in data.materials
                // (from a different type stack), append a type suffix
                if(data.materials.find(flat_name) != data.materials.end()) {
                    // The existing entry came from a different type stack.
                    // Keep both by giving the new one a distinct name.
                    flat_name = flat_name + "__mat";
                }
                mat.name = flat_name;
                data.materials[flat_name] = mat;
            }
        }
    };
    flattenStack(isotope_stacks);
    flattenStack(element_stacks);
    flattenStack(material_stacks);
    // Instance counts from the material stack only (for volume ambiguity checks)
    for(auto const & pair : material_stacks.instances) {
        data.material_instance_counts[pair.first] = (int)pair.second.size();
    }
}


// Parse all <solids> sections: convert GDML solids to SIREN Geometry objects.
// Iterates all <solids> children of root_node, running the primitive parse pass
// for each section first, then performing iterative boolean resolution across
// all sections.
static void ParseAllSolids(rapidxml::xml_node<>* root_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!root_node) return;

    // Collect all <solids> section nodes
    std::vector<rapidxml::xml_node<>*> solids_sections;
    for(auto* node = root_node->first_node("solids"); node; node = node->next_sibling("solids")) {
        solids_sections.push_back(node);
    }
    if(solids_sections.empty()) return;

    // Tolerance for comparing angular values to full-rotation defaults
    static const double ANG_TOL = 1e-6;

    // Helper lambda: parse a single non-boolean solid node and return geometry
    auto parsePrimitive = [&](rapidxml::xml_node<>* node) -> std::shared_ptr<Geometry> {
        std::string tag(node->name());
        std::shared_ptr<Geometry> geo;

        const char* lunit = SafeAttrVal(node, "lunit");
        const char* aunit = SafeAttrVal(node, "aunit");
        double lscale = LengthScale(lunit);
        double ascale = AngleScale(aunit);
        std::string name = SafeAttrVal(node, "name");

        if(tag == "box") {
            // GDML <box> x,y,z are HALF-widths; SIREN Box takes FULL widths
            double hx = SafeParseDouble(SafeAttrVal(node, "x"), data.constants) * lscale;
            double hy = SafeParseDouble(SafeAttrVal(node, "y"), data.constants) * lscale;
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;
            geo = Box(hx * 2.0, hy * 2.0, hz * 2.0).create();
        }
        else if(tag == "sphere") {
            double rmin = SafeParseDouble(SafeAttrVal(node, "rmin"), data.constants) * lscale;
            double rmax = SafeParseDouble(SafeAttrVal(node, "rmax"), data.constants) * lscale;

            double startphi = 0.0;
            double deltaphi = 2.0 * M_PI;
            double starttheta = 0.0;
            double deltatheta = M_PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
            const char* st_val = SafeAttrVal(node, "starttheta");
            if(st_val[0] != '\0') starttheta = SafeParseDouble(st_val, data.constants) * ascale;
            const char* dt_val = SafeAttrVal(node, "deltatheta");
            if(dt_val[0] != '\0') deltatheta = SafeParseDouble(dt_val, data.constants) * ascale;

            geo = Sphere(rmax, rmin, startphi, deltaphi, starttheta, deltatheta).create();
        }
        else if(tag == "orb") {
            double r = SafeParseDouble(SafeAttrVal(node, "r"), data.constants) * lscale;
            geo = Sphere(r, 0).create();
        }
        else if(tag == "tube" || tag == "tubs") {
            double rmin = SafeParseDouble(SafeAttrVal(node, "rmin"), data.constants) * lscale;
            double rmax = SafeParseDouble(SafeAttrVal(node, "rmax"), data.constants) * lscale;
            // GDML <tube> z is HALF-height; SIREN Cylinder takes full height
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;

            double startphi = 0.0;
            double deltaphi = 2.0 * M_PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;

            geo = Cylinder(rmax, rmin, hz * 2.0, startphi, deltaphi).create();
        }
        else if(tag == "cone") {
            double rmin1 = SafeParseDouble(SafeAttrVal(node, "rmin1"), data.constants) * lscale;
            double rmax1 = SafeParseDouble(SafeAttrVal(node, "rmax1"), data.constants) * lscale;
            double rmin2 = SafeParseDouble(SafeAttrVal(node, "rmin2"), data.constants) * lscale;
            double rmax2 = SafeParseDouble(SafeAttrVal(node, "rmax2"), data.constants) * lscale;
            // GDML <cone> z is HALF-height; SIREN Cone takes full height
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;

            double startphi = 0.0;
            double deltaphi = 2.0 * M_PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;

            geo = Cone(rmin1, rmax1, rmin2, rmax2, hz * 2.0, startphi, deltaphi).create();
        }
        else if(tag == "trd") {
            // GDML <trd> uses half-widths; SIREN Trd also uses half-widths
            double dx1 = SafeParseDouble(SafeAttrVal(node, "x1"), data.constants) * lscale;
            double dx2 = SafeParseDouble(SafeAttrVal(node, "x2"), data.constants) * lscale;
            double dy1 = SafeParseDouble(SafeAttrVal(node, "y1"), data.constants) * lscale;
            double dy2 = SafeParseDouble(SafeAttrVal(node, "y2"), data.constants) * lscale;
            double dz  = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;
            geo = Trd(dx1, dx2, dy1, dy2, dz).create();
        }
        else if(tag == "polycone") {
            double startphi = 0.0;
            double deltaphi = 2.0 * M_PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;

            std::vector<double> z_planes;
            std::vector<double> rmin_vec;
            std::vector<double> rmax_vec;

            for(auto* zp = node->first_node("zplane"); zp; zp = zp->next_sibling("zplane")) {
                double z = SafeParseDouble(SafeAttrVal(zp, "z"), data.constants) * lscale;
                double rmin = SafeParseDouble(SafeAttrVal(zp, "rmin"), data.constants) * lscale;
                double rmax = SafeParseDouble(SafeAttrVal(zp, "rmax"), data.constants) * lscale;
                z_planes.push_back(z);
                rmin_vec.push_back(rmin);
                rmax_vec.push_back(rmax);
            }

            if(z_planes.size() >= 2) {
                bool needs_sort = false;
                for(size_t i = 1; i < z_planes.size(); ++i) {
                    if(z_planes[i] < z_planes[i-1]) { needs_sort = true; break; }
                }
                if(needs_sort) {
                    std::vector<size_t> idx(z_planes.size());
                    for(size_t i = 0; i < idx.size(); ++i) idx[i] = i;
                    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
                        return z_planes[a] < z_planes[b];
                    });
                    std::vector<double> zs(z_planes.size()), rns(z_planes.size()), rxs(z_planes.size());
                    for(size_t i = 0; i < idx.size(); ++i) {
                        zs[i] = z_planes[idx[i]];
                        rns[i] = rmin_vec[idx[i]];
                        rxs[i] = rmax_vec[idx[i]];
                    }
                    z_planes = std::move(zs);
                    rmin_vec = std::move(rns);
                    rmax_vec = std::move(rxs);
                }
                bool z_valid = true;
                for(size_t i = 1; i < z_planes.size(); ++i) {
                    if(z_planes[i] <= z_planes[i-1]) { z_valid = false; break; }
                }
                if(z_valid) {
                    geo = Polycone(z_planes, rmin_vec, rmax_vec, startphi, deltaphi).create();
                } else {
                    EmitWarning(data, options, "polycone '" + name + "' has non-monotonic z-planes after sorting; skipping");
                }
            }
        }
        else if(tag == "polyhedra") {
            int numSide = 0;
            {
                std::string ns_str(SafeAttrVal(node, "numsides"));
                if(!ns_str.empty()) {
                    try {
                        numSide = std::stoi(ns_str);
                    } catch(std::invalid_argument &) {
                        numSide = 0;
                    } catch(std::out_of_range &) {
                        numSide = 0;
                    }
                }
            }

            if(numSide <= 0) return nullptr;

            double startphi = SafeParseDouble(SafeAttrVal(node, "startphi"), data.constants) * ascale;

            // Check deltaphi for partial polyhedra
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                double deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
                if(std::fabs(deltaphi - 2.0 * M_PI) > ANG_TOL) {
                    EmitWarning(data, options, "polyhedra '" + name + "' has partial angular extent (deltaphi=" + std::to_string(deltaphi) + "); SIREN creates full rotation");
                }
            }

            std::vector<double> z_planes;
            std::vector<double> rmin_vec;
            std::vector<double> rmax_vec;

            for(auto* zp = node->first_node("zplane"); zp; zp = zp->next_sibling("zplane")) {
                double z = SafeParseDouble(SafeAttrVal(zp, "z"), data.constants) * lscale;
                double rmin = SafeParseDouble(SafeAttrVal(zp, "rmin"), data.constants) * lscale;
                double rmax = SafeParseDouble(SafeAttrVal(zp, "rmax"), data.constants) * lscale;
                z_planes.push_back(z);
                rmin_vec.push_back(rmin);
                rmax_vec.push_back(rmax);
            }

            if(z_planes.size() >= 2 && numSide >= 3) {
                bool needs_sort = false;
                for(size_t i = 1; i < z_planes.size(); ++i) {
                    if(z_planes[i] < z_planes[i-1]) { needs_sort = true; break; }
                }
                if(needs_sort) {
                    std::vector<size_t> idx(z_planes.size());
                    for(size_t i = 0; i < idx.size(); ++i) idx[i] = i;
                    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
                        return z_planes[a] < z_planes[b];
                    });
                    std::vector<double> zs(z_planes.size()), rns(z_planes.size()), rxs(z_planes.size());
                    for(size_t i = 0; i < idx.size(); ++i) {
                        zs[i] = z_planes[idx[i]];
                        rns[i] = rmin_vec[idx[i]];
                        rxs[i] = rmax_vec[idx[i]];
                    }
                    z_planes = std::move(zs);
                    rmin_vec = std::move(rns);
                    rmax_vec = std::move(rxs);
                }
                bool z_valid = true;
                for(size_t i = 1; i < z_planes.size(); ++i) {
                    if(z_planes[i] <= z_planes[i-1]) { z_valid = false; break; }
                }
                if(z_valid) {
                    geo = Polyhedra(numSide, startphi, z_planes, rmin_vec, rmax_vec).create();
                } else {
                    EmitWarning(data, options, "polyhedra '" + name + "' has non-monotonic z-planes after sorting; skipping");
                }
            }
        }
        else if(tag == "xtru") {
            std::vector<std::vector<double>> polygon;
            std::vector<ExtrPoly::ZSection> zsections;

            for(auto* vtx = node->first_node("twoDimVertex"); vtx; vtx = vtx->next_sibling("twoDimVertex")) {
                double vx = SafeParseDouble(SafeAttrVal(vtx, "x"), data.constants) * lscale;
                double vy = SafeParseDouble(SafeAttrVal(vtx, "y"), data.constants) * lscale;
                std::vector<double> vert;
                vert.push_back(vx);
                vert.push_back(vy);
                polygon.push_back(vert);
            }

            for(auto* sec = node->first_node("section"); sec; sec = sec->next_sibling("section")) {
                double zpos = SafeParseDouble(SafeAttrVal(sec, "zPosition"), data.constants) * lscale;
                double xoff = SafeParseDouble(SafeAttrVal(sec, "xOffset"), data.constants) * lscale;
                double yoff = SafeParseDouble(SafeAttrVal(sec, "yOffset"), data.constants) * lscale;
                double scale = SafeParseDouble(SafeAttrVal(sec, "scalingFactor"), data.constants);
                if(scale == 0.0) scale = 1.0;
                double offset[2] = {xoff, yoff};
                zsections.push_back(ExtrPoly::ZSection(zpos, offset, scale));
            }

            if(!polygon.empty() && !zsections.empty()) {
                geo = ExtrPoly(polygon, zsections).create();
            }
        }
        else if(tag == "torus") {
            double rmin_t = SafeParseDouble(SafeAttrVal(node, "rmin"), data.constants) * lscale;
            double rmax_t = SafeParseDouble(SafeAttrVal(node, "rmax"), data.constants) * lscale;
            double rtor = SafeParseDouble(SafeAttrVal(node, "rtor"), data.constants) * lscale;

            double startphi = 0.0;
            double deltaphi = 2.0 * M_PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;

            if(rtor > 0 && rmax_t > 0) {
                geo = Torus(rtor, rmax_t, rmin_t, startphi, deltaphi).create();
            }
        }
        else if(tag == "eltube") {
            double dx = SafeParseDouble(SafeAttrVal(node, "dx"), data.constants) * lscale;
            double dy = SafeParseDouble(SafeAttrVal(node, "dy"), data.constants) * lscale;
            double dz = SafeParseDouble(SafeAttrVal(node, "dz"), data.constants) * lscale;
            if(dx > 0 && dy > 0 && dz > 0) {
                geo = EllipticalTube(dx, dy, dz).create();
            }
        }
        else if(tag == "cutTube") {
            double rmin = SafeParseDouble(SafeAttrVal(node, "rmin"), data.constants) * lscale;
            double rmax = SafeParseDouble(SafeAttrVal(node, "rmax"), data.constants) * lscale;
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;

            double startphi = 0.0;
            double deltaphi = 2.0 * M_PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;

            double lowX = SafeParseDouble(SafeAttrVal(node, "lowX"), data.constants);
            double lowY = SafeParseDouble(SafeAttrVal(node, "lowY"), data.constants);
            double lowZ = SafeParseDouble(SafeAttrVal(node, "lowZ"), data.constants);
            double highX = SafeParseDouble(SafeAttrVal(node, "highX"), data.constants);
            double highY = SafeParseDouble(SafeAttrVal(node, "highY"), data.constants);
            double highZ = SafeParseDouble(SafeAttrVal(node, "highZ"), data.constants);

            if(rmax > 0 && hz > 0) {
                Vector3D low_norm(lowX, lowY, lowZ);
                Vector3D high_norm(highX, highY, highZ);
                geo = CutTube(rmin, rmax, hz, low_norm, high_norm, startphi, deltaphi).create();
            }
        }
        else if(tag == "trap") {
            double dz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;
            double theta = SafeParseDouble(SafeAttrVal(node, "theta"), data.constants) * ascale;
            double phi = SafeParseDouble(SafeAttrVal(node, "phi"), data.constants) * ascale;
            double dy1 = SafeParseDouble(SafeAttrVal(node, "y1"), data.constants) * lscale;
            double dx1 = SafeParseDouble(SafeAttrVal(node, "x1"), data.constants) * lscale;
            double dx2 = SafeParseDouble(SafeAttrVal(node, "x2"), data.constants) * lscale;
            double alpha1 = SafeParseDouble(SafeAttrVal(node, "alpha1"), data.constants) * ascale;
            double dy2 = SafeParseDouble(SafeAttrVal(node, "y2"), data.constants) * lscale;
            double dx3 = SafeParseDouble(SafeAttrVal(node, "x3"), data.constants) * lscale;
            double dx4 = SafeParseDouble(SafeAttrVal(node, "x4"), data.constants) * lscale;
            double alpha2 = SafeParseDouble(SafeAttrVal(node, "alpha2"), data.constants) * ascale;
            geo = Trap(dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2).create();
        }
        else if(tag == "ellipsoid") {
            double ax = SafeParseDouble(SafeAttrVal(node, "ax"), data.constants) * lscale;
            double by = SafeParseDouble(SafeAttrVal(node, "by"), data.constants) * lscale;
            double cz = SafeParseDouble(SafeAttrVal(node, "cz"), data.constants) * lscale;
            double zcut1 = -cz;
            double zcut2 = cz;
            const char* z1_val = SafeAttrVal(node, "zcut1");
            if(z1_val[0] != '\0') zcut1 = SafeParseDouble(z1_val, data.constants) * lscale;
            const char* z2_val = SafeAttrVal(node, "zcut2");
            if(z2_val[0] != '\0') zcut2 = SafeParseDouble(z2_val, data.constants) * lscale;
            if(ax > 0 && by > 0 && cz > 0) {
                geo = Ellipsoid(ax, by, cz, zcut1, zcut2).create();
            }
        }
        else if(tag == "para") {
            double dx = SafeParseDouble(SafeAttrVal(node, "x"), data.constants) * lscale;
            double dy = SafeParseDouble(SafeAttrVal(node, "y"), data.constants) * lscale;
            double dz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;
            double alpha = SafeParseDouble(SafeAttrVal(node, "alpha"), data.constants) * ascale;
            double theta = SafeParseDouble(SafeAttrVal(node, "theta"), data.constants) * ascale;
            double phi = SafeParseDouble(SafeAttrVal(node, "phi"), data.constants) * ascale;
            if(dx > 0 && dy > 0 && dz > 0) {
                geo = Para(dx, dy, dz, alpha, theta, phi).create();
            }
        }
        else if(tag == "genericPolycone") {
            double startphi = 0.0;
            double deltaphi = 2.0 * M_PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') {
                startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            }
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
            }

            std::vector<double> r_vec, z_vec;
            for(auto* rz = node->first_node("rzpoint"); rz; rz = rz->next_sibling("rzpoint")) {
                double r = SafeParseDouble(SafeAttrVal(rz, "r"), data.constants) * lscale;
                double z = SafeParseDouble(SafeAttrVal(rz, "z"), data.constants) * lscale;
                r_vec.push_back(r);
                z_vec.push_back(z);
            }

            if(r_vec.size() >= 3) {
                bool has_phi_cut = (std::fabs(startphi) > ANG_TOL) ||
                                   (std::fabs(deltaphi - 2.0 * M_PI) > ANG_TOL);
                if(has_phi_cut) {
                    geo = GenericPolycone(r_vec, z_vec, startphi, deltaphi).create();
                } else {
                    geo = GenericPolycone(r_vec, z_vec).create();
                }
            }
        }
        else if(tag == "tessellated") {
            std::vector<std::array<math::Vector3D, 3>> triangles;

            for(auto* facet = node->first_node(); facet; facet = facet->next_sibling()) {
                std::string ftag(facet->name());
                if(ftag == "triangular") {
                    auto resolveVertex = [&](const char* attr) -> math::Vector3D {
                        std::string ref = SafeAttrVal(facet, attr);
                        auto it = data.positions.find(ref);
                        if(it != data.positions.end()) return it->second;
                        return math::Vector3D(0, 0, 0);
                    };
                    math::Vector3D v1 = resolveVertex("vertex1");
                    math::Vector3D v2 = resolveVertex("vertex2");
                    math::Vector3D v3 = resolveVertex("vertex3");
                    triangles.push_back({{v1, v2, v3}});
                } else if(ftag == "quadrangular") {
                    auto resolveVertex = [&](const char* attr) -> math::Vector3D {
                        std::string ref = SafeAttrVal(facet, attr);
                        auto it = data.positions.find(ref);
                        if(it != data.positions.end()) return it->second;
                        return math::Vector3D(0, 0, 0);
                    };
                    math::Vector3D v1 = resolveVertex("vertex1");
                    math::Vector3D v2 = resolveVertex("vertex2");
                    math::Vector3D v3 = resolveVertex("vertex3");
                    math::Vector3D v4 = resolveVertex("vertex4");
                    triangles.push_back({{v1, v2, v3}});
                    triangles.push_back({{v1, v3, v4}});
                }
            }

            if(!triangles.empty()) {
                geo = TriangularMesh(triangles).create();
            }
        }
        else if(tag == "arb8") {
            double dz = SafeParseDouble(SafeAttrVal(node, "dz"), data.constants) * lscale;
            double v1x = SafeParseDouble(SafeAttrVal(node, "v1x"), data.constants) * lscale;
            double v1y = SafeParseDouble(SafeAttrVal(node, "v1y"), data.constants) * lscale;
            double v2x = SafeParseDouble(SafeAttrVal(node, "v2x"), data.constants) * lscale;
            double v2y = SafeParseDouble(SafeAttrVal(node, "v2y"), data.constants) * lscale;
            double v3x = SafeParseDouble(SafeAttrVal(node, "v3x"), data.constants) * lscale;
            double v3y = SafeParseDouble(SafeAttrVal(node, "v3y"), data.constants) * lscale;
            double v4x = SafeParseDouble(SafeAttrVal(node, "v4x"), data.constants) * lscale;
            double v4y = SafeParseDouble(SafeAttrVal(node, "v4y"), data.constants) * lscale;
            double v5x = SafeParseDouble(SafeAttrVal(node, "v5x"), data.constants) * lscale;
            double v5y = SafeParseDouble(SafeAttrVal(node, "v5y"), data.constants) * lscale;
            double v6x = SafeParseDouble(SafeAttrVal(node, "v6x"), data.constants) * lscale;
            double v6y = SafeParseDouble(SafeAttrVal(node, "v6y"), data.constants) * lscale;
            double v7x = SafeParseDouble(SafeAttrVal(node, "v7x"), data.constants) * lscale;
            double v7y = SafeParseDouble(SafeAttrVal(node, "v7y"), data.constants) * lscale;
            double v8x = SafeParseDouble(SafeAttrVal(node, "v8x"), data.constants) * lscale;
            double v8y = SafeParseDouble(SafeAttrVal(node, "v8y"), data.constants) * lscale;

            // Vertices 1-4 at z=-dz, 5-8 at z=+dz
            math::Vector3D verts[8] = {
                math::Vector3D(v1x, v1y, -dz),
                math::Vector3D(v2x, v2y, -dz),
                math::Vector3D(v3x, v3y, -dz),
                math::Vector3D(v4x, v4y, -dz),
                math::Vector3D(v5x, v5y,  dz),
                math::Vector3D(v6x, v6y,  dz),
                math::Vector3D(v7x, v7y,  dz),
                math::Vector3D(v8x, v8y,  dz)
            };

            // Tessellate 6 faces into triangles.
            // Bottom face (v1-v4, at -dz): wound CCW when viewed from -z
            // Top face (v5-v8, at +dz): wound CCW when viewed from +z
            // Side faces connect bottom[i] to top[i]
            std::vector<std::array<math::Vector3D, 3>> triangles;
            triangles.reserve(12);

            auto addQuad = [&](math::Vector3D const & a, math::Vector3D const & b,
                              math::Vector3D const & c, math::Vector3D const & d) {
                triangles.push_back({{a, b, c}});
                triangles.push_back({{a, c, d}});
            };

            // Bottom face (outward normal in -z): v4, v3, v2, v1
            addQuad(verts[3], verts[2], verts[1], verts[0]);
            // Top face (outward normal in +z): v5, v6, v7, v8
            addQuad(verts[4], verts[5], verts[6], verts[7]);
            // Side faces (outward normals face outward)
            addQuad(verts[0], verts[1], verts[5], verts[4]);
            addQuad(verts[1], verts[2], verts[6], verts[5]);
            addQuad(verts[2], verts[3], verts[7], verts[6]);
            addQuad(verts[3], verts[0], verts[4], verts[7]);

            if(dz > 0) {
                geo = TriangularMesh(triangles).create();
            }
        }
        else if(!tag.empty()) {
            EmitWarning(data, options, "unsupported solid type '" + tag + "' (name='" + name + "'), skipping");
        }

        return geo;
    };

    // Instance-aware storage: each name can have multiple definitions.
    // solids_instances[name] = vector of shared_ptr (one per definition)
    // solid_instance_count[name] = how many definitions seen so far
    std::map<std::string, std::vector<std::shared_ptr<Geometry>>> solids_instances;
    std::map<std::string, int> solid_instance_count;

    // Helper: store a solid as a new instance.
    // If the new geometry is equal to the current (most recent) instance,
    // skip it — identical re-declarations don't create a new instance.
    auto storeSolid = [&](std::string const & name, std::shared_ptr<Geometry> geo) {
        auto & vec = solids_instances[name];
        if(!vec.empty() && *vec.back() == *geo) {
            // Identical to current instance; no-op (keeps same instance count)
            return;
        }
        vec.push_back(geo);
        solid_instance_count[name] = (int)vec.size();
    };

    // Helper: look up a solid operand by name using instance tracking.
    // If the name has been defined, returns the latest (current) instance.
    // If not defined yet, returns nullptr (forward reference).
    auto lookupSolid = [&](std::string const & name) -> std::shared_ptr<Geometry> {
        auto it = solid_instance_count.find(name);
        if(it != solid_instance_count.end()) {
            return solids_instances[name][it->second - 1];
        }
        return nullptr;
    };

    // Modified boolean parser that uses instance-aware lookup
    auto parseBooleanInstanced = [&](rapidxml::xml_node<>* node) -> std::shared_ptr<Geometry> {
        std::string tag(node->name());
        std::string name = SafeAttrVal(node, "name");

        BooleanOperation op;
        if(tag == "subtraction")   op = BooleanOperation::SUBTRACTION;
        else if(tag == "union")    op = BooleanOperation::UNION;
        else                       op = BooleanOperation::INTERSECTION;

        auto* first_node = node->first_node("first");
        auto* second_node = node->first_node("second");

        if(!first_node || !second_node) return nullptr;

        std::string first_ref = SafeAttrVal(first_node, "ref");
        std::string second_ref = SafeAttrVal(second_node, "ref");

        auto left = lookupSolid(first_ref);
        auto right_base = lookupSolid(second_ref);

        if(!left || !right_base) return nullptr;

        // The first solid may have a position and rotation (firstposition/firstrotation)
        Vector3D first_pos(0, 0, 0);
        Quaternion first_rot;

        auto* fpos_node = node->first_node("firstposition");
        if(fpos_node) {
            const char* punit = SafeAttrVal(fpos_node, "unit");
            double px = ParseLength(SafeAttrVal(fpos_node, "x"), punit, data.constants);
            double py = ParseLength(SafeAttrVal(fpos_node, "y"), punit, data.constants);
            double pz = ParseLength(SafeAttrVal(fpos_node, "z"), punit, data.constants);
            first_pos = Vector3D(px, py, pz);
        }

        auto* fposref_node = node->first_node("firstpositionref");
        if(fposref_node) {
            std::string ref = SafeAttrVal(fposref_node, "ref");
            auto it = data.positions.find(ref);
            if(it != data.positions.end()) {
                first_pos = it->second;
            }
        }

        auto* frot_node = node->first_node("firstrotation");
        if(frot_node) {
            const char* runit = SafeAttrVal(frot_node, "unit");
            double rx = ParseAngle(SafeAttrVal(frot_node, "x"), runit, data.constants);
            double ry = ParseAngle(SafeAttrVal(frot_node, "y"), runit, data.constants);
            double rz = ParseAngle(SafeAttrVal(frot_node, "z"), runit, data.constants);
            first_rot = QuatFromGDMLRotation(rx, ry, rz);
        }

        auto* frotref_node = node->first_node("firstrotationref");
        if(frotref_node) {
            std::string ref = SafeAttrVal(frotref_node, "ref");
            auto it = data.rotations.find(ref);
            if(it != data.rotations.end()) {
                first_rot = it->second;
            }
        }

        // The second solid may have a position and rotation relative to the first
        Vector3D rel_pos(0, 0, 0);
        Quaternion rel_rot;

        auto* pos_node = node->first_node("position");
        if(pos_node) {
            const char* punit = SafeAttrVal(pos_node, "unit");
            double px = ParseLength(SafeAttrVal(pos_node, "x"), punit, data.constants);
            double py = ParseLength(SafeAttrVal(pos_node, "y"), punit, data.constants);
            double pz = ParseLength(SafeAttrVal(pos_node, "z"), punit, data.constants);
            rel_pos = Vector3D(px, py, pz);
        }

        auto* posref_node = node->first_node("positionref");
        if(posref_node) {
            std::string ref = SafeAttrVal(posref_node, "ref");
            auto it = data.positions.find(ref);
            if(it != data.positions.end()) {
                rel_pos = it->second;
            }
        }

        auto* rot_node = node->first_node("rotation");
        if(rot_node) {
            const char* runit = SafeAttrVal(rot_node, "unit");
            double rx = ParseAngle(SafeAttrVal(rot_node, "x"), runit, data.constants);
            double ry = ParseAngle(SafeAttrVal(rot_node, "y"), runit, data.constants);
            double rz = ParseAngle(SafeAttrVal(rot_node, "z"), runit, data.constants);
            rel_rot = QuatFromGDMLRotation(rx, ry, rz);
        }

        auto* rotref_node = node->first_node("rotationref");
        if(rotref_node) {
            std::string ref = SafeAttrVal(rotref_node, "ref");
            auto it = data.rotations.find(ref);
            if(it != data.rotations.end()) {
                rel_rot = it->second;
            }
        }

        // Apply first-operand placement if specified
        bool has_first_placement = (first_pos.magnitude() > 0 || first_rot != Quaternion());
        std::shared_ptr<Geometry> left_placed;
        if(has_first_placement) {
            left_placed = left->create();
            Placement first_placement(first_pos, first_rot);
            left_placed->SetPlacement(first_placement);
        }

        auto right = right_base->create();
        Placement rel_placement(rel_pos, rel_rot);
        right->SetPlacement(rel_placement);

        auto left_final = has_first_placement
            ? std::const_pointer_cast<const Geometry>(left_placed)
            : std::const_pointer_cast<const Geometry>(left);

        return std::make_shared<BooleanGeometry>(op, left_final,
            std::const_pointer_cast<const Geometry>(right));
    };

    // Deferred boolean info: records a boolean node that could not be resolved
    // in the first pass due to forward references. Records which instance index
    // to use for each operand when resolved later.
    struct DeferredBoolean {
        rapidxml::xml_node<>* node;
        std::string name;
        std::string first_ref;
        std::string second_ref;
        // Instance index to use for forward-referenced operands (0 = first definition).
        // -1 means "was already defined, use latest at resolution time"
        int first_instance;
        int second_instance;
    };

    // First pass: process ALL nodes (primitives and booleans) in document order.
    // Booleans that reference not-yet-defined operands are deferred.
    std::vector<DeferredBoolean> deferred;

    for(auto* solids_node : solids_sections) {
        for(auto* node = solids_node->first_node(); node; node = node->next_sibling()) {
            std::string tag(node->name());
            std::string name = SafeAttrVal(node, "name");
            if(name.empty()) continue;

            bool is_boolean = (tag == "subtraction" || tag == "union" || tag == "intersection");

            if(!is_boolean) {
                auto geo = parsePrimitive(node);
                if(geo) {
                    storeSolid(name, geo);
                }
            } else {
                // Try to resolve the boolean using current instance tracking
                auto* first_child = node->first_node("first");
                auto* second_child = node->first_node("second");
                if(!first_child || !second_child) continue;

                std::string first_ref = SafeAttrVal(first_child, "ref");
                std::string second_ref = SafeAttrVal(second_child, "ref");

                auto left = lookupSolid(first_ref);
                auto right = lookupSolid(second_ref);

                if(left && right) {
                    // Both operands are defined; resolve now
                    auto geo = parseBooleanInstanced(node);
                    if(geo) {
                        storeSolid(name, geo);
                    }
                } else {
                    // At least one operand is a forward reference; defer.
                    DeferredBoolean db;
                    db.node = node;
                    db.name = name;
                    db.first_ref = first_ref;
                    db.second_ref = second_ref;
                    // For already-defined operands, record current instance index
                    // For forward-referenced operands, use instance 0
                    db.first_instance = left ? (solid_instance_count[first_ref] - 1) : 0;
                    db.second_instance = right ? (solid_instance_count[second_ref] - 1) : 0;
                    deferred.push_back(db);
                }
            }
        }
    }

    // Second pass: iteratively resolve deferred booleans.
    // For forward-referenced operands, use the instance index recorded at deferral time.
    bool made_progress = true;
    while(made_progress) {
        made_progress = false;
        for(auto & db : deferred) {
            // Skip if already resolved
            if(solid_instance_count.find(db.name) != solid_instance_count.end()) {
                // Check if this specific boolean was already stored
                // (name might have multiple instances; check if we already processed this deferred entry)
                if(db.node == nullptr) continue;
            }

            // Look up operands using recorded instance indices
            auto first_it = solids_instances.find(db.first_ref);
            auto second_it = solids_instances.find(db.second_ref);

            if(first_it == solids_instances.end() || second_it == solids_instances.end()) continue;
            if(db.first_instance >= (int)first_it->second.size()) continue;
            if(db.second_instance >= (int)second_it->second.size()) continue;

            auto left = first_it->second[db.first_instance];
            auto right_base = second_it->second[db.second_instance];
            if(!left || !right_base) continue;

            // Parse position/rotation from the boolean node
            std::string tag(db.node->name());
            BooleanOperation op;
            if(tag == "subtraction")   op = BooleanOperation::SUBTRACTION;
            else if(tag == "union")    op = BooleanOperation::UNION;
            else                       op = BooleanOperation::INTERSECTION;

            // First-operand placement
            Vector3D first_pos(0, 0, 0);
            Quaternion first_rot;

            auto* fpos_node = db.node->first_node("firstposition");
            if(fpos_node) {
                const char* punit = SafeAttrVal(fpos_node, "unit");
                double px = ParseLength(SafeAttrVal(fpos_node, "x"), punit, data.constants);
                double py = ParseLength(SafeAttrVal(fpos_node, "y"), punit, data.constants);
                double pz = ParseLength(SafeAttrVal(fpos_node, "z"), punit, data.constants);
                first_pos = Vector3D(px, py, pz);
            }

            auto* fposref_node = db.node->first_node("firstpositionref");
            if(fposref_node) {
                std::string ref = SafeAttrVal(fposref_node, "ref");
                auto it = data.positions.find(ref);
                if(it != data.positions.end()) {
                    first_pos = it->second;
                }
            }

            auto* frot_node = db.node->first_node("firstrotation");
            if(frot_node) {
                const char* runit = SafeAttrVal(frot_node, "unit");
                double rx = ParseAngle(SafeAttrVal(frot_node, "x"), runit, data.constants);
                double ry = ParseAngle(SafeAttrVal(frot_node, "y"), runit, data.constants);
                double rz = ParseAngle(SafeAttrVal(frot_node, "z"), runit, data.constants);
                first_rot = QuatFromGDMLRotation(rx, ry, rz);
            }

            auto* frotref_node = db.node->first_node("firstrotationref");
            if(frotref_node) {
                std::string ref = SafeAttrVal(frotref_node, "ref");
                auto it = data.rotations.find(ref);
                if(it != data.rotations.end()) {
                    first_rot = it->second;
                }
            }

            // Second-operand placement
            Vector3D rel_pos(0, 0, 0);
            Quaternion rel_rot;

            auto* pos_node = db.node->first_node("position");
            if(pos_node) {
                const char* punit = SafeAttrVal(pos_node, "unit");
                double px = ParseLength(SafeAttrVal(pos_node, "x"), punit, data.constants);
                double py = ParseLength(SafeAttrVal(pos_node, "y"), punit, data.constants);
                double pz = ParseLength(SafeAttrVal(pos_node, "z"), punit, data.constants);
                rel_pos = Vector3D(px, py, pz);
            }

            auto* posref_node = db.node->first_node("positionref");
            if(posref_node) {
                std::string ref = SafeAttrVal(posref_node, "ref");
                auto it = data.positions.find(ref);
                if(it != data.positions.end()) {
                    rel_pos = it->second;
                }
            }

            auto* rot_node = db.node->first_node("rotation");
            if(rot_node) {
                const char* runit = SafeAttrVal(rot_node, "unit");
                double rx = ParseAngle(SafeAttrVal(rot_node, "x"), runit, data.constants);
                double ry = ParseAngle(SafeAttrVal(rot_node, "y"), runit, data.constants);
                double rz = ParseAngle(SafeAttrVal(rot_node, "z"), runit, data.constants);
                rel_rot = QuatFromGDMLRotation(rx, ry, rz);
            }

            auto* rotref_node = db.node->first_node("rotationref");
            if(rotref_node) {
                std::string ref = SafeAttrVal(rotref_node, "ref");
                auto it = data.rotations.find(ref);
                if(it != data.rotations.end()) {
                    rel_rot = it->second;
                }
            }

            // Apply first-operand placement if specified
            bool has_first_placement = (first_pos.magnitude() > 0 || first_rot != Quaternion());
            std::shared_ptr<Geometry> left_placed;
            if(has_first_placement) {
                left_placed = left->create();
                Placement first_placement(first_pos, first_rot);
                left_placed->SetPlacement(first_placement);
            }

            auto right = right_base->create();
            Placement rel_placement(rel_pos, rel_rot);
            right->SetPlacement(rel_placement);

            auto left_final = has_first_placement
                ? std::const_pointer_cast<const Geometry>(left_placed)
                : std::const_pointer_cast<const Geometry>(left);

            auto geo = std::make_shared<BooleanGeometry>(op, left_final,
                std::const_pointer_cast<const Geometry>(right));

            storeSolid(db.name, geo);
            db.node = nullptr; // Mark as resolved
            made_progress = true;
        }
    }

    // Flatten multi-instance solids into data.solids with unique names.
    // For names with 1 instance: keep the original name.
    // For names with N > 1 instances: name (instance 0), name__2 (instance 1), etc.
    for(auto const & pair : solids_instances) {
        std::string const & name = pair.first;
        std::vector<std::shared_ptr<Geometry>> const & instances = pair.second;
        data.solid_instance_counts[name] = (int)instances.size();
        for(int i = 0; i < (int)instances.size(); ++i) {
            std::string flat_name;
            if(i == 0) {
                flat_name = name;
            } else {
                flat_name = name + "__" + std::to_string(i + 1);
            }
            data.solids[flat_name] = instances[i];
        }
    }
}


// Parse the <structure> section: volumes, assemblies, and physical volumes.
// Checks for ambiguous solid/material references (names with multiple instances).
static void ParseStructure(rapidxml::xml_node<>* structure_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!structure_node) return;

    // Shared physvol child parser used by both <volume> and <assembly>
    auto parsePhysVols = [&](rapidxml::xml_node<>* parent_node, GDMLVolume & vol) {
        for(auto* pv = parent_node->first_node("physvol"); pv; pv = pv->next_sibling("physvol")) {
            GDMLPhysVol physvol;
            physvol.position = Vector3D(0, 0, 0);
            physvol.rotation = Quaternion(); // identity

            auto* volref = pv->first_node("volumeref");
            if(volref) {
                physvol.volume_ref = SafeAttrVal(volref, "ref");
            }

            auto* pos = pv->first_node("position");
            if(pos) {
                const char* punit = SafeAttrVal(pos, "unit");
                double px = ParseLength(SafeAttrVal(pos, "x"), punit, data.constants);
                double py = ParseLength(SafeAttrVal(pos, "y"), punit, data.constants);
                double pz = ParseLength(SafeAttrVal(pos, "z"), punit, data.constants);
                physvol.position = Vector3D(px, py, pz);
            }
            auto* posref = pv->first_node("positionref");
            if(posref) {
                std::string ref = SafeAttrVal(posref, "ref");
                auto it = data.positions.find(ref);
                if(it != data.positions.end()) {
                    physvol.position = it->second;
                }
            }

            auto* rot = pv->first_node("rotation");
            if(rot) {
                const char* runit = SafeAttrVal(rot, "unit");
                double rx = ParseAngle(SafeAttrVal(rot, "x"), runit, data.constants);
                double ry = ParseAngle(SafeAttrVal(rot, "y"), runit, data.constants);
                double rz = ParseAngle(SafeAttrVal(rot, "z"), runit, data.constants);
                physvol.rotation = QuatFromGDMLRotation(rx, ry, rz);
            }
            auto* rotref = pv->first_node("rotationref");
            if(rotref) {
                std::string ref = SafeAttrVal(rotref, "ref");
                auto it = data.rotations.find(ref);
                if(it != data.rotations.end()) {
                    physvol.rotation = it->second;
                }
            }

            vol.children.push_back(physvol);
        }
    };

    for(auto* vol_node = structure_node->first_node(); vol_node;
        vol_node = vol_node->next_sibling()) {

        std::string tag(vol_node->name());
        if(tag != "volume" && tag != "assembly") continue;

        GDMLVolume vol;
        vol.name = SafeAttrVal(vol_node, "name");
        vol.is_assembly = (tag == "assembly");

        if(!vol.is_assembly) {
            // Get <materialref>
            auto* matref = vol_node->first_node("materialref");
            if(matref) {
                vol.material_ref = SafeAttrVal(matref, "ref");
                if(!vol.material_ref.empty()) {
                    auto it = data.material_instance_counts.find(vol.material_ref);
                    if(it != data.material_instance_counts.end() && it->second > 1) {
                        throw std::runtime_error(
                            "GDML error: volume '" + vol.name
                            + "' references ambiguous material '" + vol.material_ref
                            + "' which has " + std::to_string(it->second) + " instances");
                    }
                }
            }

            // Get <solidref>
            auto* solidref = vol_node->first_node("solidref");
            if(solidref) {
                vol.solid_ref = SafeAttrVal(solidref, "ref");
                if(!vol.solid_ref.empty()) {
                    auto it = data.solid_instance_counts.find(vol.solid_ref);
                    if(it != data.solid_instance_counts.end() && it->second > 1) {
                        throw std::runtime_error(
                            "GDML error: volume '" + vol.name
                            + "' references ambiguous solid '" + vol.solid_ref
                            + "' which has " + std::to_string(it->second) + " instances");
                    }
                }
            }
        }

        parsePhysVols(vol_node, vol);

        for(auto* rv = vol_node->first_node("replicavol"); rv; rv = rv->next_sibling("replicavol")) {
            int ncopies = 0;
            {
                std::string ns(SafeAttrVal(rv, "number"));
                if(!ns.empty()) {
                    try { ncopies = std::stoi(ns); } catch(...) { ncopies = 0; }
                }
            }
            if(ncopies <= 0) continue;

            std::string child_vol_ref;
            auto* vref = rv->first_node("volumeref");
            if(vref) child_vol_ref = SafeAttrVal(vref, "ref");
            if(child_vol_ref.empty()) continue;

            auto* raa = rv->first_node("replicate_along_axis");
            if(!raa) continue;

            double dir_x = 0, dir_y = 0, dir_z = 0;
            auto* dir = raa->first_node("direction");
            if(dir) {
                dir_x = SafeParseDouble(SafeAttrVal(dir, "x"), data.constants);
                dir_y = SafeParseDouble(SafeAttrVal(dir, "y"), data.constants);
                dir_z = SafeParseDouble(SafeAttrVal(dir, "z"), data.constants);
            }
            double dir_len = std::sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z);
            if(dir_len == 0) {
                auto* rho_attr = dir ? dir->first_attribute("rho") : nullptr;
                auto* phi_attr = dir ? dir->first_attribute("phi") : nullptr;
                if(rho_attr || phi_attr) {
                    EmitWarning(data, options, "replicavol in volume '" + vol.name + "': cylindrical axis replication not supported; skipping");
                } else {
                    EmitWarning(data, options, "replicavol in volume '" + vol.name + "': no direction specified; skipping");
                }
                continue;
            }
            dir_x /= dir_len;
            dir_y /= dir_len;
            dir_z /= dir_len;

            double width = 0;
            auto* w_node = raa->first_node("width");
            if(w_node) {
                width = ParseLength(SafeAttrVal(w_node, "value"),
                                    SafeAttrVal(w_node, "unit"), data.constants);
            }

            double offset = 0;
            auto* o_node = raa->first_node("offset");
            if(o_node) {
                offset = ParseLength(SafeAttrVal(o_node, "value"),
                                     SafeAttrVal(o_node, "unit"), data.constants);
            }

            for(int i = 0; i < ncopies; ++i) {
                double t = offset + (i + 0.5 - ncopies / 2.0) * width;
                GDMLPhysVol pv;
                pv.volume_ref = child_vol_ref;
                pv.position = Vector3D(t * dir_x, t * dir_y, t * dir_z);
                pv.rotation = Quaternion();
                vol.children.push_back(pv);
            }
        }

        if(!vol.name.empty()) {
            if(data.volumes.find(vol.name) != data.volumes.end()) {
                EmitWarning(data, options, "duplicate volume name '" + vol.name + "', overwriting");
            }
            data.volumes[vol.name] = vol;
        }
    }
}


// Parse the <setup> section to find the world volume
static void ParseSetup(rapidxml::xml_node<>* setup_node, GDMLData & data, GDMLParseOptions const & options) {
    (void)options;
    if(!setup_node) return;

    auto* world_node = setup_node->first_node("world");
    if(world_node) {
        data.world_volume = SafeAttrVal(world_node, "ref");
    }
}


// Main parser function
GDMLData ParseGDML(std::string const & filename, GDMLParseOptions const & options) {
    // Read the file into a mutable buffer (rapidxml modifies in-place)
    std::ifstream file(filename.c_str(), std::ios::binary);
    if(!file.is_open()) {
        throw std::runtime_error("Cannot open GDML file: " + filename);
    }

    // Read entire file into string
    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    file.close();

    // Determine base directory for entity resolution
    std::string base_dir;
    {
        size_t last_slash = filename.rfind('/');
        if(last_slash != std::string::npos) {
            base_dir = filename.substr(0, last_slash);
        } else {
            base_dir = ".";
        }
    }

    // Expand ENTITY references before XML parsing
    if(content.find("<!ENTITY") != std::string::npos || content.find("<!entity") != std::string::npos) {
        content = ExpandEntities(content, base_dir);
    }

    // Handle xi:include (XInclude) - expand included files inline
    if(content.find("xi:include") != std::string::npos) {
        // Simple xi:include expansion: replace <xi:include href="file"/> with file content
        std::string result;
        size_t pos = 0;
        while(pos < content.size()) {
            size_t inc_start = content.find("<xi:include", pos);
            if(inc_start == std::string::npos) {
                result.append(content, pos, content.size() - pos);
                break;
            }
            result.append(content, pos, inc_start - pos);
            size_t inc_end = content.find("/>", inc_start);
            if(inc_end == std::string::npos) {
                inc_end = content.find("</xi:include>", inc_start);
                if(inc_end != std::string::npos) inc_end += 12;
            } else {
                inc_end += 2;
            }
            if(inc_end == std::string::npos) {
                result.append(content, inc_start, content.size() - inc_start);
                break;
            }

            std::string inc_tag = content.substr(inc_start, inc_end - inc_start);
            size_t href_pos = inc_tag.find("href=");
            if(href_pos != std::string::npos) {
                href_pos += 5;
                char quote = inc_tag[href_pos];
                if(quote == '"' || quote == '\'') {
                    ++href_pos;
                    size_t href_end = inc_tag.find(quote, href_pos);
                    if(href_end != std::string::npos) {
                        std::string href = inc_tag.substr(href_pos, href_end - href_pos);
                        std::string inc_path = (href[0] == '/') ? href : base_dir + "/" + href;
                        std::ifstream inc_file(inc_path.c_str(), std::ios::binary);
                        if(inc_file.is_open()) {
                            std::string inc_content((std::istreambuf_iterator<char>(inc_file)),
                                                    std::istreambuf_iterator<char>());
                            inc_file.close();
                            size_t inc_last_slash = inc_path.rfind('/');
                            std::string inc_dir = (inc_last_slash != std::string::npos) ?
                                inc_path.substr(0, inc_last_slash) : base_dir;
                            inc_content = ExpandEntities(inc_content, inc_dir);
                            result.append(inc_content);
                        }
                    }
                }
            }
            pos = inc_end;
        }
        content = result;
    }

    // rapidxml requires a mutable, null-terminated buffer
    std::vector<char> buffer(content.begin(), content.end());
    buffer.push_back('\0');

    // Parse the XML
    rapidxml::xml_document<> doc;
    try {
        doc.parse<0>(buffer.data());
    } catch(rapidxml::parse_error & e) {
        throw std::runtime_error(std::string("GDML XML parse error: ") + e.what());
    }

    GDMLData data;

    SeedBuiltInConstants(data);

    // Find root node: accept <gdml> or any root element
    auto* gdml_node = doc.first_node("gdml");
    if(!gdml_node) {
        gdml_node = doc.first_node();
        if(!gdml_node) {
            throw std::runtime_error("GDML error: empty document: " + filename);
        }
        std::string root_tag(gdml_node->name());
        EmitWarning(data, options, "non-standard root element <" + root_tag + "> in " + filename + ", treating as <gdml>");
    }

    // Resolve all <define> sections in document order (handles loops inline).
    ResolveDefineInOrder(gdml_node, data, options);

    // Expand <loop> elements in non-define sections (structure, solids) where
    // downstream parsers need to see the expanded DOM nodes.
    ExpandLoops(doc, gdml_node, data.constants);
    // Remove loops from <define> sections (already resolved, but the DOM nodes
    // are still present; remove them so ParseDefine is not confused if called again)
    for(auto* def = gdml_node->first_node("define"); def; def = def->next_sibling("define")) {
        auto* child = def->first_node();
        while(child) {
            auto* next = child->next_sibling();
            if(std::string(child->name()) == "loop") {
                def->remove_node(child);
            }
            child = next;
        }
    }

    // Parse all <materials> sections (with instance-scoped name resolution)
    ParseAllMaterials(gdml_node, data, options);
    // Parse all <solids> sections (with instance-scoped boolean resolution)
    ParseAllSolids(gdml_node, data, options);
    // Parse all <structure> sections
    for(auto* node = gdml_node->first_node("structure"); node; node = node->next_sibling("structure")) {
        ParseStructure(node, data, options);
    }
    // Parse all <setup> sections (last wins)
    for(auto* node = gdml_node->first_node("setup"); node; node = node->next_sibling("setup")) {
        ParseSetup(node, data, options);
    }

    // Warn if no <setup> section was found (world volume undefined)
    if(data.world_volume.empty()) {
        EmitWarning(data, options, "No <setup> section found in GDML; world volume is undefined");
    }

    return data;
}

} // namespace detector
} // namespace siren
