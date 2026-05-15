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

#include <cereal/external/rapidxml/rapidxml.hpp>

namespace rapidxml = cereal::rapidxml;

#include "SIREN/math/Vector3D.h"
#include "SIREN/math/Quaternion.h"

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

using namespace siren::math;
using namespace siren::geometry;

namespace siren {
namespace detector {

namespace {

static const double PI = 3.141592653589793238462643383279502884197;

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
    return 0;
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

static std::vector<Token> Tokenize(std::string const & s, std::map<std::string, double> const & constants, std::string const & expr) {
    std::vector<Token> tokens;
    size_t i = 0;
    while(i < s.size()) {
        char c = s[i];
        if(c == ' ' || c == '\t' || c == '\n' || c == '\r') { ++i; continue; }

        if(c == '(') { tokens.push_back({TokenKind::LPAREN}); ++i; continue; }
        if(c == ')') { tokens.push_back({TokenKind::RPAREN}); ++i; continue; }
        if(c == ',') { tokens.push_back({TokenKind::COMMA}); ++i; continue; }

        if(c == '+' || c == '-' || c == '*' || c == '/') {
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

double EvalExpression(std::string const & expr, std::map<std::string, double> const & constants) {
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
            return 0.0;
        }
    }

    // Fast path: bare constant name
    {
        bool bare = true;
        for(size_t i = 0; i < s.size(); ++i) {
            if(!IsIdentChar(s[i])) { bare = false; break; }
        }
        if(bare && IsIdentStart(s[0])) {
            auto it = constants.find(s);
            if(it != constants.end()) return it->second;
            throw std::runtime_error("GDML expression error: unknown constant '" + s + "'");
        }
    }

    std::vector<Token> tokens = Tokenize(s, constants, expr);

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
            while(!op_stack.empty() && op_stack.back().kind == TokenKind::OP &&
                  OpPrecedence(op_stack.back().op) >= OpPrecedence(tok.op)) {
                output.push_back(op_stack.back());
                op_stack.pop_back();
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
double SafeParseDouble(const char* val, std::map<std::string, double> const & constants = {}) {
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
        return 0.0;
    }
    // Fall back to expression evaluation
    return EvalExpression(std::string(val), constants);
}

// Convert a GDML length value+unit to meters (SIREN base unit)
// GDML default length unit is mm
double ParseLength(const char* value, const char* unit,
                   std::map<std::string, double> const & constants = {}) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);

    if(!unit || unit[0] == '\0') {
        // GDML default is mm
        return val * 0.001;
    }

    std::string u(unit);
    if(u == "mm") return val * 0.001;
    if(u == "cm") return val * 0.01;
    if(u == "m")  return val * 1.0;

    // Unknown unit
    throw std::runtime_error("GDML error: unrecognized length unit '" + std::string(unit) + "'");
}

// Get the length scale factor for a given unit string
double LengthScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        return 0.001; // GDML default is mm
    }
    std::string u(unit);
    if(u == "mm") return 0.001;
    if(u == "cm") return 0.01;
    if(u == "m")  return 1.0;
    throw std::runtime_error("GDML error: unrecognized length unit '" + u + "'");
}

// Convert a GDML angle value+unit to radians
// GDML default angle unit is degrees
double ParseAngle(const char* value, const char* unit,
                  std::map<std::string, double> const & constants = {}) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);

    if(!unit || unit[0] == '\0') {
        // GDML default angle unit is degrees
        return val * PI / 180.0;
    }

    std::string u(unit);
    if(u == "deg") return val * PI / 180.0;
    if(u == "rad") return val;

    // Unknown unit
    throw std::runtime_error("GDML error: unrecognized angle unit '" + std::string(unit) + "'");
}

// Get the angle scale factor for a given unit string
double AngleScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        return PI / 180.0; // GDML default angle unit is degrees
    }
    std::string u(unit);
    if(u == "deg") return PI / 180.0;
    if(u == "rad") return 1.0;
    throw std::runtime_error("GDML error: unrecognized angle unit '" + u + "'");
}

// Build a quaternion from GDML rotation convention:
// Extrinsic rotations: rotate around X by rx, then Y by ry, then Z by rz
// Equivalent to intrinsic Z-Y-X, so quaternion = Qz * Qy * Qx
Quaternion QuatFromGDMLRotation(double rx, double ry, double rz) {
    // Individual axis quaternions: Quaternion(x, y, z, w)
    Quaternion qx(std::sin(rx / 2.0), 0, 0, std::cos(rx / 2.0));
    Quaternion qy(0, std::sin(ry / 2.0), 0, std::cos(ry / 2.0));
    Quaternion qz(0, 0, std::sin(rz / 2.0), std::cos(rz / 2.0));
    // Applied right-to-left: first X, then Y, then Z
    return qz * qy * qx;
}

// Record a warning. In strict mode, throws instead of continuing.
void EmitWarning(GDMLData & data, GDMLParseOptions const & options, std::string const & msg) {
    if(options.strict) {
        throw std::runtime_error("GDML error (strict mode): " + msg);
    }
    data.warnings.push_back(msg);
    std::cerr << "GDML warning: " << msg << std::endl;
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
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants);
            if(!name.empty()) {
                if(data.constants.find(name) != data.constants.end()) {
                    EmitWarning(data, options, "duplicate constant name '" + name + "', overwriting");
                }
                data.constants[name] = value;
            }
        }
        else if(tag == "quantity") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants);
            if(!name.empty()) {
                if(data.constants.find(name) != data.constants.end()) {
                    EmitWarning(data, options, "duplicate constant name '" + name + "', overwriting");
                }
                data.constants[name] = value;
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
            // Scale definitions are not directly used; skip
        }
        else if(tag == "variable") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants);
            if(!name.empty()) {
                if(data.constants.find(name) != data.constants.end()) {
                    EmitWarning(data, options, "duplicate constant name '" + name + "', overwriting");
                }
                data.constants[name] = value;
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
                    mat.density = SafeParseDouble(SafeAttrVal(d_node, "value"), data.constants);
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
            double deltaphi = 2.0 * PI;
            double starttheta = 0.0;
            double deltatheta = PI;
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

            // Check angular extent for partial tube
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') {
                double startphi = SafeParseDouble(sp_val, data.constants) * ascale;
                if(std::fabs(startphi) > ANG_TOL) {
                    EmitWarning(data, options, std::string(tag) + " '" + name + "' has partial angular extent (startphi=" + std::to_string(startphi) + "); SIREN creates full rotation");
                }
            }
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                double deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
                if(std::fabs(deltaphi - 2.0 * PI) > ANG_TOL) {
                    EmitWarning(data, options, std::string(tag) + " '" + name + "' has partial angular extent (deltaphi=" + std::to_string(deltaphi) + "); SIREN creates full rotation");
                }
            }

            geo = Cylinder(rmax, rmin, hz * 2.0).create();
        }
        else if(tag == "cone") {
            double rmin1 = SafeParseDouble(SafeAttrVal(node, "rmin1"), data.constants) * lscale;
            double rmax1 = SafeParseDouble(SafeAttrVal(node, "rmax1"), data.constants) * lscale;
            double rmin2 = SafeParseDouble(SafeAttrVal(node, "rmin2"), data.constants) * lscale;
            double rmax2 = SafeParseDouble(SafeAttrVal(node, "rmax2"), data.constants) * lscale;
            // GDML <cone> z is HALF-height; SIREN Cone takes full height
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;

            // Check angular extent for partial cone
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') {
                double startphi = SafeParseDouble(sp_val, data.constants) * ascale;
                if(std::fabs(startphi) > ANG_TOL) {
                    EmitWarning(data, options, "cone '" + name + "' has partial angular extent (startphi=" + std::to_string(startphi) + "); SIREN creates full rotation");
                }
            }
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                double deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
                if(std::fabs(deltaphi - 2.0 * PI) > ANG_TOL) {
                    EmitWarning(data, options, "cone '" + name + "' has partial angular extent (deltaphi=" + std::to_string(deltaphi) + "); SIREN creates full rotation");
                }
            }

            geo = Cone(rmin1, rmax1, rmin2, rmax2, hz * 2.0).create();
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
            // Check angular extent for partial polycone
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') {
                double startphi = SafeParseDouble(sp_val, data.constants) * ascale;
                if(std::fabs(startphi) > ANG_TOL) {
                    EmitWarning(data, options, "polycone '" + name + "' has partial angular extent (startphi=" + std::to_string(startphi) + "); SIREN creates full rotation");
                }
            }
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                double deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
                if(std::fabs(deltaphi - 2.0 * PI) > ANG_TOL) {
                    EmitWarning(data, options, "polycone '" + name + "' has partial angular extent (deltaphi=" + std::to_string(deltaphi) + "); SIREN creates full rotation");
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

            if(z_planes.size() >= 2) {
                geo = Polycone(z_planes, rmin_vec, rmax_vec).create();
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
                if(std::fabs(deltaphi - 2.0 * PI) > ANG_TOL) {
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
                geo = Polyhedra(numSide, startphi, z_planes, rmin_vec, rmax_vec).create();
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
            double deltaphi = 2.0 * PI;
            const char* sp_val = SafeAttrVal(node, "startphi");
            if(sp_val[0] != '\0') startphi = SafeParseDouble(sp_val, data.constants) * ascale;
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;

            if(rtor > 0 && rmax_t > 0) {
                geo = Torus(rtor, rmax_t, rmin_t, startphi, deltaphi).create();
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

        auto right = right_base->create();
        Placement rel_placement(rel_pos, rel_rot);
        right->SetPlacement(rel_placement);

        return std::make_shared<BooleanGeometry>(op,
            std::const_pointer_cast<const Geometry>(left),
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

            auto right = right_base->create();
            Placement rel_placement(rel_pos, rel_rot);
            right->SetPlacement(rel_placement);

            auto geo = std::make_shared<BooleanGeometry>(op,
                std::const_pointer_cast<const Geometry>(left),
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


// Parse the <structure> section: volumes and physical volumes.
// Checks for ambiguous solid/material references (names with multiple instances).
static void ParseStructure(rapidxml::xml_node<>* structure_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!structure_node) return;

    for(auto* vol_node = structure_node->first_node("volume"); vol_node;
        vol_node = vol_node->next_sibling("volume")) {

        GDMLVolume vol;
        vol.name = SafeAttrVal(vol_node, "name");

        // Get <materialref>
        auto* matref = vol_node->first_node("materialref");
        if(matref) {
            vol.material_ref = SafeAttrVal(matref, "ref");
            // Check for ambiguous material reference.
            // material_instance_counts tracks the <material> stack only
            // (not isotopes or elements), so element/material name collisions
            // don't trigger this — they are in separate namespaces.
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
            // Check for ambiguous solid reference
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

        // Parse <physvol> children
        for(auto* pv = vol_node->first_node("physvol"); pv; pv = pv->next_sibling("physvol")) {
            GDMLPhysVol physvol;
            physvol.position = Vector3D(0, 0, 0);
            physvol.rotation = Quaternion(); // identity

            // Volume reference
            auto* volref = pv->first_node("volumeref");
            if(volref) {
                physvol.volume_ref = SafeAttrVal(volref, "ref");
            }

            // Position: inline or reference
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

            // Rotation: inline or reference
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

    // Parse all <define> sections
    for(auto* node = gdml_node->first_node("define"); node; node = node->next_sibling("define")) {
        ParseDefine(node, data, options);
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

    return data;
}

} // namespace detector
} // namespace siren
