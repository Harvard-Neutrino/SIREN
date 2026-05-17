#include "GDMLParserPrivate.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace siren::math;

namespace siren {
namespace detector {
namespace gdml {
// Resolve all <define> sections in document order with inline loop handling.
// Uses an explicit stack to iterate loop bodies without recursion or DOM
// modification. Constants/variables are evaluated immediately when encountered,
// giving document-order semantics (a constant before a loop sees the pre-loop
// variable value; a constant after sees the post-loop value).
void ResolveDefineInOrder(rapidxml::xml_node<>* gdml_node,
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

} // namespace gdml
} // namespace detector
} // namespace siren
