#include "GDMLParserPrivate.h"

#include <cmath>
#include <cctype>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace siren {
namespace detector {
namespace gdml {
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
    if(entities.empty()) return content;
    std::string result;
    result.reserve(content.size());
    size_t pos = 0;
    while(pos < content.size()) {
        size_t amp = content.find('&', pos);
        if(amp == std::string::npos) {
            result.append(content, pos, content.size() - pos);
            break;
        }
        result.append(content, pos, amp - pos);
        size_t semi = content.find(';', amp + 1);
        if(semi == std::string::npos) {
            result.append(content, amp, content.size() - amp);
            break;
        }
        std::string ent_name = content.substr(amp + 1, semi - amp - 1);
        auto it = entities.find(ent_name);
        if(it != entities.end()) {
            result.append(it->second);
        } else {
            result.append(content, amp, semi - amp + 1);
        }
        pos = semi + 1;
    }
    return result;
}

// Preprocess XML content to expand ENTITY references.
// Reads the DOCTYPE for ENTITY SYSTEM declarations, loads the referenced files,
// and substitutes all &name; references with file contents.
// Uses an iterative stack to handle nested entity files without recursion.
// Each stack frame owns its own scope: entities declared in a file's DOCTYPE
// are only visible for substitution within that file's content.
std::string ExpandEntities(std::string const & content, std::string const & base_dir) {
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
            // All declarations resolved (or depth limit reached) - substitute and pop
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


static size_t FindXIncludeEnd(std::string const & content, size_t inc_start) {
    size_t self_close = content.find("/>", inc_start);
    size_t explicit_close = content.find("</xi:include>", inc_start);

    if(self_close == std::string::npos && explicit_close == std::string::npos) {
        return std::string::npos;
    }
    if(self_close != std::string::npos
       && (explicit_close == std::string::npos || self_close < explicit_close)) {
        return self_close + 2;
    }
    return explicit_close + std::string("</xi:include>").size();
}

static std::string ExtractQuotedAttribute(std::string const & tag, std::string const & attr) {
    size_t pos = 0;
    while((pos = tag.find(attr, pos)) != std::string::npos) {
        bool left_ok = (pos == 0 || !IsIdentChar(tag[pos - 1]));
        size_t after = pos + attr.size();
        bool right_ok = (after < tag.size() && tag[after] == '=');
        if(left_ok && right_ok) break;
        pos = after;
    }
    if(pos == std::string::npos) return "";

    pos += attr.size() + 1;
    if(pos >= tag.size()) return "";

    char quote = tag[pos];
    if(quote != '"' && quote != '\'') return "";

    ++pos;
    size_t end = tag.find(quote, pos);
    if(end == std::string::npos) return "";
    return tag.substr(pos, end - pos);
}

static std::string DirectoryName(std::string const & path, std::string const & fallback) {
    size_t slash = path.rfind('/');
    if(slash == std::string::npos) return fallback;
    return path.substr(0, slash);
}

static std::string ResolvePath(std::string const & href, std::string const & base_dir) {
    if(href.empty()) return "";
    if(href[0] == '/') return href;
    return base_dir + "/" + href;
}

static bool ContainsPath(std::vector<std::string> const & stack, std::string const & path) {
    for(auto const & entry : stack) {
        if(entry == path) return true;
    }
    return false;
}

static std::string ExpandXIncludesImpl(std::string const & content,
                                        std::string const & base_dir,
                                        std::vector<std::string> & include_stack,
                                        int depth) {
    static constexpr int MAX_XINCLUDE_DEPTH = 16;

    std::string result;
    size_t pos = 0;
    while(pos < content.size()) {
        size_t inc_start = content.find("<xi:include", pos);
        if(inc_start == std::string::npos) {
            result.append(content, pos, content.size() - pos);
            break;
        }
        result.append(content, pos, inc_start - pos);
        size_t inc_end = FindXIncludeEnd(content, inc_start);
        if(inc_end == std::string::npos) {
            result.append(content, inc_start, content.size() - inc_start);
            break;
        }

        std::string inc_tag = content.substr(inc_start, inc_end - inc_start);
        std::string href = ExtractQuotedAttribute(inc_tag, "href");
        std::string inc_path = ResolvePath(href, base_dir);
        if(!inc_path.empty() && depth < MAX_XINCLUDE_DEPTH && !ContainsPath(include_stack, inc_path)) {
            std::ifstream inc_file(inc_path.c_str(), std::ios::binary);
            if(inc_file.is_open()) {
                std::string inc_content((std::istreambuf_iterator<char>(inc_file)),
                                        std::istreambuf_iterator<char>());
                inc_file.close();

                std::string inc_dir = DirectoryName(inc_path, base_dir);
                inc_content = ExpandEntities(inc_content, inc_dir);
                if(inc_content.find("xi:include") != std::string::npos) {
                    include_stack.push_back(inc_path);
                    inc_content = ExpandXIncludesImpl(inc_content, inc_dir, include_stack, depth + 1);
                    include_stack.pop_back();
                }
                result.append(inc_content);
            }
        }
        pos = inc_end;
    }
    return result;
}

// Expand simple xi:include tags by replacing them with the referenced file content.
// Included files are ENTITY-expanded relative to their own directory before insertion.
std::string ExpandXIncludes(std::string const & content, std::string const & base_dir) {
    std::vector<std::string> include_stack;
    return ExpandXIncludesImpl(content, base_dir, include_stack, 0);
}

// Record a warning. In strict mode, throws instead of continuing.
void EmitWarning(GDMLData & data, GDMLParseOptions const & options, std::string const & msg) {
    if(options.strict) {
        throw std::runtime_error("GDML error (strict mode): " + msg);
    }
    data.warnings.push_back(msg);
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
    std::unordered_map<std::string, double> const & constants,
    GDMLData & data) {
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
        } catch(std::exception const & e) {
            data.errors.push_back("bracket evaluation failed for '[" + expr + "]': " + e.what());
            pos = end + 1;
        } catch(...) {
            data.errors.push_back("bracket evaluation failed for '[" + expr + "]'");
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
    std::unordered_map<std::string, double> const & constants,
    GDMLData & data) {

    auto* node = doc.allocate_node(src->type());
    node->name(src->name(), src->name_size());

    for(auto* attr = src->first_attribute(); attr; attr = attr->next_attribute()) {
        std::string val(attr->value(), attr->value_size());
        val = SubstLoopVar(val, var_name, var_value_str, constants, data);
        char* aname = doc.allocate_string(attr->name(), attr->name_size() + 1);
        char* aval = doc.allocate_string(val.c_str(), val.size() + 1);
        node->append_attribute(doc.allocate_attribute(aname, aval));
    }

    for(auto* child = src->first_node(); child; child = child->next_sibling()) {
        node->append_node(CloneNodeWithSubst(doc, child, var_name, var_value_str, constants, data));
    }
    return node;
}

// Expand all <loop> elements in a DOM subtree (iterative, stack-based).
// Clones loop body nodes for each iteration with text substitution, inserts
// them into the parent, and removes the loop node. Used for non-define
// sections (structure, solids) where the DOM must be modified for downstream
// parsers that iterate over it.
void ExpandLoops(rapidxml::xml_document<> & doc,
                        rapidxml::xml_node<>* root,
                        std::unordered_map<std::string, double> & constants,
                        GDMLData & data) {
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

        if(std::strcmp(child->name(), "loop") != 0) {
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
            if(++count > 10000) {
                throw std::runtime_error("GDML loop exceeded 10000 iterations");
            }
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
                auto* cloned = CloneNodeWithSubst(doc, body, var_name, v_str, constants, data);
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

} // namespace gdml
} // namespace detector
} // namespace siren
