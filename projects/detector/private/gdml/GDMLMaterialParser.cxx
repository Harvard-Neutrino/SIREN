#include "GDMLParserPrivate.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace siren {
namespace detector {
namespace gdml {
// Parse all <materials> sections with instance-scoped name resolution.
// Each name can have multiple definitions (instances). Composition references
// resolve to the current instance of the referenced name at definition time.
// After parsing, multi-instance names are flattened to unique names.
void ParseAllMaterials(rapidxml::xml_node<>* root_node, GDMLData & data, GDMLParseOptions const & options) {
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

                // Composite children: atom counts -> mass fractions.
                // Collect (resolved_name, atom_count) pairs, then convert
                // to mass fractions using the referenced element's atomic mass:
                //   mass_fraction_i = n_i * A_i / sum(n_j * A_j)
                std::vector<std::pair<std::string, double>> composite_entries;
                for(auto* comp = node->first_node("composite"); comp; comp = comp->next_sibling("composite")) {
                    std::string ref = SafeAttrVal(comp, "ref");
                    double n = SafeParseDouble(SafeAttrVal(comp, "n"), data.constants);
                    if(!ref.empty()) {
                        composite_entries.push_back({resolveCompRef(ref), n});
                    }
                }
                if(!composite_entries.empty()) {
                    auto lookupA = [&](std::string const & resolved_name) -> double {
                        auto iso_it = isotope_stacks.instances.find(resolved_name);
                        if(iso_it != isotope_stacks.instances.end() && !iso_it->second.empty())
                            return iso_it->second.back().A;
                        auto elem_it = element_stacks.instances.find(resolved_name);
                        if(elem_it != element_stacks.instances.end() && !elem_it->second.empty())
                            return elem_it->second.back().A;
                        auto mat_it = material_stacks.instances.find(resolved_name);
                        if(mat_it != material_stacks.instances.end() && !mat_it->second.empty())
                            return mat_it->second.back().A;
                        return 0.0;
                    };
                    double total_mass = 0.0;
                    std::vector<double> masses(composite_entries.size());
                    for(size_t ci = 0; ci < composite_entries.size(); ++ci) {
                        double A = lookupA(composite_entries[ci].first);
                        masses[ci] = composite_entries[ci].second * A;
                        total_mass += masses[ci];
                    }
                    if(total_mass > 0) {
                        for(size_t ci = 0; ci < composite_entries.size(); ++ci) {
                            mat.composition[composite_entries[ci].first] = masses[ci] / total_mass;
                        }
                    } else {
                        for(auto const & ce : composite_entries) {
                            mat.composition[ce.first] = ce.second;
                        }
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

} // namespace gdml
} // namespace detector
} // namespace siren
