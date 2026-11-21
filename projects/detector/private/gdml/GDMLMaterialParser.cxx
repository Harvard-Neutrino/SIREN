#include "GDMLParserPrivate.h"

#include <cctype>
#include <map>
#include <stdexcept>
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

    // Instance trackers for the three GDML material namespaces.
    auto mat_eq = [](GDMLMaterial const & a, GDMLMaterial const & b) { return a == b; };
    InstanceTracker<GDMLMaterial> isotope_tracker(data, options.strip_pointer_suffixes);
    InstanceTracker<GDMLMaterial> element_tracker(data, options.strip_pointer_suffixes);
    InstanceTracker<GDMLMaterial> material_tracker(data, options.strip_pointer_suffixes);

    auto storeInStack = [&](InstanceTracker<GDMLMaterial> & tracker, GDMLMaterial const & mat) {
        GDMLMaterial stored = mat;
        std::string stripped = StripPointerSuffix(mat.name, options.strip_pointer_suffixes);
        stored.name = stripped;
        tracker.Store(mat.name, stored, mat_eq);
    };

    // Resolve a composition reference across all three namespaces.
    // Isotopes first (most specific), then elements, then materials.
    auto resolveCompRef = [&](std::string const & ref_name) -> std::string {
        std::string resolved = isotope_tracker.Resolve(ref_name);
        if(resolved != StripPointerSuffix(ref_name, options.strip_pointer_suffixes)) return resolved;
        resolved = element_tracker.Resolve(ref_name);
        if(resolved != StripPointerSuffix(ref_name, options.strip_pointer_suffixes)) return resolved;
        return material_tracker.Resolve(ref_name);
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
                    storeInStack(isotope_tracker, mat);
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
                    storeInStack(element_tracker, mat);
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
                    std::string d_val_trimmed = TrimWhitespace(d_val_str);
                    auto qty_it = data.quantity_types.find(d_val_trimmed);
                    bool is_density_qty = (qty_it != data.quantity_types.end()
                                           && qty_it->second == GDMLQuantityType::DENSITY);

                    // Guard against inline-unit density expressions
                    // (e.g. value="1*gram"). SafeParseDouble would resolve
                    // 'gram' to the CLHEP-internal scale 6.24e21, and no
                    // subsequent unit-attribute conversion can undo it,
                    // so the stored density would be wrong by 20+ orders
                    // of magnitude. Typed <quantity type="density"> refs
                    // route through QuantityUnitScale and are safe.
                    if(!is_density_qty) {
                        std::string found_unit;
                        if(ExpressionContainsUnitBuiltIn(d_val_trimmed, found_unit)) {
                            throw std::runtime_error(
                                "Material '" + mat.name + "' has <D value=\"" +
                                d_val_trimmed + "\"> which embeds the unit "
                                "identifier '" + found_unit + "'. Inline-unit "
                                "expressions are not supported in density values "
                                "because SafeParseDouble would resolve the unit "
                                "to its CLHEP-internal scale. Use "
                                "<D value=\"<number>\" unit=\"g/cm3\"/> (or a "
                                "<quantity type=\"density\"> definition) instead.");
                        }
                    }

                    mat.density = SafeParseDouble(d_val_str, data.constants);

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
                if(mat.density < 0.0) {
                    throw std::runtime_error("Material '" + mat.name + "' has negative density: " + std::to_string(mat.density));
                }
                if(mat.density == 0.0) {
                    EmitWarning(data, options, "material '" + mat.name + "' has zero density, treating as vacuum (1e-25 g/cm3)");
                    mat.density = 1e-25;
                }

                auto* atom_node = node->first_node("atom");
                if(atom_node) {
                    mat.A = SafeParseDouble(SafeAttrVal(atom_node, "value"), data.constants);
                    // Z can appear on the <material> element or on <atom>.
                    // Some GDML generators put it on <atom Z="..."/>.
                    if(mat.Z == 0) {
                        mat.Z = SafeParseDouble(SafeAttrVal(atom_node, "Z"), data.constants);
                    }
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
                    // Extract original name from a flattened name like "Iron__2".
                    // The flattened format is "OriginalName__N" (double underscore + instance number).
                    // Instance 0 keeps the original name unchanged.
                    // Parse a flattened name "Base__N" into (base, index).
                    // FlatName uses 1-based suffixes, so "Iron__2" → index 1.
                    // Returns (name, 0) for bare names without a suffix.
                    auto parseFlatIdx = [](std::string const & flat) -> std::pair<std::string, int> {
                        size_t pos = flat.rfind("__");
                        if(pos == std::string::npos || pos == 0 || pos + 2 >= flat.size())
                            return {flat, 0};
                        for(size_t i = pos + 2; i < flat.size(); ++i)
                            if(!std::isdigit(static_cast<unsigned char>(flat[i])))
                                return {flat, 0};
                        int idx = std::stoi(flat.substr(pos + 2)) - 1;
                        return {flat.substr(0, pos), (idx >= 0) ? idx : 0};
                    };
                    auto lookupA = [&](std::string const & resolved_name) -> double {
                        auto [base, idx] = parseFlatIdx(resolved_name);
                        auto const & iso_inst = isotope_tracker.Instances();
                        auto iso_it = iso_inst.find(base);
                        if(iso_it != iso_inst.end() && idx < (int)iso_it->second.size())
                            return iso_it->second[idx].A;
                        auto const & elem_inst = element_tracker.Instances();
                        auto elem_it = elem_inst.find(base);
                        if(elem_it != elem_inst.end() && idx < (int)elem_it->second.size())
                            return elem_it->second[idx].A;
                        auto const & mat_inst = material_tracker.Instances();
                        auto mat_it = mat_inst.find(base);
                        if(mat_it != mat_inst.end() && idx < (int)mat_it->second.size())
                            return mat_it->second[idx].A;
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
                    storeInStack(material_tracker, mat);
                }
            }
        }
    }

    // Flatten all three trackers into GDMLData maps.
    // Flatten handles suffix-then-strip naming and alias rewriting.
    auto mat_set_name = [](GDMLMaterial & m, std::string const & n) { m.name = n; };
    auto renames_iso = isotope_tracker.Flatten(data.isotopes, nullptr, mat_set_name);
    auto renames_elem = element_tracker.Flatten(data.elements, nullptr, mat_set_name);
    auto renames_mat = material_tracker.Flatten(data.materials, &data.material_instance_counts, mat_set_name);

    // Merge all renames and rewrite composition references across all maps.
    std::unordered_map<std::string, std::string> all_renames;
    all_renames.insert(renames_iso.begin(), renames_iso.end());
    all_renames.insert(renames_elem.begin(), renames_elem.end());
    all_renames.insert(renames_mat.begin(), renames_mat.end());
    if(!all_renames.empty()) {
        auto rewriteComposition = [&](std::unordered_map<std::string, GDMLMaterial> & dest) {
            for(auto & pair : dest) {
                std::map<std::string, double> new_comp;
                for(auto const & comp : pair.second.composition) {
                    auto rit = all_renames.find(comp.first);
                    new_comp[rit != all_renames.end() ? rit->second : comp.first] = comp.second;
                }
                pair.second.composition = std::move(new_comp);
            }
        };
        rewriteComposition(data.isotopes);
        rewriteComposition(data.elements);
        rewriteComposition(data.materials);
    }
}

} // namespace gdml
} // namespace detector
} // namespace siren
