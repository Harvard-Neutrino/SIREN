#include "SIREN/detector/GDMLParser.h"

#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <cereal/external/rapidxml/rapidxml.hpp>

#include "gdml/GDMLParserPrivate.h"

namespace rapidxml = cereal::rapidxml;

namespace siren {
namespace detector {

void GDMLData::ThrowIfErrors() {
    if(errors.empty()) return;
    std::string msg = "GDML parse errors (" + std::to_string(errors.size()) + "):";
    for(auto const & e : errors) {
        msg += "\n  - " + e;
    }
    errors.clear();
    throw std::runtime_error(msg);
}

std::unordered_map<std::string, std::string> GDMLData::MergeFrom(GDMLData & other) {
    std::unordered_map<std::string, std::string> renames;

    // Find a name that doesn't collide with the parent, the sub-file,
    // or any already-chosen rename target.
    auto inAny = [](std::string const & c, GDMLData const & d) {
        return d.materials.count(c) || d.elements.count(c)
            || d.isotopes.count(c) || d.solids.count(c)
            || d.volumes.count(c) || d.positions.count(c)
            || d.rotations.count(c) || d.constants.count(c);
    };
    auto findUnique = [&](std::string const & name) -> std::string {
        std::string candidate;
        int suffix = 2;
        for(;;) {
            candidate = name + "_" + std::to_string(suffix++);
            if(inAny(candidate, *this)) continue;
            if(inAny(candidate, other)) continue;
            bool rename_hit = false;
            for(auto const & r : renames)
                if(r.second == candidate) { rename_hit = true; break; }
            if(rename_hit) continue;
            return candidate;
        }
    };

    // Detect collisions across typed maps. Identical definitions are
    // silently reused (no rename). Differing definitions get a rename.
    auto checkMaterials = [&](auto const & sub_map, auto const & parent_map) {
        for(auto const & p : sub_map) {
            auto it = parent_map.find(p.first);
            if(it != parent_map.end() && !(it->second == p.second))
                renames[p.first] = findUnique(p.first);
        }
    };
    checkMaterials(other.isotopes, isotopes);
    checkMaterials(other.elements, elements);
    checkMaterials(other.materials, materials);

    for(auto const & p : other.solids) {
        auto it = solids.find(p.first);
        if(it != solids.end() && !(*it->second == *p.second))
            renames[p.first] = findUnique(p.first);
    }

    for(auto const & p : other.volumes) {
        auto it = volumes.find(p.first);
        if(it != volumes.end()) {
            auto const & a = it->second;
            auto const & b = p.second;
            bool differ = (a.solid_ref != b.solid_ref
                || a.material_ref != b.material_ref
                || a.children.size() != b.children.size()
                || a.is_assembly != b.is_assembly);
            if(!differ) {
                for(size_t i = 0; i < a.children.size(); ++i) {
                    if(a.children[i].volume_ref != b.children[i].volume_ref
                       || !(a.children[i].position == b.children[i].position)
                       || !(a.children[i].rotation == b.children[i].rotation)) {
                        differ = true;
                        break;
                    }
                }
            }
            if(differ) renames[p.first] = findUnique(p.first);
        }
    }

    // Apply renames to the other's internal references
    auto r = [&](std::string const & name) -> std::string const & {
        auto it = renames.find(name);
        return (it != renames.end()) ? it->second : name;
    };

    if(!renames.empty()) {
        other.world_volume = r(other.world_volume);

        for(auto & vp : other.volumes) {
            vp.second.name = r(vp.second.name);
            vp.second.solid_ref = r(vp.second.solid_ref);
            vp.second.material_ref = r(vp.second.material_ref);
            for(auto & child : vp.second.children)
                child.volume_ref = r(child.volume_ref);
        }

        for(auto * map_ptr : {&other.materials, &other.elements, &other.isotopes}) {
            for(auto & mp : *map_ptr) {
                mp.second.name = r(mp.second.name);
                std::map<std::string, double> new_comp;
                for(auto & cp : mp.second.composition)
                    new_comp[r(cp.first)] = cp.second;
                mp.second.composition = std::move(new_comp);
            }
        }
    }

    // Merge with renamed keys. For identical definitions that weren't
    // renamed, emplace is a no-op (parent's version is kept).
    auto merge = [&](auto & sub_map, auto & dst_map) {
        using Map = std::remove_reference_t<decltype(dst_map)>;
        Map staged;
        for(auto & p : sub_map)
            staged.emplace(r(p.first), std::move(p.second));
        for(auto & p : staged)
            dst_map.emplace(std::move(p.first), std::move(p.second));
    };
    merge(other.solids, solids);
    merge(other.isotopes, isotopes);
    merge(other.elements, elements);
    merge(other.materials, materials);
    merge(other.volumes, volumes);
    merge(other.positions, positions);
    merge(other.rotations, rotations);
    merge(other.constants, constants);

    // Propagate sub-file warnings so the caller sees diagnostics from
    // all included files, not just the top-level parse.
    warnings.insert(warnings.end(), other.warnings.begin(), other.warnings.end());

    return renames;
}


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
        content = gdml::ExpandEntities(content, base_dir);
    }

    // Handle xi:include (XInclude) before XML parsing.
    if(content.find("xi:include") != std::string::npos) {
        content = gdml::ExpandXIncludes(content, base_dir);
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
    data.base_dir = base_dir;

    gdml::SeedBuiltInConstants(data);

    // Find root node: accept <gdml> or any root element
    auto* gdml_node = doc.first_node("gdml");
    if(!gdml_node) {
        gdml_node = doc.first_node();
        if(!gdml_node) {
            throw std::runtime_error("GDML error: empty document: " + filename);
        }
        std::string root_tag(gdml_node->name());
        gdml::EmitWarning(data, options, "non-standard root element <" + root_tag + "> in " + filename + ", treating as <gdml>");
    }

    // Resolve all <define> sections in document order (handles loops inline).
    gdml::ResolveDefineInOrder(gdml_node, data, options);

    // Expand <loop> elements in non-define sections (structure, solids) where
    // downstream parsers need to see the expanded DOM nodes.
    gdml::ExpandLoops(doc, gdml_node, data.constants, data);
    data.ThrowIfErrors();
    // Remove loops from <define> sections (already resolved, but the DOM nodes
    // are still present; remove them so ParseDefine is not confused if called again)
    for(auto* def = gdml_node->first_node("define"); def; def = def->next_sibling("define")) {
        auto* child = def->first_node();
        while(child) {
            auto* next = child->next_sibling();
            if(std::strcmp(child->name(), "loop") == 0) {
                def->remove_node(child);
            }
            child = next;
        }
    }

    // Parse all <materials> sections (with instance-scoped name resolution)
    gdml::ParseAllMaterials(gdml_node, data, options);
    // Parse all <solids> sections (with instance-scoped boolean resolution)
    gdml::ParseAllSolids(gdml_node, data, options);
    // Parse all <structure> sections
    for(auto* node = gdml_node->first_node("structure"); node; node = node->next_sibling("structure")) {
        gdml::ParseStructure(node, data, options);
    }
    // Parse all <setup> sections (last wins)
    for(auto* node = gdml_node->first_node("setup"); node; node = node->next_sibling("setup")) {
        gdml::ParseSetup(node, data, options);
    }

    // Warn if no <setup> section was found (world volume undefined)
    if(data.world_volume.empty()) {
        gdml::EmitWarning(data, options, "No <setup> section found in GDML; world volume is undefined");
    }

    return data;
}

} // namespace detector
} // namespace siren
