#include "SIREN/detector/GDMLParser.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cereal/external/rapidxml/rapidxml.hpp>

#include "gdml/GDMLParserPrivate.h"

namespace rapidxml = cereal::rapidxml;

namespace siren {
namespace detector {

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
    gdml::ExpandLoops(doc, gdml_node, data.constants);
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
