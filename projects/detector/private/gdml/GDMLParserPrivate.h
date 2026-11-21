#pragma once

#include <cctype>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "SIREN/detector/GDMLParser.h"

#include <cereal/external/rapidxml/rapidxml.hpp>

#include "SIREN/math/Quaternion.h"
#include "SIREN/math/Vector3D.h"

namespace rapidxml = cereal::rapidxml;

namespace siren {
namespace detector {
namespace gdml {

struct GDMLPlacement {
    math::Vector3D position = math::Vector3D(0, 0, 0);
    math::Quaternion rotation = math::Quaternion();
    bool specified = false;
};

struct GDMLPhiRange {
    double start = 0.0;
    double delta = 0.0;
};

struct GDMLZPlanes {
    std::vector<double> z;
    std::vector<double> rmin;
    std::vector<double> rmax;
};

const char* SafeAttrVal(rapidxml::xml_node<>* node, const char* attr_name);
std::string TrimWhitespace(std::string const & s);

// Strip trailing Geant4 pointer-address suffix (e.g. "Si0x7f3a4c002a" -> "Si").
// Requires at least 7 hex digits after "0x" to avoid false positives on
// intentional hex names (real pointer addresses are 9-12 hex digits).
// Returns the name unchanged if no valid suffix is found or if disabled.
inline std::string StripPointerSuffix(std::string const & name, bool enabled = true) {
    if(!enabled) return name;
    size_t pos = name.rfind("0x");
    if(pos == std::string::npos || pos == 0) return name;
    size_t hex_len = name.size() - (pos + 2);
    if(hex_len < 7) return name;
    for(size_t i = pos + 2; i < name.size(); ++i) {
        char c = name[i];
        bool is_hex = (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
        if(!is_hex) return name;
    }
    return name.substr(0, pos);
}

// Resolve a name through the alias map, returning the canonical name.
// If the name is not in the alias map, returns it unchanged.
inline std::string ResolveAlias(std::string const & name, GDMLData const & data) {
    auto it = data.name_aliases.find(name);
    return (it != data.name_aliases.end()) ? it->second : name;
}

inline bool IsIdentStart(char c) {
    return std::isalpha(static_cast<unsigned char>(c)) || c == '_';
}

inline bool IsIdentChar(char c) {
    return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}

bool IsBuiltInConstant(std::string const & name);
void SeedBuiltInConstants(GDMLData & data);

// Returns true if expr contains any built-in identifier that carries a
// dimensional scale (gram, cm, eV, ...). Math-only built-ins (pi, e, ...)
// are allowed. On a hit, found_name is set to the offending identifier.
// Used by the materials parser to reject <D value="N*gram"/> and similar
// inline-unit density expressions that would silently misparse.
bool ExpressionContainsUnitBuiltIn(std::string const & expr, std::string & found_name);

double EvalExpression(std::string const & expr,
                      std::unordered_map<std::string, double> const & constants,
                      std::unordered_map<std::string, std::pair<int, std::vector<double>>> const & matrices = {},
                      int depth = 0);
double SafeParseDouble(const char* val,
                       std::unordered_map<std::string, double> const & constants = {},
                       std::unordered_map<std::string, std::pair<int, std::vector<double>>> const & matrices = {});

double LengthScale(const char* unit);
double ParseLength(const char* value, const char* unit,
                   std::unordered_map<std::string, double> const & constants = {});
double AngleScale(const char* unit);
double ParseAngle(const char* value, const char* unit,
                  std::unordered_map<std::string, double> const & constants = {});
math::Quaternion QuatFromGDMLRotation(double rx, double ry, double rz);

GDMLQuantityType ResolveQuantityType(const char* type_attr, const char* unit,
                                     std::string const & name,
                                     GDMLData & data,
                                     GDMLParseOptions const & options);
double QuantityUnitScale(const char* unit, GDMLQuantityType resolved_type);

double ReadDoubleAttr(rapidxml::xml_node<>* node, const char* attr,
                      std::unordered_map<std::string, double> const & constants,
                      double scale = 1.0);
double ReadLengthAttr(rapidxml::xml_node<>* node, const char* attr,
                      double lscale,
                      std::unordered_map<std::string, double> const & constants);
double ReadAngleAttr(rapidxml::xml_node<>* node, const char* attr,
                     double ascale,
                     std::unordered_map<std::string, double> const & constants);
GDMLPhiRange ReadPhiRange(rapidxml::xml_node<>* node,
                          double ascale,
                          std::unordered_map<std::string, double> const & constants);
math::Vector3D ReadPositionXYZ(rapidxml::xml_node<>* node, GDMLData const & data);
math::Quaternion ReadRotationXYZ(rapidxml::xml_node<>* node, GDMLData const & data);
math::Vector3D ReadPositionReferenceAttr(rapidxml::xml_node<>* node,
                                         const char* attr,
                                         GDMLData const & data);
GDMLPlacement ReadPlacement(rapidxml::xml_node<>* node,
                            const char* position_tag,
                            const char* position_ref_tag,
                            const char* rotation_tag,
                            const char* rotation_ref_tag,
                            GDMLData const & data);
GDMLZPlanes ReadZPlanes(rapidxml::xml_node<>* node,
                        double lscale,
                        GDMLData const & data);
bool ValidateZPlaneMonotonicity(GDMLZPlanes const & planes);

std::string ExpandEntities(std::string const & content, std::string const & base_dir);
std::string ExpandXIncludes(std::string const & content, std::string const & base_dir);
void ExpandLoops(rapidxml::xml_document<> & doc,
                 rapidxml::xml_node<>* root,
                 std::unordered_map<std::string, double> & constants,
                 GDMLData & data);

void EmitWarning(GDMLData & data, GDMLParseOptions const & options, std::string const & msg);

// =========================================================================
// InstanceTracker: reusable instance-tracking with pointer-suffix stripping,
// dedup, alias recording, and suffix-then-strip flatten.
//
// Used by materials (GDMLMaterial), solids (shared_ptr<Geometry>), and
// volumes (GDMLVolume) to handle Geant4 exports with pointer-suffixed names
// and hand-written GDML with duplicate definitions.
// =========================================================================
template<typename T>
class InstanceTracker {
public:
    using EqualFn = std::function<bool(T const &, T const &)>;

    InstanceTracker(GDMLData & data, bool strip_suffixes)
        : data_(data), strip_(strip_suffixes) {}

    // Store an item under original_name. Strips pointer suffix, dedups
    // against the most recent instance (using eq), records alias if the
    // original name had a pointer suffix.
    void Store(std::string const & original_name, T const & item, EqualFn const & eq) {
        std::string stripped = StripPointerSuffix(original_name, strip_);
        bool had_suffix = (stripped != original_name);
        auto & vec = instances_[stripped];
        if(!vec.empty() && eq(vec.back(), item)) {
            if(had_suffix) {
                data_.name_aliases[original_name] = FlatName(stripped, (int)vec.size() - 1);
            }
            return;
        }
        vec.push_back(item);
        counts_[stripped] = (int)vec.size();
        if(had_suffix) {
            data_.name_aliases[original_name] = FlatName(stripped, (int)vec.size() - 1);
        }
    }

    // Resolve a reference name to its flattened name for the most recent
    // instance. Strips pointer suffix first.
    std::string Resolve(std::string const & ref_name) const {
        std::string name = StripPointerSuffix(ref_name, strip_);
        auto it = counts_.find(name);
        if(it != counts_.end()) {
            return FlatName(name, it->second - 1);
        }
        return name;
    }

    // Look up the most recent instance by name. Checks alias map first,
    // then falls back to stripped name.
    T const * Lookup(std::string const & original_name) const {
        auto ait = data_.name_aliases.find(original_name);
        if(ait != data_.name_aliases.end()) {
            auto parsed = ParseFlatName(ait->second);
            if(parsed.second >= 0) {
                auto sit = instances_.find(parsed.first);
                if(sit != instances_.end() && parsed.second < (int)sit->second.size()) {
                    return &sit->second[parsed.second];
                }
            }
        }
        std::string name = StripPointerSuffix(original_name, strip_);
        auto it = counts_.find(name);
        if(it != counts_.end()) {
            return &instances_.at(name)[it->second - 1];
        }
        return nullptr;
    }

    // Instance count for a stripped name.
    int Count(std::string const & name) const {
        auto it = counts_.find(name);
        return (it != counts_.end()) ? it->second : 0;
    }

    // Access raw instances map.
    std::unordered_map<std::string, std::vector<T>> const & Instances() const { return instances_; }

    using NameSetter = std::function<void(T &, std::string const &)>;

    // Flatten into dest map with suffix-then-strip naming. Populates
    // instance_counts and returns the rename map (suffixed -> bare) for
    // sole instances so callers can rewrite cross-map references.
    // If set_name is provided, it is called on each item with its final
    // map key (handles types like GDMLMaterial/GDMLVolume that have a
    // name field that must match the map key).
    std::unordered_map<std::string, std::string> Flatten(
            std::unordered_map<std::string, T> & dest,
            std::unordered_map<std::string, int> * instance_counts = nullptr,
            NameSetter set_name = nullptr) {
        for(auto const & pair : instances_) {
            std::string const & name = pair.first;
            auto const & vec = pair.second;
            if(instance_counts) {
                (*instance_counts)[name] = (int)vec.size();
            }
            for(int i = 0; i < (int)vec.size(); ++i) {
                std::string flat = FlatName(name, i);
                dest[flat] = vec[i];
                if(set_name) set_name(dest[flat], flat);
            }
        }
        // Build rename map and apply strip for sole instances
        std::unordered_map<std::string, std::string> renames;
        for(auto const & pair : instances_) {
            if(pair.second.size() == 1) {
                renames[FlatName(pair.first, 0)] = pair.first;
            }
        }
        for(auto const & r : renames) {
            auto it = dest.find(r.first);
            if(it != dest.end()) {
                T val = std::move(it->second);
                dest.erase(it);
                if(set_name) set_name(val, r.second);
                dest[r.second] = std::move(val);
            }
        }
        // Update aliases that point to renamed keys
        for(auto & alias : data_.name_aliases) {
            auto rit = renames.find(alias.second);
            if(rit != renames.end()) {
                alias.second = rit->second;
            }
        }
        return renames;
    }

    std::string FlatName(std::string const & name, int idx) const {
        return name + "__" + std::to_string(idx + 1);
    }

private:
    // Parse "name__N" into (name, N-1). Returns ("", -1) on failure.
    static std::pair<std::string, int> ParseFlatName(std::string const & flat) {
        size_t sep = flat.rfind("__");
        if(sep == std::string::npos || sep == 0 || sep + 2 >= flat.size()) {
            return {"", -1};
        }
        bool all_digits = true;
        for(size_t i = sep + 2; i < flat.size(); ++i) {
            if(!std::isdigit(static_cast<unsigned char>(flat[i]))) {
                all_digits = false;
                break;
            }
        }
        if(!all_digits) return {"", -1};
        int idx = std::stoi(flat.substr(sep + 2)) - 1;
        return {flat.substr(0, sep), idx};
    }

    GDMLData & data_;
    bool strip_;
    std::unordered_map<std::string, std::vector<T>> instances_;
    std::unordered_map<std::string, int> counts_;
};

void ResolveDefineInOrder(rapidxml::xml_node<>* gdml_node,
                          GDMLData & data,
                          GDMLParseOptions const & options);
void ParseAllMaterials(rapidxml::xml_node<>* root_node,
                       GDMLData & data,
                       GDMLParseOptions const & options);
void ParseAllSolids(rapidxml::xml_node<>* root_node,
                    GDMLData & data,
                    GDMLParseOptions const & options);
void ParseStructure(rapidxml::xml_node<>* structure_node,
                    GDMLData & data,
                    GDMLParseOptions const & options);
void ParseSetup(rapidxml::xml_node<>* setup_node,
                GDMLData & data,
                GDMLParseOptions const & options);

} // namespace gdml
} // namespace detector
} // namespace siren
