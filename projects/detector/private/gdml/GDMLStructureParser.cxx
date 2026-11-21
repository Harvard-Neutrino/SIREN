#include "GDMLParserPrivate.h"

#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>

using namespace siren::math;

namespace siren {
namespace detector {
namespace gdml {
// Parse the <structure> section: volumes, assemblies, and physical volumes.
// Checks for ambiguous solid/material references (names with multiple instances).
void ParseStructure(rapidxml::xml_node<>* structure_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!structure_node) return;

    auto vol_eq = [](GDMLVolume const & a, GDMLVolume const & b) {
        if(a.solid_ref != b.solid_ref || a.material_ref != b.material_ref
           || a.children.size() != b.children.size()
           || a.is_assembly != b.is_assembly)
            return false;
        for(size_t i = 0; i < a.children.size(); ++i) {
            if(a.children[i].volume_ref != b.children[i].volume_ref
               || !(a.children[i].position == b.children[i].position)
               || !(a.children[i].rotation == b.children[i].rotation))
                return false;
        }
        return true;
    };
    InstanceTracker<GDMLVolume> volume_tracker(data, options.strip_pointer_suffixes);

    // Scope isolation for <file> sub-modules: track which volume names
    // came from sub-files so parent volumes cannot reference them.
    std::unordered_set<std::string> file_volumes;
    std::unordered_set<std::string> file_internal_volumes;

    // Shared physvol child parser used by both <volume> and <assembly>
    auto parsePhysVols = [&](rapidxml::xml_node<>* parent_node, GDMLVolume & vol) {
        for(auto* pv = parent_node->first_node("physvol"); pv; pv = pv->next_sibling("physvol")) {
            GDMLPhysVol physvol;
            physvol.position = Vector3D(0, 0, 0);
            physvol.rotation = Quaternion(); // identity

            auto* volref = pv->first_node("volumeref");
            if(volref) {
                // Store the raw ref; alias resolution happens later in
                // BuildVolume after the volume flatten/strip pass.
                physvol.volume_ref = SafeAttrVal(volref, "ref");
            }

            // Handle <file> references: parse the external GDML as a
            // fully independent file, then merge into the parent scope
            // using GDMLData::MergeFrom (handles name collisions).
            //
            // SIREN extension attribute (not in GDML schema):
            //   as_assembly="true" -- treat the referenced volume as an
            //                         assembly (skip its own sector, place
            //                         children directly).
            auto* file_node = pv->first_node("file");
            if(file_node) {
                std::string fname = SafeAttrVal(file_node, "name");
                std::string volname = SafeAttrVal(file_node, "volname");
                std::string as_asm_str = SafeAttrVal(file_node, "as_assembly");
                bool as_assembly = (as_asm_str == "true" || as_asm_str == "1");
                if(!fname.empty()) {
                    bool is_absolute = (!fname.empty() && fname[0] == '/');
                    std::string path = (is_absolute || data.base_dir.empty()) ? fname : data.base_dir + "/" + fname;
                    GDMLData sub = ParseGDML(path, options);

                    // Collect sub-file volume names before merge (maps
                    // are partially moved-from after MergeFrom).
                    std::vector<std::string> sub_vol_keys;
                    sub_vol_keys.reserve(sub.volumes.size());
                    for(auto const & p : sub.volumes)
                        sub_vol_keys.push_back(p.first);

                    auto renames = data.MergeFrom(sub);

                    auto r = [&](std::string const & name) -> std::string const & {
                        auto it = renames.find(name);
                        return (it != renames.end()) ? it->second : name;
                    };

                    for(auto const & n : sub_vol_keys) {
                        std::string const & merged = r(n);
                        file_volumes.insert(merged);
                        file_internal_volumes.insert(merged);
                    }

                    if(!volname.empty()) {
                        physvol.volume_ref = r(volname);
                    } else if(!sub.world_volume.empty()) {
                        physvol.volume_ref = sub.world_volume;
                    }

                    // The exposed world volume is the sole gateway into
                    // the sub-file; remove it from the internal set.
                    if(!physvol.volume_ref.empty())
                        file_internal_volumes.erase(physvol.volume_ref);

                    if(as_assembly && !physvol.volume_ref.empty()) {
                        auto vit = data.volumes.find(physvol.volume_ref);
                        if(vit != data.volumes.end()) {
                            vit->second.is_assembly = true;
                        }
                    }
                }
            }

            // Validate child elements: reject unrecognized tags that almost
            // certainly indicate typos (e.g. <positoinref> instead of <positionref>).
            for(auto* child = pv->first_node(); child; child = child->next_sibling()) {
                std::string tag(child->name(), child->name_size());
                if(tag.empty() || tag == "auxiliary" || tag == "file") {
                    continue;
                } else if(tag != "volumeref" && tag != "position" && tag != "positionref"
                          && tag != "rotation" && tag != "rotationref"
                          && tag != "scale" && tag != "scaleref") {
                    throw std::runtime_error(
                        "GDML error: physvol contains unrecognized child element <" + tag + ">");
                }
            }

            GDMLPlacement placement = ReadPlacement(
                pv, "position", "positionref", "rotation", "rotationref", data);
            physvol.position = placement.position;
            physvol.rotation = placement.rotation;

            vol.children.push_back(physvol);
        }
    };

    for(auto* vol_node = structure_node->first_node(); vol_node;
        vol_node = vol_node->next_sibling()) {

        std::string tag(vol_node->name());
        if(tag != "volume" && tag != "assembly") continue;

        GDMLVolume vol;
        vol.name = StripPointerSuffix(SafeAttrVal(vol_node, "name"), options.strip_pointer_suffixes);
        vol.is_assembly = (tag == "assembly");

        if(!vol.is_assembly) {
            // Get <materialref> — resolve through alias map for pointer-suffixed names.
            // Only check for ambiguity when the raw ref name was not resolved through
            // the alias map (i.e., it was a bare name that matches multiple instances).
            auto* matref = vol_node->first_node("materialref");
            if(matref) {
                std::string raw_ref = SafeAttrVal(matref, "ref");
                vol.material_ref = ResolveAlias(raw_ref, data);
                bool was_aliased = (vol.material_ref != raw_ref);
                if(!was_aliased && !vol.material_ref.empty()) {
                    auto it = data.material_instance_counts.find(vol.material_ref);
                    if(it != data.material_instance_counts.end() && it->second > 1) {
                        throw std::runtime_error(
                            "GDML error: volume '" + vol.name
                            + "' references ambiguous material '" + vol.material_ref
                            + "' which has " + std::to_string(it->second) + " instances");
                    }
                }
            }

            // Get <solidref> — same alias + ambiguity logic
            auto* solidref = vol_node->first_node("solidref");
            if(solidref) {
                std::string raw_ref = SafeAttrVal(solidref, "ref");
                vol.solid_ref = ResolveAlias(raw_ref, data);
                bool was_aliased = (vol.solid_ref != raw_ref);
                if(!was_aliased && !vol.solid_ref.empty()) {
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
            std::string original_name(SafeAttrVal(vol_node, "name"));
            volume_tracker.Store(original_name, vol, vol_eq);
        }
    }

    auto vol_set_name = [](GDMLVolume & v, std::string const & n) { v.name = n; };
    volume_tracker.Flatten(data.volumes, nullptr, vol_set_name);

    // Parent-defined volumes shadow file-internal names: Flatten
    // overwrites sub-file entries in data.volumes, so the parent's
    // definition wins and the name is no longer file-internal.
    for(auto const & pair : volume_tracker.Instances()) {
        if(pair.second.size() == 1) {
            file_internal_volumes.erase(pair.first);
            file_volumes.erase(pair.first);
        } else {
            for(int i = 0; i < (int)pair.second.size(); ++i) {
                file_internal_volumes.erase(volume_tracker.FlatName(pair.first, i));
                file_volumes.erase(volume_tracker.FlatName(pair.first, i));
            }
        }
    }

    // Resolve volumeref and replicavol child refs now that all volumes
    // are flattened and aliases are finalized.
    for(auto & vpair : data.volumes) {
        for(auto & child : vpair.second.children) {
            std::string const & raw = child.volume_ref;
            auto ait = data.name_aliases.find(raw);
            if(ait != data.name_aliases.end()) {
                child.volume_ref = ait->second;
            } else {
                std::string stripped = StripPointerSuffix(raw, options.strip_pointer_suffixes);
                if(stripped != raw) {
                    child.volume_ref = stripped;
                }
            }
        }
    }

    // Scope isolation: parent-defined volumes must not reference
    // sub-file internal volumes. Only the exposed world volume of
    // each <file> sub-module is accessible to the parent scope.
    if(!file_internal_volumes.empty()) {
        for(auto const & vpair : data.volumes) {
            if(file_volumes.count(vpair.first)) continue;
            for(auto const & child : vpair.second.children) {
                if(file_internal_volumes.count(child.volume_ref)) {
                    data.errors.push_back(
                        "GDML scope error: volume '" + vpair.first
                        + "' references '" + child.volume_ref
                        + "' which is internal to a <file> sub-module");
                }
            }
        }
        data.ThrowIfErrors();
    }
}


// Parse the <setup> section to find the world volume
void ParseSetup(rapidxml::xml_node<>* setup_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!setup_node) return;

    auto* world_node = setup_node->first_node("world");
    if(world_node) {
        data.world_volume = StripPointerSuffix(SafeAttrVal(world_node, "ref"), options.strip_pointer_suffixes);
    }
}

} // namespace gdml
} // namespace detector
} // namespace siren
