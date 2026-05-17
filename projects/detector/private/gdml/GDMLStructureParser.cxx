#include "GDMLParserPrivate.h"

#include <cmath>
#include <stdexcept>
#include <string>

using namespace siren::math;

namespace siren {
namespace detector {
namespace gdml {
// Parse the <structure> section: volumes, assemblies, and physical volumes.
// Checks for ambiguous solid/material references (names with multiple instances).
void ParseStructure(rapidxml::xml_node<>* structure_node, GDMLData & data, GDMLParseOptions const & options) {
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
void ParseSetup(rapidxml::xml_node<>* setup_node, GDMLData & data, GDMLParseOptions const & options) {
    (void)options;
    if(!setup_node) return;

    auto* world_node = setup_node->first_node("world");
    if(world_node) {
        data.world_volume = SafeAttrVal(world_node, "ref");
    }
}

} // namespace gdml
} // namespace detector
} // namespace siren
