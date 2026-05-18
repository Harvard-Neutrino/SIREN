#pragma once

#include <map>
#include <memory>
#include <string>
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

bool IsBuiltInConstant(std::string const & name);
void SeedBuiltInConstants(GDMLData & data);

double EvalExpression(std::string const & expr,
                      std::map<std::string, double> const & constants,
                      std::map<std::string, std::pair<int, std::vector<double>>> const & matrices = {});
double SafeParseDouble(const char* val,
                       std::map<std::string, double> const & constants = {},
                       std::map<std::string, std::pair<int, std::vector<double>>> const & matrices = {});

double LengthScale(const char* unit);
double ParseLength(const char* value, const char* unit,
                   std::map<std::string, double> const & constants = {});
double AngleScale(const char* unit);
double ParseAngle(const char* value, const char* unit,
                  std::map<std::string, double> const & constants = {});
math::Quaternion QuatFromGDMLRotation(double rx, double ry, double rz);

GDMLQuantityType ResolveQuantityType(const char* type_attr, const char* unit,
                                     std::string const & name,
                                     GDMLData & data,
                                     GDMLParseOptions const & options);
double QuantityUnitScale(const char* unit, GDMLQuantityType resolved_type);

double ReadDoubleAttr(rapidxml::xml_node<>* node, const char* attr,
                      std::map<std::string, double> const & constants,
                      double scale = 1.0);
double ReadLengthAttr(rapidxml::xml_node<>* node, const char* attr,
                      double lscale,
                      std::map<std::string, double> const & constants);
double ReadAngleAttr(rapidxml::xml_node<>* node, const char* attr,
                     double ascale,
                     std::map<std::string, double> const & constants);
GDMLPhiRange ReadPhiRange(rapidxml::xml_node<>* node,
                          double ascale,
                          std::map<std::string, double> const & constants);
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
                 std::map<std::string, double> & constants);

void EmitWarning(GDMLData & data, GDMLParseOptions const & options, std::string const & msg);

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
