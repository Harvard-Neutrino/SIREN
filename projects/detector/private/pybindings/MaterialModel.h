#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/detector/MaterialModel.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"

void register_MaterialModel(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::detector;

    class_<MaterialModel, std::shared_ptr<MaterialModel>> material_model(m, "MaterialModel");
    material_model
        .def(init<>())
        .def(init<std::string>())
        .def(init<std::string, std::string>())
        .def(init<std::vector<std::string>>())
        .def(init<std::string, std::vector<std::string>>())
        .def(self == self)
        .def("SetPath", &MaterialModel::SetPath)
        .def("AddMaterial", &MaterialModel::AddMaterial)
        .def("AddModelFiles", &MaterialModel::AddModelFiles)
        .def("AddModelFile", &MaterialModel::AddModelFile)
        .def("GetMaterialTargets", &MaterialModel::GetMaterialTargets)
        .def("GetMaterialRadiationLength", &MaterialModel::GetMaterialRadiationLength)
        .def("GetMaterialName", &MaterialModel::GetMaterialName)
        .def("GetMaterialId", &MaterialModel::GetMaterialId)
        .def("HasMaterial", (bool (MaterialModel::*)(std::string const &) const)(&MaterialModel::HasMaterial))
        .def("HasMaterial", (bool (MaterialModel::*)(int) const)(&MaterialModel::HasMaterial))
        .def("GetTargetMassFraction", (double (MaterialModel::*)(int, siren::dataclasses::ParticleType) const)(&MaterialModel::GetTargetMassFraction))
        .def("GetTargetMassFraction", (std::vector<double> (MaterialModel::*)(int, std::vector<siren::dataclasses::ParticleType> const &) const)(&MaterialModel::GetTargetMassFraction))
        .def("GetTargetParticleFraction", (double (MaterialModel::*)(int, siren::dataclasses::ParticleType) const)(&MaterialModel::GetTargetParticleFraction))
        .def("GetTargetParticleFraction", (std::vector<double> (MaterialModel::*)(int, std::vector<siren::dataclasses::ParticleType> const &) const)(&MaterialModel::GetTargetParticleFraction))
        .def("GetTargetRadiationFraction", (std::vector<double> (MaterialModel::*)(int, std::vector<siren::dataclasses::ParticleType> const &) const)(&MaterialModel::GetTargetRadiationFraction))
        .def("GetMaterialConstituents", &MaterialModel::GetMaterialConstituents)
        .def_static("GetNucleonContent", &MaterialModel::GetNucleonContent)
        .def_static("GetMolarMass", &MaterialModel::GetMolarMass)
        .def_static("GetStrangeCount", &MaterialModel::GetStrangeCount)
        .def_static("GetNucleonCount", &MaterialModel::GetNucleonCount)
        .def_static("GetNeutronCount", &MaterialModel::GetNeutronCount)
        .def_static("GetProtonCount", &MaterialModel::GetProtonCount)
        .def_static("GetEmpericalNuclearBindingEnergy", &MaterialModel::GetEmpericalNuclearBindingEnergy)
        ;

    class_<MaterialModel::Component>(material_model, "Component")
        .def(init<>())
        .def(init<siren::dataclasses::ParticleType>())
        .def(self == self)
        .def_readwrite("type", &MaterialModel::Component::type)
        .def_readwrite("strange_count", &MaterialModel::Component::strange_count)
        .def_readwrite("neutron_count", &MaterialModel::Component::neutron_count)
        .def_readwrite("nucleon_count", &MaterialModel::Component::nucleon_count)
        .def_readwrite("proton_count", &MaterialModel::Component::proton_count)
        .def_readwrite("molar_mass", &MaterialModel::Component::molar_mass)
        .def_readwrite("is_atom", &MaterialModel::Component::is_atom);

    class_<MaterialModel::MaterialComponent>(material_model, "MaterialComponent")
        .def(init<>())
        .def_readwrite("component", &MaterialModel::MaterialComponent::component)
        .def_readwrite("mass_density_over_total_mass_density", &MaterialModel::MaterialComponent::mass_density_over_total_mass_density)
        .def_readwrite("particle_density_over_total_mass_density", &MaterialModel::MaterialComponent::particle_density_over_total_mass_density);
}
