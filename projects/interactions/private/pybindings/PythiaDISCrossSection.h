#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

#ifdef SIREN_HAS_PYTHIA8
#include "../../public/SIREN/interactions/PythiaDISCrossSection.h"

void register_PythiaDISCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<PythiaDISCrossSection, std::shared_ptr<PythiaDISCrossSection>, CrossSection> pythiadis(m, "PythiaDISCrossSection");

    pythiadis

        .def(init<>())
        .def(init<std::string, std::string, int, double, double,
                   std::set<siren::dataclasses::ParticleType>,
                   std::set<siren::dataclasses::ParticleType>,
                   std::string, std::string, std::string>(),
                arg("differential_filename"),
                arg("total_filename"),
                arg("interaction_type"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("pythia_data_path"),
                arg("pdf_set") = std::string("LHAPDF6:HERAPDF20_NLO_EIG"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, int, double, double,
                   std::vector<siren::dataclasses::ParticleType>,
                   std::vector<siren::dataclasses::ParticleType>,
                   std::string, std::string, std::string>(),
                arg("differential_filename"),
                arg("total_filename"),
                arg("interaction_type"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("pythia_data_path"),
                arg("pdf_set") = std::string("LHAPDF6:HERAPDF20_NLO_EIG"),
                arg("units") = std::string("cm"))
        .def(self == self)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&PythiaDISCrossSection::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double>(&PythiaDISCrossSection::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&PythiaDISCrossSection::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<double, double, double, double, double>(&PythiaDISCrossSection::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&PythiaDISCrossSection::InteractionThreshold)
        .def("FragmentationFraction",&PythiaDISCrossSection::FragmentationFraction)
        .def("GetPossibleTargets",&PythiaDISCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&PythiaDISCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&PythiaDISCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&PythiaDISCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&PythiaDISCrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&PythiaDISCrossSection::FinalStateProbability)
        .def("GetMinimumQ2",&PythiaDISCrossSection::GetMinimumQ2)
        .def("GetTargetMass",&PythiaDISCrossSection::GetTargetMass)
        .def("GetInteractionType",&PythiaDISCrossSection::GetInteractionType)
        .def_static("GeneratePythiaCharmSamples",
            [](int interaction_type, int primary_pdg, int target_pdg, double target_mass,
               std::string pdf_set, std::string pythia_data_path, double minimum_Q2,
               std::vector<double> energies, int n_events) {
                std::vector<double> sigma_mb, E, x, y;
                PythiaDISCrossSection::GeneratePythiaCharmSamples(
                    interaction_type, primary_pdg, target_pdg, target_mass,
                    pdf_set, pythia_data_path, minimum_Q2, energies, n_events,
                    sigma_mb, E, x, y);
                return pybind11::make_tuple(sigma_mb, E, x, y);
            },
            arg("interaction_type"), arg("primary_pdg"), arg("target_pdg"), arg("target_mass"),
            arg("pdf_set"), arg("pythia_data_path"), arg("minimum_Q2"),
            arg("energies"), arg("n_events"),
            "Run Pythia (init once per energy) and return (sigma_mb_per_E, E, x, y) raw samples "
            "for building charm-DIS splines. See siren.interactions.pythia_charm_splines.");
}

#else

void register_PythiaDISCrossSection(pybind11::module_ & m) {
    // Pythia8 is not available, so we do not register the PythiaDISCrossSection class.
}

#endif
