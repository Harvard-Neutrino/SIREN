#pragma once
#ifndef SIREN_pyDarkNewsCrossSection_H
#define SIREN_pyDarkNewsCrossSection_H

#include <memory>
#include <string>
#include <vector>                                       // for vector
#include <cstdint>                                      // for uint32_t
#include <stdexcept>                                    // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/interactions/DarkNewsCrossSection.h"  // for DarkNewsCrossSection

#include "SIREN/interactions/CrossSection.h"  // for CrossSection
#include "SIREN/dataclasses/Particle.h"        // for Particle
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// Trampoline class for DarkNewsCrossSection
class pyDarkNewsCrossSection : public DarkNewsCrossSection {
public:
    using DarkNewsCrossSection::DarkNewsCrossSection;
    pyDarkNewsCrossSection(DarkNewsCrossSection && parent);
    pyDarkNewsCrossSection(DarkNewsCrossSection const & parent);
    double TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const override;
    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override;
    double TotalCrossSection(siren::dataclasses::ParticleType primary, double energy, siren::dataclasses::ParticleType target) const override;
    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override;
    double DifferentialCrossSection(siren::dataclasses::ParticleType primary, siren::dataclasses::ParticleType target, double energy, double Q2) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override;
    double Q2Min(dataclasses::InteractionRecord const & interaction) const override;
    double Q2Max(dataclasses::InteractionRecord const & interaction) const override;
    double TargetMass(dataclasses::ParticleType const & target_type) const override;
    std::vector<double> SecondaryMasses(std::vector<dataclasses::ParticleType> const & secondary_types) const override;
    std::vector<double> SecondaryHelicities(dataclasses::InteractionRecord const & record) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;
    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    virtual ~pyDarkNewsCrossSection() = default;

    pybind11::object self;

    // First attempts to call a python-side override of the "get_representation" function
    // The assumption is that "get_representation" returns a python dictionary that contains the representation of the object
    // If "get_representation" is not overriden on the python side, then this function returns the contents of __dict__
    pybind11::object get_representation() {
        // First check if "self" is a valid reference to the python side of this object
        // Otherwise use "this" and search for the corresponding python side
        const DarkNewsCrossSection * ref;
        if(self) {
            ref = self.cast<DarkNewsCrossSection *>();
        } else {
            ref = dynamic_cast<const DarkNewsCrossSection *>(this);
            if (!ref) {
                int status;
                char * realname;
                const std::type_info  &ti = typeid(DarkNewsCrossSection);
                realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
                std::stringstream msg;
                msg << "Cannot cast this to " << realname;
                free(realname);
                throw std::runtime_error(msg.str());
            }
        }

        auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsCrossSection));
        pybind11::function override_func =
            tinfo ? pybind11::detail::get_type_override(static_cast<const DarkNewsCrossSection *>(ref), tinfo, "get_representation") : pybind11::function();
        if (override_func) {
            pybind11::object o = override_func();
            if(not pybind11::isinstance<pybind11::dict>(o)) {
                throw std::runtime_error("get_representation must return a dict");
            }
            return o;
        }

        pybind11::object _self;
        if(this->self) {
            _self = pybind11::reinterpret_borrow<pybind11::object>(this->self);
        } else {
            auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsCrossSection));
            pybind11::handle self_handle = get_object_handle(dynamic_cast<const DarkNewsCrossSection *>(this), tinfo);
            _self = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
        }
        pybind11::dict d;
        if (pybind11::hasattr(_self, "__dict__")) {
            d = _self.attr("__dict__");
        }
        return d;
    }

    static pybind11::tuple pickle_save(DarkNewsCrossSection & cpp_obj) {
        pybind11::object x = dynamic_cast<pyDarkNewsCrossSection *>(&cpp_obj)->get_representation();
        return pybind11::make_tuple(x);
    }

    static std::pair<std::unique_ptr<DarkNewsCrossSection>, pybind11::dict> pickle_load(const pybind11::tuple &t) {
        if (t.size() != 1) {
            std::cout << "Tuple has " << t.size() << " elements" << std::endl;
            throw std::runtime_error("Invalid state!");
        }
        auto cpp_state = std::unique_ptr<DarkNewsCrossSection>(new pyDarkNewsCrossSection);
        auto py_state = t[0].cast<pybind11::dict>();
        return std::make_pair(std::move(cpp_state), py_state);
    }

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            // Either use *self* or find the corresponsing python object for the instance of this class
            // Pass that python object (self) to pickle to get the byestream
            pybind11::object obj;
            if(this->self) {
                obj = this->self;
            } else {
                auto *tinfo = pybind11::detail::get_type_info(typeid(pyDarkNewsCrossSection));
                pybind11::handle self_handle = get_object_handle(dynamic_cast<const pyDarkNewsCrossSection *>(this), tinfo);
                obj = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
            }

			pybind11::module pkl = pybind11::module::import("pickle");
			pybind11::bytes bytes = pkl.attr("dumps")(obj);
			std::string str_repr = (std::string)(bytes.attr("hex")().cast<std::string>());

			archive(::cereal::make_nvp("PythonPickleBytesRepresentation", str_repr));

            archive(cereal::virtual_base_class<DarkNewsCrossSection>(dynamic_cast<const DarkNewsCrossSection*>(this)));

        } else {
            throw std::runtime_error("DarkNewsCrossSection only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            std::string str_repr;
			archive(::cereal::make_nvp("PythonPickleBytesRepresentation", str_repr));

            pybind11::module pkl = pybind11::module::import("pickle");

            pybind11::object bytes_module = pybind11::module::import("builtins").attr("bytes");
            pybind11::object bytes = bytes_module.attr("fromhex")(str_repr);

            pkl.attr("loads")(bytes);
            this->self = pkl.attr("loads")(bytes);

            archive(cereal::virtual_base_class<DarkNewsCrossSection>(dynamic_cast<const DarkNewsCrossSection*>(this)));

        } else {
            throw std::runtime_error("DarkNewsCrossSection only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace siren

using DarkNewsCrossSectionTrampolineType = Pybind11Trampoline<siren::interactions::DarkNewsCrossSection,siren::interactions::pyDarkNewsCrossSection>;
//RegisterTrampolineCerealMethods(siren::interactions::DarkNewsCrossSection, siren::interactions::pyDarkNewsCrossSection, DarkNewsCrossSectionTrampolineType)
CEREAL_CLASS_VERSION(siren::interactions::pyDarkNewsCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::pyDarkNewsCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::DarkNewsCrossSection, siren::interactions::pyDarkNewsCrossSection);
CEREAL_CLASS_VERSION(DarkNewsCrossSectionTrampolineType, 0);
CEREAL_REGISTER_TYPE(DarkNewsCrossSectionTrampolineType);
CEREAL_REGISTER_POLYMORPHIC_RELATION(DarkNewsCrossSectionTrampolineType, siren::interactions::pyDarkNewsCrossSection);

//CEREAL_FORCE_DYNAMIC_INIT(pyDarkNewsCrossSection);

#endif // SIREN_pyDarkNewsCrossSection_H

