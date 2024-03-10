#pragma once
#ifndef SIREN_RadialAxis1D_H
#define SIREN_RadialAxis1D_H

#include <memory>                            // for shared_ptr
#include <cstdint>                           // for uint32_t
#include <stdexcept>                         // for runtime_error

#include <cereal/cereal.hpp>                 // for CEREAL_CLASS_VERSION
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/base_class.hpp>       // for virtual_base_class
#include <cereal/types/polymorphic.hpp>      // for CEREAL_REGISTER_POLYMORP...

#include "SIREN/math/Vector3D.h"    // for Vector3D

#include "SIREN/detector/Axis1D.h"  // for Axis1D

namespace siren {
namespace detector {

class RadialAxis1D : public Axis1D {
friend cereal::access;
public:
    RadialAxis1D();
    RadialAxis1D(const math::Vector3D& fAxis, const math::Vector3D& fp0);
    RadialAxis1D(const math::Vector3D& fp0);
    ~RadialAxis1D() {};

    bool compare(const Axis1D& dens_distr) const override;

    Axis1D* clone() const override { return new RadialAxis1D(*this); };
    std::shared_ptr<Axis1D> create() const override {
        return std::shared_ptr<Axis1D>(new RadialAxis1D(*this));
    };

    double GetX(const math::Vector3D& xi) const override;
    double GetdX(const math::Vector3D& xi, const math::Vector3D& direction) const override;

    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(cereal::virtual_base_class<Axis1D>(this));
            } else {
                throw std::runtime_error("RadialAxis1D only supports version <= 0");
            }
        };
};

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::RadialAxis1D, 0);
CEREAL_REGISTER_TYPE(siren::detector::RadialAxis1D);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::detector::Axis1D, siren::detector::RadialAxis1D);

#endif // SIREN_RadialAxis1D_H
