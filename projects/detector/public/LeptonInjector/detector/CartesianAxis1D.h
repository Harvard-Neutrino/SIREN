#pragma once
#ifndef LI_CartesianAxis1D_H
#define LI_CartesianAxis1D_H
#include <memory>
#include <string>
#include <exception>
#include <functional>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/detector/Axis1D.h"

namespace LI {
namespace detector {

class CartesianAxis1D : public Axis1D {
public:
    CartesianAxis1D();
    CartesianAxis1D(const math::Vector3D& fAxis, const math::Vector3D& fp0);
    ~CartesianAxis1D() {};

    bool compare(const Axis1D& dens_distr) const override;

    Axis1D* clone() const override { return new CartesianAxis1D(*this); };
    std::shared_ptr<const Axis1D> create() const override {
        return std::shared_ptr<const Axis1D>(new CartesianAxis1D(*this));
    };

    double GetX(const math::Vector3D& xi) const override;
    double GetdX(const math::Vector3D& xi, const math::Vector3D& direction) const override;

    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(cereal::virtual_base_class<Axis1D>(this));
            } else {
                throw std::runtime_error("CartesianAxis1D only supports version <= 0");
            }
        };
};

} // namespace detector
} // namespace LI

CEREAL_CLASS_VERSION(LI::detector::CartesianAxis1D, 0);
CEREAL_REGISTER_TYPE(LI::detector::CartesianAxis1D);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::detector::Axis1D, LI::detector::CartesianAxis1D);

#endif // LI_CartesianAxis1D_H
