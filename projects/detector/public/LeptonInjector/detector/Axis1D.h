#pragma once
#ifndef LI_Axis1D_H
#define LI_Axis1D_H

#include <memory>                          // for shared_ptr
#include <cstdint>                         // for uint32_t
#include <stdexcept>                       // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "LeptonInjector/math/Vector3D.h"

namespace LI {
namespace detector {

class Axis1D {
friend cereal::access;
public:
    Axis1D();
    Axis1D(const math::Vector3D& fAxis, const math::Vector3D& fp0);
    Axis1D(const Axis1D&);

    virtual ~Axis1D() {};

    bool operator==(const Axis1D& axis) const;
    bool operator!=(const Axis1D& axis) const;
    virtual bool compare(const Axis1D& dens_distr) const = 0;

    virtual Axis1D* clone() const = 0;
    virtual std::shared_ptr<const Axis1D> create() const = 0;

    virtual double GetX(const math::Vector3D& xi) const = 0;
    virtual double GetdX(const math::Vector3D& xi, const math::Vector3D& direction) const = 0;

    math::Vector3D GetAxis() const { return fAxis_; };
    math::Vector3D GetFp0() const { return fp0_; };

    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(::cereal::make_nvp("Axis", fAxis_));
                archive(::cereal::make_nvp("Origin", fp0_));
            } else {
                throw std::runtime_error("Axis1D only supports version <= 0");
            }
        };
protected:
    math::Vector3D fAxis_;
    math::Vector3D fp0_;
};

} // namespace detector
} // namespace LI

CEREAL_CLASS_VERSION(LI::detector::Axis1D, 0);

#endif // LI_Axis1D_H
