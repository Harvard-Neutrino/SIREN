#pragma once
#ifndef LI_Placement_H
#define LI_Placement_H

#include <memory>
#include <cstdint>
#include <sstream>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/math/Quaternion.h"

namespace LI {
namespace geometry {

class Placement
{
public:
    // constructors
    Placement();
    Placement(LI::math::Vector3D const & position);
    Placement(LI::math::Quaternion const & quaternion);
    Placement(LI::math::Vector3D const & position, LI::math::Quaternion const & quaternion);
    Placement(const Placement& placement);
    Placement(Placement&& other);
    ~Placement();

    //-------------------------------------//
    // operator functions and swap
    Placement& operator=(Placement const & placement);
    Placement& operator=(Placement const && placement);
    Placement& operator=(Placement && placement);
    bool operator==(const Placement& placement) const;
    bool operator!=(const Placement& placement) const;
    bool operator<(const Placement& placement) const;
    void swap(Placement& placement);
    friend std::ostream& operator<<(std::ostream& os, Placement const& placement);

    std::shared_ptr<Placement> create() const;

    //-------------------------------------//
    // getter and setter functions
    LI::math::Vector3D GetPosition() const;
    LI::math::Quaternion GetQuaternion() const;

    void SetPosition(LI::math::Vector3D const &);
    void SetQuaternion(LI::math::Quaternion const &);

    //-------------------------------------//
    // composition function (for rotating)
    LI::math::Vector3D Rotate(LI::math::Vector3D const & p, bool inv = false) const;
    LI::math::Vector3D GlobalToLocalPosition(LI::math::Vector3D const & p) const;
    LI::math::Vector3D LocalToGlobalPosition(LI::math::Vector3D const & p) const;
    LI::math::Vector3D GlobalToLocalDirection(LI::math::Vector3D const & d) const;
    LI::math::Vector3D LocalToGlobalDirection(LI::math::Vector3D const & d) const;

    //-------------------------------------//
    // serialization
    //----------------------------------------------//
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Position", position_));
            archive(cereal::make_nvp("Quaternion", quaternion_));
        } else {
            throw std::runtime_error("Placement only supports version <= 0!");
        }
    }

private:
    LI::math::Vector3D position_;

    // Describes the "active" rotation
    // i.e. composition with the quaternion takes a position defined relative to an un-rotated
    // body (in "body coordinates") and rotates it into the global system ("global coordinates)
    // Connecting to Aerospace rotation conventions, the "active" rotation is the same as the
    // Passive Body To World" (PBTW) transformation as described in https://arxiv.org/abs/1801.07478
    LI::math::Quaternion quaternion_;
};

} // namespace geometry
} // namespace LI

CEREAL_CLASS_VERSION(LI::geometry::Placement, 0);

#endif // LI_Placement_H

