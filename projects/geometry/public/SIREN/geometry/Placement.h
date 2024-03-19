#pragma once
#ifndef SIREN_Placement_H
#define SIREN_Placement_H

#include <memory>
#include <cstdint>
#include <sstream>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include "SIREN/math/Vector3D.h"
#include "SIREN/math/Quaternion.h"

namespace siren {
namespace geometry {

class Placement
{
public:
    // constructors
    Placement();
    Placement(siren::math::Vector3D const & position);
    Placement(siren::math::Quaternion const & quaternion);
    Placement(siren::math::Vector3D const & position, siren::math::Quaternion const & quaternion);
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
    siren::math::Vector3D GetPosition() const;
    siren::math::Quaternion GetQuaternion() const;

    void SetPosition(siren::math::Vector3D const &);
    void SetQuaternion(siren::math::Quaternion const &);

    //-------------------------------------//
    // composition function (for rotating)
    siren::math::Vector3D Rotate(siren::math::Vector3D const & p, bool inv = false) const;
    siren::math::Vector3D GlobalToLocalPosition(siren::math::Vector3D const & p) const;
    siren::math::Vector3D LocalToGlobalPosition(siren::math::Vector3D const & p) const;
    siren::math::Vector3D GlobalToLocalDirection(siren::math::Vector3D const & d) const;
    siren::math::Vector3D LocalToGlobalDirection(siren::math::Vector3D const & d) const;

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
    siren::math::Vector3D position_;

    // Describes the "active" rotation
    // i.e. composition with the quaternion takes a position defined relative to an un-rotated
    // body (in "body coordinates") and rotates it into the global system ("global coordinates)
    // Connecting to Aerospace rotation conventions, the "active" rotation is the same as the
    // Passive Body To World" (PBTW) transformation as described in https://arxiv.org/abs/1801.07478
    siren::math::Quaternion quaternion_;
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Placement, 0);

#endif // SIREN_Placement_H

