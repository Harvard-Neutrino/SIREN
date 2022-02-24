#pragma once
#ifndef LI_Placement_H
#define LI_Placement_H

#include <memory>
#include <sstream>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Quaternion.h"

namespace earthmodel {

class Placement
{
public:
    // constructors
    Placement();
    Placement(Vector3D const & position);
    Placement(Quaternion const & quaternion);
    Placement(Vector3D const & position, Quaternion const & quaternion);
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

    std::shared_ptr<const Placement> create() const;

    //-------------------------------------//
    // getter and setter functions
    Vector3D GetPosition() const;
    Quaternion GetQuaternion() const;

    void SetPosition(Vector3D const &);
    void SetQuaternion(Quaternion const &);

    //-------------------------------------//
    // composition function (for rotating)
    Vector3D Rotate(Vector3D const & p, bool inv = false) const;
    Vector3D GlobalToLocalPosition(Vector3D const & p) const;
    Vector3D LocalToGlobalPosition(Vector3D const & p) const;
    Vector3D GlobalToLocalDirection(Vector3D const & d) const;
    Vector3D LocalToGlobalDirection(Vector3D const & d) const;

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
    Vector3D position_;

    // Describes the "active" rotation
    // i.e. composition with the quaternion takes a position defined relative to an un-rotated
    // body (in "body coordinates") and rotates it into the global system ("global coordinates)
    // Connecting to Aerospace rotation conventions, the "active" rotation is the same as the
    // Passive Body To World" (PBTW) transformation as described in https://arxiv.org/abs/1801.07478
    Quaternion quaternion_;
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::Placement, 0);

#endif // LI_Placement_H

