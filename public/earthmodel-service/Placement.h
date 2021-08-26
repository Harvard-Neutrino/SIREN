#ifndef LI_Placement_H
#define LI_Placement_H

#include <memory>
#include <sstream>

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Quaternion.h"

namespace earthmodel {

class Placement
{
public:
    // constructors
    Placement();
    Placement(Vector3D position, Quaternion quaternion);
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
    void swap(Placement& placement);
    friend std::ostream& operator<<(std::ostream& os, Placement const& placement);

    std::shared_ptr<const Placement> create() const;

    //-------------------------------------//
    // getter and setter functions
    Vector3D GetPosition();
    Quaternion GetQuaternion();

    //-------------------------------------//
    // composition function (for rotating)
    Vector3D Compose(Vector3D const & p, bool inv = false) const;

private:
    Vector3D position_;
    Quaternion quaternion_;
};

} // namespace earthmodel

#endif // LI_Placement_H

