#ifndef LI_Placement_H
#define LI_Placement_H

#include <sstream>
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Rotation3D.h"

namespace earthmodel {

class Placement
{
public:
    // constructors
    Placement();
    Placement(Vector3D position, Rotation3D rotation);
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

private:
    Vector3D position_;
    Rotation3D rotation_;
};

} // namespace earthmodel

#endif // LI_Placement_H

