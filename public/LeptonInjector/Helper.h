#ifndef LI_HELPER
#define LI_HELPER

#include <array>

namespace LeptonInjector {

	// these functions rotate a 3-vector about some axis
	// assumes vector is in cartesian coordinates
	std::array<double, 3> RotateY(std::array<double 3> vector, double angle);
	std::array<double, 3> RotateX(std::array<double 3> vector, double angle);
	std::array<double, 3> RotateZ(std::array<double 3> vector, double angle);

} // end namespace LeptonInjector

#endif
