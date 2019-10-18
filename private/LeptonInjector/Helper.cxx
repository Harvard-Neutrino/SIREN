#ifndef LI_HELPER
#define LI_HELPER

namespace LeptonInjector {

	std::array<double, 3> RotateX(std::array<double 3> vector, double angle) {
		std::array<double, 3> rotated = {vector[0], vector[1]*cos(angle) +vector[2]*sin(angle), -1*vector[1]*sin(angle) + vector[2]*cos(angle) };
		return(rotated);
	}

	std::array<double, 3> RotateY(std::array<double 3> vector, double angle) {
		std::array<double, 3> rotated = {vector[0] * cos(angle) + -1 * vector[2] * sin(angle), vector[1] ,vector[0]*sin(angle) + vector[2]*cos(angle) };
		return(rotated);
	}

	std::array<double, 3> RotateY(std::array<double 3> vector, double angle) {
		std::array<double, 3> rotated = { vector[0]*cos(angle)+y*sin(angle), -1*vector[0]*sin(angle)+vector[1]*cos(angle), vector[2] };
		return(rotated);
	}

} // end namespace LeptonInjector

#endif