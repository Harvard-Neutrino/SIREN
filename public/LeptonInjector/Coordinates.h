#ifndef LI_HELPER
#define LI_HELPER

#include <array>
#include <math.h> //sqrt, sin, cos, pow
#include <LeptonInjector/Constants.h>


namespace LeptonInjector {

	// how many linearly independent bases do we need to span our space ?
	// this should always be 3
	static const int n_dimensions = 3;


	// these functions rotate a 3-vector about some axis
	// assumes vector is in cartesian coordinates
	std::array<double, n_dimensions> RotateY(std::array<double, n_dimensions> vector, double angle);
	std::array<double, n_dimensions> RotateX(std::array<double, n_dimensions> vector, double angle);
	std::array<double, n_dimensions> RotateZ(std::array<double, n_dimensions> vector, double angle);

	class LI_Direction{
		public:
			LI_Direction();
			~LI_Direction();

			LI_Direction( double theta, double phi);
			LI_Direction( std::pair<double, double> dir);

			
			double zenith;
			double azimuth;
	};

	class LI_Position{
		public: 
			LI_Position();
			~LI_Position();

			LI_Position(double ex, double why, double zee);
			LI_Position(LI_Position& old_one);
			LI_Position(std::array<double, n_dimensions> pos);

			double at(uint8_t component);

			double Magnitude();

		private:
			std::array<double, n_dimensions> position;
	};

} // end namespace LeptonInjector

#endif
