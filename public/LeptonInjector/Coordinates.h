#ifndef LI_COORDS
#define LI_COORDS

#include <array> // used for many of the constructors, and for the position
#include <math.h> //sqrt, sin, cos, pow
#include <exception> //allows throwing the out_of_range exception
#include <LeptonInjector/Constants.h>

// Implements tools and classes for working with positions and directions

// Ben Smithers
// benjamin.smithers@mavs.uta.edu

namespace LeptonInjector {

	// how many linearly independent bases do we need to span our space ?
	// this should always be 3
	// but... maybe one day we want more dimensions, and don't want to have to rewrite all of this? 
	static const int n_dimensions = 3;


	// these functions rotate a 3-vector about some axis
	// assumes vector is in cartesian coordinates
	std::array<double, n_dimensions> RotateY(std::array<double, n_dimensions> vector, double angle);
	std::array<double, n_dimensions> RotateX(std::array<double, n_dimensions> vector, double angle);
	std::array<double, n_dimensions> RotateZ(std::array<double, n_dimensions> vector, double angle);

	// Creating a "LI_Direction" to mimic the I3_Direction object
	// this should behave in the same exact way, at least within the scope of this project 
	class LI_Direction{
		public:
			LI_Direction();
			virtual ~LI_Direction();

			LI_Direction( double theta, double phi);
			LI_Direction( std::array<double, 2> dir);
			LI_Direction( std::pair<double, double> dir);
			LI_Direction( LI_Direction& old_one);
			LI_Direction( LI_Position& vec);  // get the direction of a vector 

			
			double zenith;
			double azimuth;
	};

	// Creating a "LI_Position" to mimic the I3_Position object
	// this is more accurately described as a vector in cartesian coordinates
	//
	// **IMPORTANT** : member functions and operations are defined using **cartesian** coordinates.
	// 					keep this in mind when working with these coordinates! 
	class LI_Position{
		public: 
			LI_Position();
			virtual ~LI_Position();

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
