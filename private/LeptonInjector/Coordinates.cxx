#include<LeptonInjector/Coordinates.h>
#include <array>

namespace LeptonInjector {

	std::array<double, n_dimensions> RotateX(std::array<double, n_dimensions> vector, double angle) {
		std::array<double, n_dimensions> rotated = {vector[0], vector[1]*cos(angle) +vector[2]*sin(angle), -1*vector[1]*sin(angle) + vector[2]*cos(angle) };
		return(rotated);
	}

	std::array<double, n_dimensions> RotateY(std::array<double, n_dimensions> vector, double angle) {
		std::array<double, n_dimensions> rotated = {vector[0] * cos(angle) + -1 * vector[2] * sin(angle), vector[1] ,vector[0]*sin(angle) + vector[2]*cos(angle) };
		return(rotated);
	}

	std::array<double, n_dimensions> RotateZ(std::array<double, n_dimensions> vector, double angle) {
		std::array<double, n_dimensions> rotated = { vector[0]*cos(angle)+vector[1]*sin(angle), -1*vector[0]*sin(angle)+vector[1]*cos(angle), vector[2] };
		return(rotated);
	}



	// create a "direction" object to use with some of the LI dependencies 

	// default constructor, just zeroes 
	LI_Direction::LI_Direction(){
		this->zenith 	= 0.0;
		this->azimuth 	= 0.0;
	}

	// takes two doubles, {zenith, azimuth}
	LI_Direction::LI_Direction( double theta, double phi){
		this->zenith = theta;
		this->azimuth = phi;
		// if the zenith angle is too 
		while( this->azimuth >= 2*Constants::pi ){
			this->azimuth -= 2*Constants::pi;
		}
	}
	// accepts an array {zenith, azimuth}
	LI_Direction::LI_Direction(std::array<double, 2> dir ){
		this->zenith = dir[0];
		this->azimuth = dir[1];
	}

	// can also take a "pair" of doubles
	LI_Direction::LI_Direction( std::pair<double, double> dir){
		this->zenith = dir.first;
		this->azimuth = dir.second;
	}
	LI_Direction::LI_Direction(LI_Direction& old_one){
		this->zenith = old_one.zenith;
		this->azimuth = old_one.azimuth;
	}
	// construct a direction for a given vector 
	LI_Direction::LI_Direction( LI_Position& vec){
		this->azimuth = atan( vec.at(1)/ vec.at(0) );
		this->zenith  = acos( vec.at(2)/ vec.Magnitude() );
	}



	// define the LI_Position constructors and member functions 

	// default constructor sets all n components to zero
	LI_Position::LI_Position(){
		for (uint8_t iter =0; iter<n_dimensions; iter++){
			this->position[iter] = 0.0;
		}
	}
	// uses old position to make a new position
	LI_Position::LI_Position(LI_Position& old_one){
		this->position = old_one.position;
	}
	// uses an array to construct a position 
	LI_Position::LI_Position(std::array<double, n_dimensions> pos){
		for (uint8_t iter =0; iter<n_dimensions; iter++){
			this->position[iter] = pos[iter];
		}
	}
	// returns the value of a specified component 
	double LI_Position::at( uint8_t component){
		// because component is unsigned, the only achievable values are >0 and <255. So, we just check that it's 
		// 		under the dimensionality of the arrays 
		if(component >= n_dimensions){
			// throw an exception. 
			throw std::out_of_range("Invalid component requested");
		}

		return( this->position[component] );
	}

	// returns the magnitude of the position vector
	double LI_Position::Magnitude(){
		double mag = 0.0;
		// note that the pythagorean theorem trivially generalizes to n_dim>2
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			mag += pow( this->at(iter), 2); //sum the squares
		}

		return( sqrt(mag) ); //return the sum's root
	}

	// need to overload some operations on the newly formed LI_Direction and LI_Position

	// this one scales a position vector by a constant
	LI_Position operator * (LI_Position const &point, double scalar){
		std::array<double, n_dimensions> new_one; 
		for (uint8_t iter  = 0; iter<n_dimensions; iter++){
			new_one[iter] = point.at(iter) * scalar;
		}
		return(LI_Position( new_one ));
	}

	// when multiplying a direction by a scalar, you are left with a vector 
	LI_Position operator * (LI_Direction const &dir, double scalar){
		// calculate the coordinates
		double ex = scalar*cos(dir.azimuth)*sin(dir.zenith);
		double why = scalar*sin(dir.zenith)*sin(dir.azimuth);
		double zee = scalar*cos(dir.zenith);
		return( LI_Position(ex, why, zee) );
	}

	// define the dot product between a vector and a direction. The direction is first turned into a unit vector
	double operator * (LI_Position const &pos, LI_Direction const &dir){
		// construct effective position for a unit vector in the direction of dir
		std::array<double,n_dimensions> new_dir = { cos(dir.azimuth)*sin(dir.zenith), sin(dir.azimuth)*sin(dir.zenith), cos(dir.zenith)};

		// with the iterable unit vector, we now compute the inner product. 
		double projected = 0;
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			projected += pos.at(iter)*dir[iter];
		}

		return( projected );
	}

	// similar to above, this calculates the dot product of two position vectors 
	double operator * (LI_Position const &vec1, LI_Position const &vec2){
		double dot_prod = 0;
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			dot_prod += vec1.at(iter)*vec2.at(iter);
		}
		return( dot_prod );
	}

	// implement adding and subtracting positions. Like vector addition! 
	LI_Position operator + (LI_Position const &pos1, LI_Position const &pos2){
		std::array<double, n_dimensions> new_one;
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			new_one[iter] = pos1.at(iter) + pos2.at(iter);
		}
		return( LI_Position( new_one ) );
	}
	LI_Position operator - (LI_Position const &pos1, LI_Position const &pos2){
		std::array<double, n_dimensions> new_one;

		for (uint8_t iter=0; iter<n_dimensions; iter++){
			new_one[iter] = pos1.at(iter) - pos2.at(iter);
		}

		return( LI_Position( new_one ) );
	}

	LI_Position& operator += (LI_Position& one, LI_Position& two){
		one = one + two;
		return( one );
	}
	LI_Position& operator += (LI_Position& one, LI_Position& two){
		one = one - two;
		return( one );
	}

} // end namespace LeptonInjector