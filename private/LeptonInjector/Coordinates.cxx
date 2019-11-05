#include "Coordinates.h"
#include <array>

namespace LeptonInjector {

	LI_Position RotateX(LI_Position vector, double angle) {
		std::array<double, n_dimensions> rotated = {vector.at(0), vector.at(1)*cos(angle) +vector.at(2)*sin(angle), -1*vector.at(1)*sin(angle) + vector.at(2)*cos(angle) };
		return(LI_Position(rotated));
	}

	LI_Position RotateY(LI_Position vector, double angle) {
		std::array<double, n_dimensions> rotated = {vector.at(0) * cos(angle) + -1 * vector.at(2) * sin(angle), vector.at(1) ,vector.at(0)*sin(angle) + vector.at(2)*cos(angle) };
		return(LI_Position(rotated));
	}

	LI_Position RotateZ(LI_Position vector, double angle) {
		std::array<double, n_dimensions> rotated = { vector.at(0)*cos(angle)+vector.at(1)*sin(angle), -1*vector.at(0)*sin(angle)+vector.at(1)*cos(angle), vector.at(2) };
		return(LI_Position(rotated));
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
	LI_Direction::LI_Direction(LI_Direction* old_one){
		this->zenith = old_one->zenith;
		this->azimuth = old_one->azimuth;
	}

	// construct a direction for a given vector 
	LI_Direction::LI_Direction( LI_Position vec){
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
	LI_Position::LI_Position(const LI_Position&  old_one){
		this->position = old_one.position;
	}
	// uses an array to construct a position 
	LI_Position::LI_Position(std::array<double, n_dimensions> pos){
		for (uint8_t iter =0; iter<n_dimensions; iter++){
			this->position[iter] = pos[iter];
		}
	}
	LI_Position::LI_Position( double ex, double why, double z){
		this->position[0] = ex;
		this->position[1] = why;
		this->position[2] = z;
	}

	// returns the value of a specified component 
	double LI_Position::at( uint8_t component) const{
		// because component is unsigned, the only achievable values are >0 and <255. So, we just check that it's 
		// 		under the dimensionality of the arrays 
		if(component >= n_dimensions){
			// throw an exception. 
			throw std::out_of_range("Invalid component requested");
		}

		return( this->position[component] );
	}

	double LI_Position::GetX() const{return( position[0]); }
	double LI_Position::GetY() const{return( position[1]); }
	double LI_Position::GetZ() const{return( position[2]); }

	void LI_Position::SetX(double amt){this->position[0]= amt; }
	void LI_Position::SetY(double amt){this->position[1]= amt; }
	void LI_Position::SetZ(double amt){this->position[2]= amt; }


	// returns the magnitude of the position vector
	double LI_Position::Magnitude(void) const{
		double mag = 0.0;
		// note that the pythagorean theorem trivially generalizes to n_dim>2
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			mag += pow( this->at(iter), 2); //sum the squares
		}

		return( sqrt(mag) ); //return the sum's root
	}

	

	// need to overload some operations on the newly formed LI_Direction and LI_Position

	// this one scales a position vector by a constant
	LI_Position operator * ( LI_Position point, double scalar){
		std::array<double, n_dimensions> new_one; 
		for (uint8_t iter  = 0; iter<n_dimensions; iter++){
			new_one[iter] = point.at(iter) * scalar;
		}
		return(LI_Position( new_one ));
	} // and the commutation!
	LI_Position operator * (double scalar, LI_Position const point){
		return( point*scalar );
	}

	// when multiplying a direction by a scalar, you are left with a vector 
	LI_Position operator * (LI_Direction const dir, double scalar){
		// calculate the coordinates
		double ex = scalar*cos(dir.azimuth)*sin(dir.zenith);
		double why = scalar*sin(dir.zenith)*sin(dir.azimuth);
		double zee = scalar*cos(dir.zenith);
		return( LI_Position(ex, why, zee) );
	} // commutation
	LI_Position operator * (double scalar, LI_Direction const dir){
		return( dir*scalar );
	}

	// define the dot product between a vector and a direction. The direction is first turned into a unit vector
	double operator * (LI_Position  pos, LI_Direction  dir){
		// construct effective position for a unit vector in the direction of dir
		std::array<double,n_dimensions> new_dir = { cos(dir.azimuth)*sin(dir.zenith), sin(dir.azimuth)*sin(dir.zenith), cos(dir.zenith)};

		// with the iterable unit vector, we now compute the inner product. 
		double projected = 0;
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			projected += pos.at(iter)*new_dir[iter];
		}
		return( projected );
	}// commutation
	double operator * (LI_Direction  dir, LI_Position  pos){
		return( pos*dir );
	}

	// similar to above, this calculates the dot product of two position vectors 
	double operator * (LI_Position  vec1, LI_Position  vec2){
		double dot_prod = 0;
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			dot_prod += vec1.at(iter)*vec2.at(iter);
		}
		return( dot_prod );
	}//commutation is implicitly already here... 

	// implement adding and subtracting positions. Like vector addition! 
	LI_Position operator + (LI_Position pos1, LI_Position pos2){
		std::array<double, n_dimensions> new_one;
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			new_one[iter] = pos1.at(iter) + pos2.at(iter);
		}
		return( LI_Position( new_one ) );
	}// commutation is implicitly implemented
	LI_Position operator - (LI_Position pos1, LI_Position pos2){
		std::array<double, n_dimensions> new_one;

		for (uint8_t iter=0; iter<n_dimensions; iter++){
			new_one[iter] = pos1.at(iter) - pos2.at(iter);
		}
		return( LI_Position( new_one ) );
	}// implicitly blah blah blah

	LI_Position& operator += (LI_Position one, LI_Position two){
		one = one + two;
		return( one );
	}
	LI_Position& operator -= (LI_Position one, LI_Position two){
		one = one - two;
		return( one );
	}// same for those last two, too

	// turns a direction around
	LI_Direction operator - (LI_Direction obj){
		LI_Direction new_one;
		new_one.zenith = Constants::pi - obj.zenith;
		new_one.azimuth += Constants::pi; 

		if (obj.azimuth >= 2*Constants::pi){
			obj.azimuth -= 2*Constants::pi;
		}

		return(new_one);
	}


	// check if two points are identical
	// 		maybe use a different epsilon? Not sure. 
	bool operator == (LI_Position& one, LI_Position& two){
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			// check if the difference between each comonent isless than the minimum expressible distance between doubles
			//		this, as opposed to using '==' is to avoid floating point errors 
			if ( !(abs(one.at(iter)-two.at(iter)) <= std::numeric_limits<double>::epsilon()) ){
				return(false);
			}
		}
		return(true);
	}

	LI_Direction rotateRelative(LI_Direction base, double zenith, double azimuth){
		LI_Direction result;
		result.zenith += zenith*cos(azimuth);
		result.azimuth+= zenith*sin(azimuth);
		return(result);
	}

} // end namespace LeptonInjector
