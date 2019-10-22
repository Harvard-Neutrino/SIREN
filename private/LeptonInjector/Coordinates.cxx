#include<LeptonInjector/Helper.h>
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

	std::array<double, n_dimensions> RotateY(std::array<double, n_dimensions> vector, double angle) {
		std::array<double, n_dimensions> rotated = { vector[0]*cos(angle)+y*sin(angle), -1*vector[0]*sin(angle)+vector[1]*cos(angle), vector[2] };
		return(rotated);
	}



	// create a "direction" object to use with some of the LI dependencies 
	LI_Direction::LI_Direction(){
		this->zenith 	= 0.0;
		this->azimuth 	= 0.0;
	}

	LI_Direction::~LI_Direction(){
		// ... ?
	}

	LI_Direction::LI_Direction( double theta, double phi){
		this->zenith = theta;
		this->azimuth = phi;

		// if the zenith angle is too 
		while( this->azimuth >= 2*Constants::pi ){
			this->azimuth -= 2*Constants::pi;
		}

	}

	LI_Direction::LI_Direction( std::pair<double, double> dir){
		this->zenith = dir.first;
		this->azimuth = dir.second;
	}



	// create a 
	// define the LI_Position constructors 
	LI_Position::LI_Position(){
		for (uint8_t iter =0; iter<n_dimensions; iter++){
			this->position[iter] = 0.0;
		}
	}

	LI_Position::~LI_Position(){
	}

	LI_Position::LI_Position(LI_Position& old_one){
		this->position = old_one.position;
	}

	LI_Position::LI_Position(std::array<double, n_dimensions> pos){
		for (uint8_t iter =0; iter<n_dimensions; iter++){
			this->position[iter] = pos[iter];
		}
	}

	double LI_Position::at( uint8_t component){
		return( this->position[component] );
	}

	double LI_Position::Magnitude(){
		double mag = 0.0;
		for (uint8_t iter=0; iter<n_dimensions; iter++){
			mag += pow( this->at(iter), 2);
		}

		return( sqrt(mag) );
	}

	// need to overload some operations on the newly formed LI_Direction and LI_Position

	// takes a LI_Position (basically a vector), and scales each component by a constant. 
	LI_Position operator * (LI_Position const &point, double scalar){
		std::array<double, n_dimensions> new_one; 
		for (uint8_t iter  = 0; iter<n_dimensions; iter++){
			new_one[iter] = point.at(iter) * scalar;
		}
		return(LI_Position( new_one ));
	}

	// define the dot product between a vector and a direction. 
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

} // end namespace LeptonInjector