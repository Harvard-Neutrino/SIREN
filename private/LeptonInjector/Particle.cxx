#include <math.h> // adds sqrt, power functions
#include <Particle.h>
#include <assert.h>
#include <map>

namespace LeptonInjector{
    Particle::Particle(void){
        // Just sit it at the origin 
        
        // note that
        // static_cast<int32_t>(ParticleType::EMinus) == 11
        type        = ParticleType::unknown;
        
        // instantiate with minimum energy
        energy      = 0.0;
        direction   = std::make_pair( 0.0, 0.0) ;
		for (uint8_t var = 0; var<3; var++){
			position[var] = 0.0;
		}
    }
    
    // constructor for specific type
    Particle::Particle( ParticleType type ){
        this->type  = type;

        // energy will instead be built using particle's mass
        energy      = this->GetMass();
        direction   = std::make_pair(0.0, 0.0);
		for (uint8_t var = 0; var<3; var++){
			position[var] = 0.0;
		}
    }



	// returns name for particle of known type. 
	// If this code is to be expanded, this should really be modified to use the boost preprocessor libraries
	// atm, only implemented for the particles relevant to LeptonInjector 
	std::string Particle::GetTypeString(){

		// this is **BAD** and I should feel bad
		// there is a way to do this better with boost preprocessor libraries, but I think that's a little unnecessary given the scope of what LI does. 

		// this just casts the particle type to its pdg code, and uses a switch to grab the name
		switch( static_cast<int32_t>(this->type) ){
			case 0: return("Unknwon"); break;
			case 22: return("Gamma"); break;
			case 11: return("EMinus"); break;
			case -11: return("EPlus"); break;
			case 13: return("MuMinus"); break;
			case -13: return("MuPlus"); break;
			case 15: return("TauMinus"); break;
			case -15: return("TauPlus"); break;
			case 12: return("NuE"); break;
			case -12: return("NuEBar"); break;
			case 14: return("NuMu"); break;
			case -14: return("NuMuBar"); break;
			case 16: return("NuTau"); break;
			case -16: return("NuTauBar"); break;
			case -2000001006: return("Hadrons"); break;
			default: return("Unsupported"); break;
		}

	}

    bool Particle::HasMass(){
        // return the negation of the bool that (particle is massless)
        return(!( this->type == ParticleType::Gamma || 
               this->type == ParticleType::NuE   || this->type==ParticleType::NuEBar   ||
               this->type == ParticleType::NuMu  || this->type==ParticleType::NuMuBar  ||
               this->type == ParticleType::NuTau || this->type==ParticleType::NuTauBar ||
			   this->type == ParticleType::PPlus || this->type==ParticleType::Neutron) );
    }

    // only implemented for the charged leptons to stay within scope
    double Particle::GetMass(){
        switch(this->type){
            case ParticleType::EPlus:
                return( Constants::electronMass );
                break;
            case ParticleType::EMinus:
                return( Constants::electronMass );
                break;
            case ParticleType::MuPlus:
                return( Constants::muonMass );
                break;
            case ParticleType::MuMinus:
                return( Constants::muonMass );
                break;
            case ParticleType::TauPlus:
                return( Constants::tauMass );
                break;
            case ParticleType::TauMinus:
                return( Constants::tauMass );
                break;
			case ParticleType::PPlus:
				return( Constants::protonMass );
				break;
			case ParticleType::Neutron:
				return( Constants::neutronMass);
            default:
                return(0.0);
        }
    }



    // Helper functions for dealing with particle types  

    // returns true if a particle is a Lepton. False if not
    bool isLepton(Particle::ParticleType p){
		return(p==Particle::ParticleType::EMinus   || p==Particle::ParticleType::EPlus ||
			   p==Particle::ParticleType::MuMinus  || p==Particle::ParticleType::MuPlus ||
			   p==Particle::ParticleType::TauMinus || p==Particle::ParticleType::TauPlus ||
			   p==Particle::ParticleType::NuE      || p==Particle::ParticleType::NuEBar ||
			   p==Particle::ParticleType::NuMu     || p==Particle::ParticleType::NuMuBar ||
			   p==Particle::ParticleType::NuTau    || p==Particle::ParticleType::NuTauBar);
	}
	
    // returns true if the particle is either
    //        a charged lepton 
    //   (OR) a "hadrons" particle
	// If passed a disallowed particle, throws a tempter tantrum 
	bool isCharged(Particle::ParticleType p){
		if( !(isLepton(p) || p==Particle::ParticleType::Hadrons) ){
			throw "You should only be using Leptons or Hadrons!";
		}
		
		// keeps this within scope. Shouldn't be getting some other kind of charged particle
		return(p==Particle::ParticleType::EMinus   || p==Particle::ParticleType::EPlus ||
			   p==Particle::ParticleType::MuMinus  || p==Particle::ParticleType::MuPlus ||
			   p==Particle::ParticleType::TauMinus || p==Particle::ParticleType::TauPlus ||
			   p==Particle::ParticleType::Hadrons);
	}


    // returns string of particle's name
	std::string particleName(Particle::ParticleType p){
		return(Particle(p).GetTypeString());
	}
	

    // gets the mass of a particle for a given type
	double particleMass(Particle::ParticleType type){
		Particle p(type);
		if(!p.HasMass()){
            //  removed until new logging system implemented 
			//log_debug_stream("Treating particle of type " << p.GetTypeString()
			//				 << " with unknown mass as massless");
			return(0);
		}
		return(p.GetMass());
	}
	
    // Uses a particle's type (mass) and total energy to calculate kinetic energy
	double kineticEnergy(Particle::ParticleType type, double totalEnergy){
		double mass=particleMass(type);
		if(totalEnergy<mass){
            // commented out until a new logging system is implemented 
//			log_warn_stream("Negative kinetic energy (particle type = " << particleName(type) << ", mass = " << mass << ", total energy = " << totalEnergy << ')');
			return(0.);
		}
		return(sqrt(totalEnergy*totalEnergy-mass*mass));
	}
	
    // uses the particle type and kinetic energy to calculate the speed of the particle
    // relies on the constants! 
	double particleSpeed(Particle::ParticleType type, double kineticEnergy){
		Particle p=Particle(type);
		if(!p.HasMass()){
            // removing this until a new logging system is implemented... 
//			log_debug_stream("Treating particle of type " << p.GetTypeString()
//							 << " with unknown mass as massless");
			return(Constants::c);
		}
		double mass=p.GetMass();
		if(kineticEnergy<0){
            // same as always 
//			log_warn("Negative kinetic energy");
			return(0.);
		}

        // these always confuse me, so I'm leaving a comment
        // IF mass>0 THEN return mass/(stuff) ... ELSE return 0
		double r=(mass>0 ? mass/(kineticEnergy+mass) : 0.);
		return(Constants::c*sqrt(1-r*r));
	}
	
	Particle::ParticleShape decideShape(Particle::ParticleType t){
		switch(t){
			case Particle::ParticleType::MuMinus:  case Particle::ParticleType::MuPlus:
			case Particle::ParticleType::TauMinus: case Particle::ParticleType::TauPlus:
			case Particle::ParticleType::NuE:      case Particle::ParticleType::NuEBar:
			case Particle::ParticleType::NuMu:     case Particle::ParticleType::NuMuBar:
			case Particle::ParticleType::NuTau:    case Particle::ParticleType::NuTauBar:
				return(Particle::ParticleShape::MCTrack);
			case Particle::ParticleType::EMinus: case Particle::ParticleType::EPlus:
			case Particle::ParticleType::Hadrons:
				return(Particle::ParticleShape::Cascade);
			case Particle::ParticleType::unknown:
				return(Particle::ParticleShape::unknown);
			default:
                throw "BadShape"; // this replaces the previous fatal log
//				log_fatal_stream("Unable to decide shape for unexpected particle type: " << particleName(t));
		}
	}

    // This function returns the primary particle type given the final state particles
    // returns a particle type object    
	Particle::ParticleType deduceInitialType(Particle::ParticleType pType1, Particle::ParticleType pType2){
		//only accept certain particle types in general
		if(!isLepton(pType1) && pType1!=Particle::ParticleType::Hadrons)
            throw "BadParticle"; //replace log
//			log_fatal_stream("Unexpected particle type: "
//							 << particleName(pType1)
//							 << ";\nonly leptons and 'Hadrons' are supported");
		if(!isLepton(pType2) && pType2!=Particle::ParticleType::Hadrons)
            throw "BadParticle";
//			log_fatal_stream("Unexpected particle type: "
//							 << particleName(pType2)
//							 << ";\nonly leptons and 'Hadrons' are supported");
		
		bool c1=isCharged(pType1);
		bool c2=isCharged(pType2);
		bool l1=isLepton(pType1);
		bool l2=isLepton(pType2);
		
		//at least one particle should be charged
		if(!c1 && !c2)
			throw "Final state should have at least one charged particle";
//			log_fatal_stream("Final state must contain at least one charged particle\n"
//							 << "specified particles were " << particleName(pType1)
//							 << " and " << particleName(pType2)); 
		
		//first particle is charged, second is not
		if(c1 && !c2){
			//valid cases are charged lepton + matching antineutrino for GR
			if(l1){
				//!c2 => pType2 is a neutrino
				if(!((pType1==Particle::ParticleType::EMinus   && pType2==Particle::ParticleType::NuEBar) ||
					 (pType1==Particle::ParticleType::EPlus    && pType2==Particle::ParticleType::NuE) ||
					 (pType1==Particle::ParticleType::MuMinus  && pType2==Particle::ParticleType::NuMuBar) ||
					 (pType1==Particle::ParticleType::MuPlus   && pType2==Particle::ParticleType::NuMu) ||
					 (pType1==Particle::ParticleType::TauMinus && pType2==Particle::ParticleType::NuTauBar) ||
					 (pType1==Particle::ParticleType::TauPlus  && pType2==Particle::ParticleType::NuTau)))
                     throw "Final states with a charged lepton must have an anti-matching neutrino.";
//    	  		     log_fatal_stream("Final states with a charged lepton must have an anti-matching neutrino.\n"
//									 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
				//log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Glashow Resonance (leptonic)");
				return(Particle::ParticleType::NuEBar);
			}
            throw "BadFinal";
//			log_fatal_stream("Unrecognized final state type: " << particleName(pType1) << " and " << particleName(pType2));
		}
		
		//first particle is neutral, second is charged
		if(!c1 && c2){
			if(l1 && pType2==Particle::ParticleType::Hadrons){
				//particle 1 is a neutral lepton, so it must be a neutrino
//				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Neutral Current");
				return(pType1); //the incoming neutrino type is the same as the outgoing
			}
            throw "BadFinal";
//			log_fatal_stream("Unrecognized final state type: " << particleName(pType1) << " and " << particleName(pType2));
		}
		
		//have two charged particles
		if(c1 && c2){
			//no two charged lepton states
			if(l1 && l2)
                throw "BadFinal";
//				log_fatal_stream("Two charged lepton final states are not allowed.\n"
//								 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
			//lepton should be given first
			if(!l1 && l2)
                throw "BadFinal";
//				log_fatal_stream("Final states should specify charged leptons before 'Hadrons'.\n"
//								 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
			
			if(l1 && !l2){ //valid: charged lepton + Hadrons for CC
//				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Charged Current");
				switch(pType1){
					case Particle::ParticleType::EMinus: return(Particle::ParticleType::NuE);
					case Particle::ParticleType::EPlus: return(Particle::ParticleType::NuEBar);
					case Particle::ParticleType::MuMinus: return(Particle::ParticleType::NuMu);
					case Particle::ParticleType::MuPlus: return(Particle::ParticleType::NuMuBar);
					case Particle::ParticleType::TauMinus: return(Particle::ParticleType::NuTau);
					case Particle::ParticleType::TauPlus: return(Particle::ParticleType::NuTauBar);
					default: assert(false && "This point should be unreachable");
				}
			}
			if(!l1 && !l2){ //valid: two hadrons (for GR)
//				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Glashow Resonance (hadronic)");
				return(Particle::ParticleType::NuEBar);
			}
		}
        throw "You must be a wizard: this point should be unreachable";
        //log_fatal("Logic error; this point should be unreachable");
	}


    // Particle-based exceptions:

} // end namespace LI_Particle
