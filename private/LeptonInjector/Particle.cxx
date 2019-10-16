#include <math.h> // adds sqrt, power functions
#include <LeptonInjector/Particle.h>

namespace LI_Particle{
    Particle::Particle(void){
        // Just sit it at the origin 
        
        // note that
        // static_cast<int32_t>(ParticleType::EMinus) == 11
        type        = ParticleType::EMinus;
        
        // instantiate with minimum energy
        energy      = Constants::electronMass;
        direction   = {0.0, 0.0};
        position    = {0.0, 0.0, 0.0}; 
    }
    
    // constructor for specific type
    Particle::Particle( ParticleType type ){
        this->type  = type;

        // energy will instead be built using particle's mass
        energy      = this->GetMass();
        direction   = {0.0, 0.0};
        position    = {0.0, 0.0, 0.0}; 
    }

    Particle::HasMass(){
        // return the negation of the bool that (particle is massless)
        return(!( this->type == ParticleType::Gamma || 
               this->type == ParticleType::NuE   || this->type==ParticleType::NuEBar   ||
               this->type == ParticleType::NuMu  || this->type==ParticleType::NuMuBar  ||
               this->type == ParticleType::NuTau || this->type==ParticleType::NuTauBar) );
    }

    // only implemented for the charged leptons to stay within scope
    Particle::GetMass(){
        switch(this->ParticleType){
            case ParticleType::EPlus:
                return( Constants::electronMass );
                break;
            case ParticleType::EMinus:
                return( Constants::eletronMass );
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
            default:
                return(0.0);
        }
    }

    ~Particle::Particle(){}


    // Helper functions for dealing with particle types  

    // returns true if a particle is a Lepton. False if not
    bool isLepton(ParticleType p){
		return(p==ParticleType::EMinus   || p==ParticleType::EPlus ||
			   p==ParticleType::MuMinus  || p==ParticleType::MuPlus ||
			   p==ParticleType::TauMinus || p==ParticleType::TauPlus ||
			   p==ParticleType::NuE      || p==ParticleType::NuEBar ||
			   p==ParticleType::NuMu     || p==ParticleType::NuMuBar ||
			   p==ParticleType::NuTau    || p==ParticleType::NuTauBar);
	}
	
    // returns true if the particle is either
    //        a charged lepton 
    //   (OR) a "hadrons" particle
	bool isCharged(ParticleType p){
		assert(isLepton(p) || p==ParticleType::Hadrons); // keeps this within scope. Shouldn't be getting some other kind of charged particle
		return(p==ParticleType::EMinus   || p==ParticleType::EPlus ||
			   p==ParticleType::MuMinus  || p==ParticleType::MuPlus ||
			   p==ParticleType::TauMinus || p==ParticleType::TauPlus ||
			   p==ParticleType::Hadrons);
	}


    // returns string of particle's name
	std::string particleName(ParticleType p){
		return(Particle(p).GetTypeString());
	}
	

    // gets the mass of a particle for a given type
	double particleMass(ParticleType type){
		Particle p(type);
		if(!p.HasMass()){
            // VV removed until new logging system implemented VV 
			//log_debug_stream("Treating particle of type " << p.GetTypeString()
			//				 << " with unknown mass as massless");
			return(0);
		}
		return(p.GetMass());
	}
	
    // Uses a particle's type (mass) and total energy to calculate kinetic energy
	double kineticEnergy(ParticleType type, double totalEnergy){
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
	double particleSpeed(ParticleType type, double kineticEnergy){
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
	
	ParticleShape decideShape(ParticleType t){
		switch(t){
			case Particle::MuMinus:  case Particle::MuPlus:
			case Particle::TauMinus: case Particle::TauPlus:
			case Particle::NuE:      case Particle::NuEBar:
			case Particle::NuMu:     case Particle::NuMuBar:
			case Particle::NuTau:    case Particle::NuTauBar:
				return(ParticleShape::MCTrack);
			case Particle::EMinus: case Particle::EPlus:
			case Particle::Hadrons:
				return(ParticleShape::Cascade);
			default:
                throw BadShape; // this replaces the previous fatal log
//				log_fatal_stream("Unable to decide shape for unexpected particle type: " << particleName(t));
		}
	}

    // This function returns the primary particle type given the final state particles
    // returns a particle type object    
	ParticleType deduceInitialType(ParticleType pType1, ParticleType pType2){
		//only accept certain particle types in general
		if(!isLepton(pType1) && pType1!=ParticleType::Hadrons)
            throw BadParticle; //replace log
//			log_fatal_stream("Unexpected particle type: "
//							 << particleName(pType1)
//							 << ";\nonly leptons and 'Hadrons' are supported");
		if(!isLepton(pType2) && pType2!=ParticleType::Hadrons)
            throw BadParticle;
//			log_fatal_stream("Unexpected particle type: "
//							 << particleName(pType2)
//							 << ";\nonly leptons and 'Hadrons' are supported");
		
		bool c1=isCharged(pType1);
		bool c2=isCharged(pType2);
		bool l1=isLepton(pType1);
		bool l2=isLepton(pType2);
		
		//at least one particle should be charged
		if(!c1 && !c2)
			log_fatal_stream("Final state must contain at least one charged particle\n"
							 << "specified particles were " << particleName(pType1)
							 << " and " << particleName(pType2));
		
		//first particle is charged, second is not
		if(c1 && !c2){
			//valid cases are charged lepton + matching antineutrino for GR
			if(l1){
				//!c2 => pType2 is a neutrino
				if(!((pType1==ParticleType::EMinus   && pType2==ParticleType::NuEBar) ||
					 (pType1==ParticleType::EPlus    && pType2==ParticleType::NuE) ||
					 (pType1==ParticleType::MuMinus  && pType2==ParticleType::NuMuBar) ||
					 (pType1==ParticleType::MuPlus   && pType2==ParticleType::NuMu) ||
					 (pType1==ParticleType::TauMinus && pType2==ParticleType::NuTauBar) ||
					 (pType1==ParticleType::TauPlus  && pType2==ParticleType::NuTau)))
                     throw "Final states with a charged lepton must have an anti-matching neutrino.";
//    	  		     log_fatal_stream("Final states with a charged lepton must have an anti-matching neutrino.\n"
//									 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
				//log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Glashow Resonance (leptonic)");
				return(ParticleType::NuEBar);
			}
            throw BadFinal;
//			log_fatal_stream("Unrecognized final state type: " << particleName(pType1) << " and " << particleName(pType2));
		}
		
		//first particle is neutral, second is charged
		if(!c1 && c2){
			if(l1 && pType2==ParticleType::Hadrons){
				//particle 1 is a neutral lepton, so it must be a neutrino
//				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Neutral Current");
				return(pType1); //the incoming neutrino type is the same as the outgoing
			}
            throw BadFinal;
//			log_fatal_stream("Unrecognized final state type: " << particleName(pType1) << " and " << particleName(pType2));
		}
		
		//have two charged particles
		if(c1 && c2){
			//no two charged lepton states
			if(l1 && l2)
                throw BadFinal;
//				log_fatal_stream("Two charged lepton final states are not allowed.\n"
//								 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
			//lepton should be given first
			if(!l1 && l2)
                throw BadFinal;
//				log_fatal_stream("Final states should specify charged leptons before 'Hadrons'.\n"
//								 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
			
			if(l1 && !l2){ //valid: charged lepton + Hadrons for CC
//				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Charged Current");
				switch(pType1){
					case ParticleType::EMinus: return(ParticleType::NuE);
					case ParticleType::EPlus: return(ParticleType::NuEBar);
					case ParticleType::MuMinus: return(ParticleType::NuMu);
					case ParticleType::MuPlus: return(ParticleType::NuMuBar);
					case ParticleType::TauMinus: return(ParticleType::NuTau);
					case ParticleType::TauPlus: return(ParticleType::NuTauBar);
					default: assert(false && "This point should be unreachable");
				}
			}
			if(!l1 && !l2){ //valid: two hadrons (for GR)
//				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Glashow Resonance (hadronic)");
				return(ParticleType::NuEBar);
			}
		}
        throw "You must be a wizard: this point should be unreachable"i;
//log_fatal("Logic error; this point should be unreachable");
	}


} // end namespace LI_Particle
