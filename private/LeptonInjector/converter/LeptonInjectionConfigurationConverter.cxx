#include "LeptonInjectionConfigurationConverter.h"

namespace LeptonInjector{

I3TableRowDescriptionPtr EventPropertiesConverter::CreateDescription(const BasicEventProperties& p){
	I3TableRowDescriptionPtr desc(new I3TableRowDescription());
	desc->AddField<double>("totalEnergy", "GeV", "Total energy of the event");
	desc->AddField<double>("zenith", "radians", "Zenith angle of the injected primary particle");
	desc->AddField<double>("azimuth", "radians", "Azimuth angle of the injected primary particle");
	desc->AddField<double>("finalStateX", "unitless", "Feynman X of the interaction");
	desc->AddField<double>("finalStateY", "unitless", "Feynman Y of the interaction");
	MAKE_ENUM_VECTOR(type,I3Particle,I3Particle::ParticleType,I3PARTICLE_H_I3Particle_ParticleType);
	desc->AddEnumField<I3Particle::ParticleType> ("finalType1",type,"","Type of the first particle in the final state");
	desc->AddEnumField<I3Particle::ParticleType> ("finalType2",type,"","Type of the second particle in the final state");
	desc->AddEnumField<I3Particle::ParticleType> ("initialType",type,"","Type of initial state particle");
	
	desc->AddField<double>("impactParameter", "meters", "The impact parameter of the primary particle path with respect to the origin");
	desc->AddField<double>("totalColumnDepth", "g/cm^2", "The total column depth along the particle path within which the interaction was sampled");
	
	desc->AddField<double>("radius", "meters", "The sampled radial cylindrical coordinate of the interaction point");
	desc->AddField<double>("z", "meters", "The sampled vertical cylindrical coordinate of the interaction point");
	return(desc);
}

size_t EventPropertiesConverter::FillRows(const BasicEventProperties& p, I3TableRowPtr rows){
	rows->Set<double>("totalEnergy", p.totalEnergy);
	rows->Set<double>("zenith", p.zenith);
	rows->Set<double>("azimuth", p.azimuth);
	rows->Set<double>("finalStateX", p.finalStateX);
	rows->Set<double>("finalStateY", p.finalStateY);
	rows->Set<I3Particle::ParticleType>("finalType1", p.finalType1);
	rows->Set<I3Particle::ParticleType>("finalType2", p.finalType2);
	rows->Set<I3Particle::ParticleType>("initialType", p.initialType);
	
	const RangedEventProperties* rprop;
	const VolumeEventProperties* vprop;
	if((rprop=dynamic_cast<const RangedEventProperties*>(&p))){
		rows->Set<double>("impactParameter", rprop->impactParameter);
		rows->Set<double>("totalColumnDepth", rprop->totalColumnDepth);
		rows->Set<double>("radius", std::numeric_limits<double>::quiet_NaN());
		rows->Set<double>("z", std::numeric_limits<double>::quiet_NaN());
	}
	else if((vprop=dynamic_cast<const VolumeEventProperties*>(&p))){
		rows->Set<double>("radius", vprop->radius);
		rows->Set<double>("z", vprop->z);
		rows->Set<double>("impactParameter", std::numeric_limits<double>::quiet_NaN());
		rows->Set<double>("totalColumnDepth", std::numeric_limits<double>::quiet_NaN());
	}
	else
		log_error("Unrecognized sbtype of BasicEventProperties; tabulated"
		          " information is probably incomplete");
	return(1);
}

} //namespace LeptonInjector
