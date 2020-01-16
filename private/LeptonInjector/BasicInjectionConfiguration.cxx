
#include <BasicInjectionConfiguration.h>

namespace LeptonInjector{

	BasicInjectionConfiguration::BasicInjectionConfiguration():
	events(1),
	energyMinimum(10*Constants::GeV),
	energyMaximum((1e9)*Constants::GeV),
	powerlawIndex(1.0),
	azimuthMinimum(0),
	azimuthMaximum(2*Constants::pi),
	zenithMinimum(0),
	zenithMaximum(Constants::pi),
	finalType1(Particle::ParticleType::MuMinus),
	finalType2(Particle::ParticleType::Hadrons),
	injectionRadius(1200*LeptonInjector::Constants::m),
	endcapLength(1200*LeptonInjector::Constants::m),
	cylinderRadius(1200*LeptonInjector::Constants::m),
	cylinderHeight(1200*LeptonInjector::Constants::m)
	{}

} // end namespace LeptonInjector