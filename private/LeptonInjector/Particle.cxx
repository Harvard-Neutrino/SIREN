#include <LeptonInjector/Particle.h>

namespace LI_Particle{
    Particle::Particle(void){
        // Just sit it at the origin 
        
        // note that
        // static_cast<int32_t>(ParticleType::EMinus) == 11
        type        = ParticleType::EMinus;

        energy      = 0.0;
        direction   = {0.0, 0.0};
        position    = {0.0, 0.0, 0.0}; 
    }

} // end namespace LI_Particle
