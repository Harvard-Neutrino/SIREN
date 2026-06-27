#include "SIREN/interactions/DMesonELoss.h"

#include <map>                                             // for map, opera...
#include <set>                                             // for set, opera...
#include <array>                                           // for array
#include <cmath>                                           // for pow, log10
#include <tuple>                                           // for tie, opera...
#include <memory>                                          // for allocator
#include <string>                                          // for basic_string
#include <vector>                                          // for vector
#include <assert.h>                                        // for assert
#include <stddef.h>                                        // for size_t

#include <rk/rk.hh>                                        // for P4, Boost
#include <rk/geom3.hh>                                     // for Vector3

#include "SIREN/interactions/CrossSection.h"     // for CrossSection
#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/dataclasses/Particle.h"           // for Particle
#include "SIREN/utilities/Random.h"               // for SIREN_random
#include "SIREN/utilities/Constants.h"            // for electronMass
#include "SIREN/utilities/Integration.h"          // for rombergInt...
#include "SIREN/utilities/Errors.h"               // for InjectionFailure


namespace siren {
namespace interactions {

DMesonELoss::DMesonELoss() {
}


bool DMesonELoss::equal(CrossSection const & other) const {
    const DMesonELoss* x = dynamic_cast<const DMesonELoss*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(
            primary_types_,
            target_types_)
            ==
            std::tie(
            x->primary_types_,
            x->target_types_);
}



std::vector<siren::dataclasses::Particle::ParticleType> DMesonELoss::GetPossiblePrimaries() const {
    return std::vector<siren::dataclasses::Particle::ParticleType>(primary_types_.begin(), primary_types_.end());
}

// getting target should be the same for all the primary types
std::vector<siren::dataclasses::Particle::ParticleType> DMesonELoss::GetPossibleTargetsFromPrimary(siren::dataclasses::Particle::ParticleType primary_type) const {
    return std::vector<siren::dataclasses::Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<siren::dataclasses::Particle::ParticleType> DMesonELoss::GetPossibleTargets() const {
        return std::vector<siren::dataclasses::Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> DMesonELoss::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types_) {
        // hardcode the target type here, this should be fine
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParents(primary, siren::dataclasses::Particle::ParticleType::PPlus);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> DMesonELoss::GetPossibleSignaturesFromParents(siren::dataclasses::Particle::ParticleType primary_type, siren::dataclasses::Particle::ParticleType target_type) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary_type;
    signature.target_type = target_type;

    // first we deal with semileptonic decays where there are 3 final state particles
    signature.secondary_types.resize(2);
    signature.secondary_types[0] = primary_type; // same particle comes out
    signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::Hadrons; // there also is a hadronic vertex
    signatures.push_back(signature);
    return signatures;
}

// i am here

double DMesonELoss::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    siren::dataclasses::Particle::ParticleType primary_type = interaction.signature.primary_type;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    double primary_energy = interaction.primary_momentum[0];


    return TotalCrossSection(primary_type, primary_energy);
}

double DMesonELoss::TotalCrossSection(siren::dataclasses::Particle::ParticleType primary_type, double primary_energy) const {
    if(not primary_types_.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
    double log_energy = log10(primary_energy);
    double mb_to_cm2 = 1e-27;

    // current implementation uses only > 1PeV data
    double xsec = exp(1.891 + 0.205 * log_energy) - 2.157 + 1.264 * log_energy;

    return xsec * mb_to_cm2;
}

double DMesonELoss::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    double primary_energy;
    rk::P4 p1_lab;
    primary_energy = interaction.primary_momentum[0];
    p1_lab = p1;


    double final_energy = interaction.secondary_momenta[0][0];
    double z = 1 - final_energy / primary_energy;

    // The density is the truncated Gaussian in z normalized over [z_min_, z_max_].
    // Zero out-of-support records so the density support matches the sampler support
    // (the sampler rejects z outside [z_min_, z_max_] as well). This keeps closure
    // even for externally-constructed records with z < 0 (D meson gaining energy).
    if(z < z_min_ || z > z_max_) {
        return 0.0;
    }

    // now normalize the gaussian
    double total_xsec = TotalCrossSection(interaction.signature.primary_type, primary_energy);
    double z0 = 0.56;
    double sigma = 0.2;
    std::function<double(double)> integrand = [&] (double z) -> double {
            return exp(-(pow(z - z0, 2))/(2 * pow(sigma, 2)));
        };
    double unnormalized = siren::utilities::rombergIntegrate(integrand, z_min_, z_max_);
    double normalization = total_xsec / unnormalized;

    double diff_xsec = normalization * exp(-(pow(z - z0, 2))/(2 * pow(sigma, 2)));

    return diff_xsec;
}


double DMesonELoss::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    // Consider implementing DIS thershold at some point
    return 0;
}

void DMesonELoss::SampleFinalState(dataclasses::CrossSectionDistributionRecord& interaction, std::shared_ptr<siren::utilities::SIREN_random> random) const {

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    double primary_energy;
    double Dmass = interaction.primary_mass;
    rk::P4 p1_lab;
    p1_lab = p1;
    primary_energy = p1_lab.e();

    // Below the D meson mass there is no kinematically valid final state (the
    // D meson would need energy >= Dmass), and the reject loop would spin for an
    // astronomically large number of trials. Short-circuit with a recoverable
    // InjectionFailure so the Injector can retry the attempt instead of hanging.
    if (primary_energy <= Dmass) {
        throw siren::utilities::InjectionFailure(
            "DMesonELoss::SampleFinalState: primary energy below D meson mass; no valid final state");
    }

    // sample an inelasticity from gaussian using Box-Muller Transform
    double sigma = 0.2;
    double z0 = 0.56; // for mesons only, for baryons it's 0.59, but not implemented yet
    double u1, u2;
    double final_energy;
    bool accept;

    // Cap the rejection loop so a degenerate acceptance rate (primary energy just
    // above Dmass) terminates deterministically instead of hanging. Mirrors the
    // QuarkDISFromSpline reject-loop trial cap.
    const int max_trials = 1000;
    int trials = 0;

    do {
        if (++trials > max_trials) {
            throw siren::utilities::InjectionFailure(
                "DMesonELoss::SampleFinalState: failed to sample inelasticity within trial cap");
        }
        do
        {
            u1 = random->Uniform(0, 1);
        }
        while (u1 == 0);
        u2 = random->Uniform(0, 1);
        double z = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2) + z0;
        // now modify the energy of the charm hadron and the corresponding momentum
        final_energy = primary_energy * (1-z);
        // Reject z outside [z_min_, z_max_] so the realized sampling density is the
        // truncated Gaussian that DifferentialCrossSection/FinalStateProbability
        // normalize over the same interval (closure). z < z_min_ would otherwise let
        // the D meson GAIN energy (final_energy > primary_energy). The kinematic cut
        // final_energy^2 >= Dmass^2 is kept as a defensive guard.
        accept = (z >= z_min_) && (z <= z_max_) &&
                 (pow(final_energy, 2) - pow(Dmass, 2) >= 0);
    } while (!accept);
    // Guaranteed by z >= z_min_ > 0: the D meson loses energy.
    assert(final_energy <= primary_energy);
    double p3f = sqrt(pow(final_energy, 2) - pow(Dmass, 2));
    double p3i = std::sqrt(std::pow(p1.px(), 2) + std::pow(p1.py(), 2) + std::pow(p1.pz(), 2));
    double p_ratio = p3f / p3i;
    rk::P4 pf(p_ratio * geom3::Vector3(p1.px(), p1.py(), p1.pz()), Dmass);
    double E_H = primary_energy - pf.e(); // rest of the energy go into the nucleus
    double p_H = E_H; // assume collinearity, energy conservation and not momentum conservation ie massless (for now)
    double H_ratio = p_H / p3i;
    rk::P4 p4_H(H_ratio * geom3::Vector3(p1.px(), p1.py(), p1.pz()), 0);

    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = interaction.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & dmeson = secondaries[0];
    siren::dataclasses::SecondaryParticleRecord & hadron = secondaries[1];

    dmeson.SetFourMomentum({pf.e(), pf.px(), pf.py(), pf.pz()});
    dmeson.SetMass(pf.m());
    dmeson.SetHelicity(interaction.primary_helicity);

    hadron.SetFourMomentum({p4_H.e(), p4_H.px(), p4_H.py(), p4_H.pz()});
    hadron.SetMass(p4_H.m());
    hadron.SetHelicity(interaction.primary_helicity);
}

double DMesonELoss::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

// NOTE: SampleFinalState samples a single inelasticity DOF z (the D meson is set
// collinear with the parent, so there is no independent azimuth). That same z is
// fully captured by the density: DifferentialCrossSection reconstructs
// z = 1 - final_energy/primary_energy and FinalStateProbability returns the
// (truncated, normalized) Gaussian in z. Sampled DOF == density DOF, so this class
// is closure-safe in the standard unbiased configuration (the same cross-section
// object supplies both the injection and physical densities). Like the other charm
// cross sections here, BIASING the D kinematics with a separate phase-space channel
// is NOT supported and would produce incorrect weights.
std::vector<std::string> DMesonELoss::DensityVariables() const {
    return std::vector<std::string>{"Bjorken y"};
}

void DMesonELoss::InitializeSignatures() {
    return;
}

} // namespace interactions
} // namespace siren
