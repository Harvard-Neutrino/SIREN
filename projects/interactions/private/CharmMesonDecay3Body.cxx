#include "SIREN/interactions/CharmMesonDecay3Body.h"

#include <cmath>
#include <stdexcept>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle
#include "SIREN/math/Vector3D.h"                     // for Vector3D
#include "SIREN/utilities/Constants.h"               // for GeV, pi
#include "SIREN/utilities/Random.h"                  // for SIREN_random

#include "SIREN/utilities/Integration.h"          // for rombergInt...
#include "SIREN/utilities/Interpolator.h"

#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/CharmDecayKinematics.h"

namespace siren {
namespace interactions {

using charm_decay::particleMass;
using charm_decay::KStarMass;

CharmMesonDecay3Body::CharmMesonDecay3Body() {}

CharmMesonDecay3Body::CharmMesonDecay3Body(siren::dataclasses::Particle::ParticleType primary) {
  if (primary != siren::dataclasses::Particle::ParticleType::DPlus &&
      primary != siren::dataclasses::Particle::ParticleType::D0) {
    throw std::runtime_error("CharmMesonDecay3Body: only D0 and D+ are implemented. Use CharmMesonDecay, which covers D0/D+/Ds and their anti-flavors.");
  }
}

bool CharmMesonDecay3Body::equal(Decay const & other) const {
    const CharmMesonDecay3Body* x = dynamic_cast<const CharmMesonDecay3Body*>(&other);

    if(!x)
        return false;
    else
        return primary_types == x->primary_types;
}

double CharmMesonDecay3Body::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.signature.primary_type);
}

// in this implementation, should we take total decay width to be only the channels we considered?
double CharmMesonDecay3Body::TotalDecayWidth(siren::dataclasses::Particle::ParticleType primary) const {
    double total_width = 0;
    std::vector<dataclasses::InteractionSignature> possible_signatures = GetPossibleSignaturesFromParent(primary);
    for (auto sig : possible_signatures) {
      // make a fake record and full from total decay width for final state
      siren::dataclasses::InteractionRecord fake_record;
      fake_record.signature = sig;
      double this_width = TotalDecayWidthForFinalState(fake_record);
      total_width += this_width;
    }
    return total_width;
}

// current problem: in implementation we see kaons and pions both as hadrons, but they should have different branching ratios and form factors
double CharmMesonDecay3Body::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {
    // Sentinel-init: all valid branching ratios and lifetimes are strictly
    // positive, so a negative value unambiguously flags an unmatched mode.
    double branching_ratio = -1.0;
    double tau = -1.0; // total lifetime for all visible and invisible modes
    // read in the signature and types
    siren::dataclasses::Particle::ParticleType primary = record.signature.primary_type;
    std::vector<siren::dataclasses::Particle::ParticleType> secondaries_vector = record.signature.secondary_types;
    std::set<siren::dataclasses::Particle::ParticleType> secondaries = std::set<siren::dataclasses::Particle::ParticleType>(secondaries_vector.begin(), secondaries_vector.end());

    // define the decay modes
    std::set<siren::dataclasses::Particle::ParticleType> k0_eplus_nue = {siren::dataclasses::Particle::ParticleType::K0Bar,
                                                                            siren::dataclasses::Particle::ParticleType::EPlus,
                                                                            siren::dataclasses::Particle::ParticleType::NuE};
    std::set<siren::dataclasses::Particle::ParticleType> kminus_eplus_nue = {siren::dataclasses::Particle::ParticleType::KMinus,
                                                                            siren::dataclasses::Particle::ParticleType::EPlus,
                                                                            siren::dataclasses::Particle::ParticleType::NuE};
    std::set<siren::dataclasses::Particle::ParticleType> k0_muplus_numu = {siren::dataclasses::Particle::ParticleType::K0Bar,
                                                                            siren::dataclasses::Particle::ParticleType::MuPlus,
                                                                            siren::dataclasses::Particle::ParticleType::NuMu};
    std::set<siren::dataclasses::Particle::ParticleType> kminus_muplus_numu = {siren::dataclasses::Particle::ParticleType::KMinus,
                                                                            siren::dataclasses::Particle::ParticleType::MuPlus,
                                                                            siren::dataclasses::Particle::ParticleType::NuMu};
    std::set<siren::dataclasses::Particle::ParticleType> hadrons = {siren::dataclasses::Particle::ParticleType::Hadrons};
    if (primary == siren::dataclasses::Particle::ParticleType::DPlus) {
      tau = 1040 * (1e-15);
      // Exclusive K + K* semileptonic BR (PDG): K0bar l nu = 8.74%,
      // K*0bar l nu = 5.33%, summed to 14.07% (kinematic K/K* mixing, see fracK).
      // e and mu modes share the same value (lepton-mass effects sub-percent).
      if (secondaries == k0_eplus_nue) {branching_ratio = .1407;}
      else if (secondaries == k0_muplus_numu) {branching_ratio = .1407;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .1407);} // everything else
    } else if (primary == siren::dataclasses::Particle::ParticleType::D0) {
      tau = 410.1 * (1e-15);
      // Exclusive K + K* semileptonic BR (PDG): K- l nu = 3.41%,
      // K*- l nu = 2.17%, summed to 5.58% (kinematic K/K* mixing, see fracK).
      if (secondaries == kminus_eplus_nue) {branching_ratio = .0558;}
      else if (secondaries == kminus_muplus_numu) {branching_ratio = .0558;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .0558);} // everything else
    }
    else {
        // Signatures come from GetPossibleSignaturesFromParent, so an unsupported
        // primary means the lists are out of sync. Fail loudly rather than
        // return an indeterminate width.
        throw std::runtime_error("CharmMesonDecay3Body::TotalDecayWidthForFinalState: unsupported primary particle type");
    }
    // Guard the matched-primary / unmatched-secondaries case (sentinel still set).
    if (tau <= 0.0 || branching_ratio < 0.0) {
        throw std::runtime_error("CharmMesonDecay3Body::TotalDecayWidthForFinalState: no implemented decay mode matches this signature");
    }
    return branching_ratio * siren::utilities::Constants::hbar / tau * siren::utilities::Constants::GeV;
}


std::vector<dataclasses::InteractionSignature> CharmMesonDecay3Body::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> CharmMesonDecay3Body::GetPossibleSignaturesFromParent(siren::dataclasses::Particle::ParticleType primary) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    // initialize semileptonic signatures
    dataclasses::InteractionSignature semilep_signature;
    semilep_signature.primary_type = primary;
    semilep_signature.target_type = siren::dataclasses::Particle::ParticleType::Decay;
    semilep_signature.secondary_types.resize(3);
    // initialize other signatures
    dataclasses::InteractionSignature hadron_signature;
    hadron_signature.primary_type = primary;
    hadron_signature.target_type = siren::dataclasses::Particle::ParticleType::Decay;
    hadron_signature.secondary_types.resize(1);

    if (primary==siren::dataclasses::Particle::ParticleType::DPlus) {
      // semi-leptonic modes with muon and electron
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::K0Bar;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::EPlus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuE;
      signatures.push_back(semilep_signature);
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::K0Bar;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::MuPlus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuMu;
      signatures.push_back(semilep_signature);
      // all other modes implemented as one big bang
      hadron_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      signatures.push_back(hadron_signature);
    } else if (primary==siren::dataclasses::Particle::ParticleType::D0) {
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::KMinus;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::EPlus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuE;
      signatures.push_back(semilep_signature);
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::KMinus;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::MuPlus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuMu;
      signatures.push_back(semilep_signature);
      hadron_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      signatures.push_back(hadron_signature);
    }
    else {
      throw std::runtime_error("CharmMesonDecay3Body::GetPossibleSignaturesFromParent: only D0 and D+ are implemented. Use CharmMesonDecay for D0/D+/Ds and anti-flavors.");
    }
    return signatures;
}

double CharmMesonDecay3Body::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidthForFinalState(record);
}

// this is temporary implementation
void CharmMesonDecay3Body::SampleFinalStateHadronic(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & hadrons = secondaries[0];
    rk::P4 p4D_lab(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    hadrons.SetFourMomentum({p4D_lab.e(), p4D_lab.px(), p4D_lab.py(), p4D_lab.pz()});
    hadrons.SetMass(p4D_lab.m());
    hadrons.SetHelicity(record.primary_helicity);
}

void CharmMesonDecay3Body::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // Hadronic decay branch -- identical to the 2-body class.
    dataclasses::InteractionSignature signature = record.signature;
    if (signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::Hadrons) {
      SampleFinalStateHadronic(record, random);
      return;
    }

    // 3-body phase space (Pythia ParticleDecays::threeBody):
    // D (m0) -> K (m1) + lepton (m2) + neutrino (m3). Sample
    // m23 = sqrt(q^2) flat with accept-reject on the phase-space weight, then
    // apply the V-A matrix-element correction.
    //
    // K/K*(892) mixing is KINEMATIC only: a fraction (1 - fracK) of events draw
    // the K*(892) mass for the hadron, broadening the lepton spectrum. The
    // secondary's ParticleType stays the signature K species (K0bar / K-) -- we
    // neither advertise separate K* signatures nor re-type the secondary, since
    // only the mass matters for weighting the lepton/neutrino kinematics.

    double mD      = particleMass(record.signature.primary_type);
    double mK_base = particleMass(record.signature.secondary_types[0]); // K mass from signature
    double ml      = particleMass(record.signature.secondary_types[1]); // lepton mass
    double mnu     = 0.0;                                               // (massless) neutrino

    // K/K* mixing fractions from PDG semileptonic branching ratios
    // K*(892) mass read from Constants so the sampler and the
    // FinalStateProbability normalizer (KStarMass()) use the identical value
    // (required for weighting closure).
    double mKstar = KStarMass();
    double fracK;
    if (record.signature.primary_type == siren::dataclasses::Particle::ParticleType::DPlus) {
        // D+ -> Kbar0 l nu:  BR_K = 8.74%, BR_K*bar = 5.33%
        fracK = 8.74 / (8.74 + 5.33);
    } else {
        // D0 -> K- l nu:     BR_K = 3.41%, BR_K*- = 2.17%
        fracK = 3.41 / (3.41 + 2.17);
    }
    double mK = (random->Uniform(0, 1) < fracK) ? mK_base : mKstar;

    // D meson 4-momentum in lab frame (needed for final lab-frame boost)
    rk::P4 p4D_lab(geom3::Vector3(record.primary_momentum[1],
                                   record.primary_momentum[2],
                                   record.primary_momentum[3]),
                   record.primary_mass);

    // Phase-space limits for m23 = sqrt(q^2) = (lepton+nu) invariant mass
    double m23Min = ml + mnu;
    double m23Max = mD - mK;
    double mDiff  = m23Max - m23Min;

    // Kinematic envelope used as the accept-reject ceiling
    double p1Max = 0.5 * std::sqrt((mD - mK - m23Min) * (mD + mK + m23Min)
                                 * (mD + mK - m23Min) * (mD - mK + m23Min)) / mD;
    double p23Max = 0.5 * std::sqrt((m23Max - ml - mnu) * (m23Max + ml + mnu)
                                  * (m23Max + ml - mnu) * (m23Max - ml + mnu)) / m23Max;
    double wtPSmax = 0.5 * p1Max * p23Max;

    // V-A matrix element upper bound (Pythia meMode == 22)
    double wtMEmax = std::min(std::pow(mD, 4) / 16.0,
                              mD * (mD - mK - ml) * (mD - mK - mnu) * (mD - ml - mnu));

    rk::P4 p4K_Drest, p4l_Drest, p4nu_Drest;
    double wtME;

    // --- Outer accept-reject on the V-A matrix element ---
    do {
        wtME = 1.0;

        // --- Inner accept-reject on 3-body phase space ---
        double m23, p1Abs, p23Abs, wtPS;
        do {
            m23    = m23Min + random->Uniform(0, 1) * mDiff;
            p1Abs  = 0.5 * std::sqrt((mD - mK - m23) * (mD + mK + m23)
                                   * (mD + mK - m23) * (mD - mK + m23)) / mD;
            p23Abs = 0.5 * std::sqrt((m23 - ml - mnu) * (m23 + ml + mnu)
                                   * (m23 + ml - mnu) * (m23 - ml + mnu)) / m23;
            wtPS   = p1Abs * p23Abs;
        } while (wtPS < random->Uniform(0, 1) * wtPSmax);

        // Set up m23 -> lepton + neutrino, isotropic in m23 rest frame
        double cosTheta23 = random->Uniform(-1, 1);
        double sinTheta23 = std::sin(std::acos(cosTheta23));
        double phi23      = random->Uniform(0, 2 * M_PI);
        geom3::Vector3 dir23(sinTheta23 * std::cos(phi23),
                             sinTheta23 * std::sin(phi23),
                             cosTheta23);
        rk::P4 p4l_m23rest(p23Abs * dir23, ml);
        rk::P4 p4nu_m23rest(-p23Abs * dir23, mnu);

        // Set up D -> K + (m23 system), isotropic in D rest frame
        double cosTheta1 = random->Uniform(-1, 1);
        double sinTheta1 = std::sin(std::acos(cosTheta1));
        double phi1      = random->Uniform(0, 2 * M_PI);
        geom3::Vector3 dir1(sinTheta1 * std::cos(phi1),
                            sinTheta1 * std::sin(phi1),
                            cosTheta1);
        p4K_Drest = rk::P4(p1Abs * dir1, mK);
        rk::P4 p4m23_Drest(-p1Abs * dir1, m23);

        // Boost lepton/neutrino from m23 rest frame to D rest frame
        rk::Boost boost_m23_to_Drest = p4m23_Drest.labBoost();
        p4l_Drest  = p4l_m23rest.boost(boost_m23_to_Drest);
        p4nu_Drest = p4nu_m23rest.boost(boost_m23_to_Drest);

        // V-A matrix element weight:
        //   wtME = m_D * E_lepton * (p_neutrino . p_kaon) in D rest frame
        wtME = mD * p4l_Drest.e() * p4nu_Drest.dot(p4K_Drest);
    } while (wtME < random->Uniform(0, 1) * wtMEmax);

    // Boost kaon / lepton / neutrino from D rest frame to lab frame
    rk::Boost boost_D_to_lab = p4D_lab.labBoost();
    rk::P4 p4K_lab   = p4K_Drest.boost(boost_D_to_lab);
    rk::P4 p4l_lab   = p4l_Drest.boost(boost_D_to_lab);
    rk::P4 p4nu_lab  = p4nu_Drest.boost(boost_D_to_lab);

    // Write the secondaries
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & kpi      = secondaries[0];
    siren::dataclasses::SecondaryParticleRecord & lepton   = secondaries[1];
    siren::dataclasses::SecondaryParticleRecord & neutrino = secondaries[2];

    kpi.SetFourMomentum({p4K_lab.e(), p4K_lab.px(), p4K_lab.py(), p4K_lab.pz()});
    kpi.SetMass(p4K_lab.m());
    kpi.SetHelicity(record.primary_helicity);

    lepton.SetFourMomentum({p4l_lab.e(), p4l_lab.px(), p4l_lab.py(), p4l_lab.pz()});
    lepton.SetMass(p4l_lab.m());
    lepton.SetHelicity(record.primary_helicity);

    neutrino.SetFourMomentum({p4nu_lab.e(), p4nu_lab.px(), p4nu_lab.py(), p4nu_lab.pz()});
    neutrino.SetMass(p4nu_lab.m());
    neutrino.SetHelicity(record.primary_helicity);
}

// FinalStateProbability is the normalized q^2 density SampleFinalState produces
// (closure rationale in CharmDecayKinematics.h). Builds the same K/K*(892)
// mixture the sampler used, evaluates charm_decay::SampledQ2Density for the
// recorded component, and normalizes via charm_decay::SampledQ2Normalization.
double CharmMesonDecay3Body::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  dataclasses::InteractionSignature signature = record.signature;
  // Guard: finalize() does not copy the signature, so a finalized record can
  // arrive with an empty signature/secondaries -> indexing below would be UB.
  if (signature.secondary_types.empty()) {
    throw std::runtime_error("CharmMesonDecay3Body::FinalStateProbability: record has an empty signature. Set record.signature (and secondary masses/momenta) before calling; finalize() does not copy the signature.");
  }
  // Fully hadronic catch-all: deterministic final state, density 1.
  if (signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::Hadrons) {
    return 1.0;
  }

  siren::dataclasses::Particle::ParticleType primary = signature.primary_type;
  bool apply_va = true;   // 3-body class always applies the V-A matrix element

  if (record.secondary_masses.size() < 2 || record.secondary_momenta.empty()) {
    throw std::runtime_error("CharmMesonDecay3Body::FinalStateProbability: record secondaries are not populated (need >=2 masses and >=1 momentum).");
  }
  double mD = record.primary_mass;
  double ml = record.secondary_masses[1];
  double mK = record.secondary_masses[0];   // hadron mass actually sampled
  rk::P4 pD(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), mD);
  rk::P4 pK(geom3::Vector3(record.secondary_momenta[0][1], record.secondary_momenta[0][2], record.secondary_momenta[0][3]), mK);
  double q2 = (pD - pK).dot(pD - pK);

  // K/K*(892) kinematic mixture matching SampleFinalState's fracK.
  double mK_base = particleMass(signature.secondary_types[0]);
  double mKstar = KStarMass();
  double fracK;
  if (primary == siren::dataclasses::Particle::ParticleType::DPlus) {
    fracK = 8.74 / (8.74 + 5.33);
  } else {
    fracK = 3.41 / (3.41 + 2.17);
  }
  std::vector<double> masses    = {mK_base, mKstar};
  std::vector<double> fractions = {fracK, 1.0 - fracK};

  int comp = -1;
  double best = 1e30;
  for (size_t i = 0; i < masses.size(); ++i) {
    double d = std::abs(mK - masses[i]);
    if (d < best) { best = d; comp = (int)i; }
  }
  if (comp < 0) return 0.0;

  double g = charm_decay::SampledQ2Density(mD, masses[comp], ml, q2, apply_va);
  if (g <= 0.0) return 0.0;
  double norm = charm_decay::SampledQ2Normalization(mD, masses[comp], ml, apply_va, norm_cache);
  if (norm <= 0.0) return 0.0;

  // Normalized q^2 pdf summing over the K/K* mixture.
  return fractions[comp] * g / norm;
}

std::vector<std::string> CharmMesonDecay3Body::DensityVariables() const {
    return std::vector<std::string>{"Q2"};
}



} // namespace interactions
} // namespace siren

