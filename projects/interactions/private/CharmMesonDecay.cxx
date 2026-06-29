#include "SIREN/interactions/CharmMesonDecay.h"

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

CharmMesonDecay::CharmMesonDecay() {}

CharmMesonDecay::CharmMesonDecay(siren::dataclasses::Particle::ParticleType primary) {}

bool CharmMesonDecay::equal(Decay const & other) const {
    const CharmMesonDecay* x = dynamic_cast<const CharmMesonDecay*>(&other);

    if(!x)
        return false;
    else
        return primary_types == x->primary_types;
}

double CharmMesonDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.signature.primary_type);
}

// in this implementation, should we take total decay width to be only the channels we considered?
double CharmMesonDecay::TotalDecayWidth(siren::dataclasses::Particle::ParticleType primary) const {
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
double CharmMesonDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {
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
    std::set<siren::dataclasses::Particle::ParticleType> hadrons_muplus_numu = {siren::dataclasses::Particle::ParticleType::Hadrons,
                                                                            siren::dataclasses::Particle::ParticleType::MuPlus,
                                                                            siren::dataclasses::Particle::ParticleType::NuMu};
    std::set<siren::dataclasses::Particle::ParticleType> hadrons_eplus_nue = {siren::dataclasses::Particle::ParticleType::Hadrons,
                                                                            siren::dataclasses::Particle::ParticleType::EPlus,
                                                                            siren::dataclasses::Particle::ParticleType::NuE};
    std::set<siren::dataclasses::Particle::ParticleType> hadrons = {siren::dataclasses::Particle::ParticleType::Hadrons};
    // Anti-flavor (cbar) decay-mode sets, sign-conjugated from the c modes above.
    std::set<siren::dataclasses::Particle::ParticleType> k0_eminus_nuebar = {siren::dataclasses::Particle::ParticleType::K0,
                                                                            siren::dataclasses::Particle::ParticleType::EMinus,
                                                                            siren::dataclasses::Particle::ParticleType::NuEBar};
    std::set<siren::dataclasses::Particle::ParticleType> kplus_eminus_nuebar = {siren::dataclasses::Particle::ParticleType::KPlus,
                                                                            siren::dataclasses::Particle::ParticleType::EMinus,
                                                                            siren::dataclasses::Particle::ParticleType::NuEBar};
    std::set<siren::dataclasses::Particle::ParticleType> k0_muminus_numubar = {siren::dataclasses::Particle::ParticleType::K0,
                                                                            siren::dataclasses::Particle::ParticleType::MuMinus,
                                                                            siren::dataclasses::Particle::ParticleType::NuMuBar};
    std::set<siren::dataclasses::Particle::ParticleType> kplus_muminus_numubar = {siren::dataclasses::Particle::ParticleType::KPlus,
                                                                            siren::dataclasses::Particle::ParticleType::MuMinus,
                                                                            siren::dataclasses::Particle::ParticleType::NuMuBar};
    std::set<siren::dataclasses::Particle::ParticleType> hadrons_muminus_numubar = {siren::dataclasses::Particle::ParticleType::Hadrons,
                                                                            siren::dataclasses::Particle::ParticleType::MuMinus,
                                                                            siren::dataclasses::Particle::ParticleType::NuMuBar};
    std::set<siren::dataclasses::Particle::ParticleType> hadrons_eminus_nuebar = {siren::dataclasses::Particle::ParticleType::Hadrons,
                                                                            siren::dataclasses::Particle::ParticleType::EMinus,
                                                                            siren::dataclasses::Particle::ParticleType::NuEBar};
    if (primary == siren::dataclasses::Particle::ParticleType::DPlus) {
      tau = 1040 * (1e-15);
      // Exclusive K + K* semileptonic BR (PDG): K0bar l nu = 8.74%,
      // K*0bar l nu = 5.33%, summed to 14.07% because the K l nu signature
      // is sampled with kinematic K/K* mixing (see SampleFinalState fracK).
      // e and mu modes share the same value (lepton-mass effects sub-percent).
      if (secondaries == k0_eplus_nue) {branching_ratio = .1407;}
      else if (secondaries == k0_muplus_numu) {branching_ratio = .1407;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .1407);}
    } else if (primary == siren::dataclasses::Particle::ParticleType::DMinus) {
      // CP-mirror of D+ (same lifetime, same physical PDG BRs).
      tau = 1040 * (1e-15);
      if (secondaries == k0_eminus_nuebar) {branching_ratio = .1407;}
      else if (secondaries == k0_muminus_numubar) {branching_ratio = .1407;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .1407);}
    } else if (primary == siren::dataclasses::Particle::ParticleType::D0) {
      tau = 410.1 * (1e-15);
      // Exclusive K + K* semileptonic BR (PDG): K- l nu = 3.41%,
      // K*- l nu = 2.17%, summed to 5.58% (kinematic K/K* mixing, see fracK).
      if (secondaries == kminus_eplus_nue) {branching_ratio = .0558;}
      else if (secondaries == kminus_muplus_numu) {branching_ratio = .0558;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .0558);}
    } else if (primary == siren::dataclasses::Particle::ParticleType::D0Bar) {
      // CP-mirror of D0 (same lifetime, same physical PDG BRs).
      tau = 410.1 * (1e-15);
      if (secondaries == kplus_eminus_nuebar) {branching_ratio = .0558;}
      else if (secondaries == kplus_muminus_numubar) {branching_ratio = .0558;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .0558);}
    } else if (primary == siren::dataclasses::Particle::ParticleType::DsPlus) {
      tau = 504 * (1e-15);  // Ds+ lifetime: 504 fs (PDG)
      // Physical PDG BRs for Ds+ semileptonic decay. Inclusive
      // Ds -> e X = 0.0654, Ds -> mu X = 0.0654 (no tau feed-down).
      // Daughter "Hadrons" stands in for eta / eta' / phi (sampled in
      // SampleFinalState with fractions 0.46 / 0.16 / 0.38).
      if (secondaries == hadrons_eplus_nue) {branching_ratio = .0654;}
      else if (secondaries == hadrons_muplus_numu) {branching_ratio = .0654;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .0654);}
    } else if (primary == siren::dataclasses::Particle::ParticleType::DsMinus) {
      // CP-mirror of Ds+ (eta/eta'/phi are self-conjugate, only lepton/nu flip).
      tau = 504 * (1e-15);
      if (secondaries == hadrons_eminus_nuebar) {branching_ratio = .0654;}
      else if (secondaries == hadrons_muminus_numubar) {branching_ratio = .0654;}
      else if (secondaries == hadrons) {branching_ratio = (1 - 2 * .0654);}
    }
    else {
        // Signatures are produced by GetPossibleSignaturesFromParent, so an
        // unsupported primary here means the two lists are out of sync. Fail
        // loudly rather than return an indeterminate width that would silently
        // corrupt TotalDecayWidth / FinalStateProbability.
        throw std::runtime_error("CharmMesonDecay::TotalDecayWidthForFinalState: unsupported primary particle type");
    }
    // Guard the matched-primary / unmatched-secondaries case (sentinel still set).
    if (tau <= 0.0 || branching_ratio < 0.0) {
        throw std::runtime_error("CharmMesonDecay::TotalDecayWidthForFinalState: no implemented decay mode matches this signature");
    }
    return branching_ratio * siren::utilities::Constants::hbar / tau * siren::utilities::Constants::GeV;
}


std::vector<dataclasses::InteractionSignature> CharmMesonDecay::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> CharmMesonDecay::GetPossibleSignaturesFromParent(siren::dataclasses::Particle::ParticleType primary) const {
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
    } else if (primary==siren::dataclasses::Particle::ParticleType::DMinus) {
      // CP-mirror of D+: K0 (sbar d -> contains sbar), e-/mu-, nubar_e/nubar_mu.
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::K0;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::EMinus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuEBar;
      signatures.push_back(semilep_signature);
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::K0;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::MuMinus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuMuBar;
      signatures.push_back(semilep_signature);
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
    } else if (primary==siren::dataclasses::Particle::ParticleType::D0Bar) {
      // CP-mirror of D0: K+ (u sbar), e-/mu-, nubar_e/nubar_mu.
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::KPlus;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::EMinus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuEBar;
      signatures.push_back(semilep_signature);
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::KPlus;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::MuMinus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuMuBar;
      signatures.push_back(semilep_signature);
      hadron_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      signatures.push_back(hadron_signature);
    } else if (primary==siren::dataclasses::Particle::ParticleType::DsPlus) {
      // Ds: muonic + electronic semileptonic channels plus all-hadronic catch-all.
      // The Hadrons daughter stands in for eta / eta' / phi (sampled inline
      // in SampleFinalState; SIREN ParticleType has no entries for these).
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::EPlus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuE;
      signatures.push_back(semilep_signature);
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::MuPlus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuMu;
      signatures.push_back(semilep_signature);
      hadron_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      signatures.push_back(hadron_signature);
    } else if (primary==siren::dataclasses::Particle::ParticleType::DsMinus) {
      // CP-mirror of Ds+: eta/eta'/phi are self-conjugate, only lepton/nu flip.
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::EMinus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuEBar;
      signatures.push_back(semilep_signature);
      semilep_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      semilep_signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::MuMinus;
      semilep_signature.secondary_types[2] = siren::dataclasses::Particle::ParticleType::NuMuBar;
      signatures.push_back(semilep_signature);
      hadron_signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
      signatures.push_back(hadron_signature);
    }
    else {
      std::cout << "this D meson decay has not been implemented yet" << std::endl;
    }
    return signatures;
}

double CharmMesonDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidthForFinalState(record);
}

// this is temporary implementation
void CharmMesonDecay::SampleFinalStateHadronic(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & hadrons = secondaries[0];
    rk::P4 p4D_lab(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    hadrons.SetFourMomentum({p4D_lab.e(), p4D_lab.px(), p4D_lab.py(), p4D_lab.pz()});
    hadrons.SetMass(p4D_lab.m());
    hadrons.SetHelicity(record.primary_helicity);
}

void CharmMesonDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // first handle hadronic decay separately. Use size==1 so we don't
    // mis-route the Ds semileptonic mode (Hadrons + mu + nu, size 3).
    dataclasses::InteractionSignature signature = record.signature;
    if (signature.secondary_types.size() == 1 &&
        signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::Hadrons) {
      SampleFinalStateHadronic(record, random);
      return;
    }

    // 3-body phase space (Pythia ParticleDecays::threeBody):
    // D (m0) -> hadron (m1) + lepton (m2) + neutrino (m3). Sample
    // m23 = sqrt(q^2) flat with accept-reject on the phase-space weight.
    // D+/D0: V-A matrix-element correction with K vs K* mixing.
    // Ds: pure phase space, daughter sampled inline as eta/eta'/phi
    // (fractions 0.46/0.16/0.38 from BRs 2.3/0.8/1.9 %).

    siren::dataclasses::Particle::ParticleType primary = record.signature.primary_type;
    bool is_Ds = (primary == siren::dataclasses::Particle::ParticleType::DsPlus ||
                  primary == siren::dataclasses::Particle::ParticleType::DsMinus);

    double mD = particleMass(primary);                              // m0
    double ml = particleMass(record.signature.secondary_types[1]);  // m2
    double mnu = 0.0;                                               // m3

    double mK;  // m1: hadron daughter mass (chosen per-decay)
    if (is_Ds) {
        // Ds -> eta mu nu (BR 2.3%), Ds -> eta' mu nu (0.8%), Ds -> phi mu nu (1.9%)
        // Cumulative: eta=0.46, eta+eta'=0.62, all=1.0
        double mEta      = siren::utilities::Constants::EtaMass;      // 0.547862 GeV
        double mEtaPrime = siren::utilities::Constants::EtaPrimeMass; // 0.95778  GeV
        double mPhi      = siren::utilities::Constants::PhiMass;      // 1.019461 GeV
        double r = random->Uniform(0, 1);
        if (r < 2.3 / 5.0) {
            mK = mEta;
        } else if (r < 3.1 / 5.0) {
            mK = mEtaPrime;
        } else {
            mK = mPhi;
        }
    } else {
        // D+/D0: K vs K*(892) with V-A weighting
        double mK_base = particleMass(record.signature.secondary_types[0]);
        // K*(892) mass; read from Constants so the FinalStateProbability
        // normalizer (KStarMass()) and the sampler use the identical value
        // (required for weighting closure).
        double mKstar = KStarMass();
        double fracK;
        if (primary == siren::dataclasses::Particle::ParticleType::DPlus ||
            primary == siren::dataclasses::Particle::ParticleType::DMinus) {
            // D+ -> K0bar mu+ nu: BR=8.74%, D+ -> K*0bar mu+ nu: BR=5.33%
            // (D- BRs identical by CP.)
            fracK = 8.74 / (8.74 + 5.33);  // ~0.621
        } else {
            // D0 -> K- mu+ nu: BR=3.41%, D0 -> K*- mu+ nu: BR=2.17%
            // (D0bar BRs identical by CP.)
            fracK = 3.41 / (3.41 + 2.17);  // ~0.611
        }
        mK = (random->Uniform(0, 1) < fracK) ? mK_base : mKstar;
    }

    // D meson 4-momentum in lab frame
    rk::P4 p4D_lab(geom3::Vector3(record.primary_momentum[1],
                                    record.primary_momentum[2],
                                    record.primary_momentum[3]),
                    record.primary_mass);
    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();

    // Kinematic limits for m23 (lepton-neutrino invariant mass)
    double m23Min = ml + mnu;        // minimum: lepton + neutrino at rest
    double m23Max = mD - mK;         // maximum: kaon at rest
    double mDiff  = m23Max - m23Min;

    // Maximum phase space weight: p1Max * p23Max
    // p1Max = kaon momentum when m23 = m23Min
    double p1Max = 0.5 * sqrt((mD - mK - m23Min) * (mD + mK + m23Min)
                             * (mD + mK - m23Min) * (mD - mK + m23Min)) / mD;
    // p23Max = lepton momentum in m23 rest frame when m23 = m23Max
    double p23Max = 0.5 * sqrt((m23Max - ml - mnu) * (m23Max + ml + mnu)
                              * (m23Max + ml - mnu) * (m23Max - ml + mnu)) / m23Max;
    double wtPSmax = 0.5 * p1Max * p23Max;

    // V-A matrix element upper bound (from Pythia, meMode == 22).
    // For Ds we skip V-A weighting and use pure 3-body phase space -- set
    // wtMEmax = 1.0 and wtME = 1.0 inside the loop so the rejection test
    // (wtME < rand * wtMEmax) is always false and we exit after one iteration.
    double wtMEmax = is_Ds
        ? 1.0
        : std::min(std::pow(mD, 4) / 16.0,
                   mD * (mD - mK - ml) * (mD - mK - mnu) * (mD - ml - mnu));
    double wtME;

    rk::P4 p4K_Drest, p4l_Drest, p4nu_Drest;

    do {
    wtME = 1.0;

    // --- Step 1: Sample m23 flat, accept-reject on phase space weight ---
    double m23, p1Abs, p23Abs, wtPS;
    do {
        m23 = m23Min + random->Uniform(0, 1) * mDiff;

        // p1Abs = kaon momentum in D rest frame for this m23
        p1Abs = 0.5 * sqrt((mD - mK - m23) * (mD + mK + m23)
                          * (mD + mK - m23) * (mD - mK + m23)) / mD;
        // p23Abs = lepton momentum in m23 rest frame
        p23Abs = 0.5 * sqrt((m23 - ml - mnu) * (m23 + ml + mnu)
                           * (m23 + ml - mnu) * (m23 - ml + mnu)) / m23;
        wtPS = p1Abs * p23Abs;
    } while (wtPS < random->Uniform(0, 1) * wtPSmax);

    // --- Step 2: Set up m23 -> lepton + neutrino isotropically ---
    double cosTheta23 = random->Uniform(-1, 1);
    double sinTheta23 = std::sin(std::acos(cosTheta23));
    double phi23 = random->Uniform(0, 2 * M_PI);

    // Lepton and neutrino in m23 rest frame
    geom3::Vector3 dir23(sinTheta23 * std::cos(phi23),
                         sinTheta23 * std::sin(phi23),
                         cosTheta23);
    rk::P4 p4l_m23rest(p23Abs * dir23, ml);
    rk::P4 p4nu_m23rest(-p23Abs * dir23, mnu);

    // --- Step 3: Set up D -> K + (m23 system) isotropically ---
    double cosTheta1 = random->Uniform(-1, 1);
    double sinTheta1 = std::sin(std::acos(cosTheta1));
    double phi1 = random->Uniform(0, 2 * M_PI);

    geom3::Vector3 dir1(sinTheta1 * std::cos(phi1),
                        sinTheta1 * std::sin(phi1),
                        cosTheta1);
    // Kaon in D rest frame
    p4K_Drest = rk::P4(p1Abs * dir1, mK);
    // m23 system in D rest frame (opposite to kaon)
    rk::P4 p4m23_Drest(-p1Abs * dir1, m23);

    // --- Step 4: Boost lepton and neutrino from m23 rest frame to D rest frame ---
    rk::Boost boost_m23_to_Drest = p4m23_Drest.labBoost();
    p4l_Drest = p4l_m23rest.boost(boost_m23_to_Drest);
    p4nu_Drest = p4nu_m23rest.boost(boost_m23_to_Drest);

    // --- Step 5: V-A matrix element weight (skipped for Ds -- pure phase space) ---
    // wtME = m_D * E_lepton * (p_neutrino . p_kaon) in D rest frame
    // This is Lorentz invariant up to the m_D * E_l factor which equals p_D . p_l
    if (!is_Ds) {
        wtME = mD * p4l_Drest.e() * p4nu_Drest.dot(p4K_Drest);
    }
    // For Ds: wtME stays at 1.0 (set at the top of the do-loop). Combined with
    // wtMEmax = 1.0, the test (wtME < rand[0,1) * wtMEmax) is always false,
    // so we exit after one iteration with no V-A reweighting.

    } while (wtME < random->Uniform(0, 1) * wtMEmax);

    // --- Boost all particles from D rest frame to lab frame ---
    rk::Boost boost_D_to_lab = p4D_lab.labBoost();
    rk::P4 p4K_lab = p4K_Drest.boost(boost_D_to_lab);
    rk::P4 p4l_lab = p4l_Drest.boost(boost_D_to_lab);
    rk::P4 p4nu_lab = p4nu_Drest.boost(boost_D_to_lab);

    // --- Set secondary particle records ---
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & kpi = secondaries[0];
    siren::dataclasses::SecondaryParticleRecord & lepton = secondaries[1];
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
// (closure rationale in CharmDecayKinematics.h). It builds the same hadron-mass
// mixture the sampler used, evaluates charm_decay::SampledQ2Density for the
// recorded component, and normalizes via charm_decay::SampledQ2Normalization.
double CharmMesonDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  dataclasses::InteractionSignature signature = record.signature;
  // Guard: finalize() does not copy the signature, so a finalized record can
  // arrive with an empty signature/secondaries -> indexing below would be UB.
  if (signature.secondary_types.empty()) {
    throw std::runtime_error("CharmMesonDecay::FinalStateProbability: record has an empty signature. Set record.signature (and secondary masses/momenta) before calling; finalize() does not copy the signature.");
  }
  // Fully hadronic catch-all: a single Hadrons daughter carries the full
  // primary 4-momentum, so the final state is deterministic and its density is 1.
  if (signature.secondary_types.size() == 1 &&
      signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::Hadrons) {
    return 1.0;
  }

  siren::dataclasses::Particle::ParticleType primary = signature.primary_type;
  bool is_Ds = (primary == siren::dataclasses::Particle::ParticleType::DsPlus ||
                primary == siren::dataclasses::Particle::ParticleType::DsMinus);
  bool apply_va = !is_Ds;   // Ds is pure phase space (no V-A matrix element)

  if (record.secondary_masses.size() < 2 || record.secondary_momenta.empty()) {
    throw std::runtime_error("CharmMesonDecay::FinalStateProbability: record secondaries are not populated (need >=2 masses and >=1 momentum).");
  }
  // Reconstruct q^2 = (p_D - p_hadron)^2 exactly as SampleFinalState defines it.
  double mD = record.primary_mass;
  double ml = record.secondary_masses[1];
  double mK = record.secondary_masses[0];   // hadron mass actually sampled
  rk::P4 pD(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), mD);
  rk::P4 pK(geom3::Vector3(record.secondary_momenta[0][1], record.secondary_momenta[0][2], record.secondary_momenta[0][3]), mK);
  double q2 = (pD - pK).dot(pD - pK);

  // Build the hadron-mass mixture (masses + population fractions) that the
  // sampler uses, identify which component the recorded mK belongs to, and
  // normalize each component separately so the mixture weights are the true
  // population fractions.
  std::vector<double> masses;
  std::vector<double> fractions;
  if (is_Ds) {
    // Ds -> (eta / eta' / phi) l nu, fractions 0.46 / 0.16 / 0.38 (matches
    // the cumulative thresholds 2.3/5.0, 3.1/5.0 used in SampleFinalState).
    masses    = {siren::utilities::Constants::EtaMass,
                 siren::utilities::Constants::EtaPrimeMass,
                 siren::utilities::Constants::PhiMass};
    fractions = {2.3 / 5.0, 0.8 / 5.0, 1.9 / 5.0};
  } else {
    // D+/D0 -> K l nu with kinematic K/K*(892) mixing. fracK matches the
    // value used in SampleFinalState (PDG K vs K* semileptonic BRs).
    double mK_base = particleMass(signature.secondary_types[0]);
    double mKstar = KStarMass();
    double fracK;
    if (primary == siren::dataclasses::Particle::ParticleType::DPlus ||
        primary == siren::dataclasses::Particle::ParticleType::DMinus) {
      fracK = 8.74 / (8.74 + 5.33);
    } else {
      fracK = 3.41 / (3.41 + 2.17);
    }
    masses    = {mK_base, mKstar};
    fractions = {fracK, 1.0 - fracK};
  }

  // Identify the component matching the recorded hadron mass.
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

  // FinalStateProbability = mixture_fraction * (component density / component norm),
  // a true normalized q^2 pdf summing over the hadron-mass mixture.
  return fractions[comp] * g / norm;
}

std::vector<std::string> CharmMesonDecay::DensityVariables() const {
    return std::vector<std::string>{"Q2"};
}



} // namespace interactions
} // namespace siren

