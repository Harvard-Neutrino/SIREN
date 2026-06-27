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

namespace siren {
namespace interactions {

CharmMesonDecay3Body::CharmMesonDecay3Body() {
  // this is the default initialization but should never be used

  // we need to compute cdf here b/c otherwise SampleFinalState becomes non constant
  // in this case, we will need to compute all the possible dGamma's here
  // maybe add a map here, but for now we hard code first
  std::vector<double> constants;
  constants.resize(3);
  // check the primary and secondaries of the signature
  constants[0] = 0.725; // this is f^+(0)|V_cs| for charged D
  constants[1] = 0.44; // this is alpha, same for all K final states
  constants[2] = 2.01027; // this is excited charged D meson

  double mD = particleMass(siren::dataclasses::Particle::ParticleType::DPlus);
  double mK = particleMass(siren::dataclasses::Particle::ParticleType::K0Bar);

  computeDiffGammaCDF(constants, mD, mK);
}

CharmMesonDecay3Body::CharmMesonDecay3Body(siren::dataclasses::Particle::ParticleType primary) {

  //standard stuff, constant across primary types
  std::vector<double> constants;
  constants.resize(3);
  double mD;
  double mK;

  if (primary == siren::dataclasses::Particle::ParticleType::DPlus) {
    constants[0] = 0.725; // this is f^+(0)|V_cs| for charged D
    constants[1] = 0.44; // this is alpha, same for all K final states
    constants[2] = 2.01027; // this is excited charged D meson

    mD = particleMass(siren::dataclasses::Particle::ParticleType::DPlus);
    mK = particleMass(siren::dataclasses::Particle::ParticleType::K0Bar);

  } else if (primary == siren::dataclasses::Particle::ParticleType::D0) {
    constants[0] = 0.719; // this is f^+(0)|V_cs| for charged D
    constants[1] = 0.50; // this is alpha, same for all K final states
    constants[2] = 2.00697; // this is excited charged D meson

    mD = particleMass(siren::dataclasses::Particle::ParticleType::D0);
    mK = particleMass(siren::dataclasses::Particle::ParticleType::KMinus);
  }

  computeDiffGammaCDF(constants, mD, mK);

}

bool CharmMesonDecay3Body::equal(Decay const & other) const {
    const CharmMesonDecay3Body* x = dynamic_cast<const CharmMesonDecay3Body*>(&other);

    if(!x)
        return false;
    else
        return primary_types == x->primary_types;
}


double CharmMesonDecay3Body::particleMass(siren::dataclasses::ParticleType particle) {
    switch(particle){
			case siren::dataclasses::ParticleType::D0:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::D0Bar:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::DPlus:
				return( siren::utilities::Constants::DPlusMass);
			case siren::dataclasses::ParticleType::DMinus:
				return( siren::utilities::Constants::DPlusMass);
			case siren::dataclasses::ParticleType::K0:
				return( siren::utilities::Constants::K0Mass);
			case siren::dataclasses::ParticleType::K0Bar:
				return( siren::utilities::Constants::K0Mass);
			case siren::dataclasses::ParticleType::KPlus:
				return( siren::utilities::Constants::KPlusMass);
			case siren::dataclasses::ParticleType::KMinus:
				return( siren::utilities::Constants::KMinusMass);
      case siren::dataclasses::ParticleType::EPlus:
        return( siren::utilities::Constants::electronMass );
      case siren::dataclasses::ParticleType::EMinus:
        return( siren::utilities::Constants::electronMass );
      case siren::dataclasses::ParticleType::MuPlus:
        return( siren::utilities::Constants::muonMass );
      case siren::dataclasses::ParticleType::MuMinus:
        return( siren::utilities::Constants::muonMass );
      case siren::dataclasses::ParticleType::TauPlus:
        return( siren::utilities::Constants::tauMass );
      case siren::dataclasses::ParticleType::TauMinus:
        return( siren::utilities::Constants::tauMass );
      default:
        return(0.0);
    }
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
      std::cout << "this D meson decay has not been implemented yet" << std::endl;
    }
    return signatures;
}

std::vector<double> CharmMesonDecay3Body::FormFactorFromRecord(dataclasses::CrossSectionDistributionRecord const & record) const {
  dataclasses::InteractionSignature signature = record.signature;
  std::vector<double> constants;
  constants.resize(3);
  // check the primary and secondaries of the signature
  if (signature.primary_type == dataclasses::Particle::ParticleType::DPlus && signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::K0Bar) {
    constants[0] = 0.725; // this is f^+(0)|V_cs| for charged D
    constants[1] = 0.44; // this is alpha, same for all K final states
    constants[2] = 2.01027; // this is excited charged D meson
  } else if (signature.primary_type == dataclasses::Particle::ParticleType::D0 && signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::KMinus) {
    constants[0] = 0.719; // this is f^+(0)|V_cs| for neutral D
    constants[1] = 0.50; // this is alpha, same for all K final states
    constants[2] = 2.00697; // this is excited neutral D meson
  }
  return constants;
}

double CharmMesonDecay3Body::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    // first let the fully hadronic state be handled separately
    dataclasses::InteractionSignature signature = record.signature;
    if (signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::Hadrons) {
      return TotalDecayWidthForFinalState(record);
    }
    // get the form factor constants
    std::vector<double> constants = FormFactorFromRecord(record);
    // calculate the q^2
    rk::P4 pD(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::P4 pKPi(geom3::Vector3(record.secondary_momenta[0][1], record.secondary_momenta[0][2], record.secondary_momenta[0][3]), record.secondary_masses[0]);
    double Q2 = (pD - pKPi).dot(pD - pKPi);
    // primary and secondary masses are also needed
    double mD = record.primary_mass;
    double mK = record.secondary_masses[0];
    return DifferentialDecayWidth(constants, Q2, mD, mK);
}

double CharmMesonDecay3Body::DifferentialDecayWidth(std::vector<double> constants, double Q2, double mD, double mK) const {
    // get the numerical constants from the vector
    double F0CKM = constants[0];
    double alpha = constants[1];
    double ms = constants[2];
    double Q2tilde = Q2 / (ms * ms);
    // compute the 3-momentum as a function of Q2
    // double EK = 0.5 * (Q2 - pow(mD, 2) + pow(mK, 2)) / mD; // energy of Kaon
    double EK = 0.5 * (pow(mD, 2) + pow(mK, 2) - Q2) / mD; // energy of Kaon
    if (EK * EK < mK * mK) return 0.0;

    double PK = pow(pow(EK, 2) - pow(mK, 2), 0.5);
    // plug in the constants
    double dGamma = pow(siren::utilities::Constants::FermiConstant,2) / (24 * pow(siren::utilities::Constants::pi,3)) * pow(F0CKM,2) *
                    pow((1/((1-Q2tilde) * (1 - alpha * Q2tilde))),2) * pow(PK,3);
    return dGamma;
}

void CharmMesonDecay3Body::computeDiffGammaCDF(std::vector<double> constants, double mD, double mK) {

  // returns a 1D interpolater table for dGamma cdf
  // define the pdf with only Q2 as the input
  std::function<double(double)> pdf = [&] (double x) -> double {
            return DifferentialDecayWidth(constants, x, mD, mK);
        };
  // first normalize the integral
  double min = 0;
  double max = 1.4; // these set the min and max of the Q2 considered
  double normalization = siren::utilities::rombergIntegrate(pdf, min, max);
  std::function<double(double)> normed_pdf = [&] (double x) -> double {
            return DifferentialDecayWidth(constants, x, mD, mK) / normalization;
        };
  // now create the spline and compute the CDF

  // set the Q2 nodes (use 100 nodes)
  std::vector<double> Q2spline;
  for (int i = 0; i < 100; ++i) {
      Q2spline.push_back(0.01 + i * (max-min) / 100 );
  }

  // declare the cdf vectors
  std::vector<double> cdf_vector;
  std::vector<double> cdf_Q2_nodes;
  std::vector<double> pdf_vector;

  cdf_Q2_nodes.push_back(0);
  cdf_vector.push_back(0);
  pdf_vector.push_back(0);

  // compute the spline table
  for (int i = 0; i < Q2spline.size(); ++i) {
      if (i == 0) {
          double cur_Q2 = Q2spline[i];
          double cur_pdf = normed_pdf(cur_Q2);
          double area = cur_Q2 * cur_pdf * 0.5;
          pdf_vector.push_back(cur_pdf);
          cdf_vector.push_back(area);
          cdf_Q2_nodes.push_back(cur_Q2);
          continue;
      }
      double cur_Q2 = Q2spline[i];
      double cur_pdf = normed_pdf(cur_Q2);
      double area = 0.5 * (pdf_vector[i - 1] + cur_pdf) * (Q2spline[i] - Q2spline[i - 1]);
      pdf_vector.push_back(cur_pdf);
      cdf_Q2_nodes.push_back(cur_Q2);
      cdf_vector.push_back(area + cdf_vector.back());
  }

  cdf_Q2_nodes.push_back(max);
  cdf_vector.push_back(1);
  pdf_vector.push_back(0);

  // set the spline table
  siren::utilities::TableData1D<double> inverse_cdf_data;
  inverse_cdf_data.x = cdf_vector;
  inverse_cdf_data.f = cdf_Q2_nodes;

  inverseCdf = siren::utilities::Interpolator1D<double>(inverse_cdf_data);
  return;

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

    // =========================================================================
    // 3-body phase space sampling following Pythia's approach
    // (ParticleDecays::threeBody in ParticleDecays.cc)
    //
    // D (m0) -> K (m1) + lepton (m2) + neutrino (m3)
    //
    // Phase space: sample m23 (lepton-neutrino invariant mass = sqrt(q^2))
    // flat in allowed range, accept-reject on phase space weight.
    // Then apply V-A matrix element correction.
    //
    // K / K*(892) mixing: a fraction (1 - fracK) of events use mK*=0.892 GeV
    // as the hadron mass instead of the pseudoscalar K mass, which broadens
    // the lepton spectrum through phase-space and matrix-element effects.
    //
    // NOTE (design decision): this is *kinematic* K/K* mixing only. The
    // secondary particle's ParticleType is always left at the K species from
    // the signature (K0bar / K-), even when the K* mass was drawn. We do not
    // advertise separate K* signatures in GetPossibleSignaturesFromParent()
    // and we do not re-type the secondary. This is intentional -- for our use
    // case (weighting lepton/neutrino kinematics correctly in the presence of
    // the resonant K* contribution) the mass treatment is what matters, and
    // downstream propagation treats the secondary as a pseudoscalar K.
    // =========================================================================

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

// ---------------------------------------------------------------------------
// Closure-correct FinalStateProbability machinery (see CharmMesonDecay.cxx for
// the full derivation). FinalStateProbability MUST be the normalized q^2
// density that SampleFinalState produces: the sampler draws m23 = sqrt(q^2)
// flat with accept-reject on p1Abs*p23Abs, then accept-rejects on the V-A
// weight wtME = mD*E_l*(p_nu . p_K). The accepted m23 density is therefore
// proportional to p1Abs*p23Abs*<max(0,wtME)>_angle, and the q^2 density carries
// the Jacobian 1/(2 m23). Both the numerator and the per-component normalizer
// use the same SampledQ2Density, so Sample == Density by construction. The
// form-factor DifferentialDecayWidth is a separate quantity the sampler never
// used and is NOT the closure density.
// ---------------------------------------------------------------------------

double CharmMesonDecay3Body::KStarMass() {
  // Single source of truth for the K*(892) mass shared by the sampler and the
  // FinalStateProbability normalizer.
  return siren::utilities::Constants::KPrimePlusMass;
}

namespace {
// Integral over [lo, hi] of max(0, a2*c^2 + a1*c + a0). Closed form with a
// linear/constant fallback for a2 ~ 0 and numerically stable quadratic roots.
double integratePositivePartQuadratic(double a2, double a1, double a0,
                                      double lo, double hi) {
    if (hi <= lo) return 0.0;
    auto F = [&](double c) { return a2 * c * c * c / 3.0 + a1 * c * c / 2.0 + a0 * c; };
    auto seg = [&](double x0, double x1) -> double {
        if (x1 <= x0) return 0.0;
        return F(x1) - F(x0);
    };
    double scale = std::abs(a1) + std::abs(a0) + 1.0;
    if (std::abs(a2) < 1e-12 * scale) {
        if (std::abs(a1) < 1e-300) return (a0 > 0.0) ? seg(lo, hi) : 0.0;
        double root = -a0 / a1;
        if (a1 > 0.0) return seg(std::max(lo, root), hi);
        return seg(lo, std::min(hi, root));
    }
    double disc = a1 * a1 - 4.0 * a2 * a0;
    if (disc <= 0.0) {
        return (a2 > 0.0) ? seg(lo, hi) : 0.0;
    }
    double sq = std::sqrt(disc);
    double qq = -0.5 * (a1 + ((a1 >= 0.0) ? sq : -sq));
    double r1 = qq / a2;
    double r2 = a0 / qq;
    double rlo = std::min(r1, r2), rhi = std::max(r1, r2);
    if (a2 > 0.0) {
        return seg(lo, std::min(hi, rlo)) + seg(std::max(lo, rhi), hi);
    }
    return seg(std::max(lo, rlo), std::min(hi, rhi));
}
} // namespace

double CharmMesonDecay3Body::VAWeightAngleAverage(double mD, double mK, double ml, double m23) const {
    // Analytic angle-average of the sampler's ACCEPTED V-A weight,
    //   (1/2) * integral_{-1}^{1} clamp(wtME(c), 0, wtMEmax) dc.
    // wtME(c) = mD * E_l(c) * (p_nu(c) . p_K) is an exact quadratic in
    // c = cos(theta_lepton) after the boost to the D rest frame, so the positive
    // and clipped parts integrate in closed form via
    // clamp(q, 0, M) = max(0, q) - max(0, q - M). Matches SampleFinalState's
    // accepted density exactly (closure) with no quadrature error; the same
    // integral is also evaluated numerically in the unit tests as a cross-check.
    double mnu = 0.0;
    double p1Abs = 0.5 * std::sqrt((mD - mK - m23) * (mD + mK + m23)
                                 * (mD + mK - m23) * (mD - mK + m23)) / mD;
    double p23Abs = 0.5 * std::sqrt((m23 - ml - mnu) * (m23 + ml + mnu)
                                  * (m23 + ml - mnu) * (m23 - ml + mnu)) / m23;
    if (p1Abs <= 0.0 || p23Abs <= 0.0) return 0.0;
    double E23 = std::sqrt(p1Abs * p1Abs + m23 * m23);
    double bz = -p1Abs / E23;
    double gamma = E23 / m23;
    double Elrest = std::sqrt(p23Abs * p23Abs + ml * ml);
    double EK = std::sqrt(p1Abs * p1Abs + mK * mK);
    double C = Elrest;
    double D = bz * p23Abs;
    double A = EK - p1Abs * bz;
    double B = p1Abs - EK * bz;
    double pref = mD * gamma * gamma * p23Abs;
    if (pref <= 0.0) return 0.0;
    double a0 = C * A;
    double a1 = C * B + D * A;
    double a2 = D * B;
    double wtMEmax = std::min(std::pow(mD, 4) / 16.0,
                              mD * (mD - mK - ml) * (mD - mK - mnu) * (mD - ml - mnu));
    double qmax = wtMEmax / pref;
    double integral = integratePositivePartQuadratic(a2, a1, a0, -1.0, 1.0)
                    - integratePositivePartQuadratic(a2, a1, a0 - qmax, -1.0, 1.0);
    return 0.5 * pref * integral;
}

double CharmMesonDecay3Body::SampledQ2Density(double mD, double mK, double ml, double q2, bool apply_va) const {
  double m23 = std::sqrt(q2);
  double m23Max = mD - mK;
  if (m23 <= ml || m23 >= m23Max) return 0.0;
  double mnu = 0.0;
  double p1Abs = 0.5 * std::sqrt((mD - mK - m23) * (mD + mK + m23)
                               * (mD + mK - m23) * (mD - mK + m23)) / mD;
  double p23Abs = 0.5 * std::sqrt((m23 - ml - mnu) * (m23 + ml + mnu)
                                * (m23 + ml - mnu) * (m23 - ml + mnu)) / m23;
  if (p1Abs <= 0.0 || p23Abs <= 0.0) return 0.0;
  // Phase-space weight WITH the sampler's rejection ceiling. SampleFinalState
  // accept-rejects m23 against wtPSmax = 0.5*p1Max*p23Max, which p1Abs*p23Abs can
  // exceed; there the sampler saturates (always accepts), so the accepted
  // density is flat at wtPSmax. Reproduce the clip for exact closure.
  double m23Min = ml + mnu;
  double p1Max = 0.5 * std::sqrt((mD - mK - m23Min) * (mD + mK + m23Min)
                               * (mD + mK - m23Min) * (mD - mK + m23Min)) / mD;
  double p23Max = 0.5 * std::sqrt((m23Max - ml - mnu) * (m23Max + ml + mnu)
                                * (m23Max + ml - mnu) * (m23Max - ml + mnu)) / m23Max;
  double wtPSmax = 0.5 * p1Max * p23Max;
  double wtPS = p1Abs * p23Abs;
  if (wtPS > wtPSmax) wtPS = wtPSmax;
  double me = apply_va ? VAWeightAngleAverage(mD, mK, ml, m23) : 1.0;
  return wtPS * me / (2.0 * m23);
}

double CharmMesonDecay3Body::SampledQ2Normalization(double mD, double mK, double ml, bool apply_va) const {
  // Cached per (mD, mK, ml, apply_va) so per-event FinalStateProbability does
  // not re-integrate the fixed-kinematics normalizer.
  double key = ((mD * 1e3 + mK) * 1e3 + ml) * 2.0 + (apply_va ? 1.0 : 0.0);
  long ikey = (long)std::llround(key * 1e3);
  auto it = norm_cache.find(ikey);
  if (it != norm_cache.end()) return it->second;
  double q2min = ml * ml;
  double q2max = (mD - mK) * (mD - mK);
  std::function<double(double)> integrand = [&](double q2) -> double {
    return SampledQ2Density(mD, mK, ml, q2, apply_va);
  };
  double n = siren::utilities::rombergIntegrate(integrand, q2min, q2max);
  norm_cache[ikey] = n;
  return n;
}

double CharmMesonDecay3Body::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  dataclasses::InteractionSignature signature = record.signature;
  // Fully hadronic catch-all: deterministic final state, density 1.
  if (signature.secondary_types[0] == siren::dataclasses::Particle::ParticleType::Hadrons) {
    return 1.0;
  }

  siren::dataclasses::Particle::ParticleType primary = signature.primary_type;
  bool apply_va = true;   // 3-body class always applies the V-A matrix element

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

  double g = SampledQ2Density(mD, masses[comp], ml, q2, apply_va);
  if (g <= 0.0) return 0.0;
  double norm = SampledQ2Normalization(mD, masses[comp], ml, apply_va);
  if (norm <= 0.0) return 0.0;

  // Normalized q^2 pdf summing over the K/K* mixture.
  return fractions[comp] * g / norm;
}

std::vector<std::string> CharmMesonDecay3Body::DensityVariables() const {
    return std::vector<std::string>{"Q2"};
}



} // namespace interactions
} // namespace siren

