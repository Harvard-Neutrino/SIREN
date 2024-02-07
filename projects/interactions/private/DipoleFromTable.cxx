#include "LeptonInjector/interactions/DipoleFromTable.h"

#include <set>                                                // for set
#include <array>                                              // for array
#include <cmath>                                              // for pow, sqrt
#include <tuple>                                              // for tie
#include <cstdio>                                             // for snprintf
#include <memory>                                             // for allocator
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <fstream>
#include <iomanip>                                            // for operator<<
#include <utility>                                            // for pair
#include <assert.h>                                           // for assert
#include <iostream>                                           // for basic_i...
#include <stdint.h>                                           // for int32_t
#include <stdlib.h>                                           // for abs
#include <iterator>                                           // for back_in...
#include <algorithm>                                          // for find, max, min, set_int...
#include <functional>                                         // for function

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "LeptonInjector/interactions/CrossSection.h"        // for CrossSe...
#include "LeptonInjector/dataclasses/InteractionRecord.h"     // for Interac...
#include "LeptonInjector/dataclasses/InteractionSignature.h"  // for Interac...
#include "LeptonInjector/dataclasses/Particle.h"              // for Particle
#include "LeptonInjector/detector/MaterialModel.h"            // for Materia...
#include "LeptonInjector/utilities/Constants.h"               // for invGeVs...
#include "LeptonInjector/utilities/Errors.h"                  // for Injecti...
#include "LeptonInjector/utilities/Interpolator.h"            // for Interpo...
#include "LeptonInjector/utilities/Random.h"                  // for LI_random

namespace LI {
namespace interactions {

namespace {
bool fexists(const std::string filename) {
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}
}

bool DipoleFromTable::equal(CrossSection const & other) const {
    const DipoleFromTable* x = dynamic_cast<const DipoleFromTable*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(
                    z_samp,
                    primary_types,
                    hnl_mass,
                    channel,
                    differential,
                    total)
            ==
            std::tie(
                    x->z_samp,
                    x->primary_types,
                    x->hnl_mass,
                    x->channel,
                    x->differential,
                    x->total);
}

double DipoleFromTable::DipoleyMin(double Enu, double mHNL, double target_mass) {
    double yMin = 0;
    double target_mass2 = target_mass * target_mass;
    double mHNL2 = mHNL * mHNL;

    double s = 2 * Enu * target_mass + target_mass2;
    double s2 = s * s;

    double r2 = mHNL2 / s;
    double r4 = mHNL2 * mHNL2 / s2;
    double m2 = target_mass2 / s;
    double m4 = target_mass2 * target_mass2 / s2;
    double m2sub1sq = std::pow(m2 - 1, 2);
    bool small_r = r2 < 1e-6;

    // require cos(theta)<=1
    double root_term = 4*Enu*Enu*target_mass2 - 4*Enu*target_mass*mHNL2 - 4*target_mass2*mHNL2 + mHNL2*mHNL2;
    double costh_bound = (1./(2*s)) * (2*Enu*target_mass - mHNL2 - target_mass*mHNL2/Enu - std::sqrt(root_term));

    if(small_r)
        yMin = (s * m2 * r4 / m2sub1sq) / (2 * Enu * target_mass);
    else {
        double root = std::sqrt(m2sub1sq + (r4 - 2 * (1 + m2) * r2));
        yMin = 0.5 * (1 + m4 - r2 - root + m2 * (-2 - r2 + root)) * s / (2 * Enu * target_mass);
    }
    return std::max(yMin,costh_bound);
}


double DipoleFromTable::DipoleyMax(double Enu, double mHNL, double target_mass) {
    double target_mass2 = target_mass * target_mass;
    double target_mass4 = target_mass2 * target_mass2;
    double mHNL2 = mHNL * mHNL;

    double s = 2 * Enu * target_mass + target_mass2;
    double s2 = s * s;

    // require cos(theta)>=-1
    double root_term = 4*Enu*Enu*target_mass2 - 4*Enu*target_mass*mHNL2 - 4*target_mass2*mHNL2 + mHNL2*mHNL2;
    double costh_bound = (1./(2*s)) * (2*Enu*target_mass - mHNL2 - target_mass*mHNL2/Enu) + std::sqrt(root_term);

    return std::min(costh_bound,0.5 * (target_mass4 - mHNL2*s + s2 - target_mass2*(mHNL2+2*s) + (s - target_mass2) * std::sqrt(target_mass4 + std::pow(mHNL2 - s, 2) - 2*target_mass2*(mHNL2+s))) / (s * (2 * Enu * target_mass)));
}


double DipoleFromTable::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    LI::dataclasses::Particle::ParticleType const & primary_type = interaction.signature.primary_type;
    LI::dataclasses::Particle::ParticleType const & target_type = interaction.signature.target_type;
    std::array<double, 4> const & primary_momentum = interaction.primary_momentum;
    double primary_mass = interaction.primary_mass;
    rk::P4 p1(geom3::Vector3(primary_momentum[1], primary_momentum[2], primary_momentum[3]), primary_mass);
    double const & primary_energy = primary_momentum[0];
    // if we are below threshold, return 0
    if(primary_energy < InteractionThreshold(interaction))
        return 0;
    return TotalCrossSection(primary_type, primary_energy, target_type);
}

double DipoleFromTable::TotalCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, LI::dataclasses::Particle::ParticleType target_type) const {
    if(not primary_types.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }

    if(total.find(target_type) == total.end()) {
        std::cout << "Faulty target: " << target_type << std::endl;
        throw std::runtime_error("Supplied target not supported by cross section!");
    }

    LI::utilities::Interpolator1D<double> const & interp = total.at(target_type);

    if(primary_energy < interp.MinX() or primary_energy > interp.MaxX()) {
        throw std::runtime_error("Interaction energy ("+ std::to_string(primary_energy) +
                ") out of cross section table range: ["
                + std::to_string(interp.MinX()) + " GeV,"
                + std::to_string(interp.MaxX()) + " GeV]");
    }

    LI::utilities::Interpolator1D<double> const & interp_proton = total.at(LI::dataclasses::Particle::ParticleType::HNucleus);
    int nprotons = LI::detector::MaterialModel::GetProtonCount(target_type);
    if(!inelastic || target_type==LI::dataclasses::Particle::ParticleType::HNucleus) {
        nprotons = 0;
    }
    double proton_inelastic_xsec = 0;
    if(primary_energy > interp_proton.MinX() and primary_energy < interp_proton.MaxX()) {
        proton_inelastic_xsec = interp_proton(primary_energy);
    }
    double xsec = interp(primary_energy) + nprotons*proton_inelastic_xsec;

    if(in_invGeV)
        return std::pow(dipole_coupling, 2) * xsec / LI::utilities::Constants::invGeVsq_per_cmsq;
    else
        return std::pow(dipole_coupling, 2) * xsec;
}

double DipoleFromTable::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    LI::dataclasses::Particle::ParticleType primary_type = interaction.signature.primary_type;
    LI::dataclasses::Particle::ParticleType target_type = interaction.signature.target_type;
    std::array<double, 4> const & primary_momentum = interaction.primary_momentum;
    double primary_mass = interaction.primary_mass;
    double target_mass = interaction.primary_mass;
    rk::P4 p1(geom3::Vector3(primary_momentum[1], primary_momentum[2], primary_momentum[3]), primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), target_mass);
    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    primary_energy = primary_momentum[0];
    p1_lab = p1;
    p2_lab = p2;
    std::vector<LI::dataclasses::Particle::ParticleType> const & secondary_types = interaction.signature.secondary_types;
    assert(secondary_types.size() == 2);
    assert(secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4 or secondary_types[1] == LI::dataclasses::Particle::ParticleType::NuF4 or secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4Bar or secondary_types[1] == LI::dataclasses::Particle::ParticleType::NuF4Bar);
    unsigned int lepton_index = (secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4 or secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4Bar) ? 0 : 1;
    unsigned int other_index = 1 - lepton_index;

    std::array<double, 4> const & mom3 = interaction.secondary_momenta.at(lepton_index);
    std::array<double, 4> const & mom4 = interaction.secondary_momenta.at(other_index);
    rk::P4 p3(geom3::Vector3(mom3[1], mom3[2], mom3[3]), interaction.secondary_masses.at(lepton_index));
    rk::P4 p4(geom3::Vector3(mom4[1], mom4[2], mom4[3]), interaction.secondary_masses.at(other_index));

    double y = 1.0 - p2.dot(p3) / p2.dot(p1);

    double thresh = InteractionThreshold(interaction);

    return DifferentialCrossSection(primary_type, primary_energy, target_type, target_mass, y, thresh);
}

double DipoleFromTable::DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, LI::dataclasses::Particle::ParticleType target_type, double target_mass, double y) const {
    // Assume threshold is first entry of table
    LI::utilities::Interpolator2D<double> const & interp = differential.at(target_type);
    double thresh = interp.MinX();
    return DifferentialCrossSection(primary_type, primary_energy, target_type, target_mass, y, thresh);
}

double DipoleFromTable::DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, LI::dataclasses::Particle::ParticleType target_type, double target_mass, double y, double thresh) const {
    if(not primary_types.count(primary_type))
        return 0.0;

    if(total.find(target_type) == total.end())
        return 0.0;

    LI::utilities::Interpolator2D<double> const & interp = differential.at(target_type);

    LI::utilities::Interpolator2D<double> const & interp_proton = differential.at(LI::dataclasses::Particle::ParticleType::HNucleus);

    int nprotons = LI::detector::MaterialModel::GetProtonCount(target_type);
    if(!inelastic || target_type==LI::dataclasses::Particle::ParticleType::HNucleus) {
        nprotons = 0;
    }

    if(primary_energy < thresh or primary_energy > interp.MaxX())
        return 0.0;


    double yMin = DipoleyMin(primary_energy, GetHNLMass(), target_mass);
    double yMax = DipoleyMax(primary_energy, GetHNLMass(), target_mass);
    double differential_cross_section = 0.0;
    if(y < yMin or y > yMax)
        return 0.0;
    if(z_samp) {
        double z = (y - yMin) / (yMax - yMin);
        if(z < interp.MinY() or z > interp.MaxY())
            return 0.0;
        differential_cross_section = interp(primary_energy, z) + nprotons*interp_proton(primary_energy,z);
    } else {
        if(y < interp.MinY() or y > interp.MaxY())
            return 0.0;
        differential_cross_section = interp(primary_energy, y) + nprotons*interp_proton(primary_energy,y);
    }
    if(in_invGeV)
        differential_cross_section /= LI::utilities::Constants::invGeVsq_per_cmsq;
    return std::pow(dipole_coupling, 2) * differential_cross_section;
}

double DipoleFromTable::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    return hnl_mass + (hnl_mass*hnl_mass)/(2*interaction.target_mass);
}

void DipoleFromTable::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<LI::utilities::LI_random> random) const {
    LI::utilities::Interpolator2D<double> const & diff_table = differential.at(record.GetTargetType());
    LI::utilities::Interpolator2D<double> const & diff_table_proton = differential.at(LI::dataclasses::Particle::ParticleType::HNucleus);
    int nprotons = LI::detector::MaterialModel::GetProtonCount(record.GetTargetType());
    // Avoid double counting for true H nuclei
    if(!inelastic || record.signature.target_type == LI::dataclasses::Particle::ParticleType::HNucleus) {
        nprotons = 0;
    }

    // Uses Metropolis-Hastings Algorithm!
    // useful for cases where we don't know the supremum of our distribution, and the distribution is multi-dimensional

    std::array<double, 4> const & primary_momentum = record.GetPrimaryMomentum();
    double primary_mass = record.GetPrimaryMass();
    double target_mass = record.GetTargetMass();

    rk::P4 p1(geom3::Vector3(primary_momentum[1], primary_momentum[2], primary_momentum[3]), primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * target_mass_ + 2 * target_mass_ * primary_energy;
    double s = std::pow(rk::invMass(p1, p2), 2);

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    p1_lab = p1;
    p2_lab = p2;
    primary_energy = p1_lab.e();

    std::vector<LI::dataclasses::Particle::ParticleType> const & secondary_types = record.record.signature.secondary_types;

    unsigned int lepton_index = (secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4 or secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4Bar) ? 0 : 1;
    unsigned int other_index = 1 - lepton_index;
    double m = hnl_mass;
    double thresh = InteractionThreshold(record.record);
    if(primary_energy < thresh) {
        throw(LI::utilities::InjectionFailure("Primary is below interaction threshold!"));
    }

    // double m2 = p2_lab | p2_lab;
    double m2 = target_mass;
    double m3 = m;
    double E1_lab = p1_lab.e();
    double E2_lab = p2_lab.e();

    double yMin = DipoleyMin(E1_lab, GetHNLMass(), target_mass);
    double yMax = DipoleyMax(E1_lab, GetHNLMass(), target_mass);
    assert(yMin > 0);
    double log_yMax = log10(yMax);
    double log_yMin = log10(yMin);
    double min_Q2 = yMin * (s - m2 * m2);
    double z; // placeholder for z sampling

    bool accept;

    // kin_vars and its twin are 2-vectors containing [nu-energy, Bjorken Y]
    std::array<double,2> kin_vars, test_kin_vars;

    // values of cross_section from the table
    double cross_section, test_cross_section;

    // No matter what, we're evaluating at this specific energy.
    kin_vars[0] = test_kin_vars[0] = primary_energy;

    // check preconditions
    if(kin_vars[0] < thresh
            || kin_vars[0] > diff_table.MaxX())
        throw std::runtime_error("Sample: Interaction energy out of differential cross section table range: ["
                + std::to_string(diff_table.MinX()) + " GeV,"
                + std::to_string(diff_table.MaxX()) + " GeV]");

    // sample an intial point
    do {
        // rejection sample a point which is kinematically allowed by calculation limits
        double trialQ;
        do {
            kin_vars[1] = std::pow(10.0, random->Uniform(log_yMin, log_yMax));
            trialQ = (2 * E1_lab * E2_lab) * kin_vars[1];
        } while(trialQ < min_Q2);

        accept = true;
        //sanity check: demand that the sampled point be within the table extents
        z = (kin_vars[1]-yMin)/(yMax-yMin);
        if((!z_samp && (kin_vars[1] < diff_table.MinY() || kin_vars[1] > diff_table.MaxY()))
                || (z_samp && (z < diff_table.MinY() || z > diff_table.MaxY()))){
            accept = false;
        }
    } while(!accept);

    // Bx * By * xs(E, x, y)
    // evalutates the differential spline at that point
    if(z_samp) test_cross_section = diff_table(kin_vars[0], z);
    else test_cross_section = diff_table(kin_vars[0], kin_vars[1]);
    if(kin_vars[0] > diff_table_proton.MinX() and kin_vars[0] < diff_table_proton.MaxX()) {
        if(z_samp) {
            if(z > diff_table_proton.MinY() and z < diff_table_proton.MaxY()) {
                test_cross_section += nprotons*diff_table_proton(kin_vars[0],z);
            }
        }
        else {
            if(kin_vars[1] > diff_table_proton.MinY() and kin_vars[1] < diff_table_proton.MaxY()) {
                test_cross_section += nprotons*diff_table_proton(kin_vars[0],kin_vars[1]);
            }
        }
    }

    cross_section = test_cross_section;

    // this is the magic part. Metropolis Hastings Algorithm.
    // MCMC method!
    const size_t burnin = 40; // converges to the correct distribution over multiple samplings.
                              // big number means more accurate, but slower
    for(size_t j = 0; j <= burnin; j++) {
        // repeat the sampling from above to get a new valid point
        double trialQ;
        do {
            test_kin_vars[1] = std::pow(10.0, random->Uniform(log_yMin, log_yMax));
            trialQ = (2 * E1_lab * E2_lab) * test_kin_vars[1];
        } while(trialQ < min_Q2);

        accept = true;
        //sanity check: demand that the sampled point be within the table extents
        z = (test_kin_vars[1]-yMin)/(yMax-yMin);
        if((!z_samp && (test_kin_vars[1] < diff_table.MinY() || test_kin_vars[1] > diff_table.MaxY()))
                || (z_samp && (z < diff_table.MinY() || z > diff_table.MaxY()))){
            accept = false;
        }
        if(!accept)
            continue;

        // Load the differential cross section depending on sampling variable
        if(z_samp) test_cross_section = diff_table(test_kin_vars[0], z);
        else test_cross_section = diff_table(test_kin_vars[0], test_kin_vars[1]);
        if(test_kin_vars[0] > diff_table_proton.MinX() and test_kin_vars[0] < diff_table_proton.MaxX()) {
            if(z_samp) {
                if(z > diff_table_proton.MinY() and z < diff_table_proton.MaxY()) {
                    test_cross_section += nprotons*diff_table_proton(test_kin_vars[0],z);
                }
            }
            else {
                if(test_kin_vars[1] > diff_table_proton.MinY() and test_kin_vars[1] < diff_table_proton.MaxY()) {
                    test_cross_section += nprotons*diff_table_proton(test_kin_vars[0],test_kin_vars[1]);
                }
            }
        }
        if(std::isnan(test_cross_section) or test_cross_section <= 0)
            continue;

        //double odds = ((test_kin_vars[1]*test_cross_section) / (kin_vars[1]*cross_section));
        double odds = (test_cross_section / cross_section); // this gives a better match to the y distribution
        accept = (cross_section == 0 || (odds > 1.) || random->Uniform(0, 1) < odds);

        if(accept) {
            kin_vars = test_kin_vars;
            cross_section = test_cross_section;
        }
    }
    double final_y = kin_vars[1] + 1e-16; // to account for machine epsilon when adding to O(1) numbers

    std::map<std::string, double> params;
    params["energy"] = E1_lab;
    params["bjorken_y"] = final_y;
    record.SetInteractionParameters(params);

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 p1_mom = p1_lab.momentum();
    geom3::UnitVector3 p1_lab_dir = p1_mom.direction();
    geom3::Rotation3 x_to_p1_lab_rot = geom3::rotationBetween(x_dir, p1_lab_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(p1_lab_dir, phi);

    double E3_lab = E1_lab * (1.0 - final_y);
    //double p1x_lab = p1_mom.length();
    double p3_lab_sq = E3_lab * E3_lab - m3 * m3;
    //double p3x_lab_frac = (p1x_lab * p1x_lab - m3 * m3 + E1_lab * E1_lab * (1.0 - 2.0 * final_y)) / (2.0 * p1x_lab * E3_lab);
    //double p3x_lab = p3x_lab_frac * sqrt(p3_lab_sq);
    double p3x_lab = E3_lab - m2*final_y - m3*m3/(2*E1_lab);
    double p3y_lab = sqrt(p3_lab_sq - p3x_lab * p3x_lab);

    rk::P4 p3_lab(geom3::Vector3(p3x_lab, p3y_lab, 0), m3);
    p3_lab.rotate(x_to_p1_lab_rot);
    p3_lab.rotate(rand_rot);
    rk::P4 p4_lab = p2_lab + (p1_lab - p3_lab);

    rk::P4 p3;
    rk::P4 p4;
    p3 = p3_lab;
    p4 = p4_lab;

    double helicity_mul = 0.0;
    if(channel == Conserving)
        helicity_mul = 1.0;
    else if(channel == Flipping)
        helicity_mul = -1.0;

    std::vector<LI::dataclasses::Particle> secondaries = record.GetSecondaryParticles();
    secondaries[lepton_index].type = secondary_types[lepton_index];
    secondaries[lepton_index].mass = p3.m();
    secondaries[lepton_index].momentum = {p3.e(), p3.px(), p3.py(), p3.pz()};
    secondaries[lepton_index].length = 0;
    secondaries[lepton_index].helicity = std::copysign(0.5, record.GetPrimaryHelicity() * helicity_mul);

    secondaries[other_index].type = secondary_types[other_index];
    secondaries[other_index].mass = p4.m();
    secondaries[other_index].momentum = {p4.e(), p4.px(), p4.py(), p4.pz()};
    secondaries[other_index].length = 0;
    secondaries[other_index].helicity = std::copysign(0.5, record.GetPrimaryHelicity() * helicity_mul);

    record.SetSecondaryParticles(secondaries);
}

double DipoleFromTable::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else if (txs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<LI::dataclasses::Particle::ParticleType> DipoleFromTable::GetPossibleTargets() const {
    std::set<LI::dataclasses::Particle::ParticleType> diff_targets;
    std::set<LI::dataclasses::Particle::ParticleType> tot_targets;
    for(auto const & diff : differential)
        diff_targets.insert(diff.first);
    for(auto const & tot : total)
        tot_targets.insert(tot.first);
    std::vector<LI::dataclasses::Particle::ParticleType> res;
    std::set_intersection(diff_targets.begin(), diff_targets.end(), tot_targets.begin(), tot_targets.end(), std::back_inserter(res));
    return res;
}

std::vector<LI::dataclasses::Particle::ParticleType> DipoleFromTable::GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const {
    if(not primary_types.count(primary_type)) {
        return std::vector<LI::dataclasses::Particle::ParticleType>();
    }
    return GetPossibleTargets();
}

std::vector<LI::dataclasses::Particle::ParticleType> DipoleFromTable::GetPossiblePrimaries() const {
    return std::vector<LI::dataclasses::Particle::ParticleType>(primary_types.begin(), primary_types.end());
}

std::vector<dataclasses::InteractionSignature> DipoleFromTable::GetPossibleSignatures() const {
    std::vector<LI::dataclasses::Particle::ParticleType> targets = GetPossibleTargets();
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.secondary_types.resize(2);

    for(auto primary : primary_types) {
        signature.primary_type = primary;
        if(std::set<LI::dataclasses::Particle::ParticleType>{LI::dataclasses::Particle::ParticleType::NuE, LI::dataclasses::Particle::ParticleType::NuMu, LI::dataclasses::Particle::ParticleType::NuTau}.count(primary))
            signature.secondary_types[0] = LI::dataclasses::Particle::ParticleType::NuF4;
        else if(std::set<LI::dataclasses::Particle::ParticleType>{LI::dataclasses::Particle::ParticleType::NuEBar, LI::dataclasses::Particle::ParticleType::NuMuBar, LI::dataclasses::Particle::ParticleType::NuTauBar}.count(primary))
            signature.secondary_types[0] = LI::dataclasses::Particle::ParticleType::NuF4Bar;
        else
            throw std::runtime_error("Primary type not in primary_types!");
        for(auto target : targets) {
            signature.target_type = target;
            signature.secondary_types[1] = target;
            signatures.push_back(signature);
        }
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> DipoleFromTable::GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const {
    std::vector<LI::dataclasses::Particle::ParticleType> targets = GetPossibleTargets();
    if(primary_types.count(primary_type) > 0 and std::find(targets.begin(), targets.end(), target_type) != targets.end()) {
        dataclasses::InteractionSignature signature;
        signature.secondary_types.resize(2);
        signature.primary_type = primary_type;
        signature.target_type = target_type;
        signature.secondary_types[1] = target_type;
        if(std::set<LI::dataclasses::Particle::ParticleType>{LI::dataclasses::Particle::ParticleType::NuE, LI::dataclasses::Particle::ParticleType::NuMu, LI::dataclasses::Particle::ParticleType::NuTau}.count(primary_type))
            signature.secondary_types[0] = LI::dataclasses::Particle::ParticleType::NuF4;
        else if(std::set<LI::dataclasses::Particle::ParticleType>{LI::dataclasses::Particle::ParticleType::NuEBar, LI::dataclasses::Particle::ParticleType::NuMuBar, LI::dataclasses::Particle::ParticleType::NuTauBar}.count(primary_type))
            signature.secondary_types[0] = LI::dataclasses::Particle::ParticleType::NuF4Bar;
        else
            throw std::runtime_error("Primary type not in primary_types!");
        return std::vector<dataclasses::InteractionSignature>{signature};
    } else {
        return std::vector<dataclasses::InteractionSignature>();
    }
}

std::vector<std::string> DipoleFromTable::DensityVariables() const {
    return std::vector<std::string>{"Bjorken y"};
}

void DipoleFromTable::AddDifferentialCrossSectionFile(std::string filename, LI::dataclasses::Particle::ParticleType target) {
    std::string delimeter = "_";
    std::string end_delimeter = ".";
    std::string::size_type pos = filename.rfind("/") + 1;
    if(pos == std::string::npos) {
        pos = 0;
    }
    std::string::size_type next_pos = pos;
    std::string::size_type sub_len = 0;
    unsigned int Z = 0;
    unsigned int A = 0;
    bool bad = false;

    std::function<std::string()> next_substr = [&] () -> std::string {
        if(pos >= filename.size() or pos == std::string::npos) {
            bad = true;
            return std::string();
        }
        next_pos = filename.find(delimeter, pos);
        if(next_pos == std::string::npos) {
            next_pos = filename.rfind(end_delimeter, pos);
            if(next_pos == std::string::npos) {
                bad = true;
                return std::string();
            }
        }
        sub_len = std::max((int)(next_pos - pos), 0);
        next_pos = pos + sub_len;
        std::string sub = filename.substr(pos, sub_len);
        pos = next_pos + 1;
        return sub;
    };

    while(pos < filename.size() and pos != std::string::npos) {
        std::string sub = next_substr();
        if(bad)
            break;
        if(sub == "Z") {
            sub = next_substr();
            if(bad)
                break;
            Z = std::stoi(sub);
        } else if(sub == "A") {
            sub = next_substr();
            if(bad)
                break;
            A = std::stoi(sub);
        } else if(sub == "mHNL") {
            sub = next_substr();
            if(bad)
                break;
            double file_hnl_mass = std::stod(sub);
            if(std::abs(file_hnl_mass - hnl_mass) / std::max(std::abs(file_hnl_mass), std::abs(hnl_mass)) > 1e-6) {
                std::cout << std::setprecision(24);
                std::cout << "File HNL mass: "<< file_hnl_mass << std::endl;
                std::cout << "Specified HNL mass: "<< hnl_mass << std::endl;
                throw std::runtime_error("File HNL mass does not match specified HNL mass!");
            }
        }
    }
    //std::string pid_str = std::format("10{0:0>1d}{1:0>3d}{2:0>3d}{3:0>1d}", 0, Z, A, 0);
    unsigned int buffer_size = 1024;
    char buffer[buffer_size];
    unsigned int str_size = std::snprintf(buffer, buffer_size, "10%01d%03d%03d%01d", 0, Z, A, 0);
    if(str_size > 10) {
        throw std::runtime_error("Cannot create particle ID string!");
    }
    std::string pid_str(buffer);
    int32_t pid_int = std::stoul(pid_str);
    LI::dataclasses::Particle::ParticleType pid = (LI::dataclasses::Particle::ParticleType)pid_int;

    if(pid != target) {
        throw std::runtime_error("File target nucleus does not match supplied target type!");
    }

    if(fexists(filename)) {
        std::ifstream in(filename.c_str());
        std::string buf;

        LI::utilities::TableData2D<double> table_data;
        while(std::getline(in, buf)) {
            if((pos = buf.find('#')) != std::string::npos)
                buf.erase(pos);
            const char* whitespace=" \n\r\t\v";
            if((pos=buf.find_first_not_of(whitespace))!=0)
                buf.erase(0,pos);
            if(!buf.empty() && (pos=buf.find_last_not_of(whitespace))!=buf.size()-1)
                buf.erase(pos+1);
            if(buf.empty())
                continue;

            std::stringstream ss(buf);
            double x, y, f;
            ss >> x >> y >> f;
            table_data.x.push_back(x);
            table_data.y.push_back(y);
            table_data.f.push_back(f);
        }
        LI::utilities::Interpolator2D<double> interp(table_data);
        AddDifferentialCrossSection(pid, interp);
    } else {
        throw std::runtime_error("Failed open cross section file!");
    }
}

void DipoleFromTable::AddTotalCrossSectionFile(std::string filename, LI::dataclasses::Particle::ParticleType target) {
    std::string delimeter = "_";
    std::string end_delimeter = ".";
    std::string::size_type pos = filename.rfind("/") + 1;
    std::string::size_type next_pos = pos;
    std::string::size_type sub_len = 0;
    unsigned int Z = 0;
    unsigned int A = 0;
    bool bad = false;

    std::function<std::string()> next_substr = [&] () -> std::string {
        if(pos >= filename.size() or pos == std::string::npos) {
            bad = true;
            return std::string();
        }
        next_pos = filename.find(delimeter, pos);
        if(next_pos == std::string::npos) {
            next_pos = filename.find(end_delimeter, pos);
            if(next_pos == std::string::npos) {
                bad = true;
                return std::string();
            }
        }
        sub_len = std::max((int)(next_pos - pos), 0);
        next_pos = pos + sub_len;
        std::string sub = filename.substr(pos, sub_len);
        pos = next_pos + 1;
        return sub;
    };

    while(pos < filename.size() and pos != std::string::npos) {
        std::string sub = next_substr();
        if(bad)
            break;
        if(sub == "Z") {
            sub = next_substr();
            if(bad)
                break;
            Z = std::stoi(sub);
        } else if(sub == "A") {
            sub = next_substr();
            if(bad)
                break;
            A = std::stoi(sub);
        } else if(sub == "mHNL") {
            sub = next_substr();
            if(bad)
                break;
            double file_hnl_mass = std::stod(sub);
            if(std::abs(file_hnl_mass - hnl_mass) / std::max(std::abs(file_hnl_mass), std::abs(hnl_mass)) > 1e-6) {
                std::cout << std::setprecision(24);
                std::cout << "File HNL mass: "<< file_hnl_mass << std::endl;
                std::cout << "Specified HNL mass: "<< hnl_mass << std::endl;
                throw std::runtime_error("File HNL mass does not match specified HNL mass!");
            }
        }
    }
    //std::string pid_str = std::format("10{0:0>1d}{1:0>3d}{2:0>3d}{3:0>1d}", 0, Z, A, 0);
    unsigned int buffer_size = 1024;
    char buffer[buffer_size];
    unsigned int str_size = std::snprintf(buffer, buffer_size, "10%01d%03d%03d%01d", 0, Z, A, 0);
    if(str_size > 10) {
        throw std::runtime_error("Cannot create particle ID string!");
    }
    std::string pid_str(buffer);
    int32_t pid_int = std::stoul(pid_str);
    LI::dataclasses::Particle::ParticleType pid = (LI::dataclasses::Particle::ParticleType)pid_int;

    if(pid != target) {
        throw std::runtime_error("File target nucleus does not match supplied target type!");
    }

    if(fexists(filename)) {
        std::ifstream in(filename.c_str());
        std::string buf;

        LI::utilities::TableData1D<double> table_data;
        while(std::getline(in, buf)) {
            if((pos = buf.find('#')) != std::string::npos)
                buf.erase(pos);
            const char* whitespace=" \n\r\t\v";
            if((pos=buf.find_first_not_of(whitespace))!=0)
                buf.erase(0,pos);
            if(!buf.empty() && (pos=buf.find_last_not_of(whitespace))!=buf.size()-1)
                buf.erase(pos+1);
            if(buf.empty())
                continue;

            std::stringstream ss(buf);
            double x, f;
            ss >> x >> f;
            table_data.x.push_back(x);
            table_data.f.push_back(f);
        }
        LI::utilities::Interpolator1D<double> interp(table_data);
        AddTotalCrossSection(pid, interp);
    } else {
        throw std::runtime_error("Failed open cross section file!");
    }
}

void DipoleFromTable::AddDifferentialCrossSection(LI::dataclasses::Particle::ParticleType target, LI::utilities::Interpolator2D<double> interp) {
    differential.insert(std::make_pair(target, interp));
}

void DipoleFromTable::AddTotalCrossSection(LI::dataclasses::Particle::ParticleType target, LI::utilities::Interpolator1D<double> interp) {
    total.insert(std::make_pair(target, interp));
}

} // namespace interactions
} // namespace LI
