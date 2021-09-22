#include "phys-services/CrossSection.h"

#include <array>

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"
#include "stga3/STGA3.h"
#include "stga3/Typedefs.h"
#include "stga3/Utilities.h"

namespace LeptonInjector {

namespace {
    /*
    double particleMass(LeptonInjector::Particle::ParticleType type) {
        LeptonInjector::Particle p(type);
        if(!p.HasMass()) {
            return 0;
        }
        return p.GetMass();
    }
    */

    ///Check whether a given point in phase space is physically realizable.
    ///Based on equations 6-8 of http://dx.doi.org/10.1103/PhysRevD.66.113007
    ///S. Kretzer and M. H. Reno
    ///"Tau neutrino deep inelastic charged current interactions"
    ///Phys. Rev. D 66, 113007
    ///\param x Bjorken x of the interaction
    ///\param y Bjorken y of the interaction
    ///\param E Incoming neutrino in energy in the lab frame ($E_\nu$)
    ///\param M Mass of the target nucleon ($M_N$)
    ///\param m Mass of the secondary lepton ($m_\tau$)
    bool kinematicallyAllowed(double x, double y, double E, double M, double m) {
        if(x > 1) //Eq. 6 right inequality
            return false;
        if(x < ((m * m) / (2 * M * (E - m)))) //Eq. 6 left inequality
            return false;
        //denominator of a and b
        double d = 2 * (1 + (M * x) / (2 * E));
        //the numerator of a (or a*d)
        double ad = 1 - m * m * ((1 / (2 * M * E * x)) + (1 / (2 * E * E)));
        double term = 1 - ((m * m) / (2 * M * E * x));
        //the numerator of b (or b*d)
        double bd = sqrt(term * term - ((m * m) / (E * E)));
        return (ad - bd) <= d * y and d * y <= (ad + bd); //Eq. 7
    }

    inline double dot(std::array<double, 4> p0, std::array<double, 4> p1) {
        return p0[0] * p1[0] - (p0[1] * p1[1] + p0[2] * p1[2] + p0[3] * p1[3]);
    }
}

void CrossSectionCollection::InitializeTargetTypes() {
    target_types.clear();
    cross_sections_by_target.clear();
    for(unsigned int i=0; i<cross_sections.size(); ++i) {
        std::vector<Particle::ParticleType> xs_targets = cross_sections[i]->GetPossibleTargets();
        target_types.reserve(target_types.size() + std::distance(xs_targets.begin(), xs_targets.end()));
        target_types.insert(target_types.end(), xs_targets.begin(), xs_targets.end());

        for(unsigned int j=0; j<xs_targets.size(); ++j) {
            Particle::ParticleType target = xs_targets[j];
            std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>>::const_iterator it = cross_sections_by_target.find(target);
            if(it == cross_sections_by_target.end()) {
                cross_sections_by_target.insert(it, std::make_pair(target, std::vector<std::shared_ptr<CrossSection>>{cross_sections[i]}));
            } else {
                cross_sections_by_target[target].push_back(cross_sections[i]);
            }
        }
    }
}

CrossSectionCollection::CrossSectionCollection(Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections) : primary_type(primary_type), cross_sections(cross_sections) {
    InitializeTargetTypes();
}

std::vector<std::shared_ptr<CrossSection>> const & CrossSectionCollection::GetCrossSectionsForTarget(Particle::ParticleType p) const {
    if(cross_sections_by_target.find(p) != cross_sections_by_target.end()) {
        return std::vector<std::shared_ptr<CrossSection>>(cross_sections_by_target.at(p));
    } else {
        return std::vector<std::shared_ptr<CrossSection>>();
    }
}

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minimum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types) : primary_types_(primary_types), target_types_(target_types), minimum_Q2_(minimum_Q2), target_mass_(target_mass), interaction_type_(interaction) {
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
}

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types) : primary_types_(primary_types), target_types_(target_types) {
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
}

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minimum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types) : primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()), minimum_Q2_(minimum_Q2), target_mass_(target_mass), interaction_type_(interaction) {
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
}

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types) : primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()) {
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
}

void DISFromSpline::LoadFromFile(std::string dd_crossSectionFile, std::string total_crossSectionFile) {

    differential_cross_section_ = photospline::splinetable<>(dd_crossSectionFile.c_str());

    if(differential_cross_section_.get_ndim()!=3 && differential_cross_section_.get_ndim()!=2)
        throw("cross section spline has " + std::to_string(differential_cross_section_.get_ndim())
                + " dimensions, should have either 3 (log10(E), log10(x), log10(y)) or 2 (log10(E), log10(y))");

    total_cross_section_ = photospline::splinetable<>(total_crossSectionFile.c_str());

    if(total_cross_section_.get_ndim() != 1)
        throw("Total cross section spline has " + std::to_string(total_cross_section_.get_ndim())
                + " dimensions, should have 1, log10(E)");
}

void DISFromSpline::ReadParamsFromSplineTable() {
    // returns true if successfully read target mass
    bool mass_good = differential_cross_section_.read_key("TARGETMASS", target_mass_);
    // returns true if successfully read interaction type
    bool int_good = differential_cross_section_.read_key("INTERACTION", interaction_type_);
    // returns true if successfully read minimum Q2
    bool q2_good = differential_cross_section_.read_key("Q2MIN", minimum_Q2_);

    if(!int_good) {
        // assume DIS to preserve compatability with previous versions
        interaction_type_ = 1;
    }

    if(!q2_good) {
        // assume 1 GeV^2
        minimum_Q2_ = 1;
    }

    if(!mass_good) {
        if(int_good) {
            if(interaction_type_ == 1 or interaction_type_ == 2) {
                target_mass_ = (LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::PPlus)+
                        LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::Neutron))/2;
            } else if(interaction_type_ == 3) {
                target_mass_ = LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::EMinus);
            } else {
                throw("Logic error. Interaction type is not 1, 2, or 3!");
            }

        } else {
            if(differential_cross_section_.get_ndim() == 3) {
                target_mass_ = (LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::PPlus)+
                        LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::Neutron))/2;
            } else if(differential_cross_section_.get_ndim() == 2) {
                target_mass_ = LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::EMinus);
            } else {
                throw("Logic error. Spline dimensionality is not 2, or 3!");
            }
        }
    }
}

void DISFromSpline::InitializeSignatures() {
    signatures_.clear();
    for(auto primary_type : primary_types_) {
        InteractionSignature signature;
        signature.primary_type = primary_type;

        if(not isNeutrino(primary_type)) {
            throw std::runtime_error("This DIS implementation only supports neutrinos as primaries!");
        }

        Particle::ParticleType charged_lepton_product = Particle::ParticleType::unknown;
        Particle::ParticleType neutral_lepton_product = primary_type;

        if(primary_type == Particle::ParticleType::NuE) {
            charged_lepton_product = Particle::ParticleType::EMinus;
        } else if(primary_type == Particle::ParticleType::NuEBar) {
            charged_lepton_product = Particle::ParticleType::EPlus;
        } else if(primary_type == Particle::ParticleType::NuMu) {
            charged_lepton_product = Particle::ParticleType::MuMinus;
        } else if(primary_type == Particle::ParticleType::NuMuBar) {
            charged_lepton_product = Particle::ParticleType::MuPlus;
        } else if(primary_type == Particle::ParticleType::NuTau) {
            charged_lepton_product = Particle::ParticleType::TauMinus;
        } else if(primary_type == Particle::ParticleType::NuTauBar) {
            charged_lepton_product = Particle::ParticleType::TauPlus;
        } else {
            throw std::runtime_error("InitializeSignatures: Unkown parent neutrino type!");
        }

        if(interaction_type_ == 1) {
            signature.secondary_types.push_back(charged_lepton_product);
        } else if(interaction_type_ == 2) {
            signature.secondary_types.push_back(neutral_lepton_product);
        } else if(interaction_type_ == 3) {
            signature.secondary_types.push_back(Particle::ParticleType::Hadrons);
        } else {
            throw std::runtime_error("InitializeSignatures: Unkown interaction type!");
        }
        for(auto target_type : target_types_) {
            signature.target_type = target_type;

            signatures_.push_back(signature);

            std::pair<Particle::ParticleType, Particle::ParticleType> key(primary_type, target_type);
            signatures_by_parent_types_[key].push_back(signature);
        }
    }
}

double DISFromSpline::TotalCrossSection(InteractionRecord const & interaction) const {
    LeptonInjector::Particle::ParticleType primary_type = interaction.signature.primary_type;
    double primary_energy;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
    } else {
        throw std::runtime_error("Lorentz boost not implemented!");
    }
    return TotalCrossSection(primary_type, primary_energy);
}

double DISFromSpline::TotalCrossSection(LeptonInjector::Particle::ParticleType primary_type, double primary_energy) const {
    if(not primary_types_.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
    double target_mass = target_mass_;
    double log_energy = log10(primary_energy);

    if(log_energy < total_cross_section_.lower_extent(0)
            or log_energy > total_cross_section_.upper_extent(0)) {
        throw("Interaction energy ("+ std::to_string(primary_energy) +
                ") out of cross section table range: ["
                + std::to_string(pow(10.,total_cross_section_.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10.,total_cross_section_.upper_extent(0))) + " GeV]");
    }

    int center;
    total_cross_section_.searchcenters(&log_energy, &center);
    double log_xs = total_cross_section_.ndsplineeval(&log_energy, &center, 0);

    return std::pow(10.0, log_xs);
}

double DISFromSpline::DifferentialCrossSection(InteractionRecord const & interaction) const {
    LeptonInjector::Particle::ParticleType primary_type = interaction.signature.primary_type;
    double primary_energy;
    std::array<double, 4> p1;
    std::array<double, 4> p2;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
        p1 = interaction.primary_momentum;
        p2 = interaction.target_momentum;
    } else {
        throw std::runtime_error("Lorentz boost not implemented!");
    }
    assert(interaction.signature.secondary_types.size() == 2);
    assert(isLepton(interaction.signature.secondary_types[0]) or isLepton(interaction.signature.secondary_types[1]));
    unsigned int lepton_index = (isLepton(interaction.signature.secondary_types[0])) ? 0 : 1;
    std::array<double, 4> p3 = interaction.secondary_momenta[lepton_index];
    std::array<double, 4> q = {p1[0] - p3[0], p1[1] - p3[1], p1[2] - p3[2], p1[3] - p3[3]};
    double Q2 = -dot(q, q);
    double y = dot(p2, q) / dot(p2, p1);
    double x = Q2 / (2.0 * dot(p2, q));

    double lepton_mass = particleMass(interaction.signature.secondary_types[lepton_index]);

    return DifferentialCrossSection(primary_energy, x, y, lepton_mass);
}

double DISFromSpline::DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass) const {
    double log_energy = log10(energy);
    // check preconditions
    if(log_energy < differential_cross_section_.lower_extent(0)
            || log_energy>differential_cross_section_.upper_extent(0))
        throw("Interaction energy ("+ std::to_string(energy) +
                ") out of cross section table range: ["
                + std::to_string(pow(10., differential_cross_section_.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10., differential_cross_section_.upper_extent(0))) + " GeV]");
    if(x <= 0 || x >= 1)
        throw("Interaction x out of range: " + std::to_string(x));
    if(y <= 0 || y >= 1)
        throw("Interaction y out of range: " + std::to_string(y));

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    double Q2 = 2.0 * energy * target_mass_ * x * y;
    if(Q2 < minimum_Q2_) // cross section not calculated, assumed to be zero
        return 0;

    // cross section should be zero, but this check is missing from the original
    // CSMS calculation, so we must add it here
    if(!kinematicallyAllowed(x, y, energy, target_mass_, secondary_lepton_mass))
        return 0;

    std::array<double,3> coordinates{{log_energy, log10(x), log10(y)}};
    std::array<int,3> centers;
    if(!differential_cross_section_.searchcenters(coordinates.data(), centers.data()))
        return 0;
    double result = pow(10., differential_cross_section_.ndsplineeval(coordinates.data(), centers.data(), 0));
    assert(result >= 0);
    return result;
}


void DISFromSpline::SampleFinalState(LeptonInjector::InteractionRecord& interaction, std::shared_ptr<LeptonInjector::LI_random> random) const {
    // Uses Metropolis-Hastings Algorithm!
    // useful for cases where we don't know the supremum of our distribution, and the distribution is multi-dimensional
    if (differential_cross_section_.get_ndim() != 3) {
        throw("I expected 3 dimensions in the cross section spline, but got " + std::to_string(differential_cross_section_.get_ndim()) +". Maybe your fits file doesn't have the right 'INTERACTION' key?");
    }

    stga3::FourVector<double> p1{interaction.primary_momentum[0], interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]};
    stga3::FourVector<double> p2{interaction.target_momentum[0], interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]};

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * tinteraction.secondary_momentarget_mass_ + 2 * target_mass_ * primary_energy;
    double s = (p1 + p2) | (p1 + p2);

    double primary_energy;
    stga3::FourVector<double> p1_lab;
    stga3::FourVector<double> p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p1_lab = p1;
        p2_lab = p2;
        primary_energy = p1_lab.e0();
    } else {
        stga3::Beta<double> beta_start_to_lab = stga3::beta_to_rest_frame_of(p2);
        stga3::Boost<double> boost_start_to_lab = stga3::boost_from_beta(beta_start_to_lab);
        p1_lab = stga3::apply_boost(boost_start_to_lab, p1);
        p2_lab = stga3::apply_boost(boost_start_to_lab, p2);
        primary_energy = p1_lab.e0();
    }

    unsigned int lepton_index = (isLepton(interaction.signature.secondary_types[0])) ? 0 : 1;
    unsigned int other_index = 1 - lepton_index;
    double m = particleMass(interaction.signature.secondary_types[lepton_index]);

    double m1 = p1_lab | p1_lab;
    double m3 = m;
    double E1_lab = p1_lab.e0();
    double E2_lab = p2_lab.e0();

    // The out-going particle always gets at least enough energy for its rest mass
    double yMax = 1 - m / primary_energy;
    double logYMax = log10(yMax);

    // The minimum allowed value of y occurs when x = 1 and Q is minimized
    double yMin = minimum_Q2_ / s;
    double logYMin = log10(yMin);
    // The minimum allowed value of x occurs when y = yMax and Q is minimized
    // double xMin = minimum_Q2_ / ((s - target_mass_ * target_mass_) * yMax);
    double xMin = minimum_Q2_ / (2 * E1_lab * E2_lab * yMax);
    double logXMin = log10(xMin);

    bool accept;

    // kin_vars and its twin are 3-vectors containing [nu-energy, Bjorken X, Bjorken Y]
    std::array<double,3> kin_vars, test_kin_vars;

    // centers of the cross section spline tales.
    std::array<int,3> spline_table_center, test_spline_table_center;

    // values of cross_section from the splines.  By * Bx * Spline(E,x,y)
    double cross_section, test_cross_section;

    // No matter what, we're evaluating at this specific energy.
    kin_vars[0] = test_kin_vars[0] = log10(primary_energy);

    // check preconditions
    if(kin_vars[0] < differential_cross_section_.lower_extent(0)
            || kin_vars[0] > differential_cross_section_.upper_extent(0))
        throw("Interaction energy out of cross section table range: ["
                + std::to_string(pow(10.,differential_cross_section_.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10.,differential_cross_section_.upper_extent(0))) + " GeV]");

    // sample an intial point
    do {
        // rejection sample a point which is kinematically allowed by calculation limits
        double trialQ;
        do {
            kin_vars[1] = random->Uniform(logXMin,0);
            kin_vars[2] = random->Uniform(logYMin,logYMax);
            trialQ = (2 * E1_lab * E2_lab) * pow(10., kin_vars[1] + kin_vars[2]);
        } while(trialQ<minimum_Q2_ || !kinematicallyAllowed(pow(10., kin_vars[1]), pow(10., kin_vars[2]), primary_energy, target_mass_, m));

        accept = true;
        //sanity check: demand that the sampled point be within the table extents
        if(kin_vars[1] < differential_cross_section_.lower_extent(1)
                || kin_vars[1] > differential_cross_section_.upper_extent(1)) {
            accept = false;
        }
        if(kin_vars[2] < differential_cross_section_.lower_extent(2)
                || kin_vars[2] > differential_cross_section_.upper_extent(2)) {
            accept = false;
        }

        if(accept) {
            // finds the centers in the cross section spline table, returns true if it's successful
            // also sets the centers
            accept = differential_cross_section_.searchcenters(kin_vars.data(),spline_table_center.data());
        }
    } while(!accept);

    //TODO: better proposal distribution?
    double measure = pow(10., kin_vars[1] + kin_vars[2]); // Bx * By

    // Bx * By * xs(E, x, y)
    // evalutates the differential spline at that point
    cross_section = measure*pow(10., differential_cross_section_.ndsplineeval(kin_vars.data(), spline_table_center.data(), 0));

    // this is the magic part. Metropolis Hastings Algorithm.
    // MCMC method!
    const size_t burnin = 40; // converges to the correct distribution over multiple samplings.
    // big number means more accurate, but slower
    for(size_t j = 0; j <= burnin; j++) {
        // repeat the sampling from above to get a new valid point
        double trialQ;
        do {
            test_kin_vars[1] = random->Uniform(logXMin, 0);
            test_kin_vars[2] = random->Uniform(logYMin, logYMax);
            trialQ = (2 * E1_lab * E2_lab) * pow(10., test_kin_vars[1] + test_kin_vars[2]);
        } while(trialQ < minimum_Q2_ || !kinematicallyAllowed(pow(10., test_kin_vars[1]), pow(10., test_kin_vars[2]), primary_energy, target_mass_, m));

        accept = true;
        if(test_kin_vars[1] < differential_cross_section_.lower_extent(1)
                || test_kin_vars[1] > differential_cross_section_.upper_extent(1))
            accept = false;
        if(test_kin_vars[2] < differential_cross_section_.lower_extent(2)
                || test_kin_vars[2] > differential_cross_section_.upper_extent(2))
            accept = false;
        if(!accept)
            continue;

        accept = differential_cross_section_.searchcenters(test_kin_vars.data(), test_spline_table_center.data());
        if(!accept)
            continue;

        double measure = pow(10., test_kin_vars[1] + test_kin_vars[2]);
        double eval = differential_cross_section_.ndsplineeval(test_kin_vars.data(), test_spline_table_center.data(), 0);
        if(std::isnan(eval))
            continue;
        test_cross_section = measure * pow(10., eval);

        double odds = (test_cross_section / cross_section);
        accept = ((odds > 1.) || random->Uniform(0, 1) < odds);

        if(accept) {
            kin_vars = test_kin_vars;
            cross_section = test_cross_section;
        }
    }
    double final_x = pow(10., kin_vars[1]);
    double final_y = pow(10., kin_vars[2]);

    interaction.interaction_parameters.resize(3);
    interaction.interaction_parameters[0] = E1_lab;
    interaction.interaction_parameters[1] = final_x;
    interaction.interaction_parameters[2] = final_y;

    double Q2 = 2 * E1_lab * E2_lab * pow(10.0, kin_vars[1] + kin_vars[2]);
    double p1x_lab = std::sqrt(p1_lab.e1() * p1_lab.e1() + p1_lab.e2() * p1_lab.e2() + p1_lab.e3() * p1_lab.e3());
    double pqx_lab = (m1*m1 + m3*m3 + 2 * p1x_lab * p1x_lab + Q2 + 2 * E1_lab * E1_lab * (final_y - 1)) / (2.0 * p1x_lab);
    double momq_lab = std::sqrt(m1*m1 + p1x_lab*p1x_lab + Q2 + E1_lab * E1_lab * (final_y * final_y - 1));
    double pqy_lab = std::sqrt(momq_lab*momq_lab - pqx_lab *pqx_lab);
    double Eq_lab = E1_lab * final_y;

    stga3::ThreeVector<double> x_dir{1.0, 0.0, 0.0};
    stga3::Rotation<double> x_to_p1_lab_rot = stga3::rotation_between(x_dir, p1_lab);

    double phi = random->Uniform(0, 2.0 * M_PI);
    stga3::Rotation<double> rand_rot = stga3::rotation_about(p1_lab, phi);

    stga3::FourVector<double> pq_lab{Eq_lab, pqx_lab, pqy_lab, 0};
    pq_lab = stga3::apply_rotation(x_to_p1_lab_rot, pq_lab);
    pq_lab = stga3::apply_rotation(rand_rot, pq_lab);

    stga3::FourVector<double> p3_lab = p1_lab - pq_lab;
    stga3::FourVector<double> p4_lab = p2_lab + pq_lab;

    stga3::ThreeVector<double> p3_lab_vec{p3_lab.e1(), p3_lab.e2(), p3_lab.e3()};
    stga3::ThreeVector<double> p4_lab_vec{p4_lab.e1(), p4_lab.e2(), p4_lab.e3()};

    stga3::FourVector<double> p3;
    stga3::FourVector<double> p4;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p3 = p3_lab;
        p4 = p4_lab;
    } else {
        stga3::Beta<double> beta_start_to_lab = stga3::beta_to_rest_frame_of(p2);
        stga3::Boost<double> boost_lab_to_start = stga3::beta_to_boost(-beta_start_to_lab);
        p3 = stga3::apply_boost(boost_lab_to_start, p3_lab);
        p4 = stga3::apply_boost(boost_lab_to_start, p4_lab);
    }

    stga3::FourVector<double> pq_13 = p1 - p3;
    stga3::FourVector<double> pq_24 = p4 - p2;

    // Check that computed q2 in the start frame matches up with the specified Q2
    assert(std::abs(double(pq_13 | pq_13) + Q2) < std::abs(Q2 * 1e-3));
    assert(std::abs(double(pq_24 | pq_24) + Q2) < std::abs(Q2 * 1e-3));

    interaction.secondary_momenta.resize(2);
    interaction.secondary_momenta[lepton_index][0] = p3.e0(); // p3_energy
    interaction.secondary_momenta[lepton_index][1] = p3.e1(); // p3_x
    interaction.secondary_momenta[lepton_index][2] = p3.e2(); // p3_y
    interaction.secondary_momenta[lepton_index][3] = p3.e3(); // p3_z

    interaction.secondary_momenta[other_index][0] = p4.e0(); // p4_energy
    interaction.secondary_momenta[other_index][1] = p4.e1(); // p4_x
    interaction.secondary_momenta[other_index][2] = p4.e2(); // p4_y
    interaction.secondary_momenta[other_index][3] = p4.e3(); // p4_z
}

std::vector<Particle::ParticleType> DISFromSpline::GetPossiblePrimaries() const {
    return std::vector<Particle::ParticleType>(primary_types_.begin(), primary_types_.end());
}

std::vector<Particle::ParticleType> DISFromSpline::GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const {
    return std::vector<Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<InteractionSignature> DISFromSpline::GetPossibleSignatures() const {
    return std::vector<InteractionSignature>(signatures_.begin(), signatures_.end());
}

std::vector<Particle::ParticleType> DISFromSpline::GetPossibleTargets() const {
    return std::vector<Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<InteractionSignature> DISFromSpline::GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const {
    std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType> key(primary_type, target_type);
    if(signatures_by_parent_types_.find(key) != signatures_by_parent_types_.end()) {
        return signatures_by_parent_types_.at(key);
    } else {
        return std::vector<InteractionSignature>();
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace LeptonInjector

