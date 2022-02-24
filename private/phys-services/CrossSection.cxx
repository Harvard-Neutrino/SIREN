#include "phys-services/CrossSection.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"

#include <rk/rk.hh>
#include <rk/geom3.hh>

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

    template <class InIt>
	typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
		typedef typename std::iterator_traits<InIt>::value_type real;
		real sum = real(0);
		real running_error = real(0);
		real temp;
		real difference;

		for (; begin != end; ++begin) {
			difference = *begin;
			difference -= running_error;
			temp = sum;
			temp += difference;
			running_error = temp;
			running_error -= sum;
			running_error -= difference;
			sum = std::move(temp);
		}
		return sum;
	};

    template<typename T>
    T accumulate(std::initializer_list<T> list) {
        return accumulate(list.begin(), list.end());
    }

	template <typename T>
	class accumulator {
    private:
		T sum;
		T running_error;
		T temp;
		T difference;
    public:
        accumulator() : sum(0), running_error(0) {};
        accumulator(std::initializer_list<T> const & list) : sum(0), running_error(0) {
            for(auto const & val : list) {
                accumulate(val);
            }
        }
        void accumulate(T const & val) {
			difference = val;
			difference -= running_error;
			temp = sum;
			temp += difference;
			running_error = temp;
			running_error -= sum;
			running_error -= difference;
			sum = std::move(temp);
        }
        T result() {
            return sum;
        }
        operator T() {
            return sum;
        }
    };
}

bool InteractionSignature::operator==(InteractionSignature const & other) const {
    if(primary_type != other.primary_type or target_type != other.target_type) {
        return false;
    } else {
        std::map<LeptonInjector::Particle::ParticleType, int> m0;
        for(auto p : secondary_types) {
            auto it = m0.find(p);
            if(it == m0.end()) {
                m0.insert({p, 1});
            } else {
                it->second += 1;
            }
        }
        std::map<LeptonInjector::Particle::ParticleType, int> m1;
        for(auto p : other.secondary_types) {
            auto it = m1.find(p);
            if(it == m1.end()) {
                m1.insert({p, 1});
            } else {
                it->second += 1;
            }
        }
        return m0 == m1;
    }
}

void CrossSectionCollection::InitializeTargetTypes() {
    target_types.clear();
    cross_sections_by_target.clear();
    for(unsigned int i=0; i<cross_sections.size(); ++i) {
        // Gather target types
        std::vector<Particle::ParticleType> xs_targets = cross_sections[i]->GetPossibleTargets();
        //target_types.reserve(target_types.size() + std::distance(xs_targets.begin(), xs_targets.end()));
        for(auto xs : xs_targets)
            target_types.insert(xs);

        // Track cross sections by their target type
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

    // Remove duplicate target types
    // std::set<Particle::ParticleType> target_set(target_types.begin(), target_types.end());
    // target_types.resize(target_set.size());
    // std::copy(target_set.begin(), target_set.end(), target_types.begin());
}

const std::vector<std::shared_ptr<CrossSection>> CrossSectionCollection::empty = {};

CrossSectionCollection::CrossSectionCollection() {}

CrossSectionCollection::CrossSectionCollection(Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections) : primary_type(primary_type), cross_sections(cross_sections) {
    InitializeTargetTypes();
}

std::vector<std::shared_ptr<CrossSection>> const & CrossSectionCollection::GetCrossSectionsForTarget(Particle::ParticleType p) const {
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>>::const_iterator it = cross_sections_by_target.find(p);
    if(it != cross_sections_by_target.end()) {
        return it->second;
    } else {
        return empty;
    }
}

bool CrossSectionCollection::MatchesPrimary(InteractionRecord const & record) const {
    return primary_type == record.signature.primary_type;
}

DISFromSpline::DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minimum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types) : primary_types_(primary_types), target_types_(target_types), minimum_Q2_(minimum_Q2), target_mass_(target_mass), interaction_type_(interaction) {
    LoadFromMemory(differential_data, total_data);
    InitializeSignatures();
}

DISFromSpline::DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minimum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types) : primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()), minimum_Q2_(minimum_Q2), target_mass_(target_mass), interaction_type_(interaction) {
    LoadFromMemory(differential_data, total_data);
    InitializeSignatures();
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
        throw std::runtime_error("cross section spline has " + std::to_string(differential_cross_section_.get_ndim())
                + " dimensions, should have either 3 (log10(E), log10(x), log10(y)) or 2 (log10(E), log10(y))");

    total_cross_section_ = photospline::splinetable<>(total_crossSectionFile.c_str());

    if(total_cross_section_.get_ndim() != 1)
        throw std::runtime_error("Total cross section spline has " + std::to_string(total_cross_section_.get_ndim())
                + " dimensions, should have 1, log10(E)");
}

void DISFromSpline::LoadFromMemory(std::vector<char> & differential_data, std::vector<char> & total_data) {
    differential_cross_section_.read_fits_mem(differential_data.data(), differential_data.size());
    total_cross_section_.read_fits_mem(total_data.data(), total_data.size());
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
                throw std::runtime_error("Logic error. Interaction type is not 1, 2, or 3!");
            }

        } else {
            if(differential_cross_section_.get_ndim() == 3) {
                target_mass_ = (LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::PPlus)+
                        LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::Neutron))/2;
            } else if(differential_cross_section_.get_ndim() == 2) {
                target_mass_ = LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::EMinus);
            } else {
                throw std::runtime_error("Logic error. Spline dimensionality is not 2, or 3!");
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
        throw std::runtime_error("Interaction energy ("+ std::to_string(primary_energy) +
                ") out of cross section table range: ["
                + std::to_string(pow(10.,total_cross_section_.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10.,total_cross_section_.upper_extent(0))) + " GeV]");
    }

    int center;
    total_cross_section_.searchcenters(&log_energy, &center);
    double log_xs = total_cross_section_.ndsplineeval(&log_energy, &center, 0);

    return std::pow(10.0, log_xs);
}

// No implementation for DIS yet, just use non-target function
double DISFromSpline::TotalCrossSection(LeptonInjector::Particle::ParticleType primary_type, double primary_energy, Particle::ParticleType target_type) const {
		return DISFromSpline::TotalCrossSection(primary_type,primary_energy);
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
    double y = 1.0 - dot(p2, p3) / dot(p2, p1);
    double x = Q2 / (2.0 * dot(p2, q));

    double lepton_mass = particleMass(interaction.signature.secondary_types[lepton_index]);

    return DifferentialCrossSection(primary_energy, x, y, lepton_mass);
}

double DISFromSpline::DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass) const {
    double log_energy = log10(energy);
    // check preconditions
    if(log_energy < differential_cross_section_.lower_extent(0)
            || log_energy>differential_cross_section_.upper_extent(0))
        return 0.0;
        throw std::runtime_error("Interaction energy ("+ std::to_string(energy) +
                ") out of cross section table range: ["
                + std::to_string(pow(10., differential_cross_section_.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10., differential_cross_section_.upper_extent(0))) + " GeV]");
    if(x <= 0 || x >= 1) {
        return 0.0;
        throw std::runtime_error("Interaction x out of range: " + std::to_string(x));
    }
    if(y <= 0 || y >= 1) {
        return 0.0;
        throw std::runtime_error("Interaction y out of range: " + std::to_string(y));
    }

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

double DISFromSpline::InteractionThreshold(InteractionRecord const & interaction) const {
    // Consider implementing DIS thershold at some point
    return 0;
}

void DISFromSpline::SampleFinalState(LeptonInjector::InteractionRecord& interaction, std::shared_ptr<LeptonInjector::LI_random> random) const {
    // Uses Metropolis-Hastings Algorithm!
    // useful for cases where we don't know the supremum of our distribution, and the distribution is multi-dimensional
    if (differential_cross_section_.get_ndim() != 3) {
        throw std::runtime_error("I expected 3 dimensions in the cross section spline, but got " + std::to_string(differential_cross_section_.get_ndim()) +". Maybe your fits file doesn't have the right 'INTERACTION' key?");
    }

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * tinteraction.secondary_momentarget_mass_ + 2 * target_mass_ * primary_energy;
    double s = rk::invMass(p1, p2);

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p1_lab = p1;
        p2_lab = p2;
        primary_energy = p1_lab.e();
    } else {
        // Rest frame of p2 will be our "lab" frame
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p2_lab = boost_start_to_lab * p2;
        primary_energy = p1_lab.e();
    }

    unsigned int lepton_index = (isLepton(interaction.signature.secondary_types[0])) ? 0 : 1;
    unsigned int other_index = 1 - lepton_index;
    double m = particleMass(interaction.signature.secondary_types[lepton_index]);

    double m1 = interaction.primary_mass;
    double m3 = m;
    double E1_lab = p1_lab.e();
    double E2_lab = p2_lab.e();

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
        throw std::runtime_error("Interaction energy out of cross section table range: ["
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
        accept = (cross_section == 0 || (odds > 1.) || random->Uniform(0, 1) < odds);

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
    double p1x_lab = std::sqrt(p1_lab.px() * p1_lab.px() + p1_lab.py() * p1_lab.py() + p1_lab.pz() * p1_lab.pz());
    double pqx_lab = (m1*m1 + m3*m3 + 2 * p1x_lab * p1x_lab + Q2 + 2 * E1_lab * E1_lab * (final_y - 1)) / (2.0 * p1x_lab);
    double momq_lab = std::sqrt(m1*m1 + p1x_lab*p1x_lab + Q2 + E1_lab * E1_lab * (final_y * final_y - 1));
    double pqy_lab = std::sqrt(momq_lab*momq_lab - pqx_lab *pqx_lab);
    double Eq_lab = E1_lab * final_y;

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 p1_mom = p1_lab.momentum();
    geom3::UnitVector3 p1_lab_dir = p1_mom.direction();
    geom3::Rotation3 x_to_p1_lab_rot = geom3::rotationBetween(x_dir, p1_lab_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(p1_lab_dir, phi);

    rk::P4 pq_lab(Eq_lab, geom3::Vector3(pqx_lab, pqy_lab, 0));
    pq_lab.rotate(x_to_p1_lab_rot);
    pq_lab.rotate(rand_rot);

    rk::P4 p3_lab((p1_lab - pq_lab).momentum(), m3);
    rk::P4 p4_lab = p2_lab + pq_lab;

    rk::P4 p3;
    rk::P4 p4;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p3 = p3_lab;
        p4 = p4_lab;
    } else {
        rk::Boost boost_lab_to_start = p2.restBoost();
        p3 = boost_lab_to_start * p3_lab;
        p4 = boost_lab_to_start * p4_lab;
    }

    rk::P4 pq_13 = p1 - p3;
    rk::P4 pq_24 = p4 - p2;

    // Check that computed q2 in the start frame matches up with the specified Q2
    assert(std::abs(double(pq_13.dot(pq_13)) + Q2) < std::abs(Q2 * 1e-3));
    assert(std::abs(double(pq_24.dot(pq_24)) + Q2) < std::abs(Q2 * 1e-3));

    interaction.secondary_momenta.resize(2);
    interaction.secondary_masses.resize(2);
    interaction.secondary_helicity.resize(2);

    interaction.secondary_momenta[lepton_index][0] = p3.e(); // p3_energy
    interaction.secondary_momenta[lepton_index][1] = p3.px(); // p3_x
    interaction.secondary_momenta[lepton_index][2] = p3.py(); // p3_y
    interaction.secondary_momenta[lepton_index][3] = p3.pz(); // p3_z
    interaction.secondary_masses[lepton_index] = p3.m();

    interaction.secondary_helicity[lepton_index] = interaction.primary_helicity;

    interaction.secondary_momenta[other_index][0] = p4.e(); // p4_energy
    interaction.secondary_momenta[other_index][1] = p4.px(); // p4_x
    interaction.secondary_momenta[other_index][2] = p4.py(); // p4_y
    interaction.secondary_momenta[other_index][3] = p4.pz(); // p4_z
    interaction.secondary_masses[other_index] = p4.m();

    interaction.secondary_helicity[other_index] = interaction.target_helicity;
}

double DISFromSpline::FinalStateProbability(LeptonInjector::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
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

std::vector<std::string> DISFromSpline::DensityVariables() const {
    return std::vector<std::string>{"Bjorken x", "Bjorken y"};
}


double DipoleFromTable::DipoleyMin(double Enu, double mHNL, double target_mass) {
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

    if(small_r)
        return m2 * r4 / m2sub1sq;
    else {
        double root = std::sqrt(m2sub1sq + (r4 - 2 * (1 + m2) * r2));
        return 0.5 * (1 + m4 - r2 - root + m2 * (-2 - r2 + root));
    }
}


double DipoleFromTable::DipoleyMax(double Enu, double mHNL, double target_mass) {
    double target_mass2 = target_mass * target_mass;
    double target_mass4 = target_mass2 * target_mass2;
    double mHNL2 = mHNL * mHNL;

    double s = 2 * Enu * target_mass + target_mass2;
    double s2 = s * s;

    return 0.5 * (target_mass4 - mHNL2*s + s2 - target_mass2*(mHNL2+2*s) + (s - target_mass2) * std::sqrt(target_mass4 + std::pow(mHNL2 - s, 2) - 2*target_mass2*(mHNL2+s))) / s2;
}


double DipoleFromTable::TotalCrossSection(InteractionRecord const & interaction) const {
    LeptonInjector::Particle::ParticleType primary_type = interaction.signature.primary_type;
    LeptonInjector::Particle::ParticleType target_type = interaction.signature.target_type;
    double primary_energy;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
    } else {
        throw std::runtime_error("Lorentz boost not implemented!");
    }
    return TotalCrossSection(primary_type, primary_energy, target_type);
}

double DipoleFromTable::TotalCrossSection(LeptonInjector::Particle::ParticleType primary_type, double primary_energy, Particle::ParticleType target_type) const {
    if(not primary_types.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }

    if(total.find(target_type) == total.end()) {
				std::cout << "Faulty target: " << target_type << std::endl;
        throw std::runtime_error("Supplied target not supported by cross section!");
    }

    Interpolator1D<double> const & interp = total.at(target_type);

    if(primary_energy < interp.MinX() or primary_energy > interp.MaxX()) {
        throw std::runtime_error("Interaction energy ("+ std::to_string(primary_energy) +
                ") out of cross section table range: ["
                + std::to_string(interp.MinX()) + " GeV,"
                + std::to_string(interp.MaxX()) + " GeV]");
    }

    return interp(primary_energy);
}

double DipoleFromTable::DifferentialCrossSection(InteractionRecord const & interaction) const {
    LeptonInjector::Particle::ParticleType primary_type = interaction.signature.primary_type;
    LeptonInjector::Particle::ParticleType target_type = interaction.signature.target_type;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
        p1_lab = p1;
        p2_lab = p2;
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p2_lab = boost_start_to_lab * p2;
        primary_energy = p1_lab.e();
    }
    assert(interaction.signature.secondary_types.size() == 2);
    assert(interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4 or interaction.signature.secondary_types[1] == Particle::ParticleType::NuF4 or interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4Bar or interaction.signature.secondary_types[1] == Particle::ParticleType::NuF4Bar);
    unsigned int lepton_index = (interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4 or interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4Bar) ? 0 : 1;
    unsigned int other_index = 1 - lepton_index;

    std::array<double, 4> const & mom3 = interaction.secondary_momenta[lepton_index];
    std::array<double, 4> const & mom4 = interaction.secondary_momenta[other_index];
    rk::P4 p3(geom3::Vector3(mom3[1], mom3[2], mom3[3]), interaction.secondary_masses[lepton_index]);
    rk::P4 p4(geom3::Vector3(mom4[1], mom4[2], mom4[3]), interaction.secondary_masses[other_index]);

    double y = 1.0 - p2.dot(p3) / p2.dot(p1);

    return DifferentialCrossSection(primary_type, primary_energy, target_type, y);
}

double DipoleFromTable::DifferentialCrossSection(Particle::ParticleType primary_type, double primary_energy, Particle::ParticleType target_type, double y) const {
    if(not primary_types.count(primary_type)) {
        return 0.0;
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }

    if(total.find(target_type) == total.end()) {
        return 0.0;
        throw std::runtime_error("Supplied target not supported by cross section!");
    }

    Interpolator2D<double> const & interp = differential.at(target_type);

    if(primary_energy < interp.MinX() or primary_energy > interp.MaxX()) {
        return 0.0;
        throw std::runtime_error("Interaction energy ("+ std::to_string(primary_energy) +
                ") out of differential cross section table range: ["
                + std::to_string(interp.MinX()) + " GeV, "
                + std::to_string(interp.MaxX()) + " GeV]");
    }

    if(y < interp.MinY() or y > interp.MaxY()) {
        return 0.0;
        throw std::runtime_error("Bjorken y ("+ std::to_string(y) +
                ") out of cross section table range: ["
                + std::to_string(interp.MinY()) + ", "
                + std::to_string(interp.MaxY()) + "]");
    }

    return interp(primary_energy, y);
}

double DipoleFromTable::InteractionThreshold(InteractionRecord const & interaction) const {
    return hnl_mass + (hnl_mass*hnl_mass)/(2*interaction.target_mass);
}

void DipoleFromTable::SampleFinalState(LeptonInjector::InteractionRecord& interaction, std::shared_ptr<LeptonInjector::LI_random> random) const {
    Interpolator2D<double> const & diff_table = differential.at(interaction.signature.target_type);

    // Uses Metropolis-Hastings Algorithm!
    // useful for cases where we don't know the supremum of our distribution, and the distribution is multi-dimensional

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);                                                                                             rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * target_mass_ + 2 * target_mass_ * primary_energy;
    double s = rk::invMass(p1, p2);

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p1_lab = p1;
        p2_lab = p2;
        primary_energy = p1_lab.e();
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p2_lab = boost_start_to_lab * p2;
        primary_energy = p1_lab.e();
    }

    unsigned int lepton_index = (interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4 or interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4Bar) ? 0 : 1;
    unsigned int other_index = 1 - lepton_index;
    double m = hnl_mass;

    // double m1 = p1_lab | p1_lab;
    double m1 = interaction.primary_mass;
    // double m2 = p2_lab | p2_lab;
    double m2 = interaction.target_mass;
    double m3 = m;
    double E1_lab = p1_lab.e();
    double E2_lab = p2_lab.e();

    double yMin = DipoleyMin(E1_lab, GetHNLMass(), interaction.target_mass);
    double yMax = DipoleyMax(E1_lab, GetHNLMass(), interaction.target_mass);
    assert(yMin > 0);
    double log_yMax = log10(yMax);
    double log_yMin = log10(yMin);
    double min_Q2 = yMin * s;
    double z; // placeholder for z sampling

    bool accept;

    // kin_vars and its twin are 2-vectors containing [nu-energy, Bjorken Y]
    std::array<double,2> kin_vars, test_kin_vars;

    // values of cross_section from the table
    double cross_section, test_cross_section;

    // No matter what, we're evaluating at this specific energy.
    kin_vars[0] = test_kin_vars[0] = primary_energy;

    // check preconditions
    if(kin_vars[0] < diff_table.MinX()
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
        if(std::isnan(test_cross_section) or test_cross_section <= 0)
            continue;

        double odds = (test_cross_section / cross_section);
        accept = (cross_section == 0 || (odds > 1.) || random->Uniform(0, 1) < odds);

        if(accept) {
            kin_vars = test_kin_vars;
            cross_section = test_cross_section;
        }
    }
    double final_y = kin_vars[1];

    interaction.interaction_parameters.resize(2);
    interaction.interaction_parameters[0] = E1_lab;
    interaction.interaction_parameters[1] = final_y;

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 p1_mom = p1_lab.momentum();
    geom3::UnitVector3 p1_lab_dir = p1_mom.direction();
    geom3::Rotation3 x_to_p1_lab_rot = geom3::rotationBetween(x_dir, p1_lab_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(p1_lab_dir, phi);

    double E3_lab = E1_lab - E1_lab * final_y;
    double p1x_lab = p1_mom.length();
    double p3_lab_sq = E3_lab * E3_lab - m3 * m3;
    double p3x_lab_frac = (p1x_lab * p1x_lab - m3 * m3 + E1_lab * E1_lab * (1.0 - 2.0 * final_y)) / (2.0 * p1x_lab * E3_lab);
    double p3x_lab = p3x_lab_frac * sqrt(p3_lab_sq);
    double p3y_lab = sqrt(p3_lab_sq - p3x_lab * p3x_lab);

    rk::P4 p3_lab(geom3::Vector3(p3x_lab, p3y_lab, 0), m3);
    p3_lab.rotate(x_to_p1_lab_rot);
    p3_lab.rotate(rand_rot);
    rk::P4 p4_lab = p2_lab + (p1_lab - p3_lab);

    rk::P4 p3;
    rk::P4 p4;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p3 = p3_lab;
        p4 = p4_lab;
    } else {
        rk::Boost boost_lab_to_start = p2.restBoost();
        p3 = boost_lab_to_start * p3_lab;
        p4 = boost_lab_to_start * p4_lab;
    }

    rk::P4 pq_13 = p1 - p3;
    rk::P4 pq_24 = p4 - p2;

    // Check that computed q2 in the start frame matches up with the specified Q2
    //assert(std::abs(double(pq_13.dot(pq_13)) + Q2) < std::abs(Q2 * 1e-3));
    //assert(std::abs(double(pq_24.dot(pq_24)) + Q2) < std::abs(Q2 * 1e-3));

    interaction.secondary_momenta.resize(2);
    interaction.secondary_masses.resize(2);
    interaction.secondary_helicity.resize(2);

    interaction.secondary_momenta[lepton_index][0] = p3.e(); // p3_energy
    interaction.secondary_momenta[lepton_index][1] = p3.px(); // p3_x
    interaction.secondary_momenta[lepton_index][2] = p3.py(); // p3_y
    interaction.secondary_momenta[lepton_index][3] = p3.pz(); // p3_z
    interaction.secondary_masses[lepton_index] = p3.m();

    double helicity_mul = 0.0;
    if(channel == Conserving)
        helicity_mul = 1.0;
    else if(channel == Flipping)
        helicity_mul = -1.0;

    interaction.secondary_helicity[lepton_index] = std::copysign(0.5, interaction.primary_helicity * helicity_mul);

    interaction.secondary_momenta[other_index][0] = p4.e(); // p4_energy
    interaction.secondary_momenta[other_index][1] = p4.px(); // p4_x
    interaction.secondary_momenta[other_index][2] = p4.py(); // p4_y
    interaction.secondary_momenta[other_index][3] = p4.pz(); // p4_z
    interaction.secondary_masses[other_index] = p4.m();

    interaction.secondary_helicity[other_index] = std::copysign(interaction.target_helicity, interaction.target_helicity * helicity_mul);
}

double DipoleFromTable::FinalStateProbability(LeptonInjector::InteractionRecord const & interaction) const {
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

std::vector<Particle::ParticleType> DipoleFromTable::GetPossibleTargets() const {
    std::set<Particle::ParticleType> diff_targets;
    std::set<Particle::ParticleType> tot_targets;
    for(auto const & diff : differential)
        diff_targets.insert(diff.first);
    for(auto const & tot : total)
        tot_targets.insert(tot.first);
    std::vector<Particle::ParticleType> res;
    std::set_intersection(diff_targets.begin(), diff_targets.end(), tot_targets.begin(), tot_targets.end(), std::back_inserter(res));
    return res;
}

std::vector<Particle::ParticleType> DipoleFromTable::GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const {
    if(not primary_types.count(primary_type)) {
        return std::vector<Particle::ParticleType>();
    }
    return GetPossibleTargets();
}

std::vector<Particle::ParticleType> DipoleFromTable::GetPossiblePrimaries() const {
    return std::vector<Particle::ParticleType>(primary_types.begin(), primary_types.end());
}

std::vector<InteractionSignature> DipoleFromTable::GetPossibleSignatures() const {
    std::vector<Particle::ParticleType> targets = GetPossibleTargets();
    std::vector<InteractionSignature> signatures;
    InteractionSignature signature;
    signature.secondary_types.resize(2);

    for(auto primary : primary_types) {
        signature.primary_type = primary;
        if(std::set<Particle::ParticleType>{Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau}.count(primary))
            signature.secondary_types[0] = Particle::ParticleType::NuF4;
        else if(std::set<Particle::ParticleType>{Particle::ParticleType::NuEBar, Particle::ParticleType::NuMuBar, Particle::ParticleType::NuTauBar}.count(primary))
            signature.secondary_types[0] = Particle::ParticleType::NuF4Bar;
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

std::vector<InteractionSignature> DipoleFromTable::GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const {
    std::vector<Particle::ParticleType> targets = GetPossibleTargets();
    if(primary_types.count(primary_type) > 0 and std::find(targets.begin(), targets.end(), target_type) != targets.end()) {
        InteractionSignature signature;
        signature.secondary_types.resize(2);
        signature.primary_type = primary_type;
        signature.target_type = target_type;
        signature.secondary_types[1] = target_type;
        if(std::set<Particle::ParticleType>{Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau}.count(primary_type))
            signature.secondary_types[0] = Particle::ParticleType::NuF4;
        else if(std::set<Particle::ParticleType>{Particle::ParticleType::NuEBar, Particle::ParticleType::NuMuBar, Particle::ParticleType::NuTauBar}.count(primary_type))
            signature.secondary_types[0] = Particle::ParticleType::NuF4Bar;
        else
            throw std::runtime_error("Primary type not in primary_types!");
        return std::vector<InteractionSignature>{signature};
    } else {
        return std::vector<InteractionSignature>();
    }
}

std::vector<std::string> DipoleFromTable::DensityVariables() const {
    return std::vector<std::string>{"Bjorken y"};
}

namespace {
bool fexists(const char *filename)
{
    std::ifstream ifile(filename);
    return (bool)ifile;
}
bool fexists(const std::string filename)
{
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}
}

void DipoleFromTable::AddDifferentialCrossSectionFile(std::string filename, Particle::ParticleType target) {
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
    Particle::ParticleType pid = (Particle::ParticleType)pid_int;

    if(fexists(filename)) {
        std::ifstream in(filename.c_str());
        std::string buf;

        TableData2D<double> table_data;
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
        Interpolator2D<double> interp(table_data);
        AddDifferentialCrossSection(pid, interp);
    } else {
        throw std::runtime_error("Failed open cross section file!");
    }
}

void DipoleFromTable::AddTotalCrossSectionFile(std::string filename, Particle::ParticleType target) {
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
    Particle::ParticleType pid = (Particle::ParticleType)pid_int;

    if(pid != target) {
        throw std::runtime_error("File target nucleus does not match supplied target type!");
    }

    if(fexists(filename)) {
        std::ifstream in(filename.c_str());
        std::string buf;

        TableData1D<double> table_data;
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
        Interpolator1D<double> interp(table_data);
        AddTotalCrossSection(pid, interp);
    } else {
        throw std::runtime_error("Failed open cross section file!");
    }
}

void DipoleFromTable::AddDifferentialCrossSection(Particle::ParticleType target, Interpolator2D<double> interp) {
    differential.insert(std::make_pair(target, interp));
}

void DipoleFromTable::AddTotalCrossSection(Particle::ParticleType target, Interpolator1D<double> interp) {
    total.insert(std::make_pair(target, interp));
}

////////////////////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& os, InteractionSignature const& signature)
{
    std::stringstream ss;
    ss << "InteractionSignature (" << &signature << ") ";
    os << ss.str() << '\n';


    os << "PrimaryType: " << signature.primary_type << "\n";
    os << "TargetType: " << signature.target_type << "\n";
    os << "SecondaryTypes:";
    for(auto secondary: signature.secondary_types) {
        os << " " << secondary;
    }
    os << std::endl;

    return os;
}

std::ostream& operator<<(std::ostream& os, InteractionRecord const& record)
{
    std::stringstream ss;
    ss << "InteractionRecord (" << &record << ") ";
    os << ss.str() << '\n';
    os << "Signature(" << &record.signature << "): " << record.signature.primary_type << " + " << record.signature.target_type << " ->";
    for(auto secondary: record.signature.secondary_types) {
        os << " " << secondary;
    }
    os << "\n";

    os << "InteractionVertex: " << record.interaction_vertex[0] << " " << record.interaction_vertex[1] << " " << record.interaction_vertex[2] << "\n";
    os << "PrimaryMass: " << record.primary_mass << "\n";
    os << "PrimaryMomentum: " << record.primary_momentum[0] << " " << record.primary_momentum[1] << " " << record.primary_momentum[2] << " " << record.primary_momentum[3] << "\n";
    os << "TargetMass: " << record.target_mass << "\n";
    os << "TargetMomentum: " << record.target_momentum[0] << " " << record.target_momentum[1] << " " << record.target_momentum[2] << " " << record.target_momentum[3] << "\n";
    os << "SecondaryMomenta:\n";
    for(auto const & secondary: record.secondary_momenta) {
        os << "\t" << secondary[0] << " " << secondary[1] << " " << secondary[2] << " " << secondary[3] << "\n";
    }
    os << "SecondaryMasses:\n";
    for(auto const & secondary: record.secondary_masses) {
        os << "\t" << secondary << "\n";
    }
    os << "InteractionParameters:";
    for(auto param: record.interaction_parameters) {
        os << " " << param;
    }
    os << std::endl;

    return os;
}

} // namespace LeptonInjector

