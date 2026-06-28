#include "SIREN/interactions/QuarkDISFromSpline.h"

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
#include <map>
#include <limits>
#include <cmath>
#include <string>
#include <stdexcept>
#include <iostream>

#include <rk/rk.hh>                                        // for P4, Boost
#include <rk/geom3.hh>                                     // for Vector3

#include <photospline/splinetable.h>                       // for splinetable
//#include <photospline/cinter/splinetable.h>

#include "SIREN/interactions/CrossSection.h"     // for CrossSection
#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/dataclasses/Particle.h"           // for Particle
#include "SIREN/utilities/Random.h"               // for SIREN_random
#include "SIREN/utilities/Constants.h"            // for electronMass
#include "SIREN/utilities/Errors.h"                  // for PythonImplementationError


namespace siren {
namespace interactions {

namespace {
// Slow-rescaling helpers (lab frame, stationary nucleon target).
//   E   = primary neutrino lab energy
//   M   = target nucleon mass
//   mc  = charm-quark mass
inline double slowRescalingQ2(double xi, double y, double E, double M, double mc) {
    return 2.0 * M * E * y * xi - mc * mc;
}
inline double xiToBjorkenX(double xi, double y, double E, double M, double mc) {
    return xi - (mc * mc) / (2.0 * M * E * y);
}
inline double slowRescalingW2(double xi, double y, double E, double M, double mc) {
    // W^2 = M^2 + 2 M E y (1 - xi) + m_c^2  =  M^2 + (Q^2 + m_c^2)/xi - Q^2
    return M * M + 2.0 * M * E * y * (1.0 - xi) + mc * mc;
}

///\brief Slow-rescaling kinematic check.
///
/// Acceptance cuts:
///   xi in [1e-9, 1], y in [1e-9, 1 - m_lep/E],
///   Q^2 = 2 M E y xi - m_c^2 > 0  (charm threshold),
///   W^2 = M^2 + 2 M E y (1-xi) + m_c^2 > (M + M_D0)^2.
bool kinematicallyAllowed(double xi, double y, double E, double M, double m_lep) {
    if (xi < 1e-9 || xi > 1.0)              return false;
    if (y  < 1e-9)                          return false;
    if (y  > 1.0 - m_lep / E)               return false;

    const double mc  = siren::utilities::Constants::charmMass;
    const double Mch = siren::utilities::Constants::D0Mass;

    const double Q2 = slowRescalingQ2(xi, y, E, M, mc);
    if (Q2 <= 0.0)                          return false;

    const double W2 = slowRescalingW2(xi, y, E, M, mc);
    if (W2 <= (M + Mch) * (M + Mch))        return false;

    // Transverse-momentum balance: the exchanged q must have a real transverse
    // component pqy in the lab frame. This is the same q-decomposition the sampler
    // uses (SampleFinalState), so kinematicallyAllowed is the single predicate
    // shared by the sampler proposal loops and by
    // DifferentialCrossSection/FinalStateProbability. Without this check the
    // density (dxs/txs) would be nonzero on points the sampler rejects (pqy^2 < 0),
    // breaking Sample==Density closure and silently biasing low-Bjorken-x events.
    //
    // Primary is a neutrino here (InitializeSignatures enforces isNeutrino), so
    // m1 = 0 and |p1| = E. Using P1 = E (massless primary) to match the sampler.
    // pqy^2 = momq^2 - pqx^2 with:
    //   pqx  = (m_lep^2 + 2 P1^2 + Q2 + 2 E^2 (y-1)) / (2 P1)
    //   momq^2 = P1^2 + Q2 + E^2 (y^2 - 1)
    const double P1 = E;
    const double pqx = (m_lep * m_lep + 2.0 * P1 * P1 + Q2 + 2.0 * E * E * (y - 1.0)) / (2.0 * P1);
    const double momq2 = P1 * P1 + Q2 + E * E * (y * y - 1.0);
    const double pqy2 = momq2 - pqx * pqx;
    if (pqy2 < 0.0)                          return false;

    return true;
}
}

QuarkDISFromSpline::QuarkDISFromSpline() {
    // initialize the pdf normalization and cdf table for the hadronization process
    normalize_pdf();
    compute_cdf();
}

QuarkDISFromSpline::QuarkDISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minimum_Q2, std::set<siren::dataclasses::ParticleType> primary_types, std::set<siren::dataclasses::ParticleType> target_types, std::string units) : primary_types_(primary_types), target_types_(target_types), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    normalize_pdf();
    compute_cdf();
    LoadFromMemory(differential_data, total_data);
    SetInteractionType(interaction);
    InitializeSignatures();
    SetUnits(units);
}

QuarkDISFromSpline::QuarkDISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minimum_Q2, std::vector<siren::dataclasses::ParticleType> primary_types, std::vector<siren::dataclasses::ParticleType> target_types, std::string units) : primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    normalize_pdf();
    compute_cdf();
    LoadFromMemory(differential_data, total_data);
    SetInteractionType(interaction);
    InitializeSignatures();
    SetUnits(units);
}

QuarkDISFromSpline::QuarkDISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minimum_Q2, std::set<siren::dataclasses::ParticleType> primary_types, std::set<siren::dataclasses::ParticleType> target_types, std::string units) : primary_types_(primary_types), target_types_(target_types), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    normalize_pdf();
    compute_cdf();
    LoadFromFile(differential_filename, total_filename);
    SetInteractionType(interaction);
    InitializeSignatures();
    SetUnits(units);
}

QuarkDISFromSpline::QuarkDISFromSpline(std::string differential_filename, std::string total_filename, std::set<siren::dataclasses::ParticleType> primary_types, std::set<siren::dataclasses::ParticleType> target_types, std::string units) : primary_types_(primary_types), target_types_(target_types) {
    normalize_pdf();
    compute_cdf();
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
    SetUnits(units);
}

QuarkDISFromSpline::QuarkDISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minimum_Q2, std::vector<siren::dataclasses::ParticleType> primary_types, std::vector<siren::dataclasses::ParticleType> target_types, std::string units) : primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()), interaction_type_(interaction), target_mass_(target_mass), minimum_Q2_(minimum_Q2) {
    normalize_pdf();
    compute_cdf();
    LoadFromFile(differential_filename, total_filename);
    SetInteractionType(interaction);
    InitializeSignatures();
    SetUnits(units);
}

QuarkDISFromSpline::QuarkDISFromSpline(std::string differential_filename, std::string total_filename, std::vector<siren::dataclasses::ParticleType> primary_types, std::vector<siren::dataclasses::ParticleType> target_types, std::string units) : primary_types_(primary_types.begin(), primary_types.end()), target_types_(target_types.begin(), target_types.end()) {
    normalize_pdf();
    compute_cdf();
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
    SetUnits(units);
}

void QuarkDISFromSpline::SetUnits(std::string units) {
    std::transform(units.begin(), units.end(), units.begin(),
        [](unsigned char c){ return std::tolower(c); });
    if(units == "cm") {
        unit = 1.0;
    } else if(units == "m") {
        unit = 10000.0;
    } else {
        throw std::runtime_error("Cross section units not supported!");
    }
}

void QuarkDISFromSpline::SetInteractionType(int interaction) {
    interaction_type_ = interaction;
}

bool QuarkDISFromSpline::equal(CrossSection const & other) const {
    const QuarkDISFromSpline* x = dynamic_cast<const QuarkDISFromSpline*>(&other);
    // to do: include more features in the hadronization side to check equivalence
    if(!x)
        return false;
    else
        return
            std::tie(
            interaction_type_,
            target_mass_,
            minimum_Q2_,
            signatures_,
            primary_types_,
            target_types_,
            differential_cross_section_,
            total_cross_section_)
            ==
            std::tie(
            x->interaction_type_,
            x->target_mass_,
            x->minimum_Q2_,
            x->signatures_,
            x->primary_types_,
            x->target_types_,
            x->differential_cross_section_,
            x->total_cross_section_);
}

void QuarkDISFromSpline::LoadFromFile(std::string dd_crossSectionFile, std::string total_crossSectionFile) {

    differential_cross_section_ = photospline::splinetable<>(dd_crossSectionFile.c_str());

    if(differential_cross_section_.get_ndim()!=3 && differential_cross_section_.get_ndim()!=2)
        throw std::runtime_error("cross section spline has " + std::to_string(differential_cross_section_.get_ndim())
                + " dimensions, should have either 3 (log10(E), log10(x), log10(y)) or 2 (log10(E), log10(y))");

    total_cross_section_ = photospline::splinetable<>(total_crossSectionFile.c_str());

    if(total_cross_section_.get_ndim() != 1)
        throw std::runtime_error("Total cross section spline has " + std::to_string(total_cross_section_.get_ndim())
                + " dimensions, should have 1, log10(E)");
}

void QuarkDISFromSpline::LoadFromMemory(std::vector<char> & differential_data, std::vector<char> & total_data) {
    differential_cross_section_.read_fits_mem(differential_data.data(), differential_data.size());
    total_cross_section_.read_fits_mem(total_data.data(), total_data.size());
}

double QuarkDISFromSpline::GetLeptonMass(siren::dataclasses::ParticleType lepton_type) {
    int32_t lepton_number = std::abs(static_cast<int32_t>(lepton_type));
    double lepton_mass;
    switch(lepton_number) {
        case 11:
            lepton_mass = siren::utilities::Constants::electronMass;
            break;
        case 13:
            lepton_mass = siren::utilities::Constants::muonMass;
            break;
        case 15:
            lepton_mass = siren::utilities::Constants::tauMass;
            break;
        case 12:
            lepton_mass = 0;
            break;
        case 14:
            lepton_mass = 0;
            break;
        case 16:
            lepton_mass = 0;
            break;
        default:
            throw std::runtime_error("Unknown lepton type!");
    }
    return lepton_mass;
}

double QuarkDISFromSpline::getHadronMass(siren::dataclasses::ParticleType hadron_type) {
    switch(hadron_type){
			case siren::dataclasses::ParticleType::D0:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::D0Bar:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::DPlus:
				return( siren::utilities::Constants::DPlusMass);
			case siren::dataclasses::ParticleType::DMinus:
				return( siren::utilities::Constants::DPlusMass);
			case siren::dataclasses::ParticleType::DsPlus:
				return( siren::utilities::Constants::DsPlusMass);
			case siren::dataclasses::ParticleType::DsMinus:
				return( siren::utilities::Constants::DsMinusMass);
            default:
                return(0.0);
        }
}


std::map<std::string, int> QuarkDISFromSpline::getIndices(siren::dataclasses::InteractionSignature signature) {
    // Initialize to -1 so an unmatched slot is detectable instead of being read
    // back as an uninitialized index into the secondary arrays (UB / OOB).
    int lepton_id = -1, hadron_id = -1, meson_id = -1;
    for (size_t i = 0; i < signature.secondary_types.size(); i++){
        if (isLepton(signature.secondary_types[i])) {
            lepton_id = i;
            continue;
        } else if (isD(signature.secondary_types[i])) {
            meson_id = i;
            continue;
        } else {
            hadron_id = i;
            continue;
        }
    }
    // A valid charm-DIS signature is (lepton, hadron, D-meson); a missing slot
    // (e.g. the malformed {Hadrons,Hadrons,D} signature with no lepton) would
    // otherwise leave an index at -1 and corrupt downstream array access.
    if (lepton_id < 0 || hadron_id < 0 || meson_id < 0) {
        throw std::runtime_error("QuarkDISFromSpline::getIndices: signature does not contain the expected (lepton, hadron, D-meson) triplet");
    }
    return {{"lepton", lepton_id}, {"hadron", hadron_id}, {"meson", meson_id}};
}


void QuarkDISFromSpline::ReadParamsFromSplineTable() {
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
                target_mass_ = siren::utilities::Constants::isoscalarMass;
            } else if(interaction_type_ == 3) {
                target_mass_ = siren::utilities::Constants::electronMass;
            } else {
                throw std::runtime_error("Logic error. Interaction type is not 1, 2, or 3!");
            }

        } else {
            if(differential_cross_section_.get_ndim() == 3) {
                target_mass_ = siren::utilities::Constants::isoscalarMass;
            } else if(differential_cross_section_.get_ndim() == 2) {
                target_mass_ = siren::utilities::Constants::electronMass;
            } else {
                throw std::runtime_error("Logic error. Spline dimensionality is not 2, or 3!");
            }
        }
    }
}

std::set<siren::dataclasses::ParticleType> QuarkDISFromSpline::DTypesForPrimary(siren::dataclasses::ParticleType primary) {
    if (primary == siren::dataclasses::ParticleType::NuE
        || primary == siren::dataclasses::ParticleType::NuMu
        || primary == siren::dataclasses::ParticleType::NuTau) {
        return {siren::dataclasses::Particle::ParticleType::D0,
                siren::dataclasses::Particle::ParticleType::DPlus,
                siren::dataclasses::Particle::ParticleType::DsPlus};
    } else if (primary == siren::dataclasses::ParticleType::NuEBar
        || primary == siren::dataclasses::ParticleType::NuMuBar
        || primary == siren::dataclasses::ParticleType::NuTauBar) {
        return {siren::dataclasses::Particle::ParticleType::D0Bar,
                siren::dataclasses::Particle::ParticleType::DMinus,
                siren::dataclasses::Particle::ParticleType::DsMinus};
    } else {
        throw std::runtime_error("DTypesForPrimary: Unknown neutrino primary type!");
    }
}

void QuarkDISFromSpline::InitializeSignatures() {
    signatures_.clear();
    for(auto primary_type : primary_types_) {
        dataclasses::InteractionSignature signature;
        signature.primary_type = primary_type;
        if(not isNeutrino(primary_type)) {
            throw std::runtime_error("This DIS implementation only supports neutrinos as primaries!");
        }
        // first push back the charged lepton product
        siren::dataclasses::ParticleType charged_lepton_product = siren::dataclasses::ParticleType::unknown;
        siren::dataclasses::ParticleType neutral_lepton_product = primary_type;
        if(primary_type == siren::dataclasses::ParticleType::NuE) {
            charged_lepton_product = siren::dataclasses::ParticleType::EMinus;
        } else if(primary_type == siren::dataclasses::ParticleType::NuEBar) {
            charged_lepton_product = siren::dataclasses::ParticleType::EPlus;
        } else if(primary_type == siren::dataclasses::ParticleType::NuMu) {
            charged_lepton_product = siren::dataclasses::ParticleType::MuMinus;
        } else if(primary_type == siren::dataclasses::ParticleType::NuMuBar) {
            charged_lepton_product = siren::dataclasses::ParticleType::MuPlus;
        } else if(primary_type == siren::dataclasses::ParticleType::NuTau) {
            charged_lepton_product = siren::dataclasses::ParticleType::TauMinus;
        } else if(primary_type == siren::dataclasses::ParticleType::NuTauBar) {
            charged_lepton_product = siren::dataclasses::ParticleType::TauPlus;
        } else {
            throw std::runtime_error("InitializeSignatures: Unkown parent neutrino type!");
        }
        if(interaction_type_ == 1) {
            signature.secondary_types.push_back(charged_lepton_product);
        } else if(interaction_type_ == 2) {
            signature.secondary_types.push_back(neutral_lepton_product);
        } else if(interaction_type_ == 3) {
            // interaction_type 3 (electron / e-target) has no outgoing lepton in
            // this charm convention and would emit a malformed {Hadrons,Hadrons,D}
            // signature with no lepton index. The differential cross section and
            // sampler both require a lepton to form q = p1 - p3, so reject it at
            // construction rather than create an unusable signature.
            throw std::runtime_error("QuarkDISFromSpline: interaction_type 3 (e-target) not supported for charm D-meson production; use 1 (CC) or 2 (NC)");
        } else {
            throw std::runtime_error("InitializeSignatures: Unkown interaction type!");
        }
        // now push back the hadron product
        signature.secondary_types.push_back(siren::dataclasses::ParticleType::Hadrons);
        std::set<siren::dataclasses::ParticleType> d_types_local = DTypesForPrimary(primary_type);
        for (auto meson_type : d_types_local) {
            dataclasses::InteractionSignature full_signature = signature;
            full_signature.secondary_types.push_back(meson_type);
            // and finally set the target type and push back the entire signature as well as sig by target
            for(auto target_type : target_types_) {
                full_signature.target_type = target_type;

                signatures_.push_back(full_signature);

                std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType> key(primary_type, target_type);
                signatures_by_parent_types_[key].push_back(full_signature);
            }
        }
    }
}

void QuarkDISFromSpline::normalize_pdf() {
    if (fragmentation_integral == 0){
         std::function<double(double)> integrand = [&] (double x) -> double {
            return (0.8 / x ) / (std::pow(1 - (1 / x) - (0.2 / (1 - x)), 2));
        };
        fragmentation_integral = siren::utilities::rombergIntegrate(integrand, 0.001, 0.999);
    } else {
        return;
    }
}

void QuarkDISFromSpline::compute_cdf() {
    // first set the z nodes
    std::vector<double> zspline;
    for (int i = 0; i < 100; ++i) {
        zspline.push_back(0.01 + i * (0.99-0.01) / 100 );
    }

    // declare the cdf vectors
    std::vector<double> cdf_vector;
    std::vector<double> cdf_z_nodes;
    std::vector<double> pdf_vector;

    cdf_z_nodes.push_back(0);
    cdf_vector.push_back(0);
    pdf_vector.push_back(0);

    // compute the spline table
    for (int i = 0; i < zspline.size(); ++i) {
        if (i == 0) {
            double cur_z = zspline[i];
            double cur_pdf = sample_pdf(cur_z);
            double area = cur_z * cur_pdf * 0.5;
            pdf_vector.push_back(cur_pdf);
            cdf_vector.push_back(area);
            cdf_z_nodes.push_back(cur_z);
            continue;
        }
        double cur_z = zspline[i];
        double cur_pdf = sample_pdf(cur_z);
        double area = 0.5 * (pdf_vector[i - 1] + cur_pdf) * (zspline[i] - zspline[i - 1]);
        pdf_vector.push_back(cur_pdf);
        cdf_z_nodes.push_back(cur_z);
        cdf_vector.push_back(area + cdf_vector.back());
    }

    cdf_z_nodes.push_back(1);
    cdf_vector.push_back(1);
    pdf_vector.push_back(0);


    // set the spline table
    siren::utilities::TableData1D<double> inverse_cdf_data;
    inverse_cdf_data.x = cdf_vector;
    inverse_cdf_data.f = cdf_z_nodes;

    inverseCdfTable = siren::utilities::Interpolator1D<double>(inverse_cdf_data);

    return;
}

double QuarkDISFromSpline::sample_pdf(double x) const {
    return (0.8 / x ) / (std::pow(1 - (1 / x) - (0.2 / (1 - x)), 2)) / fragmentation_integral;
}

double QuarkDISFromSpline::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    siren::dataclasses::ParticleType primary_type = interaction.signature.primary_type;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    double primary_energy;
    primary_energy = interaction.primary_momentum[0];
    // if we are below threshold, return 0
    if(primary_energy < InteractionThreshold(interaction)) {
        return 0;
    }
    double total_xs = TotalCrossSection(primary_type, primary_energy);
    // Apply fragmentation fraction for the specific D meson in this signature
    // so that summing over signatures (D0 + D+) gives the correct total
    for (auto const & sec_type : interaction.signature.secondary_types) {
        if (siren::dataclasses::isD(sec_type)) {
            total_xs *= FragmentationFraction(sec_type);
            break;
        }
    }
    return total_xs;
}

double QuarkDISFromSpline::TotalCrossSection(siren::dataclasses::ParticleType primary_type, double primary_energy) const {
    if(not primary_types_.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
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
    if (std::pow(10.0, log_xs) == 0) {
    }

    return unit * std::pow(10.0, log_xs);
}

double QuarkDISFromSpline::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), interaction.target_mass);
    double primary_energy;
    primary_energy = interaction.primary_momentum[0];
    assert(interaction.signature.secondary_types.size() == 3);
    std::map<std::string, int> secondaries = getIndices(interaction.signature);
    unsigned int lepton_index = secondaries["lepton"];
    unsigned int hadron_index = secondaries["hadron"];
    unsigned int meson_index = secondaries["meson"];

    std::array<double, 4> const & mom3 = interaction.secondary_momenta[lepton_index];
    std::array<double, 4> const & mom_x = interaction.secondary_momenta[hadron_index];
    std::array<double, 4> const & mom_d = interaction.secondary_momenta[meson_index];

    rk::P4 p3(geom3::Vector3(mom3[1], mom3[2], mom3[3]), interaction.secondary_masses[lepton_index]);
    rk::P4 p_x(geom3::Vector3(mom_x[1], mom_x[2], mom_x[3]), interaction.secondary_masses[hadron_index]);
    rk::P4 p_d(geom3::Vector3(mom_d[1], mom_d[2], mom_d[3]), interaction.secondary_masses[meson_index]);
    rk::P4 p4 = p_x + p_d; // this assume that we are working in a good frame where the hadronization vertex has 4-momentum conserved
    rk::P4 q = p1 - p3;
    // however p4 is not used in computation here so we should be fine...

    double Q2 = -q.dot(q);
    double lepton_mass = GetLeptonMass(interaction.signature.secondary_types[lepton_index]);

    const double mc = siren::utilities::Constants::charmMass;
    double y  = 1.0 - p2.dot(p3) / p2.dot(p1);
    // xi from inverting Q^2 = 2 M E y xi - m_c^2; guard against y <= 0:
    double xi = (y > 0.0)
                ? (Q2 + mc * mc) / (2.0 * primary_energy * target_mass_ * y)
                : 0.0;
    double log_energy = log10(primary_energy);
    std::array<int,3> centers;

    bool use_sample_kinematics = (xi > 0.0 && y > 0.0 && Q2 >= minimum_Q2_);
    if (use_sample_kinematics) {
        std::array<double,3> coordinates{{log_energy, log10(xi), log10(y)}};
        use_sample_kinematics =
            kinematicallyAllowed(xi, y, primary_energy, target_mass_, lepton_mass)
            && differential_cross_section_.searchcenters(coordinates.data(), centers.data());
    }

    if (!use_sample_kinematics) {
        double E1_lab = interaction.interaction_parameters.at("energy");
        // Fallback: trust stored xi (records sampled by this class always have bjorken_xi).
        // If absent, std::out_of_range is intentional: this class is xi-y only.
        xi = interaction.interaction_parameters.at("bjorken_xi");
        y  = interaction.interaction_parameters.at("bjorken_y");
        Q2 = slowRescalingQ2(xi, y, E1_lab, target_mass_, mc);
    }
    return DifferentialCrossSection(primary_energy, xi, y, lepton_mass, Q2);
}

double QuarkDISFromSpline::DifferentialCrossSection(double energy, double xi, double y, double secondary_lepton_mass, double Q2) const {
    double log_energy = log10(energy);
    // Out of spline SUPPORT -> raise, never silently return 0: a silent zero on a
    // genuinely sampled event would bias that event's physical density (and hence
    // its weight) to zero.
    if (log_energy < differential_cross_section_.lower_extent(0)
            || log_energy > differential_cross_section_.upper_extent(0)) {
        throw std::runtime_error("QuarkDISFromSpline: energy " + std::to_string(energy)
            + " GeV is outside the differential spline energy range ["
            + std::to_string(std::pow(10.0, differential_cross_section_.lower_extent(0))) + ", "
            + std::to_string(std::pow(10.0, differential_cross_section_.upper_extent(0)))
            + "] GeV.");
    }
    if (xi <= 0 || xi >= 1 || y <= 0 || y >= 1) {
        throw std::runtime_error("QuarkDISFromSpline: unphysical (xi="
            + std::to_string(xi) + ", y=" + std::to_string(y) + ") outside (0, 1).");
    }

    if (std::isnan(Q2)) {
        Q2 = slowRescalingQ2(xi, y, energy, target_mass_,
                             siren::utilities::Constants::charmMass);
    }
    if (Q2 < minimum_Q2_) {
        return 0;
    }
    if (!kinematicallyAllowed(xi, y, energy, target_mass_, secondary_lepton_mass)) {
        return 0;
    }
    std::array<double,3> coordinates{{log_energy, log10(xi), log10(y)}};
    std::array<int,3> centers;
    if (!differential_cross_section_.searchcenters(coordinates.data(), centers.data())) {
        throw std::runtime_error("QuarkDISFromSpline: (xi=" + std::to_string(xi)
            + ", y=" + std::to_string(y) + ") at E=" + std::to_string(energy)
            + " GeV is outside the differential spline (xi, y) grid.");
    }
    double result = pow(10., differential_cross_section_.ndsplineeval(coordinates.data(), centers.data(), 0));
    assert(result >= 0);
    if (std::isinf(result)) {
    }
    return unit * result;
}

double QuarkDISFromSpline::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    // Consider implementing DIS thershold at some point
    return 0;
}

void QuarkDISFromSpline::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // first obtain the indices from secondaries
    std::map<std::string, int> secondary_indices = getIndices(record.signature);
    unsigned int lepton_index = secondary_indices["lepton"];
    unsigned int hadron_index = secondary_indices["hadron"];
    unsigned int meson_index = secondary_indices["meson"];

    // Uses Metropolis-Hastings Algorithm!
    // useful for cases where we don't know the supremum of our distribution, and the distribution is multi-dimensional
    if (differential_cross_section_.get_ndim() != 3) {
        throw std::runtime_error("I expected 3 dimensions in the cross section spline, but got " + std::to_string(differential_cross_section_.get_ndim()) +". Maybe your fits file doesn't have the right 'INTERACTION' key?");
    }
    rk::P4 p1(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), target_mass_);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * trecord.secondary_momentarget_mass_ + 2 * target_mass_ * primary_energy;
    // double s = std::pow(rk::invMass(p1, p2), 2);

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    p1_lab = p1;
    p2_lab = p2;
    primary_energy = p1_lab.e();

    // correctly assign lepton, hadron and meson index
    double m = GetLeptonMass(record.signature.secondary_types[lepton_index]);

    double m1 = record.primary_mass;
    double m3 = m;
    double E1_lab = p1_lab.e();
    double E2_lab = p2_lab.e();

    const double mc      = siren::utilities::Constants::charmMass;
    const double Mch     = siren::utilities::Constants::D0Mass;
    const double M_targ  = target_mass_;
    const double E_nu    = primary_energy;
    const double Q2_min  = minimum_Q2_;

    const double yMax    = 1.0 - m / E_nu;
    const double W2_thr  = (M_targ + Mch) * (M_targ + Mch);
    const double yMin    = (W2_thr - M_targ * M_targ + Q2_min) / (2.0 * M_targ * E_nu);
    const double xiMin   = (mc * mc + Q2_min) / (2.0 * M_targ * E_nu * yMax);

    if (xiMin <= 0.0 || yMin <= 0.0 || yMax <= 0.0) {
        throw std::runtime_error(
            "QuarkDISFromSpline: non-positive sampling-bound (xiMin="
            + std::to_string(xiMin) + ", yMin=" + std::to_string(yMin)
            + ", yMax=" + std::to_string(yMax) + ")");
    }
    if (xiMin >= 1.0 || yMin >= yMax) {
        throw std::runtime_error(
            "QuarkDISFromSpline: primary energy below slow-rescaling charm threshold "
            "(xiMin=" + std::to_string(xiMin) + ", yMin=" + std::to_string(yMin) +
            ", yMax=" + std::to_string(yMax) + ")");
    }

    const double logXiMin = std::log10(xiMin);
    const double logYMin  = std::log10(yMin);
    const double logYMax  = std::log10(yMax);

    bool accept;

    // kin_vars and its twin are 3-vectors containing [nu-energy, xi, Bjorken Y]
    std::array<double,3> kin_vars, test_kin_vars;

    // centers of the cross section spline tales.
    std::array<int,3> spline_table_center, test_spline_table_center;

    // values of cross_section from the splines.  By * Bxi * Spline(E,xi,y)
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
        double trials = 0;
        double xi_trial = 0.0, y_trial = 0.0;
        do {
            // Per-event recoverable failure: drop this event (InjectionFailure is
            // caught by the Injector) rather than aborting the whole run.
            if (trials >= 100) throw siren::utilities::InjectionFailure("QuarkDISFromSpline: initial rejection sampling failed to find a kinematically allowed point in 100 trials");
            trials += 1;
            kin_vars[1] = random->Uniform(logXiMin, 0.0);
            kin_vars[2] = random->Uniform(logYMin, logYMax);
            xi_trial = std::pow(10., kin_vars[1]);
            y_trial  = std::pow(10., kin_vars[2]);
            trialQ = slowRescalingQ2(xi_trial, y_trial, E_nu, M_targ, mc);
        } while (trialQ < minimum_Q2_
              || !kinematicallyAllowed(xi_trial, y_trial, E_nu, M_targ, m));

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
    double measure = pow(10., kin_vars[1] + kin_vars[2]); // Bxi * By

    // Bxi * By * xs(E, xi, y)
    // evalutates the differential spline at that point
    cross_section = measure*pow(10., differential_cross_section_.ndsplineeval(kin_vars.data(), spline_table_center.data(), 0));

    // this is the magic part. Metropolis Hastings Algorithm.
    // MCMC method!
    const size_t burnin = 40; // converges to the correct distribution over multiple samplings.
    // big number means more accurate, but slower
    for(size_t j = 0; j <= burnin; j++) {
        // repeat the sampling from above to get a new valid point
        double trialQ;
        double xi_trial = 0.0, y_trial = 0.0;
        int burnin_trials = 0;
        do {
            // Per-event recoverable failure: drop this event (InjectionFailure is
            // caught by the Injector) rather than aborting the whole run.
            if (++burnin_trials >= 100)
                throw siren::utilities::InjectionFailure("QuarkDISFromSpline: burn-in proposal failed to find allowed point in 100 trials");
            test_kin_vars[1] = random->Uniform(logXiMin, 0.0);
            test_kin_vars[2] = random->Uniform(logYMin, logYMax);
            xi_trial = std::pow(10., test_kin_vars[1]);
            y_trial  = std::pow(10., test_kin_vars[2]);
            trialQ = slowRescalingQ2(xi_trial, y_trial, E_nu, M_targ, mc);
        } while (trialQ < minimum_Q2_
              || !kinematicallyAllowed(xi_trial, y_trial, E_nu, M_targ, m));

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

    // scaling down to handle numerical issues
    double final_xi = std::pow(10., kin_vars[1]);
    double final_y = pow(10., kin_vars[2]);
    record.interaction_parameters.clear();
    record.interaction_parameters["energy"] = E1_lab;
    record.interaction_parameters["bjorken_xi"] = final_xi;
    record.interaction_parameters["bjorken_y"]  = final_y;
    record.interaction_parameters["bjorken_x"]  =
        xiToBjorkenX(final_xi, final_y, E1_lab, target_mass_,
                     siren::utilities::Constants::charmMass);

    double Q2 = slowRescalingQ2(final_xi, final_y, E1_lab, target_mass_,
                                siren::utilities::Constants::charmMass);
    // Closed form for the exchanged q decomposition. A uniform momentum rescaling
    // cannot be used to recompute Q2 here: slowRescalingQ2 holds the target/charm
    // masses fixed, so Q2 is not homogeneous of degree 2 under the scaling.
    //
    // p1x_lab is the 3-momentum MAGNITUDE P1 = |p1_lab| (not an x-component).
    // The naive expressions
    //   pqx  = (m1^2 + m3^2 + 2 P1^2 + Q2 + 2 E1^2 (y-1)) / (2 P1)
    //   momq^2 = m1^2 + P1^2 + Q2 + E1^2 (y^2 - 1)
    // lose precision because both carry dominant ~E1^2 terms that nearly cancel
    // in pqy^2 = momq^2 - pqx^2. Substituting P1^2 = E1^2 - m1^2 cancels those
    // terms analytically before the numeric evaluation:
    //   pqx    = (m3^2 - m1^2 + Q2 + 2 E1^2 y) / (2 P1)
    //   momq^2 = Q2 + E1^2 y^2
    double p1x_lab = std::sqrt(p1_lab.px() * p1_lab.px() + p1_lab.py() * p1_lab.py() + p1_lab.pz() * p1_lab.pz());
    double pqx_lab = (m3 * m3 - m1 * m1 + Q2 + 2.0 * E1_lab * E1_lab * final_y) / (2.0 * p1x_lab);
    double momq2_lab = Q2 + E1_lab * E1_lab * final_y * final_y;
    double pqy2_lab = momq2_lab - pqx_lab * pqx_lab;
    // Clamp tiny-negative pqy^2 from edge roundoff to 0; genuinely-forbidden
    // points (pqy^2 < 0 in exact arithmetic) are already excluded upstream by
    // kinematicallyAllowed, so a large-magnitude negative here signals a
    // sampler/predicate inconsistency and is a hard error.
    if (pqy2_lab < 0.0) {
        if (pqy2_lab > -1e-9 * std::max(momq2_lab, 1.0)) {
            pqy2_lab = 0.0;
        } else {
            throw(siren::utilities::InjectionFailure(
                "QuarkDISFromSpline::SampleFinalState: pqy^2 < 0 for a "
                "kinematicallyAllowed point (sampler/predicate inconsistency)"));
        }
    }
    double pqy_lab = std::sqrt(pqy2_lab);
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



    // #############################################
    // New hadronization scheme: includes partonic cross section sampling and slow rescaling for charm mass effects
    // ##############################################

    const double xi = final_xi;  // sampled directly in slow-rescaling
    if (xi >= 1.0) {
        throw(siren::utilities::InjectionFailure("xi >= 1.0; sampled slow-rescaling xi past unity"));
    }
    rk::P4 p_parton(geom3::Vector3(0, 0, 0), xi * target_mass_);  // parton at rest: (xi*M, 0, 0, 0)
    rk::P4 p4_lab = p_parton + pq_lab;                              // struck charm = xi*p2 + q
    rk::P4 p_spectator((1.0 - xi) * target_mass_, geom3::Vector3(0, 0, 0));  // spectator: ((1-xi)*M, 0,0, 0)

    rk::P4 p3;
    rk::P4 p4;
    p3 = p3_lab; // now we have our lepton momentum set, which should not be modified from here on
    p4 = p4_lab; // momentum of the virtual charm

    // D meson fragmentation: D gets fraction z of charm quark momentum
    double mCH = getHadronMass(record.signature.secondary_types[meson_index]);
    double p_charm_mag = std::sqrt(p4_lab.px()*p4_lab.px() + p4_lab.py()*p4_lab.py() + p4_lab.pz()*p4_lab.pz());
    geom3::Vector3 p4_mom = p4_lab.momentum();
    geom3::UnitVector3 p4_dir = p4_mom.direction();

    rk::P4 p4CH, p4X;
    int max_sampling = 500;
    int sampling = 0;
    do {
        sampling += 1;
        if (sampling > max_sampling) {
            throw(siren::utilities::InjectionFailure("Failed to sample hadronization!"));
        }
        double z = inverseCdfTable(random->Uniform(0, 1));
        double pCH_mag = z * p_charm_mag;
        p4CH = rk::P4(geom3::Vector3(p4_dir * pCH_mag), mCH);
        p4X = p_spectator + p4_lab - p4CH;
    } while (p4X.dot(p4X) < 0);

    // Save final state kinematics.
    //
    // NOTE: the sampled momenta are written into the record's
    // SecondaryParticleRecord vector (record.GetSecondaryParticleRecords()), NOT
    // directly into an InteractionRecord's secondary_momenta. To obtain a finalized
    // InteractionRecord with populated secondary_momenta, the caller must run
    // CrossSectionDistributionRecord::Finalize (pybind: cdr.finalize(ir)) into an
    // output record whose signature is set. Building an InteractionRecord by hand
    // with empty/zero secondary_momenta and feeding it back to
    // DifferentialCrossSection makes the primary-momentum Q2 path compute Q2 <= 0,
    // forcing the stored-(xi,y) fallback branch.
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & lepton = secondaries[lepton_index];
    siren::dataclasses::SecondaryParticleRecord & hadron = secondaries[hadron_index];
    siren::dataclasses::SecondaryParticleRecord & meson = secondaries[meson_index];

    lepton.SetFourMomentum({p3.e(), p3.px(), p3.py(), p3.pz()});
    lepton.SetMass(p3.m());
    lepton.SetHelicity(record.primary_helicity);
    hadron.SetFourMomentum({p4X.e(), p4X.px(), p4X.py(), p4X.pz()});
    // Use the already-validated dot() result instead of P4::m(), which
    // recomputes msq from e^2 - p_.lengthSquared() and can disagree with
    // dot() by FP roundoff (~1e-15) -- the assert in m() would then fire
    // even though the do-while above accepted p4X.
    const double p4X_msq = p4X.dot(p4X);
    hadron.SetMass(p4X_msq > 0.0 ? std::sqrt(p4X_msq) : 0.0);
    hadron.SetHelicity(record.target_helicity);
    meson.SetFourMomentum({p4CH.e(), p4CH.px(), p4CH.py(), p4CH.pz()});


    // #############################################
    // Included for posterity: original hadronization scheme
    // ##############################################

    /*

    rk::P4 p4_lab = p2_lab + pq_lab;

    rk::P4 p3;
    rk::P4 p4;
    p3 = p3_lab; // now we have our lepton momentum set, which should not be modified from here on
    p4 = p4_lab; // momentum of the virtual charm


    // compute the energy and 3-momentum of the virtual charm
    double p3c = std::sqrt(std::pow(p4.px(), 2) + std::pow(p4.py(), 2) + std::pow(p4.pz(), 2));
    double Ec = p4.e(); //energy of primary charm
    double mCH = getHadronMass(record.signature.secondary_types[meson_index]); // obtain charmed hadron mass

    // accept-reject sampling for a valid momentum fragmentation
    bool frag_accept;
    double randValue;
    double z;
    double ECH;

    // add a maximum number of trials in the while loop
    int max_sampling = 500;
    int sampling = 0;

    // sample again if this eenrgy is not kinematically allowed
    // this samples in the lab frame the energy of the D-meson such that mass is real
    do {
        sampling += 1;
        if (sampling > max_sampling) {
            // throw(siren::utilities::InjectionFailure("Failed to sample hadronization!"));
            break;
        }
        randValue = random->Uniform(0,1);
        z = inverseCdfTable(randValue);
        ECH = z * Ec;
        if (std::pow(ECH, 2) - std::pow(mCH, 2) <= 0) {
            frag_accept = false;
        } else {
            frag_accept = true;
        }
    } while (!frag_accept);
    // new attempt of using the isoscalar mass as the remnant hadronic shower mass
    double mX = target_mass_;
    double Mc = p4.m();
    //compute the energies in the charm rest frame
    double E_CH_c = (std::pow(Mc, 2) - std::pow(mX, 2) + std::pow(mCH, 2)) / (2 * Mc);
    double p_c = std::sqrt((std::pow(Mc, 2) - std::pow(mCH + mX, 2)) * (std::pow(Mc, 2) - std::pow(mCH - mX, 2))) / (2 * Mc);
    // compute the lorentz boost parameters
    double gamma = p4.gamma();
    double beta = p4.beta();
    // using the lab frame fragmented energy and the
    double cosTheta = std::max(std::min(((ECH - gamma * E_CH_c)/(gamma * beta * p_c)), 1.), -1.);
    // now compute the momentum vectors in the rest frame
    double sinTheta = std::sin(std::acos(cosTheta));
    rk::P4 p4CH_c(p_c * geom3::Vector3(cosTheta, sinTheta, 0), mCH);
    rk::P4 p4X_c(p_c * geom3::Vector3(-cosTheta, -sinTheta, 0), mX);
    // these all assume boost direction is charm direction. Now we should rotate back to charm lab momentum direction
    geom3::Vector3 pc_lab_momentum = p4.momentum();
    geom3::UnitVector3 pc_lab_dir = pc_lab_momentum.direction();
    geom3::Rotation3 x_to_pc_lab_rot = geom3::rotationBetween(x_dir, pc_lab_dir);
    p4X_c.rotate(x_to_pc_lab_rot);
    p4CH_c.rotate(x_to_pc_lab_rot);

    // finally, we perform a random azimuthal rotation
    double c_phi = random->Uniform(0, 2 * M_PI);
    geom3::Rotation3 azimuth_rand_rot(pc_lab_dir, c_phi);
    p4X_c.rotate(azimuth_rand_rot);
    p4CH_c.rotate(azimuth_rand_rot);

    // and boost them back to the lab frame
    rk::Boost boost_from_crest_to_lab = p4.labBoost();
    rk::P4 p4X = p4X_c.boost(boost_from_crest_to_lab);
    rk::P4 p4CH = p4CH_c.boost(boost_from_crest_to_lab);



    // now we proceed to saving the final state kinematics
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & lepton = secondaries[lepton_index];
    siren::dataclasses::SecondaryParticleRecord & hadron = secondaries[hadron_index];
    siren::dataclasses::SecondaryParticleRecord & meson = secondaries[meson_index];

    lepton.SetFourMomentum({p3.e(), p3.px(), p3.py(), p3.pz()});
    lepton.SetMass(p3.m());
    lepton.SetHelicity(record.primary_helicity);
    hadron.SetFourMomentum({p4X.e(), p4X.px(), p4X.py(), p4X.pz()});
    hadron.SetMass(p4X.m());
    hadron.SetHelicity(record.target_helicity);
    meson.SetFourMomentum({p4CH.e(), p4CH.px(), p4CH.py(), p4CH.pz()});
    meson.SetMass(p4CH.m());
    meson.SetHelicity(record.target_helicity); // this needs working on
    */
}

double QuarkDISFromSpline::FragmentationFraction(siren::dataclasses::Particle::ParticleType secondary) const {
    // D0:D+/-:Ds = 0.60:0.23:0.15, renormalized to sum to 1.0. The Lambda_c
    // channel (~0.02) is not modeled, so its fraction is redistributed into the
    // implemented D species; otherwise summing TotalCrossSection over the three
    // registered D signatures would recover only 0.98 * sigma_inclusive and
    // under-count the charm rate by ~2%.
    if (secondary == siren::dataclasses::Particle::ParticleType::D0 || secondary == siren::dataclasses::Particle::ParticleType::D0Bar) {
        return 0.6 / 0.98;
    } else if (secondary == siren::dataclasses::Particle::ParticleType::DPlus || secondary == siren::dataclasses::Particle::ParticleType::DMinus) {
        return 0.23 / 0.98;
    } else if (secondary == siren::dataclasses::Particle::ParticleType::DsPlus || secondary == siren::dataclasses::Particle::ParticleType::DsMinus) {
        return 0.15 / 0.98;
    }
    return 0;
}

// UNBIASED-ONLY CONTRACT: SampleFinalState samples (xi,y) AND an independent
// fragmentation z (the inverse-CDF draw) and uniform azimuth phi that set the
// D-meson momentum. FinalStateProbability / DifferentialCrossSection account for
// (xi,y) only; the z and phi factors are NOT included here. They cancel exactly
// in the weight ratio ONLY when the same cross-section object provides both the
// injection and physical densities and no biased phase-space channel is installed
// on the D kinematics. Biasing the D kinematics is NOT supported and would produce
// incorrect weights.
//
// NORMALIZATION CONTRACT: FinalStateProbability = dxs/txs is a normalized
// kinematic density ONLY if the external 1-D total-xs spline (txs) equals the
// integral of the differential spline (dxs) over the SAME truncated slow-rescaling
// domain: xi in [xiMin(E),1], y in [yMin(E),yMax(E)] with identical charm-threshold
// (Q2 = 2 M E xi y - m_c^2 > 0), W2 > (M + M_D0)^2, Q2 >= minimum_Q2_ cuts and the
// same TARGETMASS. If the upstream total spline integrates a different domain the
// density is mis-normalized. The fragmentation fraction is applied inside
// TotalCrossSection(record) so per-species txs carries the D-species branching.
double QuarkDISFromSpline::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    // first compute the differential and total cross section
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    // fragmentation fraction is now applied inside TotalCrossSection
    if(dxs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<siren::dataclasses::ParticleType> QuarkDISFromSpline::GetPossiblePrimaries() const {
    return std::vector<siren::dataclasses::ParticleType>(primary_types_.begin(), primary_types_.end());
}

std::vector<siren::dataclasses::ParticleType> QuarkDISFromSpline::GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const {
    return std::vector<siren::dataclasses::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> QuarkDISFromSpline::GetPossibleSignatures() const {
    return std::vector<dataclasses::InteractionSignature>(signatures_.begin(), signatures_.end());
}

std::vector<siren::dataclasses::ParticleType> QuarkDISFromSpline::GetPossibleTargets() const {
    return std::vector<siren::dataclasses::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> QuarkDISFromSpline::GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const {
    std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType> key(primary_type, target_type);
    if(signatures_by_parent_types_.find(key) != signatures_by_parent_types_.end()) {
        return signatures_by_parent_types_.at(key);
    } else {
        return std::vector<dataclasses::InteractionSignature>();
    }
}

// UNBIASED-ONLY CONTRACT: the density covers only (xi,y); the independently-sampled
// fragmentation z and azimuth phi (set in SampleFinalState) are omitted. They cancel
// in the weight ratio only in the unbiased configuration (same cross-section object on
// both sides, no biased D-kinematics channel). Biasing D kinematics is NOT supported.
std::vector<std::string> QuarkDISFromSpline::DensityVariables() const {
    return std::vector<std::string>{"Bjorken xi", "Bjorken y"};
}

} // namespace interactions
} // namespace siren
