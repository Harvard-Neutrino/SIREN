#include "SIREN/interactions/PythiaDISCrossSection.h"

#include <map>
#include <set>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <cstdlib>
#include <cctype>
#include <sys/stat.h>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include <photospline/splinetable.h>

#include <Pythia8/Pythia.h>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Errors.h"

namespace siren {
namespace interactions {

// --- SIRENRndm: bridges SIREN RNG into Pythia ---

class SIRENRndm : public Pythia8::RndmEngine {
public:
    mutable std::shared_ptr<siren::utilities::SIREN_random> rng_;
    double flat() override { return rng_->Uniform(0.0, 1.0); }
};

namespace {
// Charm-DIS Pythia configuration shared by SampleFinalState (InitializePythia)
// and the spline generator, so generated total/differential splines correspond
// to exactly the events SampleFinalState produces. Route both through here to
// keep them in lockstep.
void ApplyPythiaCharmConfig(Pythia8::Pythia & pythia, int interaction_type,
                            int beam_id, int target_pdg, double E_nu,
                            std::string const & pdf_set, double minimum_Q2) {
    if (interaction_type == 1)      pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
    else if (interaction_type == 2) pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
    pythia.readString("Beams:frameType = 2");
    pythia.readString("Beams:idA = " + std::to_string(beam_id));
    pythia.readString("Beams:idB = " + std::to_string(target_pdg));
    pythia.readString("Beams:eA = " + std::to_string(E_nu));
    pythia.readString("Beams:eB = 0.");
    if (interaction_type == 1) {
        // Charm-only CC: zero non-charm CKM so every event produces charm.
        pythia.settings.forceParm("StandardModel:Vud", 0.0);
        pythia.settings.forceParm("StandardModel:Vus", 0.0);
        pythia.settings.forceParm("StandardModel:Vub", 0.0);
        pythia.settings.forceParm("StandardModel:Vcb", 0.0);
    }
    pythia.readString("PDF:pSet = " + pdf_set);
    pythia.readString("PDF:useHard = on");
    pythia.readString("PDF:pHardSet = " + pdf_set);
    pythia.readString("HadronLevel:Hadronize = on");
    pythia.readString("HadronLevel:Decay = off");
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("PhaseSpace:mHatMin = 0.0");
    pythia.readString("PhaseSpace:Q2Min = " + std::to_string(minimum_Q2));
    pythia.readString("Print:quiet = on");
}
} // anonymous namespace

// --- Constructors ---

PythiaDISCrossSection::PythiaDISCrossSection() {}

PythiaDISCrossSection::~PythiaDISCrossSection() = default;

PythiaDISCrossSection::PythiaDISCrossSection(
    std::string differential_filename,
    std::string total_filename,
    int interaction_type,
    double target_mass,
    double minimum_Q2,
    std::set<siren::dataclasses::ParticleType> primary_types,
    std::set<siren::dataclasses::ParticleType> target_types,
    std::string pythia_data_path,
    std::string pdf_set,
    std::string units)
    : primary_types_(primary_types),
      target_types_(target_types),
      interaction_type_(interaction_type),
      target_mass_(target_mass),
      minimum_Q2_(minimum_Q2),
      pdf_set_(pdf_set),
      pythia_data_path_(pythia_data_path)
{
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
    SetUnits(units);
}

PythiaDISCrossSection::PythiaDISCrossSection(
    std::string differential_filename,
    std::string total_filename,
    int interaction_type,
    double target_mass,
    double minimum_Q2,
    std::vector<siren::dataclasses::ParticleType> primary_types,
    std::vector<siren::dataclasses::ParticleType> target_types,
    std::string pythia_data_path,
    std::string pdf_set,
    std::string units)
    : primary_types_(primary_types.begin(), primary_types.end()),
      target_types_(target_types.begin(), target_types.end()),
      interaction_type_(interaction_type),
      target_mass_(target_mass),
      minimum_Q2_(minimum_Q2),
      pdf_set_(pdf_set),
      pythia_data_path_(pythia_data_path)
{
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
    SetUnits(units);
}

// --- File I/O ---

void PythiaDISCrossSection::SetUnits(std::string units) {
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

void PythiaDISCrossSection::LoadFromFile(std::string dd_crossSectionFile, std::string total_crossSectionFile) {
    // Total spline is mandatory (drives interaction depth / position / survival).
    total_cross_section_ = photospline::splinetable<>(total_crossSectionFile.c_str());
    if(total_cross_section_.get_ndim() != 1)
        throw std::runtime_error("Total cross section spline has " + std::to_string(total_cross_section_.get_ndim())
                + " dimensions, should have 1");
    // Differential spline is OPTIONAL: an empty filename means "no differential",
    // in which case FinalStateProbability returns a constant (unbiased-only).
    if(dd_crossSectionFile.empty()) {
        has_differential_ = false;
    } else {
        differential_cross_section_ = photospline::splinetable<>(dd_crossSectionFile.c_str());
        if(differential_cross_section_.get_ndim() != 3 && differential_cross_section_.get_ndim() != 2)
            throw std::runtime_error("cross section spline has " + std::to_string(differential_cross_section_.get_ndim())
                    + " dimensions, should have either 3 or 2");
        has_differential_ = true;
    }
}

void PythiaDISCrossSection::LoadFromMemory(std::vector<char> & differential_data, std::vector<char> & total_data) {
    total_cross_section_.read_fits_mem(total_data.data(), total_data.size());
    if(!differential_data.empty()) {
        differential_cross_section_.read_fits_mem(differential_data.data(), differential_data.size());
        has_differential_ = true;
    } else {
        has_differential_ = false;
    }
}

void PythiaDISCrossSection::ReadParamsFromSplineTable() {
    bool mass_good = differential_cross_section_.read_key("TARGETMASS", target_mass_);
    bool int_good = differential_cross_section_.read_key("INTERACTION", interaction_type_);
    bool q2_good = differential_cross_section_.read_key("Q2MIN", minimum_Q2_);
    if(!int_good) interaction_type_ = 1;
    if(!q2_good) minimum_Q2_ = 1;
    if(!mass_good) target_mass_ = siren::utilities::Constants::isoscalarMass;
}

// --- Equality ---

bool PythiaDISCrossSection::equal(CrossSection const & other) const {
    const PythiaDISCrossSection* x = dynamic_cast<const PythiaDISCrossSection*>(&other);
    if(!x) return false;
    return std::tie(interaction_type_, target_mass_, minimum_Q2_, signatures_, primary_types_, target_types_)
        == std::tie(x->interaction_type_, x->target_mass_, x->minimum_Q2_, x->signatures_, x->primary_types_, x->target_types_);
}

// --- Particle ID helpers ---

bool PythiaDISCrossSection::IsCharmedHadron(int pdgId) {
    int abs_id = std::abs(pdgId);
    // Match D0 (421), D+/- (411), Ds (431). Lambda_c (4122) still routed to hadronic remnant.
    return (abs_id == 411 || abs_id == 421 || abs_id == 431);
}

siren::dataclasses::ParticleType PythiaDISCrossSection::PdgToParticleType(int pdgId) {
    // Direct cast -- SIREN ParticleType enum values match PDG codes
    return static_cast<siren::dataclasses::ParticleType>(pdgId);
}

double PythiaDISCrossSection::GetLeptonMass(siren::dataclasses::ParticleType lepton_type) {
    int32_t lepton_number = std::abs(static_cast<int32_t>(lepton_type));
    switch(lepton_number) {
        case 11: return siren::utilities::Constants::electronMass;
        case 13: return siren::utilities::Constants::muonMass;
        case 15: return siren::utilities::Constants::tauMass;
        case 12: case 14: case 16: return 0;
        default: throw std::runtime_error("Unknown lepton type!");
    }
}

std::map<std::string, int> PythiaDISCrossSection::getIndices(siren::dataclasses::InteractionSignature signature) {
    // Identify meson by elimination (not via isD(), which only covers D0/D+/-).
    // First pass: claim lepton and Hadrons. Second pass: whatever remains is the meson.
    int lepton_id = -1, hadron_id = -1, meson_id = -1;
    for (size_t i = 0; i < signature.secondary_types.size(); i++) {
        if (siren::dataclasses::isLepton(signature.secondary_types[i])) {
            lepton_id = i;
        } else if (signature.secondary_types[i] == siren::dataclasses::ParticleType::Hadrons) {
            hadron_id = i;
        }
    }
    for (size_t i = 0; i < signature.secondary_types.size(); i++) {
        if ((int)i == lepton_id || (int)i == hadron_id) continue;
        meson_id = i;
        break;
    }
    return {{"lepton", lepton_id}, {"hadron", hadron_id}, {"meson", meson_id}};
}

// --- Signatures ---

void PythiaDISCrossSection::InitializeSignatures() {
    // PythiaDISCrossSection forces charm only in charged-current (the non-charm
    // CKM elements are zeroed in ApplyPythiaCharmConfig). Z exchange has no such
    // handle, so NC charm cannot be forced and SampleFinalState would exhaust its
    // attempt budget and throw. Reject NC at construction with an actionable
    // pointer to the tool that does support NC charm (QuarkDISFromSpline).
    if(interaction_type_ != 1) {
        if(interaction_type_ == 2)
            throw std::runtime_error("PythiaDISCrossSection supports only charged-current charm production (interaction_type=1). Neutral-current charm is not forced by Z exchange in Pythia, so per-event sampling would fail; use QuarkDISFromSpline for NC charm.");
        throw std::runtime_error("PythiaDISCrossSection: unsupported interaction_type=" + std::to_string(interaction_type_) + " (expected 1 for charged current).");
    }
    signatures_.clear();
    for(auto primary_type : primary_types_) {
        dataclasses::InteractionSignature signature;
        signature.primary_type = primary_type;
        if(not siren::dataclasses::isNeutrino(primary_type)) {
            throw std::runtime_error("PythiaDISCrossSection only supports neutrinos as primaries!");
        }

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
            throw std::runtime_error("InitializeSignatures: Unknown parent neutrino type!");
        }

        if(interaction_type_ == 1) {
            signature.secondary_types.push_back(charged_lepton_product);
        } else if(interaction_type_ == 2) {
            signature.secondary_types.push_back(neutral_lepton_product);
        } else {
            throw std::runtime_error("InitializeSignatures: Unknown interaction type!");
        }

        // Hadron remnant
        signature.secondary_types.push_back(siren::dataclasses::ParticleType::Hadrons);

        // Charmed meson types. For nu the c quark fragments to D0/D+/Ds+; for nubar
        // the cbar quark fragments to Dbar0/D-/Ds-. SampleFinalState writes Pythia's
        // actual produced PID into the signature's meson slot, so the registered set
        // must include the correct charge to keep weighter signature lookups in range
        // (otherwise event_weight would be NaN).
        // TODO: Add Lambda_c (4122) support.
        bool is_antineutrino =
            (primary_type == siren::dataclasses::ParticleType::NuEBar ||
             primary_type == siren::dataclasses::ParticleType::NuMuBar ||
             primary_type == siren::dataclasses::ParticleType::NuTauBar);
        std::set<siren::dataclasses::ParticleType> D_types_local;
        if (is_antineutrino) {
            D_types_local = {siren::dataclasses::ParticleType::D0Bar,
                        siren::dataclasses::ParticleType::DMinus,
                        siren::dataclasses::ParticleType::DsMinus};
        } else {
            D_types_local = {siren::dataclasses::ParticleType::D0,
                        siren::dataclasses::ParticleType::DPlus,
                        siren::dataclasses::ParticleType::DsPlus};
        }

        for (auto meson_type : D_types_local) {
            dataclasses::InteractionSignature full_signature = signature;
            full_signature.secondary_types.push_back(meson_type);
            for(auto target_type : target_types_) {
                full_signature.target_type = target_type;
                signatures_.push_back(full_signature);
                std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType> key(primary_type, target_type);
                signatures_by_parent_types_[key].push_back(full_signature);
            }
        }
    }
}

// --- Cross sections (from splines, same as QuarkDISFromSpline) ---

double PythiaDISCrossSection::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    siren::dataclasses::ParticleType primary_type = interaction.signature.primary_type;
    double primary_energy = interaction.primary_momentum[0];
    if(primary_energy < InteractionThreshold(interaction)) return 0;
    double total_xs = TotalCrossSection(primary_type, primary_energy);
    // Partition the inclusive charm cross section across D species by fragmentation
    // fraction for the specific D meson in this signature. The total spline holds the
    // single inclusive charm cross section, so without this each of the registered
    // D-type signatures (D0 + D+ + Ds) would return the full value and summing over
    // them -- as the base-class TotalCrossSectionAllFinalStates and the Weighter do --
    // would triple-count charm production. Mirrors QuarkDISFromSpline::TotalCrossSection.
    for(auto const & sec_type : interaction.signature.secondary_types) {
        if(siren::dataclasses::isD(sec_type)) {
            total_xs *= FragmentationFraction(sec_type);
            break;
        }
    }
    return total_xs;
}

double PythiaDISCrossSection::TotalCrossSection(siren::dataclasses::ParticleType primary_type, double primary_energy) const {
    if(not primary_types_.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
    double log_energy = log10(primary_energy);
    if(log_energy < total_cross_section_.lower_extent(0)
            or log_energy > total_cross_section_.upper_extent(0)) {
        throw std::runtime_error("PythiaDISCrossSection: primary energy "
            + std::to_string(primary_energy) + " GeV is outside the total cross-section spline range ["
            + std::to_string(std::pow(10.0, total_cross_section_.lower_extent(0))) + ", "
            + std::to_string(std::pow(10.0, total_cross_section_.upper_extent(0)))
            + "] GeV; regenerate the total spline (siren.pythia_charm_splines.generate_total_spline) to cover this energy.");
    }
    int center;
    total_cross_section_.searchcenters(&log_energy, &center);
    double log_xs = total_cross_section_.ndsplineeval(&log_energy, &center, 0);
    return unit * std::pow(10.0, log_xs);
}


double PythiaDISCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    if(!has_differential_)
        throw std::runtime_error("PythiaDISCrossSection::DifferentialCrossSection called but no differential spline was loaded (the unbiased default uses a constant FinalStateProbability).");
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), interaction.target_mass);
    double primary_energy = interaction.primary_momentum[0];

    std::map<std::string, int> secondaries = getIndices(interaction.signature);
    int lepton_index_i = secondaries["lepton"];
    // Guard against a malformed/empty signature: getIndices returns -1 when no
    // lepton is found, and using that as an unsigned index reads out of bounds.
    if(lepton_index_i < 0 || static_cast<size_t>(lepton_index_i) >= interaction.secondary_momenta.size()) {
        return 0.0;
    }
    unsigned int lepton_index = static_cast<unsigned int>(lepton_index_i);

    std::array<double, 4> const & mom3 = interaction.secondary_momenta[lepton_index];
    rk::P4 p3(geom3::Vector3(mom3[1], mom3[2], mom3[3]), interaction.secondary_masses[lepton_index]);
    rk::P4 q = p1 - p3;

    double Q2 = -q.dot(q);
    double lepton_mass = GetLeptonMass(interaction.signature.secondary_types[lepton_index]);
    double y = 1.0 - p2.dot(p3) / p2.dot(p1);
    double p2q = p2.dot(q);
    double x = (p2q != 0.0) ? Q2 / (2.0 * p2q) : -1.0;   // muon-reconstructed Bjorken x

    // Use the reconstructed (x,y) when they fall inside the spline domain; the
    // spline's support is the validity region (no separate analytic charm-DIS
    // kinematic gate, which assumes a model the Pythia-derived spline does not).
    double log_energy = log10(primary_energy);
    std::array<int,3> centers;
    bool ok = (x > 0.0 && x < 1.0 && y > 0.0 && y < 1.0 && Q2 >= minimum_Q2_);
    if(ok) {
        std::array<double,3> coordinates{{log_energy, log10(x), log10(y)}};
        ok = differential_cross_section_.searchcenters(coordinates.data(), centers.data());
    }
    if(!ok) {
        // Fall back to the stored (x,y) -- the SAME muon-reconstructed Bjorken
        // variables SampleFinalState records. Non-throwing: if absent (e.g. a
        // fake record built for a rate query), the differential is undefined -> 0.
        auto itx = interaction.interaction_parameters.find("bjorken_x");
        auto ity = interaction.interaction_parameters.find("bjorken_y");
        if(itx == interaction.interaction_parameters.end() || ity == interaction.interaction_parameters.end())
            return 0.0;
        x = itx->second;
        y = ity->second;
        Q2 = 2.0 * primary_energy * p2.e() * x * y;
    }
    return DifferentialCrossSection(primary_energy, x, y, lepton_mass, Q2);
}

double PythiaDISCrossSection::DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass, double Q2) const {
    if(!has_differential_)
        throw std::runtime_error("PythiaDISCrossSection::DifferentialCrossSection called but no differential spline was loaded (the unbiased default uses a constant FinalStateProbability).");
    double log_energy = log10(energy);
    // Out of spline SUPPORT -> raise, never silently return 0: a silent zero on a
    // genuinely sampled event would bias that event's physical density (and hence
    // its weight) to zero. The energy extent and the (x, y) grid define validity.
    if(log_energy < differential_cross_section_.lower_extent(0)
            || log_energy > differential_cross_section_.upper_extent(0))
        throw std::runtime_error("PythiaDISCrossSection: energy " + std::to_string(energy)
            + " GeV is outside the differential spline energy range ["
            + std::to_string(std::pow(10.0, differential_cross_section_.lower_extent(0))) + ", "
            + std::to_string(std::pow(10.0, differential_cross_section_.upper_extent(0)))
            + "] GeV; regenerate the differential spline to cover this energy.");
    if(x <= 0 || x >= 1 || y <= 0 || y >= 1)
        throw std::runtime_error("PythiaDISCrossSection: unphysical Bjorken (x="
            + std::to_string(x) + ", y=" + std::to_string(y) + ") outside (0, 1).");
    if(std::isnan(Q2)) Q2 = 2.0 * energy * target_mass_ * x * y;
    if(Q2 < minimum_Q2_) return 0;   // below the configured Q^2 cut: a modeling threshold, not a spline-range miss
    (void)secondary_lepton_mass;     // analytic charm-DIS gate removed (see above)

    std::array<double,3> coordinates{{log_energy, log10(x), log10(y)}};
    std::array<int,3> centers;
    if(!differential_cross_section_.searchcenters(coordinates.data(), centers.data()))
        throw std::runtime_error("PythiaDISCrossSection: Bjorken (x=" + std::to_string(x)
            + ", y=" + std::to_string(y) + ") at E=" + std::to_string(energy)
            + " GeV is outside the differential spline (x, y) grid; widen the spline's logx/logy range.");
    double result = pow(10., differential_cross_section_.ndsplineeval(coordinates.data(), centers.data(), 0));
    return unit * result;
}

double PythiaDISCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    return 0;
}

// --- Fragmentation fractions ---

double PythiaDISCrossSection::FragmentationFraction(siren::dataclasses::Particle::ParticleType secondary) const {
    // Approximate fractions from Pythia (charm hadronization), renormalized to
    // sum to 1.0 over the implemented D species. Raw fractions D0:D+/-:Ds =
    // 0.60:0.23:0.15 sum to 0.98 because the Lambda_c channel is not modeled; the
    // unmodeled Lambda_c fraction is redistributed by dividing each by 0.98 so the
    // partitioned signatures exactly recover the inclusive charm cross section.
    // Values kept in lockstep with QuarkDISFromSpline::FragmentationFraction.
    if (secondary == siren::dataclasses::ParticleType::D0 || secondary == siren::dataclasses::ParticleType::D0Bar) {
        return 0.6 / 0.98;
    } else if (secondary == siren::dataclasses::ParticleType::DPlus || secondary == siren::dataclasses::ParticleType::DMinus) {
        return 0.23 / 0.98;
    } else if (secondary == siren::dataclasses::ParticleType::DsPlus || secondary == siren::dataclasses::ParticleType::DsMinus) {
        return 0.15 / 0.98;
    }
    // Lambda_c (~0.09) not yet implemented; its fraction is folded into the above.
    return 0;
}

double PythiaDISCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    // The final state is sampled by Pythia, whose per-event density is not
    // analytically available. With NO differential spline we return a constant:
    // in the standard unbiased configuration the same cross-section object
    // supplies both the injection and physical densities, so this factor cancels
    // in the weight ratio and only TotalCrossSection (interaction depth /
    // position / survival) matters. Biasing the final-state kinematics is not
    // supported in this mode.
    if(!has_differential_) return 1.0;

    // A Pythia-derived differential spline was supplied: report the true density
    // dsigma/sigma so weights remain correct under reweighting. The fragmentation
    // fraction in TotalCrossSection cancels per-signature in CrossSectionProbability.
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if (!std::isfinite(dxs) || !std::isfinite(txs) || dxs <= 0 || txs <= 0) return 0.0;
    double result = dxs / txs;
    if (!std::isfinite(result)) return 0.0;
    return result;
}

// --- Signature accessors ---

std::vector<siren::dataclasses::ParticleType> PythiaDISCrossSection::GetPossiblePrimaries() const {
    return std::vector<siren::dataclasses::ParticleType>(primary_types_.begin(), primary_types_.end());
}

std::vector<siren::dataclasses::ParticleType> PythiaDISCrossSection::GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const {
    return std::vector<siren::dataclasses::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> PythiaDISCrossSection::GetPossibleSignatures() const {
    return std::vector<dataclasses::InteractionSignature>(signatures_.begin(), signatures_.end());
}

std::vector<siren::dataclasses::ParticleType> PythiaDISCrossSection::GetPossibleTargets() const {
    return std::vector<siren::dataclasses::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> PythiaDISCrossSection::GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const {
    std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType> key(primary_type, target_type);
    if(signatures_by_parent_types_.find(key) != signatures_by_parent_types_.end()) {
        return signatures_by_parent_types_.at(key);
    }
    return std::vector<dataclasses::InteractionSignature>();
}

std::vector<std::string> PythiaDISCrossSection::DensityVariables() const {
    return std::vector<std::string>{"Bjorken x", "Bjorken y"};
}

// ======================================================================
// Pythia initialization and SampleFinalState
// ======================================================================

void PythiaDISCrossSection::InitializePythia(double E_nu, int target_pdg) const {
    // Ensure LHAPDF can find PDF sets -- derive data path from the LHAPDF library location
    // The pdf_set_ is e.g. "LHAPDF6:HERAPDF20_NLO_EIG", and LHAPDF needs LHAPDF_DATA_PATH set
    const char* lhapdf_path = std::getenv("LHAPDF_DATA_PATH");
    if (!lhapdf_path) {
        throw std::runtime_error("LHAPDF_DATA_PATH is not set; set it to your LHAPDF data directory (the parent of the PDF set folders).");
    }

    // Pre-flight: verify the requested PDF set is actually installed under
    // LHAPDF_DATA_PATH. A missing set otherwise surfaces only as Pythia's opaque
    // "initialization failed" much later; name the set and where we looked.
    {
        std::string set_name = pdf_set_;
        auto colon = set_name.find(':');                 // strip an "LHAPDF6:" prefix
        if (colon != std::string::npos) set_name = set_name.substr(colon + 1);
        std::string set_dir = std::string(lhapdf_path) + "/" + set_name;
        struct stat st;
        if (stat(set_dir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) {
            throw std::runtime_error("PythiaDISCrossSection: LHAPDF PDF set '" + set_name
                + "' not found under LHAPDF_DATA_PATH=" + std::string(lhapdf_path)
                + " (looked for " + set_dir + "). Install it (e.g. 'lhapdf install " + set_name
                + "') or construct with a pdf_set that is installed.");
        }
    }

    pythia_ = std::make_unique<Pythia8::Pythia>(pythia_data_path_, false);

    // Determine beam particle from primary types
    int beam_id = 14; // default nu_mu
    for (auto pt : primary_types_) {
        beam_id = static_cast<int>(pt);
        break; // use the first primary type
    }

    ApplyPythiaCharmConfig(*pythia_, interaction_type_, beam_id, target_pdg, E_nu, pdf_set_, minimum_Q2_);

    if (!pythia_->init()) {
        throw std::runtime_error("PythiaDISCrossSection: Pythia initialization failed");
    }

    // Bridge SIREN RNG into Pythia for reproducibility.
    // Must be done after init() -- init uses Pythia's internal RNG for setup.
    // The SIREN RNG is connected here and updated per-event in SampleFinalState.
    siren_rndm_ = std::make_shared<SIRENRndm>();
    pythia_->rndm.rndmEnginePtr(siren_rndm_);
}

void PythiaDISCrossSection::GeneratePythiaCharmSamples(
    int interaction_type, int primary_pdg, int target_pdg, double target_mass,
    std::string pdf_set, std::string pythia_data_path, double minimum_Q2,
    std::vector<double> const & energies, int n_events,
    std::vector<double> & out_sigma_mb,
    std::vector<double> & out_E, std::vector<double> & out_x, std::vector<double> & out_y)
{
    // Run Pythia once per energy (fast: ~ms/event after init) using the SAME
    // configuration SampleFinalState uses, recording the muon-reconstructed
    // Bjorken (x, y) per charm event and Pythia's generated cross section per
    // energy. The python generate_*_spline helpers fit these into FITS splines.
    out_sigma_mb.clear(); out_E.clear(); out_x.clear(); out_y.clear();
    if (!std::getenv("LHAPDF_DATA_PATH"))
        throw std::runtime_error("GeneratePythiaCharmSamples: LHAPDF_DATA_PATH is not set");
    bool is_cc = (interaction_type == 1);
    for (double E : energies) {
        Pythia8::Pythia pythia(pythia_data_path, false);
        ApplyPythiaCharmConfig(pythia, interaction_type, primary_pdg, target_pdg, E, pdf_set, minimum_Q2);
        if (!pythia.init())
            throw std::runtime_error("GeneratePythiaCharmSamples: Pythia init failed");
        for (int i = 0; i < n_events; ++i) {
            if (!pythia.next()) continue;
            int i_lep = -1;
            for (int j = 0; j < pythia.event.size(); ++j) {
                if (!pythia.event[j].isFinal()) continue;
                int ap = std::abs(pythia.event[j].id());
                bool is_lep = is_cc ? (ap == 11 || ap == 13 || ap == 15)
                                    : (ap == 12 || ap == 14 || ap == 16);
                if (is_lep) { i_lep = j; break; }
            }
            if (i_lep < 0) continue;
            Pythia8::Vec4 pl = pythia.event[i_lep].p();
            double y = 1.0 - pl.e() / E;
            if (y <= 0.0 || y >= 1.0) continue;
            Pythia8::Vec4 pnu(0.0, 0.0, E, E);   // (px,py,pz,e), nu along +z
            Pythia8::Vec4 q4 = pnu - pl;
            double Q2 = -q4.m2Calc();
            if (Q2 <= 0.0) continue;
            double x = Q2 / (2.0 * target_mass * E * y);
            if (x <= 0.0 || x >= 1.0) continue;
            out_E.push_back(E); out_x.push_back(x); out_y.push_back(y);
        }
        out_sigma_mb.push_back(pythia.info.sigmaGen());
    }
}

void PythiaDISCrossSection::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // Get indices for lepton, hadron, meson in the secondary list
    std::map<std::string, int> secondary_indices = getIndices(record.signature);
    unsigned int lepton_index = secondary_indices["lepton"];
    unsigned int hadron_index = secondary_indices["hadron"];
    unsigned int meson_index = secondary_indices["meson"];

    // Get primary neutrino 4-momentum
    rk::P4 p1(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    double E_nu = p1.e();
    geom3::UnitVector3 nu_dir = p1.momentum().direction();

    // Sample target nucleon PDG: 10/18 proton (2212), 8/18 neutron (2112)
    // to match H2O composition of ice (10 protons, 8 neutrons per molecule).
    // Pythia requires a concrete hadron beam; SIREN's "Nucleon" abstraction
    // is resolved here, per event, via the SIREN RNG.
    int target_pdg = (random->Uniform(0.0, 1.0) < (10.0 / 18.0)) ? 2212 : 2112;

    // Re-initialize Pythia every event so Beams:eA tracks this event's E_nu.
    // Pythia 8's variable-energy mode is not supported for WeakBosonExchange
    // processes, so setKinematics cannot be used here. Rebuilding the Pythia
    // object is expensive (~1 s/event) but is the only way to keep the
    // sampled muon kinematics consistent with the event's true beam energy.
    InitializePythia(E_nu, target_pdg);

    // Bridge the SIREN RNG for this event
    siren_rndm_->rng_ = random;

    // The outgoing primary lepton is fixed by the signature: the charged lepton
    // for CC (e/mu/tau, signed), or the same-flavor neutrino for NC. Match Pythia's
    // final state against that exact PDG instead of hardcoding the muon, so NC and
    // nu_e / nu_tau work, not only nu_mu CC.
    int expected_lepton_pdg = static_cast<int>(record.signature.secondary_types[lepton_index]);

    // Generate a Pythia event
    int max_attempts = 100;
    bool found_charm = false;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        if (!pythia_->next()) {
            continue;
        }

        // Find the outgoing lepton and charmed hadron in the final state
        int i_lep = -1;
        int i_charm = -1;
        Pythia8::Vec4 p_remnant(0., 0., 0., 0.);

        for (int i = 0; i < pythia_->event.size(); ++i) {
            if (!pythia_->event[i].isFinal()) continue;

            int pid = pythia_->event[i].id();

            if (pid == expected_lepton_pdg && i_lep < 0) {
                // Primary outgoing lepton matching the signature (CC charged
                // lepton or NC neutrino, with the correct charge).
                i_lep = i;
            } else if (IsCharmedHadron(pid) && i_charm < 0) {
                // First charmed hadron
                i_charm = i;
            } else {
                // Everything else goes into the remnant
                p_remnant += pythia_->event[i].p();
            }
        }

        if (i_lep >= 0 && i_charm >= 0) {
            // Trust Pythia: accept whatever charm meson it produced and
            // overwrite the signature's pre-chosen (uniform) meson slot with
            // Pythia's actual PID. Natural Lund-string fragmentation then
            // carries the physical D-type distribution without rejection.
            int pythia_pid = pythia_->event[i_charm].id();
            const_cast<dataclasses::InteractionSignature &>(record.signature)
                .secondary_types[meson_index] = PdgToParticleType(pythia_pid);

            found_charm = true;

            // Extract 4-momenta from Pythia (in Pythia's frame: nu along +z)
            Pythia8::Vec4 p_lep_pythia = pythia_->event[i_lep].p();
            Pythia8::Vec4 p_D_pythia = pythia_->event[i_charm].p();

            // DIS kinematics for weighting, reconstructed from the outgoing lepton
            // exactly as DifferentialCrossSection does (muon-reconstructed Bjorken
            // x = Q2/(2 M E y)), so an optional differential spline fit in this
            // convention closes against the sampler. NOT info.x2() (the parton
            // momentum fraction), which is a different variable.
            double E_lep = p_lep_pythia.e();
            double pythia_y = 1.0 - E_lep / E_nu;
            Pythia8::Vec4 p_nu_pythia(0.0, 0.0, E_nu, E_nu);   // (px,py,pz,e), nu along +z
            Pythia8::Vec4 q4 = p_nu_pythia - p_lep_pythia;
            double Q2_lep = -q4.m2Calc();
            double pythia_x = (pythia_y > 0.0)
                ? Q2_lep / (2.0 * target_mass_ * E_nu * pythia_y)
                : 0.0;

            // Store in interaction parameters for weighting
            record.interaction_parameters.clear();
            record.interaction_parameters["energy"] = E_nu;
            record.interaction_parameters["bjorken_x"] = pythia_x;
            record.interaction_parameters["bjorken_y"] = pythia_y;
            record.interaction_parameters["target_pdg"] = static_cast<double>(target_pdg);

            // Rotate from Pythia frame (+z) to SIREN's neutrino direction
            geom3::UnitVector3 z_dir = geom3::UnitVector3::zAxis();
            geom3::Rotation3 rot = geom3::rotationBetween(z_dir, nu_dir);

            // Helper lambda to convert Pythia Vec4 -> rotated rk::P4
            auto pythia_to_siren = [&](const Pythia8::Vec4 & pv, double mass) -> rk::P4 {
                geom3::Vector3 mom(pv.px(), pv.py(), pv.pz());
                mom = rot * mom;
                return rk::P4(mom, mass);
            };

            double lepton_mass = pythia_->event[i_lep].m();
            double D_mass = pythia_->event[i_charm].m();
            double remnant_mass = std::sqrt(std::max(0.0,
                p_remnant.e() * p_remnant.e() -
                p_remnant.px() * p_remnant.px() -
                p_remnant.py() * p_remnant.py() -
                p_remnant.pz() * p_remnant.pz()));

            rk::P4 p3 = pythia_to_siren(p_lep_pythia, lepton_mass);
            rk::P4 p_D = pythia_to_siren(p_D_pythia, D_mass);

            // Remnant
            geom3::Vector3 rem_mom(p_remnant.px(), p_remnant.py(), p_remnant.pz());
            rem_mom = rot * rem_mom;
            rk::P4 p_rem(rem_mom, remnant_mass);

            // Set secondaries
            std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
            siren::dataclasses::SecondaryParticleRecord & lepton = secondaries[lepton_index];
            siren::dataclasses::SecondaryParticleRecord & hadron = secondaries[hadron_index];
            siren::dataclasses::SecondaryParticleRecord & meson = secondaries[meson_index];

            lepton.SetFourMomentum({p3.e(), p3.px(), p3.py(), p3.pz()});
            lepton.SetMass(lepton_mass);
            lepton.SetHelicity(record.primary_helicity);

            hadron.SetFourMomentum({p_rem.e(), p_rem.px(), p_rem.py(), p_rem.pz()});
            hadron.SetMass(remnant_mass);
            hadron.SetHelicity(record.target_helicity);

            meson.SetFourMomentum({p_D.e(), p_D.px(), p_D.py(), p_D.pz()});
            meson.SetMass(D_mass);
            meson.SetHelicity(record.target_helicity);

            break;
        }
    }

    if (!found_charm) {
        // Recoverable per-event failure: throw InjectionFailure so Injector::
        // GenerateEvent catches it and drops/retries the event instead of a plain
        // std::runtime_error, which would abort the entire generation run.
        throw siren::utilities::InjectionFailure("PythiaDISCrossSection::SampleFinalState: Failed to generate charm event after " + std::to_string(max_attempts) + " attempts");
    }
}

} // namespace interactions
} // namespace siren
