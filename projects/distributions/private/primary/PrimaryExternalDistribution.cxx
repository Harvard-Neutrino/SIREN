#include "SIREN/distributions/primary/PrimaryExternalDistribution.h"

#include <algorithm>                                       // for min
#include <array>                                           // for array
#include <cmath>                                           // for isfinite, llround
#include <fstream>                                         // for ifstream
#include <sstream>                                         // for stringstream
#include <string>                                          // for basic_string
#include <stdexcept>                                       // for runtime_error
#include <set>                                             // for set
#include <tuple>                                           // for tie

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/DistributionVariable.h"
#include "SIREN/utilities/Errors.h"               // for WeightCalculationError
#include "SIREN/utilities/Random.h"               // for SIREN_random
#include "SIREN/math/Vector3D.h"

namespace siren {
namespace distributions {

//---------------
// class PrimaryExternalDistribution : PrimaryExternalDistribution
//---------------

static std::string trim(std::string const & s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

void PrimaryExternalDistribution::LoadInputFile(std::string const & _filename) {
    filename = _filename;
    keys.clear();
    input_data.clear();
    init_pos_set = false;
    vertex_set = false;
    mom_set = false;

    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        throw std::runtime_error("error: file open failed " + filename);
    }

    std::string line;
    std::getline(input_file, line);

    std::stringstream ss(line);
    std::string key;
    bool has_x0 = false, has_y0 = false, has_z0 = false;
    bool has_x = false, has_y = false, has_z = false;
    bool has_px = false, has_py = false, has_pz = false;
    while (std::getline(ss, key, ',')) {
        key = trim(key);
        keys.push_back(key);
        if (key == "x0") has_x0 = true;
        else if (key == "y0") has_y0 = true;
        else if (key == "z0") has_z0 = true;
        else if (key == "x") has_x = true;
        else if (key == "y") has_y = true;
        else if (key == "z") has_z = true;
        else if (key == "px") has_px = true;
        else if (key == "py") has_py = true;
        else if (key == "pz") has_pz = true;
    }
    init_pos_set = has_x0 && has_y0 && has_z0;
    vertex_set = has_x && has_y && has_z;
    mom_set = has_px && has_py && has_pz;

    if((has_x0 || has_y0 || has_z0) && !init_pos_set) {
        throw std::runtime_error("Incomplete initial position header in " + filename + ": need all of x0, y0, z0");
    }
    if((has_x || has_y || has_z) && !vertex_set) {
        throw std::runtime_error("Incomplete vertex header in " + filename + ": need all of x, y, z");
    }
    if((has_px || has_py || has_pz) && !mom_set) {
        throw std::runtime_error("Incomplete momentum header in " + filename + ": need all of px, py, pz");
    }

    std::string value;
    while (std::getline(input_file, line)) {
        std::string trimmed = trim(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        std::vector<double> tmp_data;
        std::stringstream _ss(line);
        size_t ikey = 0;
        bool passed = true;
        while (std::getline(_ss, value, ',')) {
            if (ikey >= keys.size()) {
                throw std::runtime_error("CSV row has more columns than header in " + filename);
            }
            if (keys[ikey] == "E") {
                if (stod(value) < emin) passed = false;
            }
            ++ikey;
            tmp_data.push_back(stod(value));
        }
        if (ikey != keys.size()) {
            throw std::runtime_error("CSV row has fewer columns than header in " + filename);
        }
        if (passed) input_data.push_back(tmp_data);
    }

    if (input_data.empty()) {
        throw std::runtime_error("No valid data rows in " + filename);
    }

    ComputeSetVariables();
    DerivePhysicalRowWeights();
}

void PrimaryExternalDistribution::ComputeSetVariables() {
    // Compute set_variables_ from the column headers
    set_variables_.clear();
    for (auto const & k : keys) {
        if (k == "E") set_variables_.insert(DistributionVariable::PrimaryEnergy);
        else if (k == "m") set_variables_.insert(DistributionVariable::PrimaryMass);
        else if (k == "px" || k == "py" || k == "pz") {
            set_variables_.insert(DistributionVariable::PrimaryDirection);
            set_variables_.insert(DistributionVariable::PrimaryEnergy);
        }
        else if (k == "x" || k == "y" || k == "z") set_variables_.insert(DistributionVariable::InteractionVertex);
        else if (k == "x0" || k == "y0" || k == "z0") set_variables_.insert(DistributionVariable::InitialPosition);
        else set_variables_.insert(DistributionVariable::InteractionParameters);
    }
}

// A column named "weight" declares per-row physical weights (see the class
// docstring). They never influence how rows are SAMPLED -- the supplied list
// is the intended injection ensemble, drawn uniformly or by an explicit bias.
// The weights only encode the return to the physical ensemble: the per-row
// physical density is served by PhysicalDensity and the weight total becomes
// the distribution's physical normalization. Must be re-run whenever
// input_data changes (construction, emin filtering, deserialization).
void PrimaryExternalDistribution::DerivePhysicalRowWeights() {
    physical_weights_.clear();
    physical_weights_sum_ = 0;
    auto weight_it = std::find(keys.begin(), keys.end(), "weight");
    if (weight_it == keys.end()) return;
    size_t weight_index = static_cast<size_t>(
        std::distance(keys.begin(), weight_it));
    physical_weights_.reserve(input_data.size());
    for (auto const & row : input_data) {
        double w = row[weight_index];
        if (!(w >= 0) || !std::isfinite(w)) {
            throw std::runtime_error(
                "PrimaryExternalDistribution: the \"weight\" column declares "
                "physical row weights, which must be finite and non-negative");
        }
        physical_weights_.push_back(w);
        physical_weights_sum_ += w;
    }
    if (!(physical_weights_sum_ > 0)) {
        throw std::runtime_error(
            "PrimaryExternalDistribution: physical row weights (the \"weight\" "
            "column) must have a positive sum");
    }
    SetNormalization(physical_weights_sum_);
}

void PrimaryExternalDistribution::BuildSamplingCDF() {
    if (sampling_weights_.empty()) {
        sampling_weights_sum_ = 0;
        sampling_cdf_.clear();
        return;
    }
    sampling_weights_sum_ = 0;
    for (double w : sampling_weights_) {
        if (w < 0) throw std::runtime_error("Sampling weights must be non-negative");
        sampling_weights_sum_ += w;
    }
    if (sampling_weights_sum_ <= 0) {
        throw std::runtime_error("Sum of sampling weights must be positive");
    }
    sampling_cdf_.resize(sampling_weights_.size());
    double cumsum = 0;
    for (size_t i = 0; i < sampling_weights_.size(); ++i) {
        cumsum += sampling_weights_[i];
        sampling_cdf_[i] = cumsum / sampling_weights_sum_;
    }
}

PrimaryExternalDistribution::PrimaryExternalDistribution(std::string _filename) : emin(0)
{
    LoadInputFile(_filename);
}

PrimaryExternalDistribution::PrimaryExternalDistribution(std::string _filename, double emin) : emin(emin)
{
    LoadInputFile(_filename);
}

PrimaryExternalDistribution::PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data)
    : emin(0)
{
    keys = std::move(_keys);
    input_data = std::move(_data);
    filename = "<in-memory>";

    // Require all three components for each flag
    bool has_x0 = false, has_y0 = false, has_z0 = false;
    bool has_x = false, has_y = false, has_z = false;
    bool has_px = false, has_py = false, has_pz = false;
    for (auto const & k : keys) {
        if (k == "x0") has_x0 = true;
        else if (k == "y0") has_y0 = true;
        else if (k == "z0") has_z0 = true;
        else if (k == "x") has_x = true;
        else if (k == "y") has_y = true;
        else if (k == "z") has_z = true;
        else if (k == "px") has_px = true;
        else if (k == "py") has_py = true;
        else if (k == "pz") has_pz = true;
    }
    init_pos_set = has_x0 && has_y0 && has_z0;
    vertex_set = has_x && has_y && has_z;
    mom_set = has_px && has_py && has_pz;

    // Compute set_variables_ from the column headers
    set_variables_.clear();
    for (auto const & k : keys) {
        if (k == "E") set_variables_.insert(DistributionVariable::PrimaryEnergy);
        else if (k == "m") set_variables_.insert(DistributionVariable::PrimaryMass);
        else if (k == "px" || k == "py" || k == "pz") {
            set_variables_.insert(DistributionVariable::PrimaryDirection);
            set_variables_.insert(DistributionVariable::PrimaryEnergy);
        }
        else if (k == "x" || k == "y" || k == "z") set_variables_.insert(DistributionVariable::InteractionVertex);
        else if (k == "x0" || k == "y0" || k == "z0") set_variables_.insert(DistributionVariable::InitialPosition);
        else set_variables_.insert(DistributionVariable::InteractionParameters);
    }
    DerivePhysicalRowWeights();
}

PrimaryExternalDistribution::PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data, double emin)
    : PrimaryExternalDistribution(std::move(_keys), std::move(_data))
{
    this->emin = emin;
    auto energy_it = std::find(keys.begin(), keys.end(), "E");
    if (energy_it != keys.end()) {
        size_t energy_index = static_cast<size_t>(
            std::distance(keys.begin(), energy_it));
        std::vector<std::vector<double>> filtered;
        for (auto const & row : input_data) {
            if (energy_index < row.size() && row[energy_index] >= emin) {
                filtered.push_back(row);
            }
        }
        input_data = std::move(filtered);
        if (input_data.empty()) {
            throw std::runtime_error("No valid in-memory PrimaryExternalDistribution rows");
        }
        // The physical row weights and normalization derived by the delegated
        // constructor refer to the unfiltered table; rebuild on the survivors.
        DerivePhysicalRowWeights();
    }
}

PrimaryExternalDistribution::PrimaryExternalDistribution(
    std::vector<std::string> _keys,
    std::vector<std::vector<double>> _data,
    std::vector<double> _sampling_weights)
    : PrimaryExternalDistribution(std::move(_keys), std::move(_data))
{
    if (!_sampling_weights.empty()) {
        if (_sampling_weights.size() != input_data.size()) {
            throw std::runtime_error(
                "sampling_weights length (" + std::to_string(_sampling_weights.size()) +
                ") must match data length (" + std::to_string(input_data.size()) + ")");
        }
        // Explicit sampling weights bias the ROW SELECTION only. They are
        // independent of the physical row weights from a "weight" column:
        // GenerationProbability reports this biased sampling density while
        // PhysicalDensity keeps reporting the physical row density, so one
        // shared instance de-biases exactly.
        sampling_weights_ = std::move(_sampling_weights);
        BuildSamplingCDF();
    }
}

PrimaryExternalDistribution::PrimaryExternalDistribution(
    std::vector<std::string> _keys,
    std::vector<std::vector<double>> _data,
    std::vector<double> _sampling_weights,
    double emin)
    : PrimaryExternalDistribution(std::move(_keys), std::move(_data), std::move(_sampling_weights))
{
    this->emin = emin;
    auto energy_it = std::find(keys.begin(), keys.end(), "E");
    if (energy_it != keys.end()) {
        size_t energy_index = static_cast<size_t>(
            std::distance(keys.begin(), energy_it));
        std::vector<std::vector<double>> filtered_data;
        std::vector<double> filtered_weights;
        for (size_t j = 0; j < input_data.size(); ++j) {
            if (energy_index < input_data[j].size() && input_data[j][energy_index] >= emin) {
                filtered_data.push_back(input_data[j]);
                if (!sampling_weights_.empty()) {
                    filtered_weights.push_back(sampling_weights_[j]);
                }
            }
        }
        input_data = std::move(filtered_data);
        if (input_data.empty()) {
            throw std::runtime_error("No valid in-memory PrimaryExternalDistribution rows");
        }
        if (!sampling_weights_.empty()) {
            sampling_weights_ = std::move(filtered_weights);
            BuildSamplingCDF();
        }
        // Physical row weights and normalization rebuild on the survivors.
        DerivePhysicalRowWeights();
    }
}

// Accounts for events above threshold only!
size_t PrimaryExternalDistribution::GetPhysicalNumEvents() const
{
    return input_data.size();
}

void PrimaryExternalDistribution::Sample(
        std::shared_ptr<siren::utilities::SIREN_random> rand,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::PrimaryDistributionRecord & record) const {
    size_t i;
    if (!sampling_cdf_.empty()) {
        double u = rand->Uniform();
        auto it = std::lower_bound(sampling_cdf_.begin(), sampling_cdf_.end(), u);
        i = static_cast<size_t>(std::distance(sampling_cdf_.begin(), it));
        if (i >= input_data.size()) i = input_data.size() - 1;
    } else {
        i = std::min(size_t(rand->Uniform() * input_data.size()), input_data.size() - 1);
    }
    std::array<double, 3> _initial_position;
    std::array<double, 3> _vertex;
    std::array<double, 3> _momentum;
    for(size_t i_key = 0; i_key < keys.size(); ++i_key) {
        double value = input_data[i][i_key];
        if (keys[i_key] == "x0") {
            _initial_position[0] = value;
        }
        else if (keys[i_key] == "y0") {
            _initial_position[1] = value;
        }
        else if (keys[i_key] == "z0") {
            _initial_position[2] = value;
        }
        else if (keys[i_key] == "x") {
            _vertex[0] = value;
        }
        else if (keys[i_key] == "y") {
            _vertex[1] = value;
        }
        else if (keys[i_key] == "z") {
            _vertex[2] = value;
        }
        else if (keys[i_key] == "px") {
            _momentum[0] = value;
        }
        else if (keys[i_key] == "py") {
            _momentum[1] = value;
        }
        else if (keys[i_key] == "pz") {
            _momentum[2] = value;
        }
        else if (keys[i_key] == "t0") {
            record.SetInitialTime(value);
        }
        else if (keys[i_key] == "E") {
            record.SetEnergy(value);
        }
        else if (keys[i_key] == "m") {
            record.SetMass(value);
        }
        else {
            record.SetInteractionParameter(keys[i_key], value);
        }
    }
    // Cache the sampled row so any evaluating instance over the same table
    // can report its own density for this row (injection and physical
    // instances may weight the rows differently, e.g. a biased sampler paired
    // with a physically weighted evaluator).
    record.SetInteractionParameter("PrimaryExternalDistribution_row",
                                   static_cast<double>(i));
    if (!sampling_weights_.empty()) {
        double gen_prob = static_cast<double>(input_data.size()) * sampling_weights_[i] / sampling_weights_sum_;
        record.SetInteractionParameter("PrimaryExternalDistribution_gen_prob", gen_prob);
    }
    if(mom_set) record.SetThreeMomentum(_momentum);
    if(init_pos_set) {
        record.SetInitialPosition(_initial_position);
        if(!vertex_set) record.SetInteractionVertex(_initial_position);
        _cached_position = _initial_position;
    }
    if(vertex_set) {
        record.SetInteractionVertex(_vertex);
        if(!init_pos_set) record.SetInitialPosition(_vertex);
        _cached_position = _vertex;
    }
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> PrimaryExternalDistribution::SamplePosition(
        std::shared_ptr<siren::utilities::SIREN_random> rand,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::PrimaryDistributionRecord & record) const {
    siren::math::Vector3D pos(_cached_position[0], _cached_position[1], _cached_position[2]);
    siren::math::Vector3D dir(0, 0, 1);
    return std::make_tuple(pos, dir);
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> PrimaryExternalDistribution::InjectionBounds(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::InteractionRecord const & interaction) const {
    siren::math::Vector3D pos(interaction.interaction_vertex[0],
                               interaction.interaction_vertex[1],
                               interaction.interaction_vertex[2]);
    siren::math::Vector3D end = pos;
    return std::make_tuple(pos, end);
}

std::vector<std::string> PrimaryExternalDistribution::DensityVariables() const {
    return std::vector<std::string>{"External"};
}

std::set<DistributionVariable> PrimaryExternalDistribution::SetVariables() const {
    return set_variables_;
}

std::set<DistributionVariable> PrimaryExternalDistribution::RequiredVariables() const {
    return {};
}

std::string PrimaryExternalDistribution::Name() const {
    return "PrimaryExternalDistribution";
}

double PrimaryExternalDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model,
                                                          std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
                                                          siren::dataclasses::InteractionRecord const & record) const {
    double energy = record.primary_momentum[0];
    // The SAMPLING density over the supplied rows, relative to the
    // uniform-over-rows base measure: flat 1 without explicit sampling
    // weights, N * s_i / sum(s) with them. Which row the record came from is
    // read from the row index cached at sampling time. Physical row weights
    // (the "weight" column) never appear here; they enter through
    // PhysicalDensity on the physical side of the weight ratio. The legacy
    // cached density is honored for records from before the row cache
    // existed.
    if (!sampling_weights_.empty()) {
        auto row_it = record.interaction_parameters.find("PrimaryExternalDistribution_row");
        if (row_it != record.interaction_parameters.end()) {
            size_t row = static_cast<size_t>(std::llround(row_it->second));
            if (row < sampling_weights_.size()) {
                if (energy < emin) return 0;
                return static_cast<double>(input_data.size())
                    * sampling_weights_[row] / sampling_weights_sum_;
            }
            throw siren::utilities::WeightCalculationError(
                "PrimaryExternalDistribution: cached row index "
                + std::to_string(row) + " is out of range for this table ("
                + std::to_string(sampling_weights_.size()) + " rows); the "
                "record was sampled from a different table");
        }
        auto gp_it = record.interaction_parameters.find("PrimaryExternalDistribution_gen_prob");
        if (gp_it != record.interaction_parameters.end()) {
            if (energy >= emin) return gp_it->second;
            return 0;
        }
        // Records produced by this distribution's own Sample always carry the
        // cached row index (and, when sampling weights are in use, the cached
        // generation probability); their absence means the record did not
        // originate from this table, so the weighted density cannot be
        // recovered.
        throw siren::utilities::WeightCalculationError(
            "PrimaryExternalDistribution: missing cached interaction parameter "
            "\"PrimaryExternalDistribution_row\" (or legacy "
            "\"PrimaryExternalDistribution_gen_prob\")");
    }
    if (energy >= emin) return 1;
    return 0;
}

// The PHYSICAL row density N * w_i / sum(w) relative to the uniform-over-rows
// base measure, from the "weight" column: the importance weights encode how
// to return from the supplied (intentionally biased) parent ensemble to the
// physical one, given the sampling density GenerationProbability reports.
// Together with the sum(w) normalization this makes the event weight carry
// exactly w_i / (sampling pdf of row i). Tables without a weight column fall
// back to the sampling density, restoring the cancel-with-injection default.
double PrimaryExternalDistribution::PhysicalDensity(std::shared_ptr<siren::detector::DetectorModel const> detector_model,
                                                    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
                                                    siren::dataclasses::InteractionRecord const & record) const {
    if (physical_weights_.empty()) {
        return GenerationProbability(detector_model, interactions, record);
    }
    double energy = record.primary_momentum[0];
    auto row_it = record.interaction_parameters.find("PrimaryExternalDistribution_row");
    if (row_it == record.interaction_parameters.end()) {
        // Records produced by any instance over this table always carry the
        // cached row index; its absence means the record did not originate
        // from an external table, so the physical row weight cannot be
        // recovered.
        throw siren::utilities::WeightCalculationError(
            "PrimaryExternalDistribution: missing cached interaction parameter "
            "\"PrimaryExternalDistribution_row\"; the physical row density "
            "requires a record sampled from an external table");
    }
    size_t row = static_cast<size_t>(std::llround(row_it->second));
    if (row >= physical_weights_.size()) {
        throw siren::utilities::WeightCalculationError(
            "PrimaryExternalDistribution: cached row index "
            + std::to_string(row) + " is out of range for this table ("
            + std::to_string(physical_weights_.size()) + " rows); the "
            "record was sampled from a different table");
    }
    if (energy < emin) return 0;
    return static_cast<double>(input_data.size())
        * physical_weights_[row] / physical_weights_sum_;
}

bool PrimaryExternalDistribution::PhysicalDensityDiffers() const {
    return !physical_weights_.empty();
}

std::shared_ptr<PrimaryInjectionDistribution> PrimaryExternalDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new PrimaryExternalDistribution(*this));
}

bool PrimaryExternalDistribution::equal(WeightableDistribution const & other) const {
    const PrimaryExternalDistribution* x = dynamic_cast<const PrimaryExternalDistribution*>(&other);
    if(!x)
        return false;
    return emin == x->emin
        && keys == x->keys
        && input_data == x->input_data
        && sampling_weights_ == x->sampling_weights_;
}

bool PrimaryExternalDistribution::less(WeightableDistribution const & other) const {
    const PrimaryExternalDistribution* x = dynamic_cast<const PrimaryExternalDistribution*>(&other);
    if(!x)
        return false;
    return std::tie(emin, keys, input_data, sampling_weights_)
        < std::tie(x->emin, x->keys, x->input_data, x->sampling_weights_);
}


} // namespace distributions
} // namespace siren
