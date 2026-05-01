#include "SIREN/distributions/primary/PrimaryExternalDistribution.h"

#include <array>                                           // for array
#include <fstream>                                         // for ifstream
#include <sstream>                                         // for stringstream
#include <string>                                          // for basic_string
#include <stdexcept>                                       // for runtime_error
#include <tuple>                                           // for tie

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren {
namespace distributions {

//---------------
// class PrimaryExternalDistribution : PrimaryExternalDistribution
//---------------

void PrimaryExternalDistribution::LoadInputFile(std::string _filename)
{
    filename = _filename;
    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        throw std::runtime_error("error: file open failed " + filename);
    }

    std::string line;
    std::getline(input_file, line);

    std::stringstream ss(line);
    std::string key;
    bool has_x0 = false, has_y0 = false, has_z0 = false;
    bool has_px = false, has_py = false, has_pz = false;
    while (std::getline(ss, key, ',')) {
        keys.push_back(key);
        if (key == "x0") has_x0 = true;
        else if (key == "y0") has_y0 = true;
        else if (key == "z0") has_z0 = true;
        else if (key == "px") has_px = true;
        else if (key == "py") has_py = true;
        else if (key == "pz") has_pz = true;
    }
    init_pos_set = has_x0 && has_y0 && has_z0;
    mom_set = has_px && has_py && has_pz;

    std::string value;
    while (std::getline(input_file, line)) {
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
}

PrimaryExternalDistribution::PrimaryExternalDistribution(std::string _filename) : emin(0)
{
    LoadInputFile(_filename);
}

PrimaryExternalDistribution::PrimaryExternalDistribution(std::string _filename, double emin) : emin(emin)
{
    LoadInputFile(_filename);
}

// Accounts for events above threshold only!
int PrimaryExternalDistribution::GetPhysicalNumEvents() const
{
    return input_data.size();
}

void PrimaryExternalDistribution::Sample(
        std::shared_ptr<siren::utilities::SIREN_random> rand,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::PrimaryDistributionRecord & record) const {

    size_t max_tries = 1000;
    size_t num_tries = 0;
    bool success = false;
    while(!success && num_tries < max_tries) {
        ++num_tries;
        int i = int(rand->Uniform() * input_data.size());
        int i_key = 0;
        std::array<double, 3> _initial_position;
        std::array<double, 3> _momentum;
        for (auto value : input_data[i]) {
            if (keys[i_key] == "x0") {
                _initial_position[0] = value;
            }
            else if (keys[i_key] == "y0") {
                _initial_position[1] = value;
            }
            else if (keys[i_key] == "z0") {
                _initial_position[2] = value;
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
            else if (keys[i_key] == "E") {
                if (value > emin) success=true;
                else success=false;
                record.SetEnergy(value);
            }
            else if (keys[i_key] == "m") {
                record.SetMass(value);
            }
            else {
                record.SetInteractionParameter(keys[i_key],value);
            }
            ++i_key;
        }
        if(mom_set) record.SetThreeMomentum(_momentum);
        if(init_pos_set) record.SetInitialPosition(_initial_position);
    }
    if(!success) {
        throw std::runtime_error("Failed to generate a physical primary after " + std::to_string(max_tries) + " attempts");
    }
}

std::vector<std::string> PrimaryExternalDistribution::DensityVariables() const {
    return std::vector<std::string>{"External"};
}

std::string PrimaryExternalDistribution::Name() const {
    return "PrimaryExternalDistribution";
}

double PrimaryExternalDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model,
                                                          std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
                                                          siren::dataclasses::InteractionRecord const & record) const {
    double energy = record.primary_momentum[0];
    if (energy > emin) return 1;
    return 0;
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
        && input_data == x->input_data;
}

bool PrimaryExternalDistribution::less(WeightableDistribution const & other) const {
    const PrimaryExternalDistribution* x = dynamic_cast<const PrimaryExternalDistribution*>(&other);
    return std::tie(emin, keys, input_data)
        < std::tie(x->emin, x->keys, x->input_data);
}


} // namespace distributions
} // namespace siren
