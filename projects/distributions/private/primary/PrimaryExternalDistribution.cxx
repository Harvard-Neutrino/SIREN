#include "SIREN/distributions/primary/PrimaryExternalDistribution.h"

#include <array>                                           // for array
#include <string>                                          // for basic_string

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren {
namespace distributions {

//---------------
// class PrimaryExternalDistribution : PrimaryExternalDistribution
//---------------
PrimaryExternalDistribution::PrimaryExternalDistribution(std::string _filename)
{
    filename = _filename;
    input_file.open(filename);
    init_pos_set = false;
    mom_set = false;

    std::string line;
    if (!input_file.is_open()) {
        std::cerr << "error: file open failed " << filename << ".\n";
        exit(0);
    }
    std::getline(input_file,line);

    std::stringstream ss(line);
    std::string key;
    while (std::getline(ss, key, ',')) {
        keys.push_back(key);
        if (key == "x0" ||
            key == "y0" ||
            key == "z0") init_pos_set = true;
        if (key == "px" ||
            key == "py" ||
            key == "pz") mom_set = true;
    }

    std::string value;
    // fill input data
    while (std::getline(input_file,line)) {
        std::vector<double> tmp_data;
        ss.str(line);
        while (std::getline(ss, value, ',')) {
            tmp_data.push_back(stod(value));
        }
        input_data.push_back(tmp_data);
    }

}

PrimaryExternalDistribution::PrimaryExternalDistribution(PrimaryExternalDistribution const & other)
{
    PrimaryExternalDistribution(other.filename);
    // input_file.open(filename);
    // init_pos_set = false;
    // mom_set = false;

    // std::string line;
    // if (!input_file.is_open()) {
    //     std::cerr << "error: file open failed " << filename << ".\n";
    //     exit(0);
    // }
    // std::getline(input_file,line);

    // std::stringstream ss(line);
    // std::string key;
    // while (std::getline(ss, key, ',')) {
    //     keys.push_back(key);
    //     if (key == "x0" ||
    //         key == "y0" ||
    //         key == "z0") init_pos_set = true;
    //     if (key == "px" ||
    //         key == "py" ||
    //         key == "pz") mom_set = true;
    // }

    // std::string value;
    // // fill input data
    // while (std::getline(input_file,line)) {
    //     std::vector<double> tmp_data;
    //     ss.str(line);
    //     while (std::getline(ss, value, ',')) {
    //         tmp_data.push_back(stod(value));
    //     }
    //     input_data.push_back(tmp_data);
    // }

}

void PrimaryExternalDistribution::Sample(
        std::shared_ptr<siren::utilities::SIREN_random> rand,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::PrimaryDistributionRecord & record) const {

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

std::vector<std::string> PrimaryExternalDistribution::DensityVariables() const {
    return std::vector<std::string>{"External"};
}

std::string PrimaryExternalDistribution::Name() const {
    return "PrimaryExternalDistribution";
}

double PrimaryExternalDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model,
                                                          std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
                                                          siren::dataclasses::InteractionRecord const & record) const {
    return 1;
}

std::shared_ptr<PrimaryInjectionDistribution> PrimaryExternalDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new PrimaryExternalDistribution(*this));
}

bool PrimaryExternalDistribution::equal(WeightableDistribution const & other) const {
    const PrimaryExternalDistribution* x = dynamic_cast<const PrimaryExternalDistribution*>(&other);

    if(!x)
        return false;
    else
        return filename == x->filename;
}

bool PrimaryExternalDistribution::less(WeightableDistribution const & other) const {
    const PrimaryExternalDistribution* x = dynamic_cast<const PrimaryExternalDistribution*>(&other);
    return filename != x->filename;
}


} // namespace distributions
} // namespace sirenREN