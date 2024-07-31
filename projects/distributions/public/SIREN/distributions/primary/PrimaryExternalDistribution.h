#pragma once
#ifndef SIREN_PrimaryExternalDistribution_H
#define SIREN_PrimaryExternalDistribution_H

#include <memory>                                        // for shared_ptr
#include <string>                                        // for string
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/Distributions.h"  // for WeightableDi...

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }
namespace cereal { class access; }

namespace siren {
namespace distributions {

class PrimaryExternalDistribution : virtual public PrimaryInjectionDistribution {
friend cereal::access;
protected:
    PrimaryExternalDistribution() {};
public:
    ~PrimaryExternalDistribution() {
        input_file.close();
    };
    std::string filename;
private:
    std::ifstream input_file;
    std::vector<std::vector<double>> input_data;
    std::vector<std::string> keys;
    std::vector<std::string> possible_keys = {"x0","y0","z0","m","E","px","py","pz"};
    bool init_pos_set;
    bool mom_set;
public:
    PrimaryExternalDistribution(std::string _filename);
    PrimaryExternalDistribution(PrimaryExternalDistribution const & other);
    void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryExternalDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryExternalDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PrimaryExternalDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::PrimaryExternalDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryInjectionDistribution, siren::distributions::PrimaryExternalDistribution);

#endif // SIREN_PrimaryExternalDistribution_H