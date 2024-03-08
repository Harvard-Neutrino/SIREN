#pragma once
#ifndef LI_ConstantDistribution1D_H
#define LI_ConstantDistribution1D_H

#include <memory>                                    // for shared_ptr
#include <cstdint>                                   // for uint32_t
#include <stdexcept>                                 // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/detector/Distribution1D.h"

namespace SI {
namespace detector {

class ConstantDistribution1D : public Distribution1D {
public:
    ConstantDistribution1D();
    ConstantDistribution1D(const ConstantDistribution1D&);
    ConstantDistribution1D(double val);
    bool compare(const Distribution1D& dist) const override;
    Distribution1D* clone() const override { return new ConstantDistribution1D(*this); };
    std::shared_ptr<Distribution1D> create() const override {
        return std::shared_ptr<Distribution1D>(new ConstantDistribution1D(*this));
    };
    double Derivative(double x) const override;
    double AntiDerivative(double x) const override;
    double Evaluate(double x) const override;
    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(::cereal::make_nvp("Value", val_));
                archive(cereal::virtual_base_class<Distribution1D>(this));
            } else {
                throw std::runtime_error("ConstantDistribution1D only supports version <= 0");
            }
        };
protected:
    double val_;
};

} // namespace detector
} // namespace SI

CEREAL_CLASS_VERSION(SI::detector::ConstantDistribution1D, 0);
CEREAL_REGISTER_TYPE(SI::detector::ConstantDistribution1D);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::detector::Distribution1D, SI::detector::ConstantDistribution1D);

#endif // LI_ConstantDistribution1D_H
