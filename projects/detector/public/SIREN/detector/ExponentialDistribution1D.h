#pragma once
#ifndef SIREN_detector_ExponentialDistribution1D_H
#define SIREN_detector_ExponentialDistribution1D_H

#include <memory>                                    // for shared_ptr
#include <cstdint>                                   // for uint32_t
#include <stdexcept>                                 // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/detector/Distribution1D.h"

namespace siren {
namespace detector {

class ExponentialDistribution1D : public Distribution1D {
friend cereal::access;
public:
    ExponentialDistribution1D();
    ExponentialDistribution1D(const ExponentialDistribution1D&);
    ExponentialDistribution1D(double sigma);
    // Full form: rho(x) = amplitude * exp(sigma * (x - x0)).
    // The x0 offset lets the exponent be measured from a reference point near
    // the region of interest (e.g. a shell radius), avoiding overflow/underflow
    // when the raw axis coordinate is large (e.g. a radial atmosphere at ~6.37e6 m).
    ExponentialDistribution1D(double sigma, double amplitude, double x0);
    bool compare(const Distribution1D& dist) const override;
    Distribution1D* clone() const override { return new ExponentialDistribution1D(*this); };
    std::shared_ptr<Distribution1D> create() const override {
        return std::shared_ptr<Distribution1D>(new ExponentialDistribution1D(*this));
    };
    double Derivative(double x) const override;
    double AntiDerivative(double x) const override;
    double Evaluate(double x) const override;
    double GetSigma() const { return sigma_; };
    double GetAmplitude() const { return amplitude_; };
    double GetX0() const { return x0_; };
    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                // Legacy: sigma only; amplitude=1, x0=0 (identical Evaluate).
                archive(::cereal::make_nvp("Sigma", sigma_));
                amplitude_ = 1.0;
                x0_ = 0.0;
                archive(cereal::virtual_base_class<Distribution1D>(this));
            } else if(version == 1) {
                archive(::cereal::make_nvp("Sigma", sigma_));
                archive(::cereal::make_nvp("Amplitude", amplitude_));
                archive(::cereal::make_nvp("X0", x0_));
                archive(cereal::virtual_base_class<Distribution1D>(this));
            } else {
                throw std::runtime_error("ExponentialDistribution1D only supports version <= 1");
            }
        };
protected:
    double sigma_;
    double amplitude_;
    double x0_;
};

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::ExponentialDistribution1D, 1);
CEREAL_REGISTER_TYPE(siren::detector::ExponentialDistribution1D);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::detector::Distribution1D, siren::detector::ExponentialDistribution1D);

CEREAL_FORCE_DYNAMIC_INIT(siren_ExponentialDistribution1D);

#endif // SIREN_detector_ExponentialDistribution1D_H
