#pragma once
#ifndef SIREN_PolynomialDistribution1D_H
#define SIREN_PolynomialDistribution1D_H

#include <memory>                                    // for shared_ptr
#include <vector>                                    // for vector
#include <cstdint>                                   // for uint32_t
#include <stdexcept>                                 // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/math/Polynomial.h"          // for Polynom

#include "SIREN/detector/Distribution1D.h"  // for Distribution1D

namespace siren {
namespace detector {

class PolynomialDistribution1D : public Distribution1D {
friend cereal::access;
protected:
public:
    PolynomialDistribution1D();
    PolynomialDistribution1D(const PolynomialDistribution1D&);
    PolynomialDistribution1D(const math::Polynom&);
    PolynomialDistribution1D(const std::vector<double>&);
    bool compare(const Distribution1D& dist) const override;
    Distribution1D* clone() const override { return new PolynomialDistribution1D(*this); };
    std::shared_ptr<Distribution1D> create() const override {
        return std::shared_ptr<Distribution1D>(new PolynomialDistribution1D(*this));
    };
    double Derivative(double x) const override;
    double AntiDerivative(double x) const override;
    double Evaluate(double x) const override;
    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(::cereal::make_nvp("Polynomial", polynom_));
                archive(::cereal::make_nvp("PolynomialIntegral", Ipolynom_));
                archive(::cereal::make_nvp("PolynomialDerivative", dpolynom_));
                archive(cereal::virtual_base_class<Distribution1D>(this));
            } else {
                throw std::runtime_error("PolynomialDistribution1D only supports version <= 0");
            }
        };
protected:
    math::Polynom polynom_;
    math::Polynom Ipolynom_;
    math::Polynom dpolynom_;
};

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::PolynomialDistribution1D, 0);
CEREAL_REGISTER_TYPE(siren::detector::PolynomialDistribution1D);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::detector::Distribution1D, siren::detector::PolynomialDistribution1D);

#endif // SIREN_PolynomialDistribution1D_H
