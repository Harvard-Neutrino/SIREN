#pragma once
#ifndef SIREN_Distribution1D_H
#define SIREN_Distribution1D_H
#include <memory>                 // for shared_ptr
#include <cstdint>                // for uint32_t

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

namespace siren {
namespace detector {

class Distribution1D {
friend cereal::access;
public:
    virtual ~Distribution1D() = default;
    bool operator==(const Distribution1D& dist) const;
    bool operator!=(const Distribution1D& dist) const;
    virtual bool compare(const Distribution1D& dist) const = 0;
    virtual Distribution1D* clone() const = 0;
    virtual std::shared_ptr<Distribution1D> create() const = 0;
    virtual double Derivative(double x) const = 0;
    virtual double AntiDerivative(double x) const = 0;
    virtual double Evaluate(double x) const = 0;

    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {};
};

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::Distribution1D, 0);

#endif // SIREN_Distribution1D_H
