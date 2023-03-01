#pragma once
#ifndef LI_Distribution1D_H
#define LI_Distribution1D_H
#include <memory>
#include <string>
#include <exception>
#include <functional>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

namespace LI {
namespace detector {

class Distribution1D {
friend cereal::access;
public:
    virtual ~Distribution1D() = default;
    bool operator==(const Distribution1D& dist) const;
    bool operator!=(const Distribution1D& dist) const;
    virtual bool compare(const Distribution1D& dist) const = 0;
    virtual Distribution1D* clone() const = 0;
    virtual std::shared_ptr<const Distribution1D> create() const = 0;
    virtual double Derivative(double x) const = 0;
    virtual double AntiDerivative(double x) const = 0;
    virtual double Evaluate(double x) const = 0;

    template<class Archive>
        void serialize(Archive & archive, std::uint32_t const version) {};
};

} // namespace detector
} // namespace LI

CEREAL_CLASS_VERSION(LI::detector::Distribution1D, 0);

#endif // LI_Distribution1D_H
