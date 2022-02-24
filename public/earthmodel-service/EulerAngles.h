#pragma once
#ifndef LI_EulerAngles_H
#define LI_EulerAngles_H

#include <sstream>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include "earthmodel-service/Matrix3D.h"

namespace earthmodel {

enum EulerFrame {
    Static = 0,
    Rotating = 1
};

constexpr EulerFrame GetEulerFrame(unsigned int order) {
    return EulerFrame(order & 1);
};

enum EulerRepetition {
    No = 0,
    Yes = 1
};

constexpr EulerRepetition GetEulerRepetition(unsigned int order) {
    return EulerRepetition((order >> 1) & 1);
};

enum EulerParity {
    Even = 0,
    Odd = 1
};

constexpr EulerParity GetEulerParity(unsigned int order) {
    return EulerParity((order >> 2) & 1);
};

enum EulerAxis {
    X = 0,
    Y = 1,
    Z = 2
};

constexpr EulerAxis GetEulerAxis(unsigned int order) {
    return EulerAxis(((order)>>3)&3);
}

constexpr unsigned int GetEulerOrder(EulerAxis axis, EulerParity parity, EulerRepetition repetition, EulerFrame frame) {
    return (((((((axis)<<1)+(parity))<<1)+(repetition))<<1)+(frame));
};

enum EulerOrder {
    XYZs = GetEulerOrder(EulerAxis::X, EulerParity::Even, EulerRepetition::No,  EulerFrame::Static),
    XYXs = GetEulerOrder(EulerAxis::X, EulerParity::Even, EulerRepetition::Yes, EulerFrame::Static),
    XZYs = GetEulerOrder(EulerAxis::X, EulerParity::Odd,  EulerRepetition::No,  EulerFrame::Static),
    XZXs = GetEulerOrder(EulerAxis::X, EulerParity::Odd,  EulerRepetition::Yes, EulerFrame::Static),
    YZXs = GetEulerOrder(EulerAxis::Y, EulerParity::Even, EulerRepetition::No,  EulerFrame::Static),
    YZYs = GetEulerOrder(EulerAxis::Y, EulerParity::Even, EulerRepetition::Yes, EulerFrame::Static),
    YXZs = GetEulerOrder(EulerAxis::Y, EulerParity::Odd,  EulerRepetition::No,  EulerFrame::Static),
    YXYs = GetEulerOrder(EulerAxis::Y, EulerParity::Odd,  EulerRepetition::Yes, EulerFrame::Static),
    ZXYs = GetEulerOrder(EulerAxis::Z, EulerParity::Even, EulerRepetition::No,  EulerFrame::Static),
    ZXZs = GetEulerOrder(EulerAxis::Z, EulerParity::Even, EulerRepetition::Yes, EulerFrame::Static),
    ZYXs = GetEulerOrder(EulerAxis::Z, EulerParity::Odd,  EulerRepetition::No,  EulerFrame::Static),
    ZYZs = GetEulerOrder(EulerAxis::Z, EulerParity::Odd,  EulerRepetition::Yes, EulerFrame::Static),
    ZYXr = GetEulerOrder(EulerAxis::X, EulerParity::Even, EulerRepetition::No,  EulerFrame::Rotating),
    XYXr = GetEulerOrder(EulerAxis::X, EulerParity::Even, EulerRepetition::Yes, EulerFrame::Rotating),
    YZXr = GetEulerOrder(EulerAxis::X, EulerParity::Odd,  EulerRepetition::No,  EulerFrame::Rotating),
    XZXr = GetEulerOrder(EulerAxis::X, EulerParity::Odd,  EulerRepetition::Yes, EulerFrame::Rotating),
    XZYr = GetEulerOrder(EulerAxis::Y, EulerParity::Even, EulerRepetition::No,  EulerFrame::Rotating),
    YZYr = GetEulerOrder(EulerAxis::Y, EulerParity::Even, EulerRepetition::Yes, EulerFrame::Rotating),
    ZXYr = GetEulerOrder(EulerAxis::Y, EulerParity::Odd,  EulerRepetition::No,  EulerFrame::Rotating),
    YXYr = GetEulerOrder(EulerAxis::Y, EulerParity::Odd,  EulerRepetition::Yes, EulerFrame::Rotating),
    YXZr = GetEulerOrder(EulerAxis::Z, EulerParity::Even, EulerRepetition::No,  EulerFrame::Rotating),
    ZXZr = GetEulerOrder(EulerAxis::Z, EulerParity::Even, EulerRepetition::Yes, EulerFrame::Rotating),
    XYZr = GetEulerOrder(EulerAxis::Z, EulerParity::Odd,  EulerRepetition::No,  EulerFrame::Rotating),
    ZYZr = GetEulerOrder(EulerAxis::Z, EulerParity::Odd,  EulerRepetition::Yes, EulerFrame::Rotating)
};

constexpr unsigned int EulerSafe[4] = {0, 1, 2, 0};
constexpr unsigned int EulerNext[4] = {1, 2, 0, 1};

constexpr EulerAxis GetEulerAxisI(EulerOrder order) {
    return (EulerAxis)EulerSafe[(unsigned int)GetEulerAxis(order)];
};

constexpr EulerAxis GetEulerAxisJ(EulerOrder order) {
    return (EulerAxis)EulerNext[(unsigned int)(GetEulerAxisI(order) + (GetEulerParity(order) == EulerParity::Odd))];
};

constexpr EulerAxis GetEulerAxisK(EulerOrder order) {
    return (EulerAxis)EulerNext[(unsigned int)(GetEulerAxisI(order) + (GetEulerParity(order) != EulerParity::Odd))];
};

constexpr EulerAxis GetEulerAxisH(EulerOrder order) {
    return (EulerAxis)((GetEulerRepetition(order) == EulerRepetition::No) ? GetEulerAxisK(order) : GetEulerAxisI(order));
};

class EulerAngles
{
public:
    // constructors
    EulerAngles();
    EulerAngles(const EulerOrder order, const double alpha, const double beta, const double gamma);
    EulerAngles(const EulerAngles& euler);
    EulerAngles(EulerAngles&& other);
    ~EulerAngles();

    //-------------------------------------//
    // operator functions and swap
    EulerAngles& operator=(EulerAngles const & euler);
    EulerAngles& operator=(EulerAngles const && euler);
    EulerAngles& operator=(EulerAngles && euler);
    bool operator==(const EulerAngles& euler) const;
    bool operator!=(const EulerAngles& euler) const;
    void swap(EulerAngles& euler);
    friend std::ostream& operator<<(std::ostream& os, EulerAngles const& euler);

    EulerOrder GetOrder() const {return order_;}
    double GetAlpha() const {return alpha_;}
    double GetBeta() const {return beta_;}
    double GetGamma() const {return gamma_;}

    void SetOrder(EulerOrder order) {order_ = order;}
    void SetAlpha(double alpha) {alpha_ = alpha;}
    void SetBeta(double beta) {beta_ = beta;}
    void SetGamma(double gamma) {gamma_ = gamma;}

    Matrix3D GetMatrix() const;

    //-------------------------------------//
    // serialization
    //----------------------------------------------//
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("EulerOrder", order_));
            archive(cereal::make_nvp("Alpha", alpha_));
            archive(cereal::make_nvp("Beta", beta_));
            archive(cereal::make_nvp("Gamma", gamma_));
        } else {
            throw std::runtime_error("EulerAngles only supports version <= 0!");
        }
    }
private:
    EulerOrder order_;
    double alpha_;
    double beta_;
    double gamma_;
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::EulerAngles, 0);

#endif // LI_EulerAngles_H

