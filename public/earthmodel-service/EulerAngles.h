#ifndef LI_EulerAngles_H
#define LI_EulerAngles_H

#include <sstream>

namespace earthmodel {

enum EulerAxis {
    X = 0,
    Y = 1,
    Z = 2
};

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

constexpr unsigned int EulerAxisI(EulerOrder order) {
    return 0;
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
private:
    EulerOrder order_;
    double alpha_;
    double beta_;
    double gamma_;
};

} // namespace earthmodel

#endif // LI_EulerAngles_H

