
#include "LeptonInjector/math/Interpolation.h"

namespace LI {
namespace math {

template<> class Transform<double>;

template<> class IdentityTransform<double>;

template<> class GenericTransform<double>;

template<> class LogTransform<double>;

template<> class SymLogTransform<double>;

template<> class RangeTransform<double>;

template<> class FunctionalRangeTransform<double>;

template<> class LinearInterpolator<double>;

template<> class DropLinearInterpolator<double>;

template<> std::tuple<std::shared_ptr<Transform<double>>, std::shared_ptr<Transform<double>>> DetermineInterpolationSpace1D(
        std::vector<double> const & x,
        std::vector<double> const & y,
        std::shared_ptr<LinearInterpolator<double>> interp);

} // namespace math
} // namespace LI
