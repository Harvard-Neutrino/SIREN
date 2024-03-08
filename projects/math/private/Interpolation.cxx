#include "SIREN/math/Interpolation.h"

#include <tuple>
#include <memory>

namespace SI {
namespace math {

//template<> class Transform<double>;

//template<> class IdentityTransform<double>;

//template<> class GenericTransform<double>;

//template<> class LogTransform<double>;

//template<> class SymLogTransform<double>;

//template<> class RangeTransform<double>;

//template<> struct LinearInterpolationOperator<double>;

//template<> struct DropLinearInterpolationOperator<double>;

template<> std::tuple<std::shared_ptr<Transform<double>>, std::shared_ptr<Transform<double>>> SelectInterpolationSpace1D(
        std::vector<double> const & x,
        std::vector<double> const & y,
        std::shared_ptr<LinearInterpolationOperator<double>> interp);

//template<> class RegularIndexer1D<double>;

//template<> class IrregularIndexer1D<double>;

//template<> class TransformIndexer1D<double>;

template<>
std::shared_ptr<Indexer1D<double>> SelectIndexer1D(
        std::vector<double> x,
        std::shared_ptr<Transform<double>> x_transform);

template<> class LinearInterpolator1D<double>;

template<> struct BiLinearInterpolationOperator<double>;

template<> struct DropBiLinearInterpolationOperator<double>;

using RegularGridIndexer2DDouble = RegularGridIndexer2D<double>;

using IrregularGridIndexer2DDouble = IrregularGridIndexer2D<double>;

template<> class GridLinearInterpolator2D<double>;

template<> class DelaunayIndexer2D<double>;

template<> class LinearDelaunayInterpolator2D<double>;

} // namespace math
} // namespace SI
