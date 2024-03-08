#include "SIREN/math/Matrix3D.h"

#include <utility>
#include <iostream>
#include <stdexcept>
#include <initializer_list>

#include "SIREN/math/Vector3D.h"

using namespace SI::math;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
Matrix3D::Matrix3D() :
    xx_(0.0),
    xy_(0.0),
    xz_(0.0),
    yx_(0.0),
    yy_(0.0),
    yz_(0.0),
    zx_(0.0),
    zy_(0.0),
    zz_(0.0)
{
}

Matrix3D::Matrix3D(
    const double xx,
    const double xy,
    const double xz,
    const double yx,
    const double yy,
    const double yz,
    const double zx,
    const double zy,
    const double zz
) :
    xx_(xx),
    xy_(xy),
    xz_(xz),
    yx_(yx),
    yy_(yy),
    yz_(yz),
    zx_(zx),
    zy_(zy),
    zz_(zz)
{
}

// copy constructor
Matrix3D::Matrix3D(const Matrix3D& matrix_3d) :
    xx_(matrix_3d.xx_),
    xy_(matrix_3d.xy_),
    xz_(matrix_3d.xz_),
    yx_(matrix_3d.yx_),
    yy_(matrix_3d.yy_),
    yz_(matrix_3d.yz_),
    zx_(matrix_3d.zx_),
    zy_(matrix_3d.zy_),
    zz_(matrix_3d.zz_)
{
}

Matrix3D::Matrix3D(Matrix3D&& other) :
    xx_(std::move(other.xx_)),
    xy_(std::move(other.xy_)),
    xz_(std::move(other.xz_)),
    yx_(std::move(other.yx_)),
    yy_(std::move(other.yy_)),
    yz_(std::move(other.yz_)),
    zx_(std::move(other.zx_)),
    zy_(std::move(other.zy_)),
    zz_(std::move(other.zz_))
{
}

// destructor
Matrix3D::~Matrix3D() {}

//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//

Matrix3D& Matrix3D::operator=(Matrix3D const & matrix_3d) {
    if (this != &matrix_3d)
    {
        Matrix3D tmp(matrix_3d);
        swap(tmp);
    }
    return *this;
}

Matrix3D& Matrix3D::operator=(Matrix3D && other) {
    xx_ = std::move(other.xx_);
    xy_ = std::move(other.xy_);
    xz_ = std::move(other.xz_);
    yx_ = std::move(other.yx_);
    yy_ = std::move(other.yy_);
    yz_ = std::move(other.yz_);
    zx_ = std::move(other.zx_);
    zy_ = std::move(other.zy_);
    zz_ = std::move(other.zz_);
    return *this;
}

Matrix3D& Matrix3D::operator=(Matrix3D const && other) {
    xx_ = other.xx_;
    xy_ = other.xy_;
    xz_ = other.xz_;
    yx_ = other.yx_;
    yy_ = other.yy_;
    yz_ = other.yz_;
    zx_ = other.zx_;
    zy_ = other.zy_;
    zz_ = other.zz_;
    return *this;
}

bool Matrix3D::operator==(const Matrix3D& matrix_3d) const
{
    return (this == &matrix_3d) or (
        xx_ == matrix_3d.xx_ and
        xy_ == matrix_3d.yy_ and
        xz_ == matrix_3d.zz_ and
        yx_ == matrix_3d.xx_ and
        yy_ == matrix_3d.yy_ and
        yz_ == matrix_3d.zz_ and
        zx_ == matrix_3d.xx_ and
        zy_ == matrix_3d.yy_ and
        zz_ == matrix_3d.zz_);
}

bool Matrix3D::operator!=(const Matrix3D& matrix_3d) const
{
    return !(*this == matrix_3d);
}

void Matrix3D::swap(Matrix3D& matrix_3d)
{
    using std::swap;

    swap(xx_, matrix_3d.xx_);
    swap(xy_, matrix_3d.xy_);
    swap(xz_, matrix_3d.xz_);
    swap(yx_, matrix_3d.yx_);
    swap(yy_, matrix_3d.yy_);
    swap(yz_, matrix_3d.yz_);
    swap(zx_, matrix_3d.zx_);
    swap(zy_, matrix_3d.zy_);
    swap(zz_, matrix_3d.zz_);
}

namespace SI {
namespace math {
std::ostream& operator<<(std::ostream& os, Matrix3D const& matrix_3d)
{
    std::stringstream ss;
    ss << " Matrix3D (" << &matrix_3d << ") ";
    os << ss.str() << '\n';
    return os;
}
} // namespace math
} // namespace SI

//----------------------------------------------------------------------//
//-----------------------operator basic arithmetic ---------------------//
//----------------------------------------------------------------------//

namespace SI {
namespace math {

Matrix3D operator+(const Matrix3D& mat1, const Matrix3D& mat2)
{
    Matrix3D matrix_sum;
    matrix_sum.xx_ = mat1.xx_ + mat2.xx_;
    matrix_sum.xy_ = mat1.xy_ + mat2.xy_;
    matrix_sum.xz_ = mat1.xz_ + mat2.xz_;
    matrix_sum.yx_ = mat1.yx_ + mat2.yx_;
    matrix_sum.yy_ = mat1.yy_ + mat2.yy_;
    matrix_sum.yz_ = mat1.yz_ + mat2.yz_;
    matrix_sum.zx_ = mat1.zx_ + mat2.zx_;
    matrix_sum.zy_ = mat1.zy_ + mat2.zy_;
    matrix_sum.zz_ = mat1.zz_ + mat2.zz_;
    return matrix_sum;
}

Matrix3D& operator+=(Matrix3D& mat1, const Matrix3D& mat2)
{
    mat1.xx_ += mat2.xx_;
    mat1.xy_ += mat2.xy_;
    mat1.xz_ += mat2.xz_;
    mat1.yx_ += mat2.yx_;
    mat1.yy_ += mat2.yy_;
    mat1.yz_ += mat2.yz_;
    mat1.zx_ += mat2.zx_;
    mat1.zy_ += mat2.zy_;
    mat1.zz_ += mat2.zz_;
    return mat1;
}

Matrix3D operator-(const Matrix3D& mat1, const Matrix3D& mat2)
{
    Matrix3D matrix_diff;
    matrix_diff.xx_ = mat1.xx_ - mat2.xx_;
    matrix_diff.xy_ = mat1.xy_ - mat2.xy_;
    matrix_diff.xz_ = mat1.xz_ - mat2.xz_;
    matrix_diff.yx_ = mat1.yx_ - mat2.yx_;
    matrix_diff.yy_ = mat1.yy_ - mat2.yy_;
    matrix_diff.yz_ = mat1.yz_ - mat2.yz_;
    matrix_diff.zx_ = mat1.zx_ - mat2.zx_;
    matrix_diff.zy_ = mat1.zy_ - mat2.zy_;
    matrix_diff.zz_ = mat1.zz_ - mat2.zz_;
    return matrix_diff;
}

Matrix3D operator*(const double factor1, const Matrix3D& mat1)
{
    Matrix3D product;
    product.xx_ = factor1 * mat1.xx_;
    product.xy_ = factor1 * mat1.xy_;
    product.xz_ = factor1 * mat1.xz_;
    product.yx_ = factor1 * mat1.yx_;
    product.yy_ = factor1 * mat1.yy_;
    product.yz_ = factor1 * mat1.yz_;
    product.zx_ = factor1 * mat1.zx_;
    product.zy_ = factor1 * mat1.zy_;
    product.zz_ = factor1 * mat1.zz_;
    return product;
}

Matrix3D operator*(const Matrix3D& mat1, const double factor1)
{
    Matrix3D product;
    product.xx_ = factor1 * mat1.xx_;
    product.xy_ = factor1 * mat1.xy_;
    product.xz_ = factor1 * mat1.xz_;
    product.yx_ = factor1 * mat1.yx_;
    product.yy_ = factor1 * mat1.yy_;
    product.yz_ = factor1 * mat1.yz_;
    product.zx_ = factor1 * mat1.zx_;
    product.zy_ = factor1 * mat1.zy_;
    product.zz_ = factor1 * mat1.zz_;
    return product;
}

Vector3D operator*(const Matrix3D& mat1, const Vector3D& vec1)
{
    Vector3D product;
    product.cartesian_.x_ =
        mat1.xx_ * vec1.cartesian_.x_ +
        mat1.xy_ * vec1.cartesian_.y_ +
        mat1.xz_ * vec1.cartesian_.z_;
    product.cartesian_.y_ =
        mat1.yx_ * vec1.cartesian_.x_ +
        mat1.yy_ * vec1.cartesian_.y_ +
        mat1.yz_ * vec1.cartesian_.z_;
    product.cartesian_.z_ = mat1.zx_ * vec1.cartesian_.x_ +
        mat1.zy_ * vec1.cartesian_.y_ +
        mat1.zz_ * vec1.cartesian_.z_;
    return product;
}

Matrix3D operator/(const Matrix3D& mat1, const double factor1)
{
    Matrix3D product;
    product.xx_ = mat1.xx_ / factor1;
    product.xy_ = mat1.xy_ / factor1;
    product.xz_ = mat1.xz_ / factor1;
    product.yx_ = mat1.yx_ / factor1;
    product.yy_ = mat1.yy_ / factor1;
    product.yz_ = mat1.yz_ / factor1;
    product.zx_ = mat1.zx_ / factor1;
    product.zy_ = mat1.zy_ / factor1;
    product.zz_ = mat1.zz_ / factor1;
    return product;
}

Matrix3D& operator*=(Matrix3D& mat1, const double factor1)
{
    mat1.xx_ *= factor1;
    mat1.xy_ *= factor1;
    mat1.xz_ *= factor1;
    mat1.yx_ *= factor1;
    mat1.yy_ *= factor1;
    mat1.yz_ *= factor1;
    mat1.zx_ *= factor1;
    mat1.zy_ *= factor1;
    mat1.zz_ *= factor1;
    return mat1;
}

Matrix3D& operator/=(Matrix3D& mat1, const double factor1)
{
    mat1.xx_ /= factor1;
    mat1.xy_ /= factor1;
    mat1.xz_ /= factor1;
    mat1.yx_ /= factor1;
    mat1.yy_ /= factor1;
    mat1.yz_ /= factor1;
    mat1.zx_ /= factor1;
    mat1.zy_ /= factor1;
    mat1.zz_ /= factor1;
    return mat1;
}

Matrix3D operator*(const Matrix3D& mat1, const Matrix3D& mat2)
{
    return matrix_product(mat1, mat2);
}

Matrix3D scalar_product(const Matrix3D& mat1, const Matrix3D& mat2)
{

    Matrix3D matrix_product;
    matrix_product.xx_ = mat1.xx_ * mat2.xx_;
    matrix_product.xy_ = mat1.xy_ * mat2.xy_;
    matrix_product.xz_ = mat1.xz_ * mat2.xz_;
    matrix_product.yx_ = mat1.yx_ * mat2.yx_;
    matrix_product.yy_ = mat1.yy_ * mat2.yy_;
    matrix_product.yz_ = mat1.yz_ * mat2.yz_;
    matrix_product.zx_ = mat1.zx_ * mat2.zx_;
    matrix_product.zy_ = mat1.zy_ * mat2.zy_;
    matrix_product.zz_ = mat1.zz_ * mat2.zz_;
    return matrix_product;
}

Matrix3D matrix_product(const Matrix3D& mat1, const Matrix3D& mat2)
{
    Matrix3D p;
    p.xx_ = mat1.xx_ * mat2.xx_ +
            mat1.xy_ * mat2.yx_ +
            mat1.xz_ * mat2.zx_;
    p.xx_ = mat1.xx_ * mat2.xy_ +
            mat1.xy_ * mat2.yy_ +
            mat1.xz_ * mat2.zy_;
    p.xz_ = mat1.xx_ * mat2.xz_ +
            mat1.xy_ * mat2.yz_ +
            mat1.xz_ * mat2.zz_;
    p.yx_ = mat1.yx_ * mat2.xx_ +
            mat1.yy_ * mat2.yx_ +
            mat1.yz_ * mat2.zx_;
    p.yx_ = mat1.yx_ * mat2.xy_ +
            mat1.yy_ * mat2.yy_ +
            mat1.yz_ * mat2.zy_;
    p.yz_ = mat1.yx_ * mat2.xz_ +
            mat1.yy_ * mat2.yz_ +
            mat1.yz_ * mat2.zz_;
    p.zx_ = mat1.zx_ * mat2.xx_ +
            mat1.zy_ * mat2.yx_ +
            mat1.zz_ * mat2.zx_;
    p.zx_ = mat1.zx_ * mat2.xy_ +
            mat1.zy_ * mat2.yy_ +
            mat1.zz_ * mat2.zy_;
    p.zz_ = mat1.zx_ * mat2.xz_ +
            mat1.zy_ * mat2.yz_ +
            mat1.zz_ * mat2.zz_;
    return p;
}

} // namespace math
} // namespace SI

Matrix3D Matrix3D::operator-() const
{
    Matrix3D matrix_3d;
    matrix_3d.xx_ = -xx_;
    matrix_3d.xy_ = -xy_;
    matrix_3d.xz_ = -xz_;
    matrix_3d.yx_ = -yx_;
    matrix_3d.yy_ = -yy_;
    matrix_3d.yz_ = -yz_;
    matrix_3d.zx_ = -zx_;
    matrix_3d.zy_ = -zy_;
    matrix_3d.zz_ = -zz_;
    return matrix_3d;
}

double const & Matrix3D::operator[](std::initializer_list<unsigned int> index) const {
    unsigned int row = *index.begin();
    unsigned int column = *(index.begin()+1);
    unsigned int i = ((row & 3) << 2) | (column & 3);
    switch(i) {
        case 0:
            return xx_; break;
        case 1:
            return xy_; break;
        case 2:
            return xz_; break;
        case 4:
            return yx_; break;
        case 5:
            return yy_; break;
        case 6:
            return yz_; break;
        case 8:
            return zx_; break;
        case 9:
            return zy_; break;
        case 10:
            return zz_; break;
        default:
            throw(std::runtime_error("Out of bounds!"));
    }
}

double & Matrix3D::operator[](std::initializer_list<unsigned int> index) {
    return const_cast<double &>(const_cast<const Matrix3D*>(this)->operator[](index));
}

//----------------------------------------------------------------------//
//---------------------private member function--------------------------//
//----------------------------------------------------------------------//
