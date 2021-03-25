#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>

#include "earthmodel-service/Polynomial.h"
#include "earthmodel-service/Geometry.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace earthmodel;

Polynom::Polynom(const std::vector<double>& coefficients) : N_(coefficients.size()) {
    coeff_ = coefficients;
}

Polynom::Polynom(const Polynom& poly) : N_(poly.N_), coeff_(poly.coeff_) {}

bool Polynom::operator==(const Polynom& polynom) const {
    if (N_ != polynom.N_)
        return false;
    for (int i = 0; i < N_; ++i) {
        if(coeff_[i] != polynom.coeff_[i])
            return false;
    }
    return true;
}

bool Polynom::operator!=(const Polynom& polynom) const {
    return !(*this == polynom);
}

double Polynom::evaluate(double x) const {
    if(N_ == 0) {
        return 0.0;
    }
    double aux = coeff_[N_ - 1];

    for (int i = N_ - 2; i >= 0; --i)
        aux = aux * x + coeff_[i];

    return aux;
}

void Polynom::shift(double x) {
    // Shaw and Traub method for the Taylor shift
    // https://planetcalc.com/7726/#fnref1:shaw

    std::function<int(int)> rev = [=] (int i) -> int {
        return N_-1 - i;
    };
    if (std::fabs(x) > Geometry::GEOMETRY_PRECISION) {
        int n = N_ - 1;
        double** t = new double*[N_];
        for (int count = 0; count < N_; ++count)
            t[count] = new double[N_];

        for (int i = 0; i < n; ++i) {
            t[i][0] = coeff_[n - i - 1] * std::pow(x, n - i - 1);
            t[i][i + 1] = coeff_[n] * std::pow(x, n);
        }

        for (int j = 0; j < n; ++j) {
            for (int i = j + 1; i <= n; ++i) {
                t[i][j + 1] = t[i - 1][j] + t[i - 1][j + 1];
            }
        }

        for (int i = 0; i < n; ++i) {
            coeff_[i] = t[n][i + 1] / std::pow(x, i);
        }

        for (int count = 0; count < N_; ++count)
            delete t[count];
    }
}

void Polynom::scale(double x) {
    for(int i = 0; i<N_; ++i) {
        coeff_[i] *= std::pow(x, i);
    }
}

Polynom Polynom::GetDerivative() const {
    std::vector<double> derivative_coeff_;

    for (auto i = 1; i < N_; ++i)
        derivative_coeff_.push_back(coeff_[i] * i);

    return Polynom(derivative_coeff_);
}

Polynom Polynom::GetAntiderivative(double constant) const {
    std::vector<double> derivative_coeff_{constant};

    for (auto i = 0; i < N_; ++i)
        derivative_coeff_.push_back(coeff_[i] / (i + 1));

    return Polynom(derivative_coeff_);
}

std::vector<double> Polynom::GetCoefficient() const {
    return coeff_;
}

std::function<double(double)> Polynom::GetFunction() {
    return (std::function<double(double)>)std::bind(&Polynom::evaluate, this,
                                                    std::placeholders::_1);
}

namespace earthmodel {
std::ostream& operator<<(std::ostream& os, const Polynom& p) {
    os << "p(x) =";
    for (int i = 0; i < p.N_; ++i) {
        if (p.coeff_[i] != 0) {
            if (!std::signbit(p.coeff_[i]))
                os << "+";
            os << p.coeff_[i] << "*x^{" << i << "}";
        }
    }
    return os;
}

double NewtonRaphson(std::function<double(double)> func,
                     std::function<double(double)> dfunc,
                     double x1,
                     double x2,
                     double xinit,
                     int MAX_STEPS,
                     double xacc) {
    /*
     * Method adapted from rtsafe method from:
     * Numerical Recipes in C (2Nd Ed.): The Art of Scientific Computing, 1992
     * Press, William H. and Teukolsky, Saul A. and Vetterling, William T. and
     * Flannery, Brian P.
     */

    double df, dx, dxold, f, fh, fl;
    double temp, xh, xl, rts;

    fl = func(x1);
    fh = func(x2);

    if (fl * fh > 0.0) {
        throw MathException("Root must be bracketed in NewtonRaphson method!");
        // log_error("Root must be bracketed in NewtonRaphson method!");
        // return 0;
    }

    if (fl == 0.0) {
        return x1;
    }
    if (fh == 0.0) {
        return x2;
    }

    if (fl < 0.0) {
        // swap orientation of search
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }

    rts = xinit;  // initial guess for root
    dxold = std::abs(x2 - x1);
    dx = dxold;

    f = func(rts);
    df = dfunc(rts);

    for (int j = 1; j < MAX_STEPS; j++) {
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
            (std::abs(2.0 * f) > std::abs(dxold * df))) {
            // use bisection if Newton method produces an x that is out of range
            // or if Newton method is not decreasing fast enough
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) {
                return rts;
            }
        } else {
            // use Newton step
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) {
                return rts;
            }
        }

        if (std::abs(dx) < xacc) {
            // convergence criterion reached
            return rts;
        }
        f = func(rts);
        df = dfunc(rts);
        // update bracket
        if (f < 0.0) {
            xl = rts;
        } else {
            xh = rts;
        };
    }
    //log_warn("Maximum number of iteration exeeded in NewtonRaphson");
    return rts;
}

}  // namespace earthmodel

