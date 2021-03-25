#ifndef LI_Polynomial_H
#define LI_Polynomial_H

/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#include <functional>
#include <vector>
#include <fstream>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace earthmodel {

class Polynom {
   public:
    Polynom(const std::vector<double>& coefficients);
    Polynom(const Polynom&);

    ~Polynom() {};
    Polynom& operator=(const Polynom& poly) = default;

    Polynom* clone() const { return new Polynom(*this); };

    bool operator==(const Polynom& polynom) const;
    bool operator!=(const Polynom& polynom) const;

    double evaluate(double x) const;
    void shift(double x);
    void scale(double x);

    Polynom GetDerivative() const;
    Polynom GetAntiderivative(double constant) const;
    std::vector<double> GetCoefficient() const;

    friend std::ostream& operator<<(std::ostream& os, const Polynom& p);

    std::function<double(double)> GetFunction();

   protected:
    int N_;
    std::vector<double> coeff_;
};

class MathException: public std::exception {
    public:
        MathException(char const* message) : message_(message) {};
        const char* what() const throw()
        {
            return message_.c_str();
        }
    protected:
        std::string message_;
};

/// @brief Netwon-Raphson method and bisection to find the root of the function f
/// @param f    Function to find the root. Root must exist and be inside the interval [x1, x2]
/// @param df   Derivative of f
/// @param x1   lower limit of x to find the root
/// @param x2   upper limit of x to find the root
/// @param xacc convergence criterion: if $ \frac{f(x)}{df(x)} < xacc $ than accept x as the root
/// @return root of x

double NewtonRaphson(std::function<double(double)> f, std::function<double(double)> df, double x1, double x2,
        double xinit, int MAX_STEPS = 101, double xacc = 1.e-6);

}  // namespace earthmodel

#endif // LI_Polynomial_H

