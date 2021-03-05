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

#pragma once
#include <exception>
#include <functional>
#include <string>
#include <memory>
#include <earthmodel-service/Vector3D.h>
#include <earthmodel-service/Polynomial.h>

namespace earthmodel {
class Axis {
   public:
    Axis();
    Axis(const Vector3D& fAxis, const Vector3D& fp0);
    Axis(const Axis&);

    virtual ~Axis() {};

    bool operator==(const Axis& axis) const;
    bool operator!=(const Axis& axis) const;

    virtual Axis* clone() const = 0;
    virtual std::shared_ptr<const Axis> create() const = 0;

    virtual double GetDepth(const Vector3D& xi) const = 0;
    virtual double GetEffectiveDistance(const Vector3D& xi,
                                        const Vector3D& direction) const = 0;

    Vector3D GetAxis() const { return fAxis_; };
    Vector3D GetFp0() const { return fp0_; };

   protected:
    Vector3D fAxis_;
    Vector3D fp0_;
};

class RadialAxis : public Axis {
   public:
    RadialAxis();
    RadialAxis(const Vector3D& fAxis, const Vector3D& fp0);

    Axis* clone() const override { return new RadialAxis(*this); };
    virtual std::shared_ptr<const Axis> create() const override {
        return std::shared_ptr<const Axis>(new RadialAxis(*this));
    };
    ~RadialAxis() {};

    double GetDepth(const Vector3D& xi) const override;
    double GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const override;
};

class CartesianAxis : public Axis {
   public:
    CartesianAxis();
    CartesianAxis(const Vector3D& fAxis, const Vector3D& fp0);
    ~CartesianAxis() {};

    Axis* clone() const override { return new CartesianAxis(*this); };
    virtual std::shared_ptr<const Axis> create() const override {
        return std::shared_ptr<const Axis>(new CartesianAxis(*this));
    };

    double GetDepth(const Vector3D& xi) const override;
    double GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const override;
};

class DensityException : public std::exception {
   public:
    DensityException(const char* m) : message_(m){};
    const char* what() const throw() { return message_.c_str(); };

   private:
    std::string message_;
};

class DensityDistribution {
   public:
    DensityDistribution();
    DensityDistribution(const DensityDistribution&);

    virtual bool operator==(const DensityDistribution& dens_distr) const;
    virtual bool operator!=(const DensityDistribution& dens_distr) const;
    virtual bool compare(const DensityDistribution& dens_distr) const = 0;


    virtual DensityDistribution* clone() const = 0;
    virtual std::shared_ptr<const DensityDistribution> create() const = 0;

    virtual double Derivative(const Vector3D& xi,
                              const Vector3D& direction) const = 0;
    virtual double AntiDerivative(const Vector3D& xi,
                                  const Vector3D& direction) const = 0;
    virtual double Integral(const Vector3D& xi,
                            const Vector3D& direction,
                            double distance) const = 0;
    virtual double Integral(const Vector3D& xi,
                            const Vector3D& xj) const = 0;
    virtual double InverseIntegral(const Vector3D& xi,
                                   const Vector3D& direction,
                                   double integral,
                                   double max_distance) const =0;
    virtual double Evaluate(const Vector3D& xi) const = 0;
};

class Density_homogeneous : public DensityDistribution {
   public:
    Density_homogeneous();
    Density_homogeneous(const Density_homogeneous&);
    Density_homogeneous(double correction_factor);

    // ~Density_homogeneous(){};

    bool compare(const DensityDistribution& dens_distr) const override;

    DensityDistribution* clone() const override {
        return new Density_homogeneous(*this);
    };

    std::shared_ptr<const DensityDistribution> create() const override {
        return std::shared_ptr<const DensityDistribution>(new Density_homogeneous(*this));
    };

    double Derivative(const Vector3D& xi,
                      const Vector3D& direction) const override;
    double AntiDerivative(const Vector3D& xi,
                          const Vector3D& direction) const override;
    double Integral(const Vector3D& xi,
                    const Vector3D& direction,
                    double distance) const override;
    double Integral(const Vector3D& xi,
                    const Vector3D& xj) const override;
    double InverseIntegral(const Vector3D& xi,
                           const Vector3D& direction,
                           double integral,
                           double max_distance) const override;
    double Evaluate(const Vector3D& xi) const override;

    double GetInverseIntegralionfactor() const { return correction_factor_; }

   private:
    double correction_factor_;
   protected:
    Axis* axis_;
};

class Density_polynomial : public DensityDistribution {
   public:
    Density_polynomial(const Axis&, const Polynom&);
    Density_polynomial(const Density_polynomial&);
    ~Density_polynomial();

    bool compare(const DensityDistribution& dens_distr) const override;

    DensityDistribution* clone() const override {
        return new Density_polynomial(*this);
    };
    std::shared_ptr<const DensityDistribution> create() const override { return std::shared_ptr<const DensityDistribution>( new Density_polynomial(*this) ); };

    double Derivative(const Vector3D& xi,
                      const Vector3D& direction) const override;
    double AntiDerivative(const Vector3D& xi,
                          const Vector3D& direction) const override;
    double Integral(const Vector3D& xi,
                    const Vector3D& direction,
                    double distance) const override;
    double Integral(const Vector3D& xi,
                    const Vector3D& xj) const override;
    double InverseIntegral(const Vector3D& xi,
                           const Vector3D& direction,
                           double integral,
                           double max_distance) const override;
    double Evaluate(const Vector3D& xi) const override;

    double Helper_function(const Vector3D& xi,
                           const Vector3D& direction,
                           double res,
                           double l) const;
    double helper_function(const Vector3D& xi,
                           const Vector3D& direction,
                           double res,
                           double l) const;

   protected:
    Polynom polynom_;
    Polynom Polynom_;

    std::function<double(double)> density_distribution;
    std::function<double(double)> antiderived_density_distribution;
   protected:
    Axis* axis_;
};

class Density_exponential : public DensityDistribution {
   public:
    Density_exponential(const Axis& axis, double sigma);

    ~Density_exponential(){};

    bool compare(const DensityDistribution& dens_distr) const override;

    DensityDistribution* clone() const override {
        return new Density_exponential(*this);
    };
    std::shared_ptr<const DensityDistribution> create() const override { return std::shared_ptr<const DensityDistribution>( new Density_exponential(*this) ); };

    double Derivative(const Vector3D& xi,
                      const Vector3D& direction) const override;
    double AntiDerivative(const Vector3D& xi,
                          const Vector3D& direction) const override;
    double Integral(const Vector3D& xi,
                    const Vector3D& direction,
                    double distance) const override;
    double Integral(const Vector3D& xi,
                    const Vector3D& xj) const override;
    double InverseIntegral(const Vector3D& xi,
                           const Vector3D& direction,
                           double integral,
                           double max_distance) const override;
    double Evaluate(const Vector3D& xi) const override;

    double GetDepth(const Vector3D& xi) const;
    double GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const;

   private:
    double sigma_;
   protected:
    Axis* axis_;
};

}  // namespace earthmodel
