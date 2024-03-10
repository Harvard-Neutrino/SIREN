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
#ifndef SIREN_DensityDistribution_H
#define SIREN_DensityDistribution_H

#include <string>                          // for basic_string
#include <memory>                          // for shared_ptr
#include <cstdint>                         // for uint32_t
#include <exception>                       // for exception

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/math/Vector3D.h"

namespace siren {
namespace detector {

class DensityException : public std::exception {
   public:
    DensityException(const char* m) : message_(m){};
    const char* what() const throw() { return message_.c_str(); };

   private:
    std::string message_;
};

class DensityDistribution {
friend cereal::access;
public:
    DensityDistribution();
    virtual ~DensityDistribution() = default;
    DensityDistribution(const DensityDistribution&);

    bool operator==(const DensityDistribution& dens_distr) const;
    bool operator!=(const DensityDistribution& dens_distr) const;
    virtual bool compare(const DensityDistribution& dens_distr) const = 0;


    virtual DensityDistribution* clone() const = 0;
    virtual std::shared_ptr<DensityDistribution> create() const = 0;

    virtual double Derivative(const math::Vector3D& xi,
                              const math::Vector3D& direction) const = 0;
    virtual double AntiDerivative(const math::Vector3D& xi,
                                  const math::Vector3D& direction) const = 0;
    virtual double Integral(const math::Vector3D& xi,
                            const math::Vector3D& direction,
                            double distance) const = 0;
    virtual double Integral(const math::Vector3D& xi,
                            const math::Vector3D& xj) const = 0;
    virtual double InverseIntegral(const math::Vector3D& xi,
                                   const math::Vector3D& direction,
                                   double integral,
                                   double max_distance) const = 0;
    virtual double InverseIntegral(const math::Vector3D& xi,
                                   const math::Vector3D& direction,
                                   double constant,
                                   double integral,
                                   double max_distance) const = 0;
    virtual double Evaluate(const math::Vector3D& xi) const = 0;

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {};
};

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::DensityDistribution, 0);

#endif // SIREN_DensityDistribution_H

