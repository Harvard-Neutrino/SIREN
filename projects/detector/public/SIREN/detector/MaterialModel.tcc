#pragma once
#ifndef SIREN_MaterialModel_TCC
#define SIREN_MaterialModel_TCC

#include "SIREN/detector/MaterialModel.h"

#include <cmath>

namespace siren {
namespace detector {

template<typename Iterator, class>
std::vector<double> MaterialModel::GetTargetMassFraction(int material_id, Iterator begin, Iterator end) const {
    std::vector<double> fractions;
    fractions.reserve(std::distance(begin, end));

    for(Iterator it = begin; it != end; ++it) {
        std::pair<int, siren::dataclasses::ParticleType> key(material_id, *it);
        if(material_components_by_id_.find(key) != material_components_by_id_.end())
            fractions.push_back(material_components_by_id_.at(key).mass_density_over_total_mass_density);
        else
            fractions.push_back(0.0);
    }
    return fractions;
}

template<typename Iterator, class>
std::vector<double> MaterialModel::GetTargetParticleFraction(int material_id, Iterator begin, Iterator end) const {
    std::vector<double> fractions;
    fractions.reserve(std::distance(begin, end));

    for(Iterator it = begin; it != end; ++it) {
        std::pair<int, siren::dataclasses::ParticleType> key(material_id, *it);
        if(material_components_by_id_.find(key) != material_components_by_id_.end())
            fractions.push_back(material_components_by_id_.at(key).particle_density_over_total_mass_density);
        else
            fractions.push_back(0.0);
    }
    return fractions;
}

template<typename Iterator, class>
std::vector<double> MaterialModel::GetTargetRadiationFraction(int material_id, Iterator begin, Iterator end) const {
    double X0inv = 0;
    std::vector<double> fractions;
    fractions.reserve(std::distance(begin, end));
    for(Iterator it = begin; it != end; ++it) {
        std::pair<int, siren::dataclasses::ParticleType> key(material_id, *it);
        if(material_components_by_id_.find(key) != material_components_by_id_.end()) {
            fractions.push_back(0.0);
            continue;
        }
        MaterialComponent const & component = material_components_by_id_.at(key);
        if(not component.component.is_atom) {
            fractions.push_back(0.0);
            continue;
        }
        int A = component.component.nucleon_count;
        int Z = component.component.proton_count;
        double f = component.mass_density_over_total_mass_density;
        double X0i = 716.4 * A / ( Z*(Z + 1) * std::log(287./std::sqrt(Z))); // g/cm^2, Grupen eq 1.59
        double frac = f / X0i;
        fractions.push_back(frac);
        X0inv += frac;
    }
    for(unsigned int i=0; i<fractions.size(); ++i) {
        fractions[i] /= X0inv;
    }
    return fractions;
}

} // namespace detector
} // namespace siren

#endif
