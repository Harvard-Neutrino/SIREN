#pragma once
#ifndef LI_InteractionRecord_H
#define LI_InteractionRecord_H

#include <array>                                              // for array
#include <iosfwd>                                             // for ostream
#include <vector>                                             // for vector
#include <cstdint>                                            // for uint32_t
#include <stdexcept>                                          // for runtime...

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/ParticleID.h"
#include "LeptonInjector/dataclasses/ParticleType.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"

namespace LI { namespace dataclasses { struct InteractionRecord; } }

std::ostream& operator<<(std::ostream& os, LI::dataclasses::InteractionRecord const& record);

namespace LI {
namespace dataclasses {

class InteractionRecord {
private:
    // Stateful information
    bool signature_set = false;
    InteractionSignature reference_signature;

    // Data
    InteractionSignature signature;
    ParticleID primary_id;
    std::array<double, 3> primary_initial_position = {0, 0, 0};
    double primary_mass = 0;
    std::array<double, 4> primary_momentum = {0, 0, 0, 0};
    double primary_helicity = 0;
    ParticleID target_id;
    double target_mass = 0;
    std::array<double, 4> target_momentum = {0, 0, 0, 0};
    double target_helicity = 0;
    std::array<double, 3> interaction_vertex = {0, 0, 0};
    std::vector<ParticleID> secondary_ids;
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicity;
    std::map<std::string, double> interaction_parameters;
public:
    // Constructors
    InteractionRecord() = default;
    InteractionRecord(InteractionSignature const & signature);
    InteractionRecord & operator=(InteractionRecord const & record) = default;

    // Getters
    InteractionSignature const & GetSignature() const;

    Particle GetPrimary() const;
    ParticleID const & GetPrimaryID() const;
    std::array<double, 3> const & GetPrimaryInitialPosition() const;
    double const & GetPrimaryMass() const;
    std::array<double, 4> const & GetPrimaryMomentum() const;
    double const & GetPrimaryHelicity() const;

    Particle GetTarget() const;
    ParticleID const & GetTargetID() const;
    double const & GetTargetMass() const;
    std::array<double, 4> const & GetTargetMomentum() const;
    double const & GetTargetHelicity() const;
    std::array<double, 3> const & GetInteractionVertex() const;

    std::vector<Particle> GetSecondaries() const;
    std::vector<ParticleID> const & GetSecondaryIDs() const;
    std::vector<double> const & GetSecondaryMasses() const;
    std::vector<std::array<double, 4>> const & GetSecondaryMomenta() const;
    std::vector<double> const & GetSecondaryHelicity() const;
    std::map<std::string, double> const & GetInteractionParameters() const;

    // Setters
    void SetSignature(InteractionSignature const & signature);

    ParticleID SetPrimary(Particle & primary);
    ParticleID SetPrimary(Particle const & primary);
    ParticleID SetPrimaryID(ParticleID const & primary_id);
    void SetPrimaryType(Particle::ParticleType const & primary_type);
    void SetPrimaryInitialPosition(std::array<double, 3> const & primary_initial_position);
    void SetPrimaryMass(double const & primary_mass);
    void SetPrimaryMomentum(std::array<double, 4> const & primary_momentum);
    void SetPrimaryHelicity(double const & primary_helicity);

    ParticleID SetTarget(Particle & target);
    ParticleID SetTarget(Particle const & target);
    ParticleID SetTargetID(ParticleID & target_id);
    void SetTargetType(Particle::ParticleType const & target_type);
    void SetTargetMass(double const & target_mass);
    void SetTargetMomentum(std::array<double, 4> const & target_momentum);
    void SetTargetHelicity(double const & target_helicity);
    void SetInteractionVertex(std::array<double, 3> const & interaction_vertex);

    std::vector<ParticleID> SetSecondaries(std::vector<Particle> & secondaries);
    std::vector<ParticleID> SetSecondaries(std::vector<Particle> const & secondaries);
    void SetSecondaryIDs(std::vector<ParticleID> const & secondary_ids);
    void SetSecondaryTypes(std::vector<Particle::ParticleType> const & secondary_types);
    void SetSecondaryMasses(std::vector<double> const & secondary_masses);
    void SetSecondaryMomenta(std::vector<std::array<double, 4>> const & secondary_momenta);
    void SetSecondaryHelicity(std::vector<double> const & secondary_helicity);
    void SetInteractionParameters(std::map<std::string, double> const & interaction_parameters);

    ParticleID SetSecondary(size_t const & index, Particle & secondary);
    ParticleID SetSecondary(size_t const & index, Particle const & secondary);
    void SetSecondaryID(size_t const & index, ParticleID const & secondary_id);
    void SetSecondaryType(size_t const & index, Particle::ParticleType const & secondary_type);
    void SetSecondaryMass(size_t const & index, double const & secondary_mass);
    void SetSecondaryMomentum(size_t const & index, std::array<double, 4> const & secondary_momentum);
    void SetSecondaryHelicity(size_t const & index, double const & secondary_helicity);
    void SetInteractionParameter(std::string const & key, double const & value);

    // Adders
    ParticleID AddSecondary();
    ParticleID AddSecondary(ParticleID const & secondary_id);
    ParticleID AddSecondary(Particle & secondary);
    ParticleID AddSecondary(Particle const & secondary);
    ParticleID AddSecondary(Particle::ParticleType const & secondary_type, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity);
    ParticleID AddSecondary(ParticleID const & secondary_id, Particle::ParticleType const & secondary_type, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity);
    ParticleID AddSecondary(double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity);
    ParticleID AddSecondary(ParticleID const & secondary_id, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity);

    bool CheckSignature() const;

    bool operator==(InteractionRecord const & other) const;
    bool operator<(InteractionRecord const & other) const;
    friend std::ostream& ::operator<<(std::ostream& os, InteractionRecord const& record);

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            CheckSignature();
            archive(::cereal::make_nvp("InteractionSignature", signature));
            archive(::cereal::make_nvp("PrimaryID", primary_id));
            archive(::cereal::make_nvp("PrimaryInitialPosition", primary_initial_position));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("TargetID", target_id));
            archive(::cereal::make_nvp("TargetMass", target_mass));
            archive(::cereal::make_nvp("TargetMomentum", target_momentum));
            archive(::cereal::make_nvp("TargetHelicity", target_helicity));
            archive(::cereal::make_nvp("InteractionVertex", interaction_vertex));
            archive(::cereal::make_nvp("SecondaryIDs", secondary_ids));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("InteractionParameters", interaction_parameters));
        } else {
            throw std::runtime_error("InteractionRecord only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InteractionSignature", signature));
            archive(::cereal::make_nvp("PrimaryID", primary_id));
            archive(::cereal::make_nvp("PrimaryInitialPosition", primary_initial_position));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("TargetID", target_id));
            archive(::cereal::make_nvp("TargetMass", target_mass));
            archive(::cereal::make_nvp("TargetMomentum", target_momentum));
            archive(::cereal::make_nvp("TargetHelicity", target_helicity));
            archive(::cereal::make_nvp("InteractionVertex", interaction_vertex));
            archive(::cereal::make_nvp("SecondaryIDs", secondary_ids));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("InteractionParameters", interaction_parameters));
            signature_set = true;
            reference_signature = signature;
        } else {
            throw std::runtime_error("InteractionRecord only supports version <= 0!");
        }
    }
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::InteractionRecord, 0);

#endif // LI_InteractionRecord_H
