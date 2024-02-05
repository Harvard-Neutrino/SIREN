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

namespace LI { namespace dataclasses { class InteractionRecord; } }

std::ostream& operator<<(std::ostream& os, LI::dataclasses::InteractionRecord const& record);

namespace LI {
namespace dataclasses {

class PrimaryRecord {
private:
    ParticleType type;

    bool mass_set = false;
    bool energy_set = false;
    bool kinetic_energy_set = false;
    bool direction_set = false;
    bool momentum_set = false;
    bool length_set = false;
    bool initial_position_set = false;
    bool interaction_vertex_set = false;

    double mass;
    double energy;
    double kinetic_energy;
    std::array<double, 3> direction;
    std::array<double, 3> momentum;
    double length;
    std::array<double, 3> initial_position;
    std::array<double, 3> interaction_vertex;

public:
    PrimaryRecord(ParticleType const & type);

    PrimaryRecord & operator=(PrimaryRecord const & record);

    ParticleType const & GetType() const;
    double GetMass() const;
    double GetEnergy() const;
    double GetKineticEnergy() const;
    std::array<double, 3> GetDirection() const;
    std::array<double, 3> GetThreeMomentum() const;
    std::array<double, 4> GetFourMomentum() const;
    double GetLength() const;
    std::array<double, 3> GetIntialPosition() const;
    std::array<double, 3> GetInteractionVertex() const;

    double const & GetMass();
    double const & GetEnergy();
    double const & GetKineticEnergy();
    std::array<double, 3> const & GetDirection();
    std::array<double, 3> const & GetThreeMomentum();
    std::array<double, 4> GetFourMomentum();
    double const & GetLength();
    std::array<double, 3> const & GetIntialPosition();
    std::array<double, 3> const & GetInteractionVertex();

    //void SetType(ParticleType type);
    void SetMass(double mass);
    void SetEnergy(double energy);
    void SetKineticEnergy(double kinetic_energy);
    void SetDirection(std::array<double, 3> direction);
    void SetThreeMomentum(std::array<double, 3> momentum);
    void SetFourMomentum(std::array<double, 4> momentum);
    void SetLength(double length);
    void SetIntialPosition(std::array<double, 3> initial_position);
    void SetInteractionVertex(std::array<double, 3> initial_position);

    void UpdateMass();
    void UpdateEnergy();
    void UpdateKineticEnergy();
    void UpdateDirection();
    void UpdateMomentum();
    void UpdateLength();
    void UpdateIntialPosition();
    void UpdateInteractionVertex();
};

class SecondaryRecord : protected PrimaryRecord {
public:
    SecondaryRecord(ParticleType const & type, std::array<double, 3> initial_position, double mass, std::array<double, 4> momentum);

    SecondaryRecord & operator=(SecondaryRecord const & record);

    ParticleType const & GetType() const;
    double const & GetMass();
    double const & GetEnergy();
    double const & GetKineticEnergy();
    std::array<double, 3> const & GetDirection();
    std::array<double, 3> const & GetThreeMomentum();
    std::array<double, 4> GetFourMomentum();
    double const & GetLength();
    std::array<double, 3> const & GetIntialPosition();
    std::array<double, 3> const & GetInteractionVertex();

    void SetLength(double length);
};

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
    double target_helicity = 0;
    std::array<double, 3> interaction_vertex = {0, 0, 0};
    std::vector<ParticleID> secondary_ids;
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicities;
    std::map<std::string, double> interaction_parameters;
public:
    // Constructors
    InteractionRecord() = default;
    InteractionRecord(InteractionSignature const & signature);
    InteractionRecord(InteractionSignature const & signature, PrimaryRecord const & primary);
    InteractionRecord & operator=(InteractionRecord const & record) = default;

    // Special getters
    PrimaryRecord GetPrimaryRecord() const;

    // Special setters
    void SetPrimaryRecord(PrimaryRecord const & primary);

    // Getters
    InteractionSignature const & GetSignature() const;

    Particle GetPrimary() const;
    ParticleID const & GetPrimaryID() const;
    ParticleType const & GetPrimaryType() const;
    std::array<double, 3> const & GetPrimaryInitialPosition() const;
    double const & GetPrimaryMass() const;
    std::array<double, 4> const & GetPrimaryMomentum() const;
    double const & GetPrimaryHelicity() const;

    Particle GetTarget() const;
    ParticleID const & GetTargetID() const;
    ParticleType const & GetTargetType() const;
    double const & GetTargetMass() const;
    double const & GetTargetHelicity() const;
    std::array<double, 3> const & GetInteractionVertex() const;

    std::vector<Particle> GetSecondaries() const;
    std::vector<ParticleID> const & GetSecondaryIDs() const;
    std::vector<ParticleType> const & GetSecondaryTypes() const;
    std::vector<double> const & GetSecondaryMasses() const;
    std::vector<std::array<double, 4>> const & GetSecondaryMomenta() const;
    std::vector<double> const & GetSecondaryHelicities() const;
    std::map<std::string, double> const & GetInteractionParameters() const;

    Particle GetSecondary(size_t const & index) const;
    ParticleID const & GetSecondaryID(size_t const & index) const;
    ParticleType const & GetSecondaryType(size_t const & index) const;
    double const & GetSecondaryMass(size_t const & index) const;
    std::array<double, 4> const & GetSecondaryMomentum(size_t const & index) const;
    double const & GetSecondaryHelicity(size_t const & index) const;
    double const & GetInteractionParameter(std::string const & key) const;

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
    ParticleID SetTargetID(ParticleID const & target_id);
    void SetTargetType(Particle::ParticleType const & target_type);
    void SetTargetMass(double const & target_mass);
    void SetTargetHelicity(double const & target_helicity);
    void SetInteractionVertex(std::array<double, 3> const & interaction_vertex);

    std::vector<ParticleID> SetSecondaries(std::vector<Particle> & secondaries);
    std::vector<ParticleID> SetSecondaries(std::vector<Particle> const & secondaries);
    std::vector<ParticleID> SetSecondaryIDs(std::vector<ParticleID> const & secondary_ids);
    void SetSecondaryTypes(std::vector<Particle::ParticleType> const & secondary_types);
    void SetSecondaryMasses(std::vector<double> const & secondary_masses);
    void SetSecondaryMomenta(std::vector<std::array<double, 4>> const & secondary_momenta);
    void SetSecondaryHelicities(std::vector<double> const & secondary_helicities);
    void SetInteractionParameters(std::map<std::string, double> const & interaction_parameters);

    ParticleID SetSecondary(size_t const & index, Particle & secondary);
    ParticleID SetSecondary(size_t const & index, Particle const & secondary);
    ParticleID SetSecondaryID(size_t const & index, ParticleID const & secondary_id);
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
            archive(::cereal::make_nvp("TargetHelicity", target_helicity));
            archive(::cereal::make_nvp("InteractionVertex", interaction_vertex));
            archive(::cereal::make_nvp("SecondaryIDs", secondary_ids));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicities", secondary_helicities));
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
            archive(::cereal::make_nvp("TargetHelicity", target_helicity));
            archive(::cereal::make_nvp("InteractionVertex", interaction_vertex));
            archive(::cereal::make_nvp("SecondaryIDs", secondary_ids));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicities", secondary_helicities));
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
