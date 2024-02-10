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

// Record intended to be passed to the primary particle injection distributions
class PrimaryDistributionRecord {
public:
    // The signature should be constant
    InteractionRecord const & record;
    InteractionSignature const & signature;
    ParticleID const id;
    ParticleType const & type;
private:
    // The rest of the primary particle properties can be changed
    bool mass_set = false;
    bool energy_set = false;
    bool kinetic_energy_set = false;
    bool direction_set = false;
    bool momentum_set = false;
    bool length_set = false;
    bool initial_position_set = false;
    bool interaction_vertex_set = false;
    bool helicity_set = false;

    double mass;
    double energy;
    double kinetic_energy;
    std::array<double, 3> direction;
    std::array<double, 3> momentum;
    double length;
    std::array<double, 3> initial_position;
    std::array<double, 3> interaction_vertex;
    double helicity = 0;
public:

    PrimaryDistributionRecord(InteractionRecord const & record);

    InteractionSignature const & GetSignature() const;

    Particle GetParticle();
    void SetParticle(Particle const & particle);

    ParticleID const & GetID() const;
    ParticleType const & GetType() const;

    double GetMass() const;
    double GetEnergy() const;
    double GetKineticEnergy() const;
    std::array<double, 3> GetDirection() const;
    std::array<double, 3> GetThreeMomentum() const;
    std::array<double, 4> GetFourMomentum() const;
    double GetLength() const;
    std::array<double, 3> GetInitialPosition() const;
    std::array<double, 3> GetInteractionVertex() const;

    double const & GetMass();
    double const & GetEnergy();
    double const & GetKineticEnergy();
    std::array<double, 3> const & GetDirection();
    std::array<double, 3> const & GetThreeMomentum();
    std::array<double, 4> GetFourMomentum();
    double const & GetLength();
    std::array<double, 3> const & GetInitialPosition();
    std::array<double, 3> const & GetInteractionVertex();
    double const & GetHelicity() const;

    void SetMass(double mass);
    void SetEnergy(double energy);
    void SetKineticEnergy(double kinetic_energy);
    void SetDirection(std::array<double, 3> direction);
    void SetThreeMomentum(std::array<double, 3> momentum);
    void SetFourMomentum(std::array<double, 4> momentum);
    void SetLength(double length);
    void SetInitialPosition(std::array<double, 3> initial_position);
    void SetInteractionVertex(std::array<double, 3> initial_position);
    void SetHelicity(double helicity);

    void UpdateMass();
    void UpdateEnergy();
    void UpdateKineticEnergy();
    void UpdateDirection();
    void UpdateMomentum();
    void UpdateLength();
    void UpdateInitialPosition();
    void UpdateInteractionVertex();

    void Finalize(InteractionRecord & record);
};

class SecondaryParticleRecord {
public:
    ParticleID const id;
    ParticleType const & type;
    std::array<double, 3> const & initial_position;
private:
    // The rest of the primary particle properties can be changed
    bool mass_set = false;
    bool energy_set = false;
    bool kinetic_energy_set = false;
    bool direction_set = false;
    bool momentum_set = false;
    bool helicity_set = false;

    double mass = 0;
    double energy = 0;
    double kinetic_energy = 0;
    std::array<double, 3> direction = {0, 0, 0};
    std::array<double, 3> momentum = {0, 0, 0};
    double helicity = 0;
public:

    SecondaryParticleRecord(InteractionRecord const & record, size_t secondary_index);

    InteractionSignature const & GetSignature() const;

    Particle GetParticle();
    void SetParticle(Particle const & particle);

    ParticleID const & GetID() const;

    ParticleType const & GetType() const;
    double GetMass() const;
    double GetEnergy() const;
    double GetKineticEnergy() const;
    std::array<double, 3> GetDirection() const;
    std::array<double, 3> GetThreeMomentum() const;
    std::array<double, 4> GetFourMomentum() const;
    std::array<double, 3> const & GetInitialPosition() const;
    double const & GetHelicity() const;

    double const & GetMass();
    double const & GetEnergy();
    double const & GetKineticEnergy();
    std::array<double, 3> const & GetDirection();
    std::array<double, 3> const & GetThreeMomentum();
    std::array<double, 4> GetFourMomentum();

    void SetMass(double mass);
    void SetEnergy(double energy);
    void SetKineticEnergy(double kinetic_energy);
    void SetDirection(std::array<double, 3> direction);
    void SetThreeMomentum(std::array<double, 3> momentum);
    void SetFourMomentum(std::array<double, 4> momentum);
    void SetHelicity(double helicity);

    void UpdateMass();
    void UpdateEnergy();
    void UpdateKineticEnergy();
    void UpdateDirection();
    void UpdateMomentum();

    void Finalize(InteractionRecord & record);
};

class CrossSectionDistributionRecord {
public:
    InteractionRecord const & record;
    InteractionSignature const & signature;
    ParticleID const & primary_id;
    ParticleType const & primary_type;
    std::array<double, 3> const & primary_initial_position;
    double const & primary_mass;
    std::array<double, 4> const & primary_momentum;
    double const & primary_helicity;
    std::array<double, 3> const & interaction_vertex;

    ParticleID const target_id;
    ParticleType const & target_type;

    double target_mass = 0;
    double target_helicity = 0;
    std::map<std::string, double> interaction_parameters;
private:
    std::vector<SecondaryParticleRecord> secondary_particles;
public:
    CrossSectionDistributionRecord(InteractionRecord const & record);

    InteractionSignature const & GetSignature() const;
    ParticleID const & GetPrimaryID() const;
    ParticleType const & GetPrimaryType() const;
    std::array<double, 3> const & GetPrimaryInitialPosition() const;
    double const & GetPrimaryMass() const;
    std::array<double, 4> const & GetPrimaryMomentum() const;
    double const & GetPrimaryHelicity() const;
    Particle GetPrimaryParticle() const;

    std::array<double, 3> const & GetInteractionVertex() const;
    ParticleID const & GetTargetID() const;
    ParticleType const & GetTargetType() const;

    double const & GetTargetMass() const;
    double const & GetTargetHelicity() const;
    double & GetTargetMass();
    double & GetTargetHelicity();
    std::map<std::string, double> & GetInteractionParameters();

    void SetTargetMass(double mass);
    void SetTargetHelicity(double helicity);
    void SetInteractionParameters(std::map<std::string, double> const & parameters);
    void SetInteractionParameter(std::string const & name, double value);

    Particle GetTargetParticle() const;
    void SetTargetParticle(Particle const & particle);

    SecondaryParticleRecord & GetSecondaryParticleRecord(size_t index);

    Particle GetSecondaryParticle(size_t index);
    std::vector<Particle> GetSecondaryParticles();

    void SetSecondaryParticle(size_t index, Particle const & particle);
    void SetSecondaryParticles(std::vector<Particle> const & particles);

    void Finalize(InteractionRecord & record);
};

class SecondaryDistributionRecord {
public:
    InteractionRecord const & parent_record;
    size_t const secondary_index;

    ParticleID const id;
    ParticleType const & type;
    double const & mass;
    std::array<double, 3> direction;
    std::array<double, 4> const & momentum;
    double const & helicity;
    std::array<double, 3> const & initial_position;
private:
    double length;
public:
    SecondaryDistributionRecord(InteractionRecord const & parent_record, size_t secondary_index);

    void SetLength(double const & length);
    double & GetLength();

    void Finalize(InteractionRecord & record);
};

class InteractionRecord {
public:
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

    bool operator==(InteractionRecord const & other) const;
    bool operator<(InteractionRecord const & other) const;
    friend std::ostream& ::operator<<(std::ostream& os, InteractionRecord const& record);

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
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
        } else {
            throw std::runtime_error("InteractionRecord only supports version <= 0!");
        }
    }
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::InteractionRecord, 0);

#endif // LI_InteractionRecord_H
