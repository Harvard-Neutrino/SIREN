#pragma once
#ifndef SIREN_InteractionRecord_H
#define SIREN_InteractionRecord_H

#include <array>                                              // for array
#include <iosfwd>                                             // for ostream
#include <vector>                                             // for vector
#include <cstdint>                                            // for uint32_t
#include <stdexcept>                                          // for runtime...

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/serialization/array.h"

#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/ParticleID.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/dataclasses/InteractionSignature.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class PrimaryDistributionRecord; } }
namespace siren { namespace dataclasses { class SecondaryParticleRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { class SecondaryDistributionRecord; } }

std::ostream& operator<<(std::ostream& os, siren::dataclasses::InteractionRecord const & record);
std::ostream& operator<<(std::ostream& os, siren::dataclasses::PrimaryDistributionRecord const & record);
std::ostream& operator<<(std::ostream& os, siren::dataclasses::SecondaryParticleRecord const & record);
std::ostream& operator<<(std::ostream& os, siren::dataclasses::CrossSectionDistributionRecord const & record);
std::ostream& operator<<(std::ostream& os, siren::dataclasses::SecondaryDistributionRecord const & record);

namespace siren {
namespace dataclasses {

// Record intended to be passed to the primary particle injection distributions
class PrimaryDistributionRecord {
public:
    // The id and type should be constant
    ParticleID const id;
    ParticleType const type;
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

    mutable double mass;
    mutable double energy;
    mutable double kinetic_energy;
    mutable std::array<double, 3> direction;
    mutable std::array<double, 3> momentum;
    mutable double length;
    mutable std::array<double, 3> initial_position;
    mutable std::array<double, 3> interaction_vertex;
    mutable double helicity = 0;
public:
    friend std::ostream& ::operator<<(std::ostream& os, PrimaryDistributionRecord const& record);

    PrimaryDistributionRecord(PrimaryDistributionRecord const & other) = delete;
    PrimaryDistributionRecord(PrimaryDistributionRecord && other) = default;
    PrimaryDistributionRecord & operator=(PrimaryDistributionRecord const & other) = delete;
    PrimaryDistributionRecord & operator=(PrimaryDistributionRecord && other) = delete;

    PrimaryDistributionRecord(ParticleType type);

    Particle GetParticle() const;
    void SetParticle(Particle const & particle);

    ParticleID const & GetID() const;
    ParticleType const & GetType() const;

    double const & GetMass() const;
    double const & GetEnergy() const;
    double const & GetKineticEnergy() const;
    std::array<double, 3> const & GetDirection() const;
    std::array<double, 3> const & GetThreeMomentum() const;
    std::array<double, 4> GetFourMomentum() const;
    double const & GetLength() const;
    std::array<double, 3> const & GetInitialPosition() const;
    std::array<double, 3> const & GetInteractionVertex() const;
    double const & GetHelicity() const;

    void SetMass(double mass);
    void SetEnergy(double energy);
    void SetKineticEnergy(double kinetic_energy);
    void SetDirection(std::array<double, 3> direction);
    void SetThreeMomentum(std::array<double, 3> momentum);
    void SetFourMomentum(std::array<double, 4> momentum);
    void SetLength(double length);
    void SetInitialPosition(std::array<double, 3> initial_position);
    void SetInteractionVertex(std::array<double, 3> interaction_vertex);
    void SetHelicity(double helicity);

    void UpdateMass() const;
    void UpdateEnergy() const;
    void UpdateKineticEnergy() const;
    void UpdateDirection() const;
    void UpdateMomentum() const;
    void UpdateLength() const;
    void UpdateInitialPosition() const;
    void UpdateInteractionVertex() const;

    void FinalizeAvailable(InteractionRecord & record) const;
    void Finalize(InteractionRecord & record) const;
};

class SecondaryParticleRecord {
public:
    size_t const secondary_index;
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

    mutable double mass = 0;
    mutable double energy = 0;
    mutable double kinetic_energy = 0;
    mutable std::array<double, 3> direction = {0, 0, 0};
    mutable std::array<double, 3> momentum = {0, 0, 0};
    mutable double helicity = 0;
public:
    friend std::ostream& ::operator<<(std::ostream& os, SecondaryParticleRecord const& record);

    SecondaryParticleRecord(SecondaryParticleRecord const & other) = delete;
    SecondaryParticleRecord(SecondaryParticleRecord && other) = default;
    SecondaryParticleRecord & operator=(SecondaryParticleRecord const & other) = delete;
    SecondaryParticleRecord & operator=(SecondaryParticleRecord && other) = delete;

    SecondaryParticleRecord(InteractionRecord const & record, size_t secondary_index);

    InteractionSignature const & GetSignature() const;

    Particle GetParticle() const;
    void SetParticle(Particle const & particle);

    ParticleID const & GetID() const;

    ParticleType const & GetType() const;
    double const & GetMass() const;
    double const & GetEnergy() const;
    double const & GetKineticEnergy() const;
    std::array<double, 3> const & GetDirection() const;
    std::array<double, 3> const & GetThreeMomentum() const;
    std::array<double, 4> GetFourMomentum() const;
    std::array<double, 3> const & GetInitialPosition() const;
    double const & GetHelicity() const;

    void SetMass(double mass);
    void SetEnergy(double energy);
    void SetKineticEnergy(double kinetic_energy);
    void SetDirection(std::array<double, 3> direction);
    void SetThreeMomentum(std::array<double, 3> momentum);
    void SetFourMomentum(std::array<double, 4> momentum);
    void SetHelicity(double helicity);

    void UpdateMass() const;
    void UpdateEnergy() const;
    void UpdateKineticEnergy() const;
    void UpdateDirection() const;
    void UpdateMomentum() const;

    void Finalize(InteractionRecord & record) const;
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
    friend std::ostream& ::operator<<(std::ostream& os, CrossSectionDistributionRecord const& record);

    CrossSectionDistributionRecord(CrossSectionDistributionRecord const & other) = delete;
    CrossSectionDistributionRecord(CrossSectionDistributionRecord && other) = default;
    CrossSectionDistributionRecord & operator=(CrossSectionDistributionRecord const & other) = delete;
    CrossSectionDistributionRecord & operator=(CrossSectionDistributionRecord && other) = delete;

    CrossSectionDistributionRecord(InteractionRecord const & record);

    InteractionSignature const & GetSignature() const;
    ParticleID const & GetPrimaryID() const;
    ParticleType const & GetPrimaryType() const;
    std::array<double, 3> const & GetPrimaryInitialPosition() const;
    double const & GetPrimaryMass() const;
    std::array<double, 4> const & GetPrimaryMomentum() const;
    double const & GetPrimaryHelicity() const;

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

    std::vector<SecondaryParticleRecord> & GetSecondaryParticleRecords();
    std::vector<SecondaryParticleRecord> const & GetSecondaryParticleRecords() const;
    SecondaryParticleRecord & GetSecondaryParticleRecord(size_t index);
    SecondaryParticleRecord const & GetSecondaryParticleRecord(size_t index) const;

    void Finalize(InteractionRecord & record) const;
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

class SecondaryDistributionRecord {
    size_t const secondary_index = 0;
public:
    InteractionRecord const record;

    ParticleID const id;
    ParticleType const & type;
    double const & mass;
    std::array<double, 3> const direction;
    std::array<double, 4> const & momentum;
    double const & helicity;
    std::array<double, 3> const & initial_position;
private:
    bool length_set = false;
    mutable double length;
public:
    friend std::ostream& ::operator<<(std::ostream& os, SecondaryDistributionRecord const& record);

    SecondaryDistributionRecord(SecondaryDistributionRecord const & other) = default;
    SecondaryDistributionRecord & operator=(SecondaryDistributionRecord const & other) = delete;
    SecondaryDistributionRecord & operator=(SecondaryDistributionRecord && other) = delete;

    static InteractionRecord CreateSecondaryRecord(InteractionRecord const & parent_record, size_t secondary_index);
    SecondaryDistributionRecord(InteractionRecord & record);
    SecondaryDistributionRecord(InteractionRecord const & parent_record, size_t secondary_index);

    void SetLength(double const & length);
    double const & GetLength() const;

    void Finalize(InteractionRecord & record) const;
};


} // namespace dataclasses
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::InteractionRecord, 0);

#endif // SIREN_InteractionRecord_H
