#include "LeptonInjector/injection/Process.h"

#include <tuple>

namespace LI {
namespace injection {

Process::Process(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions) : primary_type(_primary_type), interactions(_interactions) {}

Process::Process(Process const & other) : primary_type(other.primary_type), interactions(other.interactions) {}

Process::Process(Process && other) : primary_type(other.primary_type), interactions(other.interactions) {}

Process & Process::operator=(Process const & other) {
    primary_type = other.primary_type;
    interactions = other.interactions;
    return *this;
}

Process & Process::operator=(Process && other) {
    primary_type = other.primary_type;
    interactions = other.interactions;
    return *this;
}

void Process::SetInteractions(std::shared_ptr<interactions::InteractionCollection> _interactions) {
    interactions = _interactions;
}

std::shared_ptr<interactions::InteractionCollection> Process::GetInteractions() const {
    return interactions;
}

void Process::SetPrimaryType(LI::dataclasses::Particle::ParticleType _primary_type) {
    primary_type = _primary_type;
}

LI::dataclasses::Particle::ParticleType Process::GetPrimaryType() const {
    return primary_type;
}

bool Process::operator==(Process const & other) const {
    return std::tie(
        primary_type,
        interactions)
        ==
        std::tie(
        other.primary_type,
        other.interactions);
}

bool Process::MatchesHead(std::shared_ptr<Process> const & other) const {
    return std::tie(
        primary_type,
        interactions)
        ==
        std::tie(
        other->primary_type,
        other->interactions);
}

PhysicalProcess::PhysicalProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions) : Process(_primary_type, _interactions) {};

PhysicalProcess::PhysicalProcess(PhysicalProcess const & other) : Process(other), physical_distributions(other.physical_distributions) {};

PhysicalProcess::PhysicalProcess(PhysicalProcess && other) : Process(other), physical_distributions(other.physical_distributions) {};

PhysicalProcess & PhysicalProcess::operator=(PhysicalProcess const & other) {
    Process::operator=(other);
    physical_distributions = other.physical_distributions;
    return *this;
};

PhysicalProcess & PhysicalProcess::operator=(PhysicalProcess && other) {
    Process::operator=(other);
    physical_distributions = other.physical_distributions;
    return *this;
};

void PhysicalProcess::AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) {
    for(auto _dist: physical_distributions) {
        if((*_dist) == (*dist))
            throw std::runtime_error("Cannot add duplicate WeightableDistributions");
    }
    physical_distributions.push_back(dist);
}

std::vector<std::shared_ptr<distributions::WeightableDistribution>> const & PhysicalProcess::GetPhysicalDistributions() const {
    return physical_distributions;
}

InjectionProcess::InjectionProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions) : PhysicalProcess(_primary_type, _interactions) {};

InjectionProcess::InjectionProcess(InjectionProcess const & other) : PhysicalProcess(other), injection_distributions(other.injection_distributions) {};

InjectionProcess::InjectionProcess(InjectionProcess && other) : PhysicalProcess(other), injection_distributions(other.injection_distributions) {};

InjectionProcess & InjectionProcess::operator=(InjectionProcess const & other) {
    PhysicalProcess::operator=(other);
    injection_distributions = other.injection_distributions;
    return *this;
};

InjectionProcess & InjectionProcess::operator=(InjectionProcess && other) {
    PhysicalProcess::operator=(other);
    injection_distributions = other.injection_distributions;
    return *this;
};

void InjectionProcess::AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) {
    throw std::runtime_error("Cannot add a physical distribution to an InjectionProcess");
}
void InjectionProcess::AddInjectionDistribution(std::shared_ptr<distributions::InjectionDistribution> dist) {
    for(auto _dist: injection_distributions) {
        if((*_dist) == (*dist))
            throw std::runtime_error("Cannot add duplicate InjectionDistributions");
    }
    injection_distributions.push_back(dist);
    physical_distributions.push_back(std::static_pointer_cast<distributions::WeightableDistribution>(dist));
}

std::vector<std::shared_ptr<distributions::InjectionDistribution>> const & InjectionProcess::GetInjectionDistributions() const {
    return injection_distributions;
}

/////////////////////////////////////////////

SecondaryInjectionProcess::SecondaryInjectionProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions) : InjectionProcess(_primary_type, _interactions) {};

SecondaryInjectionProcess::SecondaryInjectionProcess(SecondaryInjectionProcess const & other) : InjectionProcess(other), secondary_injection_distributions(other.secondary_injection_distributions) {};

SecondaryInjectionProcess::SecondaryInjectionProcess(SecondaryInjectionProcess && other) : InjectionProcess(other), secondary_injection_distributions(other.secondary_injection_distributions) {};

SecondaryInjectionProcess & SecondaryInjectionProcess::operator=(SecondaryInjectionProcess const & other) {
    InjectionProcess::operator=(other);
    secondary_injection_distributions = other.secondary_injection_distributions;
    return *this;
};

SecondaryInjectionProcess & SecondaryInjectionProcess::operator=(SecondaryInjectionProcess && other) {
    InjectionProcess::operator=(other);
    secondary_injection_distributions = other.secondary_injection_distributions;
    return *this;
};

void SecondaryInjectionProcess::AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) {
    throw std::runtime_error("Cannot add a physical distribution to an SecondaryInjectionProcess");
}

void SecondaryInjectionProcess::AddInjectionDistribution(std::shared_ptr<distributions::InjectionDistribution> dist) {
    throw std::runtime_error("Cannot add an injection distribution to an SecondaryInjectionProcess");
}

void SecondaryInjectionProcess::AddSecondaryInjectionDistribution(std::shared_ptr<distributions::SecondaryInjectionDistribution> dist) {
    for(auto _dist: injection_distributions) {
        if((*_dist) == (*dist))
            throw std::runtime_error("Cannot add duplicate SecondaryInjectionDistributions");
    }
    physical_distributions.push_back(std::static_pointer_cast<distributions::WeightableDistribution>(dist));
    injection_distributions.push_back(std::static_pointer_cast<distributions::InjectionDistribution>(dist));
    secondary_injection_distributions.push_back(dist);
}

std::vector<std::shared_ptr<distributions::SecondaryInjectionDistribution>> const & SecondaryInjectionProcess::GetSecondaryInjectionDistributions() const {
    return secondary_injection_distributions;
}

} // namespace injection
} // namespace LI
