#include "LeptonInjector/injection/Process.h"

#include <tuple>

namespace LI {
namespace injection {

Process::Process(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<crosssections::CrossSectionCollection> _cross_sections) : primary_type(_primary_type), cross_sections(_cross_sections) {}

Process::Process(Process const & other) : primary_type(other.primary_type), cross_sections(other.cross_sections) {}

Process::Process(Process && other) : primary_type(other.primary_type), cross_sections(other.cross_sections) {}

Process & Process::operator=(Process const & other) {
    primary_type = other.primary_type;
    cross_sections = other.cross_sections;
    return *this;
}

Process & Process::operator=(Process && other) {
    primary_type = other.primary_type;
    cross_sections = other.cross_sections;
    return *this;
}

void Process::SetCrossSections(std::shared_ptr<crosssections::CrossSectionCollection> _cross_sections) {
    cross_sections = _cross_sections;
}

std::shared_ptr<crosssections::CrossSectionCollection> Process::GetCrossSections() const {
    return cross_sections;
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
        cross_sections)
        ==
        std::tie(
        other.primary_type,
        other.cross_sections);
}

bool Process::MatchesHead(std::shared_ptr<Process> const & other) const {
    return std::tie(
        primary_type,
        cross_sections)
        ==
        std::tie(
        other->primary_type,
        other->cross_sections);
}

PhysicalProcess::PhysicalProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<crosssections::CrossSectionCollection> _cross_sections) : Process(_primary_type, _cross_sections) {};

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
        if((*_dist) == (*dist)) return;
    }
    physical_distributions.push_back(dist);
}

std::vector<std::shared_ptr<distributions::WeightableDistribution>> const & PhysicalProcess::GetPhysicalDistributions() const {
    return physical_distributions;
}

InjectionProcess::InjectionProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<crosssections::CrossSectionCollection> _cross_sections) : PhysicalProcess(_primary_type, _cross_sections) {};

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
        if((*_dist) == (*dist)) return;
    }
    injection_distributions.push_back(dist);
    physical_distributions.push_back(std::static_pointer_cast<distributions::WeightableDistribution>(dist));
}

std::vector<std::shared_ptr<distributions::InjectionDistribution>> const & InjectionProcess::GetInjectionDistributions() const {
    return injection_distributions;
}

} // namespace injection
} // namespace LI
