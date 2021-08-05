#ifndef LI_CrossSection_H
#define LI_CrossSection_H

#include <map>
#include <set>
#include <array>
#include <string>

#include <photospline/splinetable.h>

#include "LeptonInjector/Particle.h"

namespace LeptonInjector {

struct InteractionSignature {
    LeptonInjector::Particle::ParticleType primary_type;
    LeptonInjector::Particle::ParticleType target_type;
    std::vector<LeptonInjector::Particle::ParticleType> secondary_types;
};

struct InteractionRecord {
    InteractionSignature signature;
    std::array<double, 4> primary_momentum = {0, 0, 0, 0};
    std::array<double, 4> target_momentum = {0, 0, 0, 0};
    std::vector<std::array<double, 4>> secondary_momenta;
};

class CrossSection {
private:
    //
public:
    CrossSection() {};
    virtual double TotalCrossSection(InteractionRecord const &) const = 0;
    virtual double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy) const = 0;
    virtual double DifferentialCrossSection(InteractionRecord const &) const = 0;
    virtual void SampleFinalState(InteractionRecord &) const = 0;

    virtual std::vector<Particle::ParticleType> GetPossibleTargets() const = 0;
    virtual std::vector<Particle::ParticleType> GetPossiblePrimaries() const = 0;
    virtual std::vector<InteractionSignature> GetPossibleSignatures() const = 0;
};

class DISFromSpline : public CrossSection {
private:
    photospline::splinetable<> differential_cross_section_;
    photospline::splinetable<> total_cross_section_;

    std::map<std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>, InteractionSignature> signatures_;
    std::set<LeptonInjector::Particle::ParticleType> primary_types_;
    std::set<LeptonInjector::Particle::ParticleType> target_types_;

    double minimum_Q2_;
    double target_mass_;
    int interaction_type_;

public:
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);

    double TotalCrossSection(InteractionRecord const &) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy) const;
    double DifferentialCrossSection(InteractionRecord const &) const;
    void SampleFinalState(InteractionRecord &) const;

    std::vector<Particle::ParticleType> GetPossibleTargets() const;
    std::vector<Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<InteractionSignature> GetPossibleSignatures() const;

    void LoadFromFile(std::string differential_filename, std::string total_filename);

private:
    void ReadParamsFromSplineTable();
    void InitializeSignatures();
};

} // namespace LeptonInjector

#endif // LI_CrossSection_H

