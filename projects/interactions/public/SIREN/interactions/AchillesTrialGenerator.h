#pragma once
#ifndef SIREN_AchillesTrialGenerator_H
#define SIREN_AchillesTrialGenerator_H

#include <memory>                                         // for unique_ptr
#include <string>                                         // for string
#include <vector>                                         // for vector
#include <utility>                                        // for pair

#include "SIREN/interactions/ExternalGenerator.h"          // for TrialGenerator
#include "SIREN/dataclasses/InteractionSignature.h"        // for InteractionSignature

namespace siren {
namespace interactions {

// The Achilles event generator as a SIREN TrialGenerator, through Achilles'
// geometry-driver interface (AchillesGen/Achilles, feature/nugeometry-interface).
// Wrap it in a BatchResampler and an ExternalGeneratorCrossSection to use it.
//
// One trial is EventGen::InjectRay for the vertex's neutrino four-momentum and
// target, then EventGen::GenerateSingleEvent. An emitted trial carries
// Achilles' partial-unweighting weight W = max(1, w / w_max); the envelope is
// EventGen::VertexEnvelope and the total cross section EventGen::TotalXSec.
// Energies cross the boundary in MeV (Achilles) and GeV (SIREN); cross
// sections are cm^2 on both sides. Achilles events are in the lab frame.
//
// A filled vertex holds the outgoing leptons first, then every hadron Achilles
// reports except the spectator nucleons of the nucleus (status 25), each with
// Achilles' own status code (1 final; 3 the struck nucleon; 26 captured; 29
// interacted in the cascade; ...). Achilles particles carry no parent links,
// so none are set. The vertex also gets the integer attributes process_id
// (Achilles' process ID) and achilles_spectators, and the double parameter
// achilles_weight (W).
//
// The run card is an Achilles run card, at its top level or under an
// "Achilles" key, with a beam of Type: Geometry and an EnergyRange in MeV.
// Achilles resolves its data paths against the working directory, and its
// models open some data files on first use. The generator therefore changes
// into data_directory (an Achilles build or install directory with data/)
// while Achilles initializes and during every trial, and changes back after.
//
// Each channel is nu + nucleus -> charged lepton + Hadrons, for a neutrino
// and nucleus the run card configures. Achilles keeps its own random stream
// (seeded by the run card), so a SIREN seed alone does not fix the events.
// Achilles is not thread-safe: use one generator per process.
class AchillesTrialGenerator : public TrialGenerator {
public:
    AchillesTrialGenerator(std::string const & run_card,
                           std::string const & data_directory,
                           std::vector<siren::dataclasses::InteractionSignature> const & channels);
    ~AchillesTrialGenerator() override;

    std::vector<siren::dataclasses::InteractionSignature> GetChannels() const override;
    double TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                             double primary_energy) const override;
    double Threshold(siren::dataclasses::InteractionSignature const & channel) const override;
    double Envelope(siren::dataclasses::InteractionSignature const & channel,
                    double primary_energy) const override;
    double Trial(siren::dataclasses::CrossSectionDistributionRecord const & record,
                 std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    void Keep() const override;
    void Fill(siren::dataclasses::CrossSectionDistributionRecord & record) const override;
    std::string SamplerKey() const override;

    // The beam EnergyRange from the run card, in GeV.
    std::pair<double, double> EnergyRange() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace interactions
} // namespace siren

#endif // SIREN_AchillesTrialGenerator_H
