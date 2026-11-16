#pragma once
#ifndef SIREN_GammaDISFromSpline_H
#define SIREN_GammaDISFromSpline_H

// DIS with the kinematics of DISFromSpline, where the hadronic system (and the
// electron in CC nu_e) is returned as GammaSubshower secondaries: gamma-shaped
// longitudinal Cherenkov profiles, sampled for each final-state hadron.
//
// Inputs:
//   shower model    npz with the gamma-mixture parameter distributions per species
//   PYTHIA library  npz with hadronic final states (E_had, PDG codes, 4-momenta)

#include <set>
#include <map>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>

#include "SIREN/dataclasses/InteractionSignature.h"

#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/interactions/DISFromSpline.h"
#include "SIREN/dataclasses/Particle.h"

namespace siren { namespace dataclasses { struct InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class GammaDISFromSpline : public DISFromSpline {
friend cereal::access;
public:
    // parameter curves of one species
    struct SpeciesModel {
        int pid = 0;
        std::vector<double> logE;              // [L] training grid in log10(E/GeV)
        std::vector<double> z;                 // [Z] depth-bin centers [cm]
        std::vector<double> pm;                // [L*Kmax] p(m) over energy (row-major)
        std::vector<double> logN;              // [L] log mean yield
        std::vector<double> logNsig;           // [L] lognormal spread of the yield
        std::vector<int>    m_avail;           // which m have fitted params
        // per m: mean [L*(3m-1)] and cov [L*(3m-1)*(3m-1)], row-major
        std::map<int, std::vector<double>> mean;
        std::map<int, std::vector<double>> cov;
        double logE_lo = 0.0, logE_hi = 0.0;
    };

    // one gamma component of a sampled shower
    struct Component { double alpha; double beta; double amplitude; };

private:
    std::string shower_model_path_;
    std::string pythia_events_path_;
    int    Kmax_ = 2;
    int    top_k_ = 10;
    int    n_bins_ = 12;

    // shower model, keyed by PDG
    std::map<int, SpeciesModel> species_;

    // PYTHIA event library (parallel arrays over events)
    std::vector<double>              pythia_Ehad_;      // [Nev]
    std::vector<std::vector<int>>    pythia_pids_;      // [Nev][<=20]
    std::vector<std::vector<double>> pythia_energies_;  // [Nev][<=20]
    // 4-momenta, used for the Lorentz transformation
    std::vector<std::vector<std::array<double,4>>> pythia_momenta_; // [Nev][<=20][E,px,py,pz]
    std::vector<std::array<double,4>> pythia_phad_;     // [Nev] total hadronic 4-mom
    bool have_momenta_ = false;                          // true if the library has 4-momenta
    std::vector<double>              bin_edges_;         // [n_bins+1] over E_had

    // lepton slot (mu, tau, NC neutrino) + kMaxSub_ GammaSubshower slots
    std::vector<dataclasses::InteractionSignature> gamma_signatures_;
    std::map<std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType>,
             std::vector<dataclasses::InteractionSignature>> gamma_signatures_by_parents_;
    int kMaxSub_ = 0;                                    // number of GammaSubshower slots

    // 1 = CC, 2 = NC, 0 = use the base class value (see the constructor)
    int interaction_type_override_ = 0;

public:
    GammaDISFromSpline();

    // interaction type, target mass and minimum Q2 are read from the spline table
    GammaDISFromSpline(
        std::string differential_filename, std::string total_filename,
        std::vector<siren::dataclasses::ParticleType> primary_types,
        std::vector<siren::dataclasses::ParticleType> target_types,
        std::string shower_model_path, std::string pythia_events_path,
        std::string units = "cm");

    GammaDISFromSpline(
        std::string differential_filename, std::string total_filename,
        int interaction, double target_mass, double minimum_Q2,
        std::vector<siren::dataclasses::ParticleType> primary_types,
        std::vector<siren::dataclasses::ParticleType> target_types,
        std::string shower_model_path, std::string pythia_events_path,
        std::string units = "cm");

    virtual ~GammaDISFromSpline() {}

    bool equal(CrossSection const & other) const override;

    // Returns the hadronic system (and the electron in CC nu_e) as GammaSubshower
    // secondaries whose 4-momentum slots hold [amplitude, alpha, beta, parent index].
    // Muons, taus and NC neutrinos are regular secondaries in slot 0.
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &,
                          std::shared_ptr<siren::utilities::SIREN_random> random) const override;

    // The base class reads the lepton from the secondaries, which are GammaSubshowers
    // here, so the lepton is taken from the lep_* interaction parameters.
    using DISFromSpline::DifferentialCrossSection;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;

    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(
        siren::dataclasses::ParticleType primary_type,
        siren::dataclasses::ParticleType target_type) const override;

    std::string Name() const;

    void LoadShowerModel(std::string const & npz_path);
    void LoadPythiaEvents(std::string const & npz_path);
    bool HasSpecies(int pid) const { return species_.count(pid) > 0; }
    double YieldMean(int pid, double E_GeV) const;

    // 1 = CC, 2 = NC
    int InteractionType() const;
    void SetInteractionType(int interaction_type);

    // components (alpha, beta, amplitude) of one particle's shower
    std::vector<Component> SampleShower(
        int pid, double E_GeV,
        std::shared_ptr<siren::utilities::SIREN_random> random) const;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("ShowerModelPath", shower_model_path_));
            archive(::cereal::make_nvp("PythiaEventsPath", pythia_events_path_));
            archive(::cereal::make_nvp("TopK", top_k_));
            archive(::cereal::make_nvp("InteractionTypeOverride", interaction_type_override_));
            archive(cereal::virtual_base_class<DISFromSpline>(this));
        } else {
            throw std::runtime_error("GammaDISFromSpline only supports version <= 0");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("ShowerModelPath", shower_model_path_));
            archive(::cereal::make_nvp("PythiaEventsPath", pythia_events_path_));
            archive(::cereal::make_nvp("TopK", top_k_));
            archive(::cereal::make_nvp("InteractionTypeOverride", interaction_type_override_));
            archive(cereal::virtual_base_class<DISFromSpline>(this));
            if(!shower_model_path_.empty())  LoadShowerModel(shower_model_path_);
            if(!pythia_events_path_.empty()) LoadPythiaEvents(pythia_events_path_);
            InitializeGammaSignatures();
        } else {
            throw std::runtime_error("GammaDISFromSpline only supports version <= 0");
        }
    }

private:
    double interp1d_(std::vector<double> const & grid, double const * y, int ny_stride, int comp, double lq) const;
    double lq_clamped_(SpeciesModel const & s, double E) const;
    int    bin_index_(double E_had) const;
    static void from_z_(std::vector<double> const & z, int m,
                        std::vector<double> & w, std::vector<double> & alpha, std::vector<double> & beta);

    void InitializeGammaSignatures();

    // (pid, 4-momentum) of the final-state hadrons for a hadronic system p4: top-K
    // hadrons and a pi0 with the rest. The 4-momenta sum to p4.
    std::vector<std::pair<int,std::array<double,4>>> BoostedBreakdown(
        std::array<double,4> const & p4,
        std::shared_ptr<siren::utilities::SIREN_random> random) const;
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::GammaDISFromSpline, 0);
CEREAL_REGISTER_TYPE(siren::interactions::GammaDISFromSpline);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::GammaDISFromSpline);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::DISFromSpline, siren::interactions::GammaDISFromSpline);

#endif // SIREN_GammaDISFromSpline_H
