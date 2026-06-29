#pragma once
#ifndef SIREN_PythiaDISCrossSection_H
#define SIREN_PythiaDISCrossSection_H

#include <set>
#include <map>
#include <memory>
#include <vector>
#include <cstdint>
#include <utility>
#include <string>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/Particle.h"

// Forward declarations for Pythia
namespace Pythia8 {
    class Pythia;
}

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// Forward declaration -- defined in .cxx
class SIRENRndm;

class PythiaDISCrossSection : public CrossSection {
friend cereal::access;
private:
    // Splines for total/differential cross section (SIREN weighting).
    // The total spline is always required (drives interaction depth / position /
    // survival). The differential spline is OPTIONAL: when absent,
    // FinalStateProbability returns a constant that cancels in the unbiased
    // weight (final state comes from Pythia, injection == physical).
    photospline::splinetable<> differential_cross_section_;
    photospline::splinetable<> total_cross_section_;
    bool has_differential_ = false;

    // Pythia instance (mutable because SampleFinalState is const)
    mutable std::unique_ptr<Pythia8::Pythia> pythia_;
    mutable std::shared_ptr<SIRENRndm> siren_rndm_;

    // Signature bookkeeping
    std::vector<dataclasses::InteractionSignature> signatures_;
    std::set<siren::dataclasses::ParticleType> primary_types_;
    std::set<siren::dataclasses::ParticleType> target_types_;
    std::map<std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType>, std::vector<dataclasses::InteractionSignature>> signatures_by_parent_types_;

    // DIS parameters
    int interaction_type_;  // 1=CC, 2=NC
    double target_mass_;
    double minimum_Q2_;
    double unit;

    // Pythia configuration
    std::string pdf_set_;
    std::string pythia_data_path_;

    // Helper methods
    void InitializePythia(double E_nu, int target_pdg) const;
    void InitializeSignatures();
    void LoadFromFile(std::string differential_filename, std::string total_filename);
    void LoadFromMemory(std::vector<char> & differential_data, std::vector<char> & total_data);
    void ReadParamsFromSplineTable();
    void SetUnits(std::string units);

    // Particle ID helpers
    static bool IsCharmedHadron(int pdgId);
    static siren::dataclasses::ParticleType PdgToParticleType(int pdgId);
    static double GetLeptonMass(siren::dataclasses::ParticleType lepton_type);
    static std::map<std::string, int> getIndices(siren::dataclasses::InteractionSignature signature);

public:
    PythiaDISCrossSection();
    ~PythiaDISCrossSection();

    // Main constructor
    PythiaDISCrossSection(
        std::string differential_filename,
        std::string total_filename,
        int interaction_type,
        double target_mass,
        double minimum_Q2,
        std::set<siren::dataclasses::ParticleType> primary_types,
        std::set<siren::dataclasses::ParticleType> target_types,
        std::string pythia_data_path,
        std::string pdf_set = "LHAPDF6:HERAPDF20_NLO_EIG",
        std::string units = "cm"
    );

    // Constructor with vectors
    PythiaDISCrossSection(
        std::string differential_filename,
        std::string total_filename,
        int interaction_type,
        double target_mass,
        double minimum_Q2,
        std::vector<siren::dataclasses::ParticleType> primary_types,
        std::vector<siren::dataclasses::ParticleType> target_types,
        std::string pythia_data_path,
        std::string pdf_set = "LHAPDF6:HERAPDF20_NLO_EIG",
        std::string units = "cm"
    );

    virtual bool equal(CrossSection const & other) const override;

    // Cross section from splines
    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(siren::dataclasses::ParticleType primary, double energy) const;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass, double Q2=std::numeric_limits<double>::quiet_NaN()) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;

    // Fragmentation fraction (from Pythia output statistics)
    double FragmentationFraction(siren::dataclasses::Particle::ParticleType secondary) const;

    // Final state sampling via Pythia
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random> random) const override;

    // Generate raw Pythia charm-DIS samples for building total/differential
    // splines (init once per energy, using the same Pythia config as
    // SampleFinalState). out_sigma_mb is per energy (Pythia's generated cross
    // section, mb); out_E/out_x/out_y are flat per-event muon-reconstructed
    // (E, Bjorken x, y). Consumed by the python generate_*_spline helpers.
    static void GeneratePythiaCharmSamples(
        int interaction_type, int primary_pdg, int target_pdg, double target_mass,
        std::string pdf_set, std::string pythia_data_path, double minimum_Q2,
        std::vector<double> const & energies, int n_events,
        std::vector<double> & out_sigma_mb,
        std::vector<double> & out_E, std::vector<double> & out_x, std::vector<double> & out_y);

    // Signature methods
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;

    // Getters
    double GetMinimumQ2() const { return minimum_Q2_; }
    double GetTargetMass() const { return target_mass_; }
    int GetInteractionType() const { return interaction_type_; }

public:
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            splinetable_buffer buf;
            archive(::cereal::make_nvp("HasDifferential", has_differential_));

            // The differential spline is optional; only serialize it when present.
            if(has_differential_) {
                buf.size = 0;
                auto diff_obj = differential_cross_section_.write_fits_mem();
                buf.data = diff_obj.first;
                buf.size = diff_obj.second;
                std::vector<char> diff_blob;
                diff_blob.resize(buf.size);
                std::copy((char*)buf.data, (char*)buf.data + buf.size, &diff_blob[0]);
                archive(::cereal::make_nvp("DifferentialCrossSectionSpline", diff_blob));
            }

            buf.size = 0;
            auto result_obj = total_cross_section_.write_fits_mem();
            buf.data = result_obj.first;
            buf.size = result_obj.second;

            std::vector<char> total_blob;
            total_blob.resize(buf.size);
            std::copy((char*)buf.data, (char*)buf.data + buf.size, &total_blob[0]);

            archive(::cereal::make_nvp("TotalCrossSectionSpline", total_blob));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(::cereal::make_nvp("TargetTypes", target_types_));
            archive(::cereal::make_nvp("InteractionType", interaction_type_));
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            archive(::cereal::make_nvp("Unit", unit));
            archive(::cereal::make_nvp("PdfSet", pdf_set_));
            archive(::cereal::make_nvp("PythiaDataPath", pythia_data_path_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("PythiaDISCrossSection only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            std::vector<char> differential_data;
            std::vector<char> total_data;
            archive(::cereal::make_nvp("HasDifferential", has_differential_));
            if(has_differential_) {
                archive(::cereal::make_nvp("DifferentialCrossSectionSpline", differential_data));
            }
            archive(::cereal::make_nvp("TotalCrossSectionSpline", total_data));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(::cereal::make_nvp("TargetTypes", target_types_));
            archive(::cereal::make_nvp("InteractionType", interaction_type_));
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            archive(::cereal::make_nvp("Unit", unit));
            archive(::cereal::make_nvp("PdfSet", pdf_set_));
            archive(::cereal::make_nvp("PythiaDataPath", pythia_data_path_));
            archive(cereal::virtual_base_class<CrossSection>(this));
            LoadFromMemory(differential_data, total_data);
            InitializeSignatures();
        } else {
            throw std::runtime_error("PythiaDISCrossSection only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::PythiaDISCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::PythiaDISCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::PythiaDISCrossSection);

#endif // SIREN_PythiaDISCrossSection_H
