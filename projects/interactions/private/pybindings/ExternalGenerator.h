#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/ExternalGenerator.h"
#ifdef SIREN_HAS_ACHILLES
#include "../../public/SIREN/interactions/AchillesTrialGenerator.h"
#endif
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

namespace siren {
namespace interactions {

// Records cannot be copied, so overrides that take one pass it to Python by
// reference (as SIREN_OVERRIDE_REF does for the CrossSection trampoline).
#define SIREN_EXTERNAL_OVERRIDE_REF(ret, Base, name, ...) \
    do { \
        pybind11::gil_scoped_acquire gil; \
        pybind11::function override = pybind11::get_override(static_cast<Base const *>(this), #name); \
        if(!override) \
            pybind11::pybind11_fail("Tried to call pure virtual function \"" #Base "::" #name "\""); \
        auto o = override.operator()<pybind11::return_value_policy::reference>(__VA_ARGS__); \
        return pybind11::detail::cast_safe<ret>(std::move(o)); \
    } while(false)

// Python subclasses of the generator interfaces. The Python object must stay
// alive while SIREN uses the generator; the constructors that take one keep
// it alive for as long as the object they build.
class pyExternalGenerator : public ExternalGenerator {
public:
    using ExternalGenerator::ExternalGenerator;
    std::vector<siren::dataclasses::InteractionSignature> GetChannels() const override {
        PYBIND11_OVERRIDE_PURE(std::vector<siren::dataclasses::InteractionSignature>,
                               ExternalGenerator, GetChannels, );
    }
    double TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                             double primary_energy) const override {
        PYBIND11_OVERRIDE_PURE(double, ExternalGenerator, TotalCrossSection, channel, primary_energy);
    }
    double Threshold(siren::dataclasses::InteractionSignature const & channel) const override {
        PYBIND11_OVERRIDE(double, ExternalGenerator, Threshold, channel);
    }
    void GenerateFinalState(siren::dataclasses::CrossSectionDistributionRecord & record,
                            std::shared_ptr<siren::utilities::SIREN_random> random) const override {
        SIREN_EXTERNAL_OVERRIDE_REF(void, ExternalGenerator, GenerateFinalState, record, random);
    }
    std::string SamplerKey() const override {
        PYBIND11_OVERRIDE_PURE(std::string, ExternalGenerator, SamplerKey, );
    }
};

class pyTrialGenerator : public TrialGenerator {
public:
    using TrialGenerator::TrialGenerator;
    std::vector<siren::dataclasses::InteractionSignature> GetChannels() const override {
        PYBIND11_OVERRIDE_PURE(std::vector<siren::dataclasses::InteractionSignature>,
                               TrialGenerator, GetChannels, );
    }
    double TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                             double primary_energy) const override {
        PYBIND11_OVERRIDE_PURE(double, TrialGenerator, TotalCrossSection, channel, primary_energy);
    }
    double Threshold(siren::dataclasses::InteractionSignature const & channel) const override {
        PYBIND11_OVERRIDE(double, TrialGenerator, Threshold, channel);
    }
    double Envelope(siren::dataclasses::InteractionSignature const & channel,
                    double primary_energy) const override {
        PYBIND11_OVERRIDE_PURE(double, TrialGenerator, Envelope, channel, primary_energy);
    }
    double Trial(siren::dataclasses::CrossSectionDistributionRecord const & record,
                 std::shared_ptr<siren::utilities::SIREN_random> random) const override {
        SIREN_EXTERNAL_OVERRIDE_REF(double, TrialGenerator, Trial, record, random);
    }
    void Keep() const override {
        PYBIND11_OVERRIDE_PURE(void, TrialGenerator, Keep, );
    }
    void Fill(siren::dataclasses::CrossSectionDistributionRecord & record) const override {
        SIREN_EXTERNAL_OVERRIDE_REF(void, TrialGenerator, Fill, record);
    }
    std::string SamplerKey() const override {
        PYBIND11_OVERRIDE_PURE(std::string, TrialGenerator, SamplerKey, );
    }
};

} // namespace interactions
} // namespace siren

void register_ExternalGenerator(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<ExternalGenerator, std::shared_ptr<ExternalGenerator>, pyExternalGenerator>(m, "ExternalGenerator")
        .def(init<>())
        .def("GetChannels", &ExternalGenerator::GetChannels)
        .def("TotalCrossSection", &ExternalGenerator::TotalCrossSection,
             arg("channel"), arg("primary_energy"))
        .def("Threshold", &ExternalGenerator::Threshold, arg("channel"))
        .def("GenerateFinalState", &ExternalGenerator::GenerateFinalState,
             arg("record"), arg("random"))
        .def("SamplerKey", &ExternalGenerator::SamplerKey);

    class_<TrialGenerator, std::shared_ptr<TrialGenerator>, pyTrialGenerator>(m, "TrialGenerator")
        .def(init<>())
        .def("GetChannels", &TrialGenerator::GetChannels)
        .def("TotalCrossSection", &TrialGenerator::TotalCrossSection,
             arg("channel"), arg("primary_energy"))
        .def("Threshold", &TrialGenerator::Threshold, arg("channel"))
        .def("Envelope", &TrialGenerator::Envelope, arg("channel"), arg("primary_energy"))
        .def("Trial", &TrialGenerator::Trial, arg("record"), arg("random"))
        .def("Keep", &TrialGenerator::Keep)
        .def("Fill", &TrialGenerator::Fill, arg("record"))
        .def("SamplerKey", &TrialGenerator::SamplerKey);

    class_<BatchResampler, std::shared_ptr<BatchResampler>, ExternalGenerator>(m, "BatchResampler")
        .def(init<std::shared_ptr<TrialGenerator>, unsigned int>(),
             arg("trials"), arg("batch_size"), keep_alive<1, 2>())
        .def("GetTrialGenerator", &BatchResampler::GetTrialGenerator)
        .def("GetBatchSize", &BatchResampler::GetBatchSize);

    class_<ExternalGeneratorCrossSection, std::shared_ptr<ExternalGeneratorCrossSection>, CrossSection>(
            m, "ExternalGeneratorCrossSection")
        .def(init<std::shared_ptr<ExternalGenerator>>(), arg("generator"), keep_alive<1, 2>())
        .def("GetGenerator", &ExternalGeneratorCrossSection::GetGenerator)
        .def("TotalCrossSection",
             overload_cast<siren::dataclasses::InteractionRecord const &>(
                 &ExternalGeneratorCrossSection::TotalCrossSection, const_))
        .def("GetPossibleSignatures", &ExternalGeneratorCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents", &ExternalGeneratorCrossSection::GetPossibleSignaturesFromParents)
        .def("GetPossibleTargets", &ExternalGeneratorCrossSection::GetPossibleTargets)
        .def("GetPossiblePrimaries", &ExternalGeneratorCrossSection::GetPossiblePrimaries)
        .def("FinalStateDensityEvaluable", &ExternalGeneratorCrossSection::FinalStateDensityEvaluable)
        .def("FinalStateSamplerKey", &ExternalGeneratorCrossSection::FinalStateSamplerKey);

#ifdef SIREN_HAS_ACHILLES
    class_<AchillesTrialGenerator, std::shared_ptr<AchillesTrialGenerator>, TrialGenerator>(
            m, "AchillesTrialGenerator")
        .def(init<std::string const &, std::string const &,
                  std::vector<siren::dataclasses::InteractionSignature> const &>(),
             arg("run_card"), arg("data_directory"), arg("channels"))
        .def("EnergyRange", &AchillesTrialGenerator::EnergyRange);
    m.attr("HAS_ACHILLES") = true;
#else
    m.attr("HAS_ACHILLES") = false;
#endif
}
