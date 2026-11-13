#include "SIREN/interactions/AchillesTrialGenerator.h"

#include <algorithm>                                      // for find
#include <cmath>                                          // for sqrt
#include <cstdint>                                        // for uint64_t
#include <filesystem>                                     // for current_path
#include <fstream>                                        // for ifstream
#include <iomanip>                                        // for setw
#include <sstream>                                        // for ostringstream
#include <stdexcept>                                      // for runtime_error

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-literal-operator"
#endif
#include "Achilles/Event.hh"
#include "Achilles/EventGen.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "SIREN/dataclasses/InteractionRecord.h"          // for CrossSectionDistributionRecord
#include "SIREN/dataclasses/ParticleType.h"                // for ParticleType

namespace siren {
namespace interactions {

namespace {

constexpr double kMeVPerGeV = 1000.0;

std::string ReadFile(std::string const & path) {
    std::ifstream in(path);
    if(!in)
        throw std::invalid_argument("AchillesTrialGenerator: cannot read run card " + path);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

// FNV-1a, to name the sampler by its run card.
std::string Digest(std::string const & text) {
    std::uint64_t h = 1469598103934665603ull;
    for(unsigned char c : text) {
        h ^= c;
        h *= 1099511628211ull;
    }
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << h;
    return out.str();
}

std::string Describe(siren::dataclasses::InteractionSignature const & s) {
    std::ostringstream oss;
    oss << s.primary_type << " + " << s.target_type << " ->";
    for(auto const & t : s.secondary_types) oss << ' ' << t;
    return oss.str();
}

// Changes the working directory for its lifetime.
class WorkingDirectory {
public:
    explicit WorkingDirectory(std::filesystem::path const & dir)
        : previous_(std::filesystem::current_path()) {
        std::filesystem::current_path(dir);
    }
    ~WorkingDirectory() {
        std::error_code ec;
        std::filesystem::current_path(previous_, ec);
    }
private:
    std::filesystem::path previous_;
};

} // namespace

struct AchillesTrialGenerator::Impl {
    std::unique_ptr<achilles::EventGen> generator;
    std::filesystem::path data_directory;
    std::vector<siren::dataclasses::InteractionSignature> channels;
    std::string key;
    double min_energy = 0.0;  // MeV
    double max_energy = 0.0;  // MeV
    achilles::Event kept;
    bool has_kept = false;

    void RequireChannel(siren::dataclasses::InteractionSignature const & channel) const {
        if(std::find(channels.begin(), channels.end(), channel) == channels.end())
            throw std::invalid_argument("AchillesTrialGenerator: " + Describe(channel)
                + " is not one of this generator's channels");
    }
};

AchillesTrialGenerator::AchillesTrialGenerator(
        std::string const & run_card,
        std::string const & data_directory,
        std::vector<siren::dataclasses::InteractionSignature> const & channels)
    : impl_(std::make_unique<Impl>())
{
    if(channels.empty())
        throw std::invalid_argument("AchillesTrialGenerator: no channels");
    for(auto const & c : channels) {
        if(c.secondary_types.size() != 2
           || c.secondary_types[1] != siren::dataclasses::ParticleType::Hadrons) {
            throw std::invalid_argument("AchillesTrialGenerator: channel " + Describe(c)
                + " must be nu + nucleus -> lepton + Hadrons");
        }
    }
    std::string const text = ReadFile(run_card);
    std::filesystem::path const card = std::filesystem::absolute(run_card);
    YAML::Node root = YAML::Load(text);
    YAML::Node node = root["Achilles"] ? root["Achilles"] : root;
    YAML::Node beam = node["Beam"];
    if(!beam || !beam["Type"] || beam["Type"].as<std::string>() != "Geometry" || !beam["EnergyRange"])
        throw std::invalid_argument("AchillesTrialGenerator: the run card needs a Beam of Type: "
            "Geometry with an EnergyRange in MeV");
    impl_->min_energy = beam["EnergyRange"][0].as<double>();
    impl_->max_energy = beam["EnergyRange"][1].as<double>();
    if(!(impl_->min_energy < impl_->max_energy))
        throw std::invalid_argument("AchillesTrialGenerator: the beam EnergyRange is empty");

    impl_->data_directory = std::filesystem::absolute(data_directory);
    {
        WorkingDirectory here(impl_->data_directory);
        impl_->generator = std::make_unique<achilles::EventGen>(
            node, std::vector<std::string>{}, card.string());
        impl_->generator->Initialize();
    }
    impl_->generator->SetSamplingMode(achilles::EventGen::SamplingMode::EnvelopeNoRetry);

    for(auto const & c : channels) {
        double const envelope = impl_->generator->VertexEnvelope(
            achilles::PID{static_cast<int>(c.primary_type)}, achilles::PID{static_cast<int>(c.target_type)});
        if(!(envelope > 0.0))
            throw std::invalid_argument("AchillesTrialGenerator: the run card has no process for "
                + Describe(c) + " (VertexEnvelope is zero)");
    }
    impl_->channels = channels;
    impl_->key = "Achilles|card=" + Digest(text);
}

AchillesTrialGenerator::~AchillesTrialGenerator() = default;

std::vector<siren::dataclasses::InteractionSignature> AchillesTrialGenerator::GetChannels() const {
    return impl_->channels;
}

double AchillesTrialGenerator::TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                                                 double primary_energy) const {
    impl_->RequireChannel(channel);
    double const e = primary_energy * kMeVPerGeV;
    // Outside the warm-up range Achilles has no model support and rejects
    // every trial, so the channel must not be selected there.
    if(e < impl_->min_energy || e > impl_->max_energy) return 0.0;
    return impl_->generator->TotalXSec(e, achilles::PID{static_cast<int>(channel.primary_type)},
                                       achilles::PID{static_cast<int>(channel.target_type)});
}

double AchillesTrialGenerator::Threshold(siren::dataclasses::InteractionSignature const &) const {
    return impl_->min_energy / kMeVPerGeV;
}

double AchillesTrialGenerator::Envelope(siren::dataclasses::InteractionSignature const & channel,
                                        double) const {
    impl_->RequireChannel(channel);
    return impl_->generator->VertexEnvelope(achilles::PID{static_cast<int>(channel.primary_type)},
                                            achilles::PID{static_cast<int>(channel.target_type)});
}

double AchillesTrialGenerator::Trial(siren::dataclasses::CrossSectionDistributionRecord const & record,
                                     std::shared_ptr<siren::utilities::SIREN_random>) const {
    auto const & signature = record.GetSignature();
    impl_->RequireChannel(signature);
    auto const & p = record.GetPrimaryMomentum();
    achilles::FourVector const ray{p[0] * kMeVPerGeV, p[1] * kMeVPerGeV,
                                   p[2] * kMeVPerGeV, p[3] * kMeVPerGeV};
    impl_->generator->InjectRay(ray, achilles::PID{static_cast<int>(signature.primary_type)},
                                achilles::PID{static_cast<int>(signature.target_type)});
    // Achilles' models open data files on first use, by paths relative to the
    // working directory, so trials run from the data directory too.
    WorkingDirectory here(impl_->data_directory);
    if(!impl_->generator->GenerateSingleEvent()) return 0.0;
    return impl_->generator->LastEvent().Weight();
}

void AchillesTrialGenerator::Keep() const {
    impl_->kept = impl_->generator->LastEvent();
    impl_->has_kept = true;
}

void AchillesTrialGenerator::Fill(siren::dataclasses::CrossSectionDistributionRecord & record) const {
    if(!impl_->has_kept)
        throw std::runtime_error("AchillesTrialGenerator::Fill: no emitted trial was kept");
    achilles::Event const & event = impl_->kept;

    std::vector<achilles::Particle const *> particles;
    for(auto const & lepton : event.Leptons()) {
        auto const s = lepton.Status();
        if(s == achilles::ParticleStatus::initial_state || s == achilles::ParticleStatus::beam) continue;
        particles.push_back(&lepton);
    }
    std::int64_t spectators = 0;
    for(auto const & hadron : event.Hadrons()) {
        if(hadron.Status() == achilles::ParticleStatus::background) {
            ++spectators;
            continue;
        }
        particles.push_back(&hadron);
    }
    if(particles.empty() || !particles.front()->Info().IsLepton())
        throw std::runtime_error("AchillesTrialGenerator::Fill: the Achilles event has no outgoing lepton");

    std::vector<siren::dataclasses::ParticleType> types;
    types.reserve(particles.size());
    for(auto const * part : particles)
        types.push_back(static_cast<siren::dataclasses::ParticleType>(static_cast<int>(part->ID())));
    record.SetSecondaryTypes(types);

    auto & secondaries = record.GetSecondaryParticleRecords();
    for(std::size_t i = 0; i < particles.size(); ++i) {
        auto const & mom = particles[i]->Momentum();
        std::array<double, 4> const p{mom.E() / kMeVPerGeV, mom.Px() / kMeVPerGeV,
                                      mom.Py() / kMeVPerGeV, mom.Pz() / kMeVPerGeV};
        double const m2 = p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3];
        secondaries[i].SetFourMomentum(p);
        secondaries[i].SetMass(std::sqrt(std::max(m2, 0.0)));
        secondaries[i].SetHelicity(0.0);
        int const status = static_cast<int>(particles[i]->Status());
        if(status != 1) record.SetSecondaryStatus(i, status);
    }
    record.interaction_int_parameters["process_id"] = event.ProcessId();
    record.interaction_int_parameters["achilles_spectators"] = spectators;
    record.interaction_parameters["achilles_weight"] = event.Weight();
}

std::string AchillesTrialGenerator::SamplerKey() const {
    return impl_->key;
}

std::pair<double, double> AchillesTrialGenerator::EnergyRange() const {
    return {impl_->min_energy / kMeVPerGeV, impl_->max_energy / kMeVPerGeV};
}

} // namespace interactions
} // namespace siren
