#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace siren {
namespace injection {

std::string PhaseSpaceConventionName(PhaseSpaceConvention convention) {
    switch (convention) {
        case PhaseSpaceConvention::RestFrameSolidAngle:
            return "RestFrameSolidAngle";
        case PhaseSpaceConvention::LabFrameSolidAngle:
            return "LabFrameSolidAngle";
        case PhaseSpaceConvention::Recursive2Body:
            return "Recursive2Body";
        case PhaseSpaceConvention::Dalitz:
            return "Dalitz";
        case PhaseSpaceConvention::HelicityAngles:
            return "HelicityAngles";
        case PhaseSpaceConvention::BjorkenXY:
            return "BjorkenXY";
        case PhaseSpaceConvention::MandelstamST:
            return "MandelstamST";
        case PhaseSpaceConvention::Custom:
            return "Custom";
    }
    return "Unknown";
}

namespace {

int ConventionPriority(PhaseSpaceConvention convention) {
    switch (convention) {
        case PhaseSpaceConvention::RestFrameSolidAngle:
            return 0;
        case PhaseSpaceConvention::Recursive2Body:
            return 1;
        case PhaseSpaceConvention::BjorkenXY:
            return 2;
        case PhaseSpaceConvention::MandelstamST:
            return 3;
        case PhaseSpaceConvention::HelicityAngles:
            return 4;
        case PhaseSpaceConvention::Dalitz:
            return 5;
        case PhaseSpaceConvention::LabFrameSolidAngle:
            return 6;
        case PhaseSpaceConvention::Custom:
            return 7;
    }
    return 8;
}

} // anonymous namespace

int MultiChannelPhaseSpace::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    WarnOnConventionMismatch();

    if (channels.empty()) {
        throw std::runtime_error("MultiChannelPhaseSpace has no channels");
    }

    // Pick a channel according to the weights
    double r = random->Uniform(0, 1);
    double cumulative = 0.0;
    int selected = static_cast<int>(channels.size()) - 1;
    for (int i = 0; i < static_cast<int>(channels.size()); ++i) {
        cumulative += weights[i];
        if (r < cumulative) {
            selected = i;
            break;
        }
    }

    channels[selected]->Sample(random, detector_model, record);
    return selected;
}

double MultiChannelPhaseSpace::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    WarnOnConventionMismatch();

    double density = 0.0;
    for (size_t i = 0; i < channels.size(); ++i) {
        density += weights[i] * channels[i]->Density(detector_model, record);
    }
    return density;
}

PhaseSpaceConvention MultiChannelPhaseSpace::CommonConvention() const
{
    if (channels.empty()) return PhaseSpaceConvention::Custom;

    std::map<PhaseSpaceConvention, int> counts;
    int n_custom = 0;
    for (auto const & channel : channels) {
        if (!channel) {
            throw std::runtime_error("MultiChannelPhaseSpace contains a null channel");
        }
        PhaseSpaceConvention convention = channel->Convention();
        counts[convention] += 1;
        if (convention == PhaseSpaceConvention::Custom) ++n_custom;
    }

    if (n_custom > 0 && n_custom != static_cast<int>(channels.size())) {
        throw std::runtime_error(
            "MultiChannelPhaseSpace cannot mix Custom phase-space conventions with non-Custom conventions");
    }
    if (n_custom == static_cast<int>(channels.size())) {
        return PhaseSpaceConvention::Custom;
    }

    PhaseSpaceConvention best = channels.front()->Convention();
    int best_count = -1;
    int best_priority = 100;
    for (auto const & kv : counts) {
        int priority = ConventionPriority(kv.first);
        if (kv.second > best_count ||
            (kv.second == best_count && priority < best_priority)) {
            best = kv.first;
            best_count = kv.second;
            best_priority = priority;
        }
    }
    return best;
}

std::vector<std::string> MultiChannelPhaseSpace::ValidateConventions() const
{
    std::vector<std::string> diagnostics;
    if (channels.empty()) return diagnostics;

    PhaseSpaceConvention common = CommonConvention();
    for (size_t i = 0; i < channels.size(); ++i) {
        PhaseSpaceConvention convention = channels[i]->Convention();
        if (convention != common) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << ") uses " << PhaseSpaceConventionName(convention)
                << " while the selected common convention is "
                << PhaseSpaceConventionName(common)
                << "; densities must be transformed before they are combined";
            diagnostics.push_back(oss.str());
        }
    }
    return diagnostics;
}

void MultiChannelPhaseSpace::WarnOnConventionMismatch() const
{
    auto diagnostics = ValidateConventions();
    if (!diagnostics.empty() && !convention_warning_emitted) {
        for (auto const & diagnostic : diagnostics) {
            std::cerr << "SIREN phase-space warning: " << diagnostic << std::endl;
        }
        convention_warning_emitted = true;
    }
}

std::vector<std::string> MultiChannelPhaseSpace::ValidateChannels(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord template_record,
    int samples_per_channel) const
{
    std::vector<std::string> diagnostics;
    auto convention_diagnostics = ValidateConventions();
    diagnostics.insert(diagnostics.end(),
                       convention_diagnostics.begin(),
                       convention_diagnostics.end());

    for (size_t i = 0; i < channels.size(); ++i) {
        int n_zero = 0;
        int n_nan = 0;
        int n_neg = 0;
        double max_ratio = 0.0;

        for (int s = 0; s < samples_per_channel; ++s) {
            siren::dataclasses::InteractionRecord record = template_record;
            channels[i]->Sample(random, detector_model, record);

            double d_self = channels[i]->Density(detector_model, record);
            if (d_self <= 0 || !std::isfinite(d_self)) continue;

            for (size_t j = 0; j < channels.size(); ++j) {
                if (i == j) continue;
                double d_other = channels[j]->Density(detector_model, record);

                if (std::isnan(d_other)) {
                    ++n_nan;
                } else if (d_other < 0) {
                    ++n_neg;
                } else if (d_other == 0) {
                    ++n_zero;
                } else {
                    double ratio = d_self / d_other;
                    if (ratio > max_ratio) max_ratio = ratio;
                }
            }
        }

        if (n_nan > 0) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): " << n_nan << "/" << samples_per_channel
                << " samples produced NaN density from other channels";
            diagnostics.push_back(oss.str());
        }
        if (n_neg > 0) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): " << n_neg << "/" << samples_per_channel
                << " samples produced negative density from other channels";
            diagnostics.push_back(oss.str());
        }
        if (n_zero > samples_per_channel / 2) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): " << n_zero << "/" << samples_per_channel
                << " samples got zero density from other channels "
                << "(possible measure mismatch or non-overlapping support)";
            diagnostics.push_back(oss.str());
        }
        if (max_ratio > 1e6) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): max density ratio " << max_ratio
                << " (extreme ratio may indicate measure mismatch)";
            diagnostics.push_back(oss.str());
        }
    }

    return diagnostics;
}

} // namespace injection
} // namespace siren
