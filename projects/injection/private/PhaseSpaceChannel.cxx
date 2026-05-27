#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace siren {
namespace injection {

// ================================================================ //
//  PhaseSpaceChannel default Convention() from Topology+Measure      //
// ================================================================ //

PhaseSpaceConvention PhaseSpaceChannel::Convention() const {
    PhaseSpaceMeasure m = Measure();
    switch (m) {
        case PhaseSpaceMeasure::SolidAngleRest:  return PhaseSpaceConvention::RestFrameSolidAngle;
        case PhaseSpaceMeasure::SolidAngleLab:   return PhaseSpaceConvention::LabFrameSolidAngle;
        case PhaseSpaceMeasure::Recursive2Body:  return PhaseSpaceConvention::Recursive2Body;
        case PhaseSpaceMeasure::DalitzPair:      return PhaseSpaceConvention::Dalitz;
        case PhaseSpaceMeasure::HelicityAngles:  return PhaseSpaceConvention::HelicityAngles;
        case PhaseSpaceMeasure::MandelstamQ2:    return PhaseSpaceConvention::MandelstamST;
        case PhaseSpaceMeasure::BjorkenXY:       return PhaseSpaceConvention::BjorkenXY;
        case PhaseSpaceMeasure::Unspecified:      return PhaseSpaceConvention::Custom;
    }
    return PhaseSpaceConvention::Custom;
}

// ================================================================ //
//  MultiChannelPhaseSpace                                            //
// ================================================================ //
//
//  All channels MUST report density in the same measure.  The
//  combiner sums densities directly:  g(x) = sum_i alpha_i g_i(x).
//  No cross-channel Jacobian conversion is performed.
//
//  This follows the MadGraph/MadEvent design: each channel defines
//  a complete mapping and evaluates its own density at any phase-space
//  point.  Different channels may use different internal
//  parameterizations (different spectator indices, different random
//  number mappings), but the density they REPORT must be with
//  respect to the same integration measure so the sum is meaningful.
//

namespace {

int MeasurePriority(PhaseSpaceMeasure m) {
    switch (m) {
        case PhaseSpaceMeasure::SolidAngleRest:  return 0;
        case PhaseSpaceMeasure::Recursive2Body:  return 1;
        case PhaseSpaceMeasure::MandelstamQ2:    return 2;
        case PhaseSpaceMeasure::BjorkenXY:       return 3;
        case PhaseSpaceMeasure::HelicityAngles:  return 4;
        case PhaseSpaceMeasure::DalitzPair:      return 5;
        case PhaseSpaceMeasure::SolidAngleLab:   return 6;
        case PhaseSpaceMeasure::Unspecified:      return 7;
    }
    return 8;
}

} // anonymous namespace

int MultiChannelPhaseSpace::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    WarnOnIncompatibility();

    if (channels.empty()) {
        throw std::runtime_error("MultiChannelPhaseSpace has no channels");
    }

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
    WarnOnIncompatibility();

    // Plain summation -- all channels must report density in the same measure.
    double density = 0.0;
    for (size_t i = 0; i < channels.size(); ++i) {
        density += weights[i] * channels[i]->Density(detector_model, record);
    }
    return density;
}

PhaseSpaceTopology MultiChannelPhaseSpace::CommonTopology() const {
    if (channels.empty()) return PhaseSpaceTopology::Unspecified;

    PhaseSpaceTopology topo = channels.front()->Topology();
    for (size_t i = 1; i < channels.size(); ++i) {
        if (channels[i]->Topology() != topo) {
            throw std::runtime_error(
                "MultiChannelPhaseSpace topology mismatch: channel 0 ("
                + channels[0]->Name() + ") is "
                + PhaseSpaceTopologyName(topo)
                + " but channel " + std::to_string(i) + " ("
                + channels[i]->Name() + ") is "
                + PhaseSpaceTopologyName(channels[i]->Topology()));
        }
    }
    return topo;
}

PhaseSpaceMeasure MultiChannelPhaseSpace::CommonMeasure() const {
    if (channels.empty()) return PhaseSpaceMeasure::Unspecified;

    std::map<PhaseSpaceMeasure, int> counts;
    for (auto const & channel : channels) {
        counts[channel->Measure()] += 1;
    }

    PhaseSpaceMeasure best = channels.front()->Measure();
    int best_count = -1;
    int best_priority = 100;
    for (auto const & kv : counts) {
        int priority = MeasurePriority(kv.first);
        if (kv.second > best_count ||
            (kv.second == best_count && priority < best_priority)) {
            best = kv.first;
            best_count = kv.second;
            best_priority = priority;
        }
    }
    return best;
}

PhaseSpaceConvention MultiChannelPhaseSpace::CommonConvention() const {
    PhaseSpaceMeasure m = CommonMeasure();
    switch (m) {
        case PhaseSpaceMeasure::SolidAngleRest:  return PhaseSpaceConvention::RestFrameSolidAngle;
        case PhaseSpaceMeasure::SolidAngleLab:   return PhaseSpaceConvention::LabFrameSolidAngle;
        case PhaseSpaceMeasure::Recursive2Body:  return PhaseSpaceConvention::Recursive2Body;
        case PhaseSpaceMeasure::DalitzPair:      return PhaseSpaceConvention::Dalitz;
        case PhaseSpaceMeasure::HelicityAngles:  return PhaseSpaceConvention::HelicityAngles;
        case PhaseSpaceMeasure::MandelstamQ2:    return PhaseSpaceConvention::MandelstamST;
        case PhaseSpaceMeasure::BjorkenXY:       return PhaseSpaceConvention::BjorkenXY;
        case PhaseSpaceMeasure::Unspecified:      return PhaseSpaceConvention::Custom;
    }
    return PhaseSpaceConvention::Custom;
}

std::vector<std::string> MultiChannelPhaseSpace::ValidateChannels() const {
    std::vector<std::string> diagnostics;
    if (channels.empty()) return diagnostics;

    // 1. Check topology agreement
    PhaseSpaceTopology topo0 = channels.front()->Topology();
    for (size_t i = 1; i < channels.size(); ++i) {
        PhaseSpaceTopology topo_i = channels[i]->Topology();
        if (topo_i != topo0) {
            std::ostringstream oss;
            oss << "Topology mismatch: channel 0 (" << channels[0]->Name()
                << ") is " << PhaseSpaceTopologyName(topo0)
                << " but channel " << i << " (" << channels[i]->Name()
                << ") is " << PhaseSpaceTopologyName(topo_i);
            diagnostics.push_back(oss.str());
        }
    }
    if (!diagnostics.empty()) return diagnostics;

    // 2. Check measure agreement (strict: all must match)
    PhaseSpaceMeasure meas0 = channels.front()->Measure();
    for (size_t i = 1; i < channels.size(); ++i) {
        PhaseSpaceMeasure meas_i = channels[i]->Measure();
        if (meas_i != meas0) {
            std::ostringstream oss;
            oss << "Measure mismatch: channel 0 (" << channels[0]->Name()
                << ") uses " << PhaseSpaceMeasureName(meas0)
                << " but channel " << i << " (" << channels[i]->Name()
                << ") uses " << PhaseSpaceMeasureName(meas_i)
                << ". All channels must report density in the same "
                << "measure; no cross-channel conversion is performed.";
            diagnostics.push_back(oss.str());
        }
    }
    return diagnostics;
}

void MultiChannelPhaseSpace::WarnOnIncompatibility() const {
    if (compatibility_warning_emitted) return;
    auto diagnostics = ValidateChannels();
    if (!diagnostics.empty()) {
        for (auto const & d : diagnostics) {
            std::cerr << "SIREN phase-space: " << d << std::endl;
        }
        compatibility_warning_emitted = true;
    }
}

std::vector<std::string> MultiChannelPhaseSpace::ValidateChannelDensities(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord template_record,
    int samples_per_channel) const
{
    std::vector<std::string> diagnostics;
    auto compat = ValidateChannels();
    diagnostics.insert(diagnostics.end(), compat.begin(), compat.end());

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
