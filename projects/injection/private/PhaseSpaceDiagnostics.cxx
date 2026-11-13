#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Errors.h"

#include <cmath>
#include <sstream>

namespace siren {
namespace injection {

std::vector<std::string> ValidateChannelDensities(
    MultiChannelPhaseSpace const & mixture,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord template_record,
    int samples_per_channel)
{
    std::vector<std::string> diagnostics;
    for (auto const & diagnostic : mixture.ValidateChannelsDetailed())
        diagnostics.push_back(diagnostic.message);

    auto const & channels = mixture.Channels();
    auto convention = mixture.CommonConvention();
    for (size_t i = 0; i < channels.size(); ++i) {
        int n_zero = 0;
        int n_nan = 0;
        int n_neg = 0;
        int n_unsampleable = 0;
        double max_ratio = 0.0;

        for (int s = 0; s < samples_per_channel; ++s) {
            siren::dataclasses::InteractionRecord record = template_record;
            try {
                channels[i]->Sample(random, detector_model, record);
            } catch (siren::utilities::InjectionFailure const &) {
                ++n_unsampleable;
                continue;
            }

            double d_self = ConvertDensity(
                channels[i]->Density(detector_model, record),
                channels[i]->Measure(), convention.measure,
                convention.topology, record);
            if (d_self <= 0 || !std::isfinite(d_self)) continue;

            for (size_t j = 0; j < channels.size(); ++j) {
                if (i == j) continue;
                double d_other = ConvertDensity(
                    channels[j]->Density(detector_model, record),
                    channels[j]->Measure(), convention.measure,
                    convention.topology, record);

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
        if (n_unsampleable == samples_per_channel) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): the probe template lies outside the channel's "
                << "support (" << n_unsampleable << "/" << samples_per_channel
                << " draws raised InjectionFailure); density probe skipped "
                << "for this channel";
            diagnostics.push_back(oss.str());
        }
    }

    return diagnostics;
}

} // namespace injection
} // namespace siren
