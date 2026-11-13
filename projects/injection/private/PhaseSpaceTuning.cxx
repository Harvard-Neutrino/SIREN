#include "PhaseSpaceTuning.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/utilities/Errors.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <set>

namespace siren { namespace injection {
namespace {
using Mixture = MultiChannelPhaseSpace;

std::vector<Mixture *> TuningTree(Mixture * root, bool recurse) {
    std::vector<Mixture *> out;
    std::set<Mixture *> seen;
    std::function<void(Mixture *)> visit = [&](Mixture * m) {
        if(!seen.insert(m).second)
            throw std::invalid_argument("Tuning requires an unshared, acyclic mixture tree");
        out.push_back(m);
        if(recurse) for(auto const & ch : m->Channels()) {
            auto nested = std::dynamic_pointer_cast<NestedMixtureChannel>(ch);
            if(nested) {
                if(!nested->mixture) throw std::invalid_argument("Null nested mixture");
                visit(nested->mixture.get());
            }
        }
    };
    visit(root);
    return out;
}

double LogAdd(double a, double b) {
    if(a == -std::numeric_limits<double>::infinity()) return b;
    if(b == -std::numeric_limits<double>::infinity()) return a;
    double hi = std::max(a, b);
    return hi + std::log1p(std::exp(std::min(a, b) - hi));
}
}

void MultiChannelPhaseSpace::CheckTuningState() const {
    if(tuning_count_ && (weights_ != tuning_weights_ || channels_ != tuning_channels_))
        throw std::invalid_argument("Mixture changed during tuning; reset accumulators first");
}

void MultiChannelPhaseSpace::Accumulate(
        std::shared_ptr<siren::detector::DetectorModel const> detector,
        siren::dataclasses::InteractionRecord const & record,
        double weight, bool recurse) {
    if(!std::isfinite(weight))
        throw siren::utilities::WeightCalculationError("Non-finite tuning weight");
    auto mixtures = TuningTree(this, recurse);
    for(auto m : mixtures) m->CheckTuningState();
    std::vector<double> bare;
    auto const conventions = RequireCompatibleConvention();
    double g = ComputeDensity(detector, record, conventions, nullptr, &bare);
    if(!std::isfinite(g) || g <= 0)
        throw siren::utilities::WeightCalculationError("Tuning point has no proposal support");
    auto topology = conventions.common.topology;
    auto measure = conventions.common.measure;
    std::vector<std::vector<double>> updates;
    for(auto m : mixtures) {
        if(m != this) {
            auto const local_conventions = m->RequireCompatibleConvention();
            if(local_conventions.common.topology != topology)
                throw siren::utilities::MeasureCompatibilityError("Nested tuning topology mismatch");
            m->ComputeDensity(detector, record, local_conventions, nullptr, &bare);
            auto local = local_conventions.common.measure;
            for(auto & d : bare) d = ConvertDensity(d, local, measure, topology, record);
        }
        std::vector<double> sums = m->tuning_log_sums_;
        if(!m->tuning_count_) sums.assign(m->Channels().size(), -std::numeric_limits<double>::infinity());
        if(m->tuning_count_ == std::numeric_limits<std::uint64_t>::max())
            throw std::overflow_error("Tuning count overflow");
        for(std::size_t i = 0; i < bare.size(); ++i) {
            if(!std::isfinite(bare[i]) || bare[i] < 0)
                throw siren::utilities::WeightCalculationError("Invalid tuning channel density");
            if(weight != 0 && bare[i] > 0) {
                // Log sums avoid overflow in w^2 and retain very small contributions.
                double term = 2 * std::log(std::abs(weight)) + std::log(bare[i]) - std::log(g);
                sums[i] = LogAdd(sums[i], term);
            }
        }
        updates.push_back(std::move(sums));
    }
    for(std::size_t i = 0; i < mixtures.size(); ++i) {
        auto m = mixtures[i];
        m->tuning_log_sums_ = std::move(updates[i]);
        m->tuning_weights_ = m->Weights();
        m->tuning_channels_ = m->Channels();
        ++m->tuning_count_;
    }
}

void MultiChannelPhaseSpace::UpdateWeights(std::string const & rule,
        double damping, double min_weight, bool recurse) {
    if(rule != "alpha_sqrt_W" && rule != "sqrt_W")
        throw std::invalid_argument("Unknown tuning rule");
    if(!std::isfinite(damping) || damping < 0 || damping > 1 ||
       !std::isfinite(min_weight) || min_weight <= 0)
        throw std::invalid_argument("Tuning requires 0 <= damping <= 1 and a positive min_weight");
    auto mixtures = TuningTree(this, recurse);
    std::vector<std::vector<double>> proposals;
    for(auto m : mixtures) {
        m->CheckTuningState();
        std::size_t n = m->Weights().size();
        if(n * min_weight > 1)
            throw std::invalid_argument("min_weight exceeds the mixture simplex");
        auto candidate = m->Weights();
        if(m->tuning_count_) {
            std::vector<double> logs(n);
            for(std::size_t i = 0; i < n; ++i) {
                logs[i] = m->tuning_log_sums_[i] / 2;
                if(rule == "alpha_sqrt_W") logs[i] += std::log(m->Weights()[i]);
            }
            double peak = *std::max_element(logs.begin(), logs.end());
            if(std::isfinite(peak)) {
                double total = 0;
                for(std::size_t i = 0; i < n; ++i) total += (candidate[i] = std::exp(logs[i] - peak));
                for(std::size_t i = 0; i < n; ++i) {
                    double target = min_weight + (1 - n * min_weight) * candidate[i] / total;
                    candidate[i] = damping * target + (1 - damping) * m->Weights()[i];
                }
            }
        }
        // The floor is a support guarantee, including zero-statistic pilots.
        double free = 0;
        for(auto & w : candidate) { w = std::max(0.0, w - min_weight); free += w; }
        for(auto & w : candidate)
            w = min_weight + (1 - n * min_weight) * (free > 0 ? w / free : 1.0 / n);
        proposals.push_back(std::move(candidate));
    }
    for(std::size_t i = 0; i < mixtures.size(); ++i) {
        mixtures[i]->weights_ = std::move(proposals[i]);
        mixtures[i]->ResetAccumulators(false);
    }
}

void CopyTunedWeights(MultiChannelPhaseSpace & destination,
                      MultiChannelPhaseSpace const & source) {
    if (source.tuning_count_ || destination.tuning_count_) {
        throw std::invalid_argument("Finish or reset mixture tuning before copying weights");
    }

    auto const source_conventions = source.RequireCompatibleConvention();
    auto const destination_conventions = destination.RequireCompatibleConvention();
    if (!(source_conventions.common == destination_conventions.common)
            || source_conventions.measures != destination_conventions.measures) {
        throw std::invalid_argument("Tuned weights require matching mixture conventions and channel counts");
    }
    source.ValidateArchivedWeights();
    auto weights = source.weights_;
    destination.weights_ = std::move(weights);
}

void MultiChannelPhaseSpace::ResetAccumulators(bool recurse) {
    for(auto m : TuningTree(this, recurse)) {
        m->tuning_count_ = 0;
        m->tuning_log_sums_.clear();
        m->tuning_weights_.clear();
        m->tuning_channels_.clear();
    }
}
}}
