#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace siren {
namespace injection {

using MType = PhaseSpaceMeasure::Type;

// ================================================================ //
//  PhaseSpaceChannel default Convention() from Topology+Measure      //
// ================================================================ //

PhaseSpaceConvention PhaseSpaceChannel::Convention() const {
    PhaseSpaceMeasure m = Measure();
    switch (m.type) {
        case MType::SolidAngleRest:  return PhaseSpaceConvention::RestFrameSolidAngle;
        case MType::SolidAngleLab:   return PhaseSpaceConvention::LabFrameSolidAngle;
        case MType::Recursive2Body:  return PhaseSpaceConvention::Recursive2Body;
        case MType::DalitzPair:      return PhaseSpaceConvention::Dalitz;
        case MType::HelicityAngles:  return PhaseSpaceConvention::HelicityAngles;
        case MType::MandelstamQ2:    return PhaseSpaceConvention::MandelstamST;
        case MType::BjorkenXY:       return PhaseSpaceConvention::BjorkenXY;
        case MType::Unspecified:     return PhaseSpaceConvention::Custom;
    }
    return PhaseSpaceConvention::Custom;
}

// ================================================================ //
//  Measure conversion (density Jacobians)                            //
// ================================================================ //

namespace {

int MeasurePriority(PhaseSpaceMeasure const & m) {
    switch (m.type) {
        case MType::SolidAngleRest:  return 0;
        case MType::Recursive2Body:  return 1;
        case MType::MandelstamQ2:    return 2;
        case MType::BjorkenXY:       return 3;
        case MType::HelicityAngles:  return 4;
        case MType::DalitzPair:      return 5;
        case MType::SolidAngleLab:   return 6;
        case MType::Unspecified:     return 7;
    }
    return 8;
}

// Signal a conversion ConvertDensity cannot perform, instead of silently
// returning the input density unchanged (former gap G6).  The from==to
// short-circuit handles "no conversion needed"; this throw covers
// "conversion not implemented" and "required inputs missing", which a
// caller forming w = f/g must never silently absorb.
[[noreturn]] void ThrowUnconvertible(
    const char * reason,
    PhaseSpaceMeasure const & from,
    PhaseSpaceMeasure const & to,
    PhaseSpaceTopology topology)
{
    std::ostringstream oss;
    oss << "ConvertDensity: " << reason << " (from "
        << PhaseSpaceMeasureName(from) << " to " << PhaseSpaceMeasureName(to)
        << " in " << PhaseSpaceTopologyName(topology) << " topology)";
    throw std::runtime_error(oss.str());
}

double ConvertDensity(
    double density,
    PhaseSpaceMeasure const & from,
    PhaseSpaceMeasure const & to,
    PhaseSpaceTopology topology,
    siren::dataclasses::InteractionRecord const & record)
{
    if (from == to || density == 0.0) return density;
    namespace J = siren::injection::phase_space_jacobian;

    // ---- Decay2Body: SolidAngleRest <-> SolidAngleLab ----

    if (topology == PhaseSpaceTopology::Decay2Body &&
        ((from.type == MType::SolidAngleRest &&
          to.type == MType::SolidAngleLab) ||
         (from.type == MType::SolidAngleLab &&
          to.type == MType::SolidAngleRest)))
    {
        if (record.secondary_masses.size() < 2)
            ThrowUnconvertible(
                "Decay2Body rest<->lab requires 2 secondary masses",
                from, to, topology);

        double M = record.primary_mass;
        double m_A = record.secondary_masses[0];
        double m_B = record.secondary_masses[1];
        double p_rest = TwoBodyRestMomentum(M, m_A, m_B);
        double E_rest = TwoBodyRestEnergy(M, m_A, m_B);

        double E_parent = record.primary_momentum[0];
        double px = record.primary_momentum[1];
        double py = record.primary_momentum[2];
        double pz = record.primary_momentum[3];
        double p_parent = std::sqrt(px*px + py*py + pz*pz);
        if (p_parent < 1e-15) return density;

        double beta_parent = p_parent / E_parent;
        double gamma_parent = E_parent / M;

        double px_A = record.secondary_momenta[0][1];
        double py_A = record.secondary_momenta[0][2];
        double pz_A = record.secondary_momenta[0][3];
        double p_A = std::sqrt(px_A*px_A + py_A*py_A + pz_A*pz_A);
        if (p_A < 1e-15) return density;

        double cos_theta_lab =
            (px_A*px + py_A*py + pz_A*pz) / (p_A * p_parent);

        double E_A_lab = record.secondary_momenta[0][0];
        double p_par_lab = p_A * cos_theta_lab;
        double p_par_rest = gamma_parent * (p_par_lab - beta_parent * E_A_lab);
        double cos_theta_rest_actual = (p_rest > 0) ? p_par_rest / p_rest : 0.0;

        auto solutions = SolveLabAngle(
            beta_parent, gamma_parent, p_rest, E_rest,
            m_A, cos_theta_lab);

        double best_J = 0.0;
        double best_dist = 1e30;
        for (auto const & sol : solutions) {
            if (!sol.valid) continue;
            double dist = std::abs(sol.cos_theta_rest - cos_theta_rest_actual);
            if (dist < best_dist) {
                best_dist = dist;
                best_J = sol.jacobian;
            }
        }
        if (best_J <= 0.0 || !std::isfinite(best_J)) return density;

        if (from.type == MType::SolidAngleRest) {
            return density / best_J;
        } else {
            return density * best_J;
        }
    }

    // ---- Scatter2to2 conversions ----

    if (topology == PhaseSpaceTopology::Scatter2to2) {
        // MandelstamQ2 <-> BjorkenXY
        if (from.type == MType::BjorkenXY && to.type == MType::MandelstamQ2) {
            auto it = record.interaction_parameters.find("bjorken_y");
            if (it == record.interaction_parameters.end())
                it = record.interaction_parameters.find("y");
            if (it == record.interaction_parameters.end())
                ThrowUnconvertible(
                    "requires 'bjorken_y' (or 'y') in interaction_parameters",
                    from, to, topology);
            return J::BjorkenXYDensityToQ2YDensity(
                density, it->second, record.target_mass,
                record.primary_momentum[0]);
        }
        if (from.type == MType::MandelstamQ2 && to.type == MType::BjorkenXY) {
            auto it = record.interaction_parameters.find("bjorken_y");
            if (it == record.interaction_parameters.end())
                it = record.interaction_parameters.find("y");
            if (it == record.interaction_parameters.end())
                ThrowUnconvertible(
                    "requires 'bjorken_y' (or 'y') in interaction_parameters",
                    from, to, topology);
            return J::Q2YDensityToBjorkenXYDensity(
                density, it->second, record.target_mass,
                record.primary_momentum[0]);
        }

        // SolidAngleRest <-> MandelstamQ2
        if (from.type == MType::SolidAngleRest && to.type == MType::MandelstamQ2) {
            double E = record.primary_momentum[0];
            double s = record.primary_mass * record.primary_mass
                     + record.target_mass * record.target_mass
                     + 2.0 * record.target_mass * E;
            return J::SolidAngleRestDensityToMandelstamQ2Density(
                density, s, record.primary_mass, record.target_mass);
        }
        if (from.type == MType::MandelstamQ2 && to.type == MType::SolidAngleRest) {
            double E = record.primary_momentum[0];
            double s = record.primary_mass * record.primary_mass
                     + record.target_mass * record.target_mass
                     + 2.0 * record.target_mass * E;
            return J::MandelstamQ2DensityToSolidAngleRestDensity(
                density, s, record.primary_mass, record.target_mass);
        }

        // Compose remaining pairs through intermediates
        if (from.type == MType::SolidAngleRest && to.type == MType::BjorkenXY) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::MandelstamQ2(), topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::MandelstamQ2(),
                to, topology, record);
        }
        if (from.type == MType::BjorkenXY && to.type == MType::SolidAngleRest) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::MandelstamQ2(), topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::MandelstamQ2(),
                to, topology, record);
        }

        // SolidAngleLab <-> others: compose through SolidAngleRest
        if (from.type == MType::SolidAngleLab || to.type == MType::SolidAngleLab) {
            if (from.type == MType::SolidAngleLab) {
                double d = ConvertDensity(density, from,
                    PhaseSpaceMeasure::SolidAngleRest(),
                    PhaseSpaceTopology::Decay2Body, record);
                return ConvertDensity(d, PhaseSpaceMeasure::SolidAngleRest(),
                    to, topology, record);
            } else {
                double d = ConvertDensity(density, from,
                    PhaseSpaceMeasure::SolidAngleRest(), topology, record);
                return ConvertDensity(d, PhaseSpaceMeasure::SolidAngleRest(),
                    to, PhaseSpaceTopology::Decay2Body, record);
            }
        }
    }

    // ---- Decay3Body / Scatter2to3: 3-body measure conversions ----
    // Uses the measure's index fields to read the correct momenta.

    if (topology == PhaseSpaceTopology::Decay3Body ||
        topology == PhaseSpaceTopology::Scatter2to3) {

        // For scattering, the parent mass in the Jacobian is the CM
        // energy sqrt(s), not the beam particle mass.
        double parent_mass = record.primary_mass;
        if (topology == PhaseSpaceTopology::Scatter2to3) {
            double E = record.primary_momentum[0];
            double m_beam = record.primary_mass;
            double m_target = record.target_mass;
            double s = m_beam * m_beam + m_target * m_target + 2.0 * m_target * E;
            if (s > 0) parent_mass = std::sqrt(s);
        }

        auto compute_s_pair = [&](PhaseSpaceMeasure const & m) -> double {
            int i1 = m.pair_first;
            int i2 = m.pair_second;
            if (i1 < 0 || i2 < 0 ||
                i1 >= static_cast<int>(record.secondary_momenta.size()) ||
                i2 >= static_cast<int>(record.secondary_momenta.size())) {
                return -1.0;
            }
            auto const & p1 = record.secondary_momenta[i1];
            auto const & p2 = record.secondary_momenta[i2];
            double E = p1[0] + p2[0];
            double px = p1[1] + p2[1];
            double py = p1[2] + p2[2];
            double pz = p1[3] + p2[3];
            return E*E - px*px - py*py - pz*pz;
        };

        // Same type, different indices: cross-factorization conversion.
        // Route through Dalitz (index-independent) as intermediate.
        if (from.type == to.type && from != to) {
            // Convert from's factorization to Dalitz, then Dalitz to to's factorization.
            PhaseSpaceMeasure dalitz_from = PhaseSpaceMeasure::DalitzPair(
                from.spectator, from.pair_first, from.pair_second);
            PhaseSpaceMeasure dalitz_to = PhaseSpaceMeasure::DalitzPair(
                to.spectator, to.pair_first, to.pair_second);

            // Step 1: convert from's Recursive2Body to Dalitz using from's indices
            double s_pair_from = compute_s_pair(from);
            if (s_pair_from < 0)
                ThrowUnconvertible(
                    "invalid 'from' pair invariant mass (bad factorization indices)",
                    from, to, topology);
            double dalitz_density = J::Recursive2BodyDensityToDalitzDensity(
                density, parent_mass,
                record.secondary_masses[from.spectator],
                record.secondary_masses[from.pair_first],
                record.secondary_masses[from.pair_second], s_pair_from);

            // Step 2: convert Dalitz to to's Recursive2Body using to's indices
            double s_pair_to = compute_s_pair(to);
            if (s_pair_to < 0)
                ThrowUnconvertible(
                    "invalid 'to' pair invariant mass (bad factorization indices)",
                    from, to, topology);
            return J::DalitzDensityToRecursive2BodyDensity(
                dalitz_density, parent_mass,
                record.secondary_masses[to.spectator],
                record.secondary_masses[to.pair_first],
                record.secondary_masses[to.pair_second], s_pair_to);
        }

        // Recursive2Body <-> HelicityAngles: trivial (same measure, different angular reference)
        if ((from.type == MType::Recursive2Body && to.type == MType::HelicityAngles) ||
            (from.type == MType::HelicityAngles && to.type == MType::Recursive2Body)) {
            return density;
        }

        // Recursive2Body <-> DalitzPair
        if (from.type == MType::Recursive2Body && to.type == MType::DalitzPair) {
            double s_pair = compute_s_pair(from);
            if (s_pair < 0)
                ThrowUnconvertible(
                    "invalid pair invariant mass (bad factorization indices)",
                    from, to, topology);
            return J::Recursive2BodyDensityToDalitzDensity(
                density, parent_mass,
                record.secondary_masses[from.spectator],
                record.secondary_masses[from.pair_first],
                record.secondary_masses[from.pair_second], s_pair);
        }
        if (from.type == MType::DalitzPair && to.type == MType::Recursive2Body) {
            double s_pair = compute_s_pair(to);
            if (s_pair < 0)
                ThrowUnconvertible(
                    "invalid pair invariant mass (bad factorization indices)",
                    from, to, topology);
            return J::DalitzDensityToRecursive2BodyDensity(
                density, parent_mass,
                record.secondary_masses[to.spectator],
                record.secondary_masses[to.pair_first],
                record.secondary_masses[to.pair_second], s_pair);
        }

        // HelicityAngles <-> DalitzPair: compose through Recursive2Body
        if (from.type == MType::HelicityAngles && to.type == MType::DalitzPair) {
            PhaseSpaceMeasure intermediate = PhaseSpaceMeasure::Recursive2Body(
                from.spectator, from.pair_first, from.pair_second);
            return ConvertDensity(density, intermediate, to, topology, record);
        }
        if (from.type == MType::DalitzPair && to.type == MType::HelicityAngles) {
            PhaseSpaceMeasure intermediate = PhaseSpaceMeasure::Recursive2Body(
                to.spectator, to.pair_first, to.pair_second);
            return ConvertDensity(density, from, intermediate, topology, record);
        }
    }

    ThrowUnconvertible(
        "no conversion implemented for this measure pair and topology",
        from, to, topology);
}

} // anonymous namespace

// ================================================================ //
//  MultiChannelPhaseSpace                                            //
// ================================================================ //

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

    PhaseSpaceTopology common_topo = CommonTopology();
    PhaseSpaceMeasure common_meas = CommonMeasure();

    double density = 0.0;
    for (size_t i = 0; i < channels.size(); ++i) {
        double d = channels[i]->Density(detector_model, record);
        PhaseSpaceMeasure ch_meas = channels[i]->Measure();
        if (ch_meas != common_meas &&
            ch_meas.type != MType::Unspecified) {
            d = ConvertDensity(d, ch_meas, common_meas, common_topo, record);
        }
        density += weights[i] * d;
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
    if (channels.empty()) return PhaseSpaceMeasure::Unspecified();

    // Group by measure type (ignoring indices) for majority vote
    std::map<MType, int> type_counts;
    for (auto const & channel : channels) {
        type_counts[channel->Measure().type] += 1;
    }

    MType best_type = channels.front()->Measure().type;
    int best_count = -1;
    int best_priority = 100;
    for (auto const & kv : type_counts) {
        PhaseSpaceMeasure tmp;
        tmp.type = kv.first;
        int priority = MeasurePriority(tmp);
        if (kv.second > best_count ||
            (kv.second == best_count && priority < best_priority)) {
            best_type = kv.first;
            best_count = kv.second;
            best_priority = priority;
        }
    }

    // Return the full measure (with indices) from the first channel
    // that has the winning type.
    for (auto const & channel : channels) {
        if (channel->Measure().type == best_type) {
            return channel->Measure();
        }
    }
    return PhaseSpaceMeasure::Unspecified();
}

PhaseSpaceConvention MultiChannelPhaseSpace::CommonConvention() const {
    PhaseSpaceMeasure m = CommonMeasure();
    switch (m.type) {
        case MType::SolidAngleRest:  return PhaseSpaceConvention::RestFrameSolidAngle;
        case MType::SolidAngleLab:   return PhaseSpaceConvention::LabFrameSolidAngle;
        case MType::Recursive2Body:  return PhaseSpaceConvention::Recursive2Body;
        case MType::DalitzPair:      return PhaseSpaceConvention::Dalitz;
        case MType::HelicityAngles:  return PhaseSpaceConvention::HelicityAngles;
        case MType::MandelstamQ2:    return PhaseSpaceConvention::MandelstamST;
        case MType::BjorkenXY:       return PhaseSpaceConvention::BjorkenXY;
        case MType::Unspecified:     return PhaseSpaceConvention::Custom;
    }
    return PhaseSpaceConvention::Custom;
}

std::vector<std::string> MultiChannelPhaseSpace::ValidateChannels() const {
    std::vector<std::string> diagnostics;
    if (channels.empty()) return diagnostics;

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

    PhaseSpaceMeasure common = CommonMeasure();
    for (size_t i = 0; i < channels.size(); ++i) {
        PhaseSpaceMeasure meas_i = channels[i]->Measure();
        if (meas_i == common) continue;
        if (!PhaseSpaceCompatible(topo0, common, topo0, meas_i)) {
            std::ostringstream oss;
            oss << "Measure incompatibility: channel " << i
                << " (" << channels[i]->Name()
                << ") uses " << PhaseSpaceMeasureName(meas_i)
                << " which is not convertible to "
                << PhaseSpaceMeasureName(common)
                << " within " << PhaseSpaceTopologyName(topo0)
                << " topology";
            diagnostics.push_back(oss.str());
        } else {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << ") uses " << PhaseSpaceMeasureName(meas_i)
                << "; will auto-convert to "
                << PhaseSpaceMeasureName(common);
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

// ================================================================ //
//  NestedMixtureChannel                                              //
// ================================================================ //

void NestedMixtureChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    if (!mixture) {
        throw std::runtime_error("NestedMixtureChannel has no inner mixture");
    }
    mixture->Sample(random, detector_model, record);
}

double NestedMixtureChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!mixture) return 0.0;
    return mixture->Density(detector_model, record);
}

std::string NestedMixtureChannel::Name() const {
    return label;
}

PhaseSpaceTopology NestedMixtureChannel::Topology() const {
    if (!mixture) return PhaseSpaceTopology::Unspecified;
    return mixture->CommonTopology();
}

PhaseSpaceMeasure NestedMixtureChannel::Measure() const {
    if (!mixture) return PhaseSpaceMeasure::Unspecified();
    return mixture->CommonMeasure();
}

bool NestedMixtureChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const {
    // The group genuinely directs iff any member does.  A member with no
    // fallback notion (physical / isotropic) returns true by default, so a group
    // containing one is always active; a group of directed channels all in their
    // isotropic fallback returns false, letting the optimizer discount it as one.
    if (!mixture) return true;
    for (auto const & channel : mixture->channels) {
        if (channel && channel->DirectingActive(record)) return true;
    }
    return false;
}

} // namespace injection
} // namespace siren
