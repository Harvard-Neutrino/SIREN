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
//  Measure conversion (density Jacobians)                            //
// ================================================================ //

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

double ConvertDensity(
    double density,
    PhaseSpaceMeasure from,
    PhaseSpaceMeasure to,
    PhaseSpaceTopology topology,
    siren::dataclasses::InteractionRecord const & record)
{
    if (from == to || density == 0.0) return density;
    namespace J = siren::injection::phase_space_jacobian;

    // ---- Decay2Body: SolidAngleRest <-> SolidAngleLab ----

    if (topology == PhaseSpaceTopology::Decay2Body &&
        ((from == PhaseSpaceMeasure::SolidAngleRest &&
          to == PhaseSpaceMeasure::SolidAngleLab) ||
         (from == PhaseSpaceMeasure::SolidAngleLab &&
          to == PhaseSpaceMeasure::SolidAngleRest)))
    {
        if (record.secondary_masses.size() < 2) return density;

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

        // Match the correct Lorentz-boost solution
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

        if (from == PhaseSpaceMeasure::SolidAngleRest) {
            return density / best_J;   // rest -> lab
        } else {
            return density * best_J;   // lab -> rest
        }
    }

    // ---- Scatter2to2: MandelstamQ2 <-> BjorkenXY ----

    if (topology == PhaseSpaceTopology::Scatter2to2) {
        if (from == PhaseSpaceMeasure::BjorkenXY &&
            to == PhaseSpaceMeasure::MandelstamQ2) {
            auto it = record.interaction_parameters.find("bjorken_y");
            if (it == record.interaction_parameters.end())
                it = record.interaction_parameters.find("y");
            if (it == record.interaction_parameters.end()) return density;
            return J::BjorkenXYDensityToQ2YDensity(
                density, it->second, record.target_mass,
                record.primary_momentum[0]);
        }
        if (from == PhaseSpaceMeasure::MandelstamQ2 &&
            to == PhaseSpaceMeasure::BjorkenXY) {
            auto it = record.interaction_parameters.find("bjorken_y");
            if (it == record.interaction_parameters.end())
                it = record.interaction_parameters.find("y");
            if (it == record.interaction_parameters.end()) return density;
            return J::Q2YDensityToBjorkenXYDensity(
                density, it->second, record.target_mass,
                record.primary_momentum[0]);
        }

        // SolidAngleRest <-> MandelstamQ2
        if (from == PhaseSpaceMeasure::SolidAngleRest &&
            to == PhaseSpaceMeasure::MandelstamQ2) {
            double E = record.primary_momentum[0];
            double s = record.primary_mass * record.primary_mass
                     + record.target_mass * record.target_mass
                     + 2.0 * record.target_mass * E;
            return J::SolidAngleRestDensityToMandelstamQ2Density(
                density, s, record.primary_mass, record.target_mass);
        }
        if (from == PhaseSpaceMeasure::MandelstamQ2 &&
            to == PhaseSpaceMeasure::SolidAngleRest) {
            double E = record.primary_momentum[0];
            double s = record.primary_mass * record.primary_mass
                     + record.target_mass * record.target_mass
                     + 2.0 * record.target_mass * E;
            return J::MandelstamQ2DensityToSolidAngleRestDensity(
                density, s, record.primary_mass, record.target_mass);
        }

        // SolidAngleRest <-> BjorkenXY: compose through MandelstamQ2
        if (from == PhaseSpaceMeasure::SolidAngleRest &&
            to == PhaseSpaceMeasure::BjorkenXY) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::MandelstamQ2, topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::MandelstamQ2,
                to, topology, record);
        }
        if (from == PhaseSpaceMeasure::BjorkenXY &&
            to == PhaseSpaceMeasure::SolidAngleRest) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::MandelstamQ2, topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::MandelstamQ2,
                to, topology, record);
        }

        // SolidAngleLab <-> SolidAngleRest for scattering:
        // same Jacobian as Decay2Body (Lorentz boost)
        if ((from == PhaseSpaceMeasure::SolidAngleRest &&
             to == PhaseSpaceMeasure::SolidAngleLab) ||
            (from == PhaseSpaceMeasure::SolidAngleLab &&
             to == PhaseSpaceMeasure::SolidAngleRest)) {
            // Delegate to the Decay2Body path which does the boost
            return ConvertDensity(density, from, to,
                PhaseSpaceTopology::Decay2Body, record);
        }

        // SolidAngleLab <-> MandelstamQ2/BjorkenXY: compose
        if (from == PhaseSpaceMeasure::SolidAngleLab &&
            (to == PhaseSpaceMeasure::MandelstamQ2 ||
             to == PhaseSpaceMeasure::BjorkenXY)) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::SolidAngleRest, topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::SolidAngleRest,
                to, topology, record);
        }
        if ((from == PhaseSpaceMeasure::MandelstamQ2 ||
             from == PhaseSpaceMeasure::BjorkenXY) &&
            to == PhaseSpaceMeasure::SolidAngleLab) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::SolidAngleRest, topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::SolidAngleRest,
                to, topology, record);
        }
    }

    // ---- Decay3Body / Scatter2to3: Recursive2Body <-> DalitzPair <-> HelicityAngles ----

    if (topology == PhaseSpaceTopology::Decay3Body ||
        topology == PhaseSpaceTopology::Scatter2to3) {

        // Recursive2Body <-> HelicityAngles: trivial
        if ((from == PhaseSpaceMeasure::Recursive2Body &&
             to == PhaseSpaceMeasure::HelicityAngles) ||
            (from == PhaseSpaceMeasure::HelicityAngles &&
             to == PhaseSpaceMeasure::Recursive2Body)) {
            return density;
        }

        // Recursive2Body <-> DalitzPair
        auto compute_s_pair = [&]() -> double {
            if (record.secondary_masses.size() < 3) return -1.0;
            auto const & p1 = record.secondary_momenta[1];
            auto const & p2 = record.secondary_momenta[2];
            double E = p1[0] + p2[0];
            double px = p1[1] + p2[1];
            double py = p1[2] + p2[2];
            double pz = p1[3] + p2[3];
            return E*E - px*px - py*py - pz*pz;
        };

        if (from == PhaseSpaceMeasure::Recursive2Body &&
            to == PhaseSpaceMeasure::DalitzPair) {
            double s_pair = compute_s_pair();
            if (s_pair < 0) return density;
            return J::Recursive2BodyDensityToDalitzDensity(
                density, record.primary_mass,
                record.secondary_masses[0],
                record.secondary_masses[1],
                record.secondary_masses[2], s_pair);
        }
        if (from == PhaseSpaceMeasure::DalitzPair &&
            to == PhaseSpaceMeasure::Recursive2Body) {
            double s_pair = compute_s_pair();
            if (s_pair < 0) return density;
            return J::DalitzDensityToRecursive2BodyDensity(
                density, record.primary_mass,
                record.secondary_masses[0],
                record.secondary_masses[1],
                record.secondary_masses[2], s_pair);
        }

        // HelicityAngles <-> DalitzPair: compose through Recursive2Body
        if (from == PhaseSpaceMeasure::HelicityAngles &&
            to == PhaseSpaceMeasure::DalitzPair) {
            return ConvertDensity(density,
                PhaseSpaceMeasure::Recursive2Body,
                PhaseSpaceMeasure::DalitzPair, topology, record);
        }
        if (from == PhaseSpaceMeasure::DalitzPair &&
            to == PhaseSpaceMeasure::HelicityAngles) {
            return ConvertDensity(density,
                PhaseSpaceMeasure::DalitzPair,
                PhaseSpaceMeasure::Recursive2Body, topology, record);
        }
    }

    return density;
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
            ch_meas != PhaseSpaceMeasure::Unspecified) {
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

    // 2. Check measure compatibility
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

} // namespace injection
} // namespace siren
