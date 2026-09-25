#include "SIREN/injection/RestFrameEnvelope2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Errors.h"
#include "InteractionRecordUtils.h"
#include "OnShellDecayKinematics.h"
#include "RestFrameEnvelope.h"

namespace siren { namespace injection {
namespace {
detail::RestFrameEnvelope Envelope(
    siren::dataclasses::InteractionRecord const & r, int index,
    siren::geometry::Geometry const & target)
{
    auto parent=detail::OnShellParent(r);
    return detail::BuildRestFrameEnvelope(parent.e,parent.p,r.primary_mass,
        TwoBodyRestEnergy(r.primary_mass,r.secondary_masses[index],r.secondary_masses[1-index]),
        TwoBodyRestMomentum(r.primary_mass,r.secondary_masses[index],r.secondary_masses[1-index]),
        detail::ReadVertex(r),target);
}
bool Allowed(siren::dataclasses::InteractionRecord const & r) {
    return detail::OnShellParentValid(r) && detail::HasSecondaryStorage(r,2) && r.primary_mass>r.secondary_masses[0]+r.secondary_masses[1]
        && r.secondary_masses[0]>=0 && r.secondary_masses[1]>=0;
}
}

RestFrameEnvelope2BodyChannel::RestFrameEnvelope2BodyChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,int index)
    :target_(std::move(target)),daughter_index_(index)
{
    if (!target_ || index<0 || index>1)
        throw std::invalid_argument("RestFrameEnvelope2BodyChannel requires a target and daughter index 0 or 1");
}

bool RestFrameEnvelope2BodyChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & r) const
{ return Allowed(r) && Envelope(r,daughter_index_,*target_).active; }

void RestFrameEnvelope2BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & r) const
{
    detail::RequireSecondaryStorage(r,2,Name().c_str());
    detail::RequireOnShellParent(r);
    if (!Allowed(r)) throw siren::utilities::InjectionFailure(
        siren::utilities::FailureReason::KinematicallyForbidden,"Rest-frame envelope has no two-body phase space");
    auto box=Envelope(r,daughter_index_,*target_);
    double draw=random->Uniform(0,1)*box.length;
    double u=box.intervals.back().second;
    for (auto const & interval : box.intervals) {
        double length=interval.second-interval.first;
        if (draw<=length) {u=interval.first+draw;break;}
        draw-=length;
    }
    u=std::clamp(u,-1.0,1.0);
    double phi=(2*random->Uniform(0,1)-1)*box.half_phi;
    double st=std::sqrt(std::max(0.0,1-u*u));
    double p=TwoBodyRestMomentum(r.primary_mass,r.secondary_masses[0],r.secondary_masses[1]);
    double e=TwoBodyRestEnergy(r.primary_mass,r.secondary_masses[daughter_index_],r.secondary_masses[1-daughter_index_]);
    auto v=(box.axis*u+box.first*(st*std::cos(phi))+box.second*(st*std::sin(phi)))*p;
    auto parent=detail::OnShellParent(r);
    auto daughter=detail::BoostOnShell(parent,r.primary_mass,{e,v});
    detail::WriteSecondary(r,daughter_index_,daughter);
    detail::WriteSecondary(r,1-daughter_index_,{parent.e-daughter.e,parent.p-daughter.p});
}

double RestFrameEnvelope2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & r) const
{
    if (!Allowed(r)) return 0;
    auto box=Envelope(r,daughter_index_,*target_);
    if (!box.active) return 1/(4*M_PI);
    auto parent=detail::OnShellParent(r), daughter=detail::ReadSecondary(r,daughter_index_);
    auto rest=detail::BoostOnShell(parent,r.primary_mass,daughter,true);
    double rest_parallel=siren::math::scalar_product(rest.p,box.axis);
    double p=TwoBodyRestMomentum(r.primary_mass,r.secondary_masses[0],r.secondary_masses[1]);
    double u=rest_parallel/p;
    double phi=std::atan2(siren::math::scalar_product(daughter.p,box.second),
                          siren::math::scalar_product(daughter.p,box.first));
    return box.Contains(u,phi) ? 1/box.Area() : 0;
}
}}
