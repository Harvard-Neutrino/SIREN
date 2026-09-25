#include "SIREN/injection/OnShellCascadeChannel.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/RestFrameEnvelope2BodyChannel.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/utilities/Random.h"
#include "InteractionRecordUtils.h"
#include "OnShellDecayKinematics.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/utilities/Errors.h"
#include <limits>
#include <algorithm>
#include <cmath>

namespace siren { namespace injection {
namespace {
using Record=siren::dataclasses::InteractionRecord;
using V=siren::math::Vector3D;
using P4=detail::FourVector;
bool Allowed(Record const & r,double mass) {
    return detail::OnShellParentValid(r) && detail::HasSecondaryStorage(r,3) && r.secondary_masses[0]>=0
        && r.secondary_masses[1]>=0 && r.secondary_masses[1]==r.secondary_masses[2]
        && mass>2*r.secondary_masses[1] && r.primary_mass>mass+r.secondary_masses[0];
}
Record Effective(Record const & r,double energy,int daughter) {
    Record e=r;
    e.primary_momentum[0]=detail::OnShellParent(r).e;
    double m=r.secondary_masses[daughter], M=r.primary_mass;
    double recoil2=M*M+m*m-2*M*energy;
    if (recoil2 < -64*std::numeric_limits<double>::epsilon()*(M*M+m*m+2*M*energy))
        throw siren::utilities::InjectionFailure(siren::utilities::FailureReason::KinematicallyForbidden,
            "OnShellCascade: invalid effective recoil mass");
    recoil2=std::max(0.0,recoil2);
    e.signature.secondary_types.resize(2);
    e.secondary_masses={m,std::sqrt(recoil2)};
    e.secondary_momenta.resize(2);
    e.secondary_helicities.resize(2);
    auto v=detail::ReadSecondary(r,daughter), p=detail::OnShellParent(r);
    detail::WriteSecondary(e,0,v);
    detail::WriteSecondary(e,1,{p.e-v.e,p.p-v.p});
    return e;
}
// Parent-frame energy of pair daughter d, rebuilt from its lab vector. The lab
// vector carries the boost error eps*gamma*EA, which inversion amplifies by up
// to 2*gamma; the other terms bound the rounding of the inversion itself.
// Within that error an energy just outside EA/2 +- PA*pd/pair_mass is clamped
// onto the kinematic range, so the effective recoil mass stays physical.
bool RestDaughterEnergy(Record const & r,P4 const & parent,int d,double pair_mass,double & energy) {
    double M=r.primary_mass, gamma=parent.e/M, eps=std::numeric_limits<double>::epsilon();
    double EA=TwoBodyRestEnergy(M,pair_mass,r.secondary_masses[0]);
    double PA=TwoBodyRestMomentum(M,pair_mass,r.secondary_masses[0]);
    double pd=TwoBodyRestMomentum(pair_mass,r.secondary_masses[1],r.secondary_masses[2]);
    double scale=PA*pd/pair_mass;
    auto lab=detail::ReadSecondary(r,d);
    energy=detail::BoostOnShell(parent,M,lab,true).e;
    double error=std::max({128*eps*std::max(gamma*std::abs(lab.e),EA),256*eps*gamma*gamma*EA,1e-12*scale});
    if (!(std::abs(energy-EA/2)<=scale+error)) return false;
    energy=std::clamp(energy,EA/2-scale,EA/2+scale);
    return true;
}
double OrientationDensity(Record const & e,std::array<double,3> const & weights,
    PhaseSpaceChannel const * envelope,PhaseSpaceChannel const * directed)
{
    double q=weights[0]/(4*M_PI);
    if (weights[1]>0) q+=weights[1]*envelope->Density(nullptr,e);
    if (weights[2]>0) q+=weights[2]*directed->Density(nullptr,e);
    return q;
}
}

OnShellCascadeChannel::OnShellCascadeChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,double mass,double kappa,
    std::array<double,3> weights,double rho,double volume)
    :target_(std::move(target)),pair_mass_(mass),kappa_(kappa),rho_(rho),weights_(weights),target_volume_(volume)
{
    Validate(false);
    double sum=weights_[0]+weights_[1]+weights_[2];
    for (auto & w:weights_) w/=sum;
    if (weights_[2]>0) target_volume_=ResolveDetectorDirectedVolume(*target_,true,volume);
    BuildComponents();
}

void OnShellCascadeChannel::BuildComponents() {
    envelope_=weights_[1]>0 ? std::make_shared<RestFrameEnvelope2BodyChannel>(target_) : nullptr;
    // The volume was checked when it was resolved or loaded; do not repeat it.
    directed_=weights_[2]>0 ? std::shared_ptr<PhaseSpaceChannel const>(new DetectorDirected2BodyChannel(
        target_,0,DetectorDirected2BodyChannel::CheckedVolume{target_volume_})) : nullptr;
}

void OnShellCascadeChannel::Validate(bool archived) const {
    if (!target_ || !(pair_mass_>0) || !std::isfinite(pair_mass_)
        || !std::isfinite(kappa_) || kappa_<-1 || !std::isfinite(rho_) || rho_<0 || rho_>1)
        throw std::invalid_argument("Invalid OnShellCascade target, mass, kappa or daughter probability");
    double sum=0;
    for (auto w:weights_) {
        if (!std::isfinite(w) || w<0) throw std::invalid_argument("Invalid cascade mixture weight");
        sum+=w;
    }
    if (!(weights_[0]>0) || !std::isfinite(sum))
        throw std::invalid_argument("Cascade orientation requires a positive isotropic support weight");
    if (archived && (std::abs(sum-1)>1e-12 || (weights_[2]>0 && !(target_volume_>0))))
        throw std::invalid_argument("Invalid archived cascade normalization or volume");
    if (archived) ValidateArchivedDetectorDirectedVolume(target_.get(),weights_[2]>0,target_volume_);
}

double OnShellCascadeChannel::InternalDensity(Record const & r) const {
    if (!Allowed(r,pair_mass_)) return 0;
    auto parent=detail::OnShellParent(r);
    auto a=detail::ReadSecondary(r,1),b=detail::ReadSecondary(r,2);
    P4 pair{a.e+b.e,a.p+b.p};
    long double ep=pair.e, px=pair.p.GetX(), py=pair.p.GetY(), pz=pair.p.GetZ();
    long double pair2=ep*ep-px*px-py*py-pz*pz;
    double EA=TwoBodyRestEnergy(r.primary_mass,pair_mass_,r.secondary_masses[0]);
    // Rounding of the invariant, plus the boost error eps*gamma*(E*+p*) of each
    // lab daughter: E*_1+E*_2=EA and |p*_1|+|p*_2|<=EA. The second term
    // dominates for pairs emitted backward from a fast parent.
    double eps=std::numeric_limits<double>::epsilon(), gamma=parent.e/r.primary_mass;
    double tolerance=64*eps*(ep*ep+px*px+py*py+pz*pz)
        +256*eps*gamma*EA*(ep+std::sqrt(px*px+py*py+pz*pz))+1e-10*pair_mass_*pair_mass_;
    if (std::abs(pair2-pair_mass_*pair_mass_)>tolerance) return 0;
    double energy;
    if (!RestDaughterEnergy(r,parent,1,pair_mass_,energy)) return 0;
    double PA=TwoBodyRestMomentum(r.primary_mass,pair_mass_,r.secondary_masses[0]);
    double pd=TwoBodyRestMomentum(pair_mass_,r.secondary_masses[1],r.secondary_masses[2]);
    double c=std::clamp((energy-EA/2)/(PA*pd/pair_mass_),-1.0,1.0);
    // f(c)/(4*pi*2*pi), in dOmega_pair dOmega_sub, with no ds_pair.
    return (1+kappa_*c*c)/(16*M_PI*M_PI*(1+kappa_/3));
}

double OnShellCascadeChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,Record const & r) const
{
    if (!detail::OnShellParentValid(r)) return 0;
    double p=InternalDensity(r);
    if (p==0) return 0;
    auto parent=detail::OnShellParent(r);
    double sum=0;
    for (int d=1;d<=2;++d) {
        double rho=d==1 ? rho_ : 1-rho_;
        if (rho==0) continue;
        double energy;
        if (!RestDaughterEnergy(r,parent,d,pair_mass_,energy)) return 0;
        sum+=rho*OrientationDensity(Effective(r,energy,d),weights_,envelope_.get(),directed_.get());
    }
    return p*4*M_PI*sum;
}

void OnShellCascadeChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector,Record & r) const
{
    detail::RequireOnShellParent(r);
    if (!Allowed(r,pair_mass_)) throw siren::utilities::InjectionFailure(
        siren::utilities::FailureReason::KinematicallyForbidden,
        "OnShellCascade requires open phase space and equal nonnegative pair masses");
    int d=random->Uniform(0,1)<rho_ ? 1 : 2;
    double c;
    do { c=random->Uniform(-1,1); }
    while(random->Uniform(0,1)*(1+std::max(0.0,kappa_))>1+kappa_*c*c);
    double EA=TwoBodyRestEnergy(r.primary_mass,pair_mass_,r.secondary_masses[0]);
    double PA=TwoBodyRestMomentum(r.primary_mass,pair_mass_,r.secondary_masses[0]);
    double md=r.secondary_masses[d], pd=TwoBodyRestMomentum(pair_mass_,md,md);
    double energy=EA/2+PA*pd*c/pair_mass_;
    auto effective=Effective(r,energy,d);
    double component=random->Uniform(0,1);
    if (component<weights_[0]) Isotropic2BodyChannel().Sample(random,detector,effective);
    else if(component<weights_[0]+weights_[1]) envelope_->Sample(random,detector,effective);
    else directed_->Sample(random,detector,effective);
    auto parent=detail::OnShellParent(r);
    auto daughter=detail::BoostOnShell(parent,r.primary_mass,detail::ReadSecondary(effective,0),true);
    V axis=daughter.p/daughter.p.magnitude();
    double p=std::sqrt(std::max(0.0,(energy-md)*(energy+md)));
    daughter={energy,axis*p};
    double ca=std::clamp((EA*energy-pair_mass_*pair_mass_/2)/(PA*p),-1.0,1.0);
    V ref=std::abs(axis.GetZ())<0.9 ? V(0,0,1) : V(1,0,0);
    V first=ref-axis*siren::math::scalar_product(ref,axis);
    first=first/first.magnitude();
    V second=siren::math::cross_product(axis,first);
    double phi=random->Uniform(-M_PI,M_PI), sa=std::sqrt(std::max(0.0,1-ca*ca));
    P4 pair{EA,(axis*ca+(first*std::cos(phi)+second*std::sin(phi))*sa)*PA};
    P4 spectator{r.primary_mass-EA,pair.p*(-1)};
    P4 other{EA-energy,pair.p-daughter.p};
    detail::WriteSecondary(r,0,detail::BoostOnShell(parent,r.primary_mass,spectator));
    detail::WriteSecondary(r,d,detail::BoostOnShell(parent,r.primary_mass,daughter));
    detail::WriteSecondary(r,3-d,detail::BoostOnShell(parent,r.primary_mass,other));
    auto lab=detail::BoostOnShell(parent,r.primary_mass,pair);
    r.interaction_parameters["cascade_pair_energy"]=lab.e;
    r.interaction_parameters["cascade_pair_px"]=lab.p.GetX();
    r.interaction_parameters["cascade_pair_py"]=lab.p.GetY();
    r.interaction_parameters["cascade_pair_pz"]=lab.p.GetZ();
    r.interaction_parameters["cascade_pair_mass"]=pair_mass_;
    r.interaction_parameters["cascade_internal_kappa"]=kappa_;
    r.interaction_parameters["cascade_target_daughter"]=d;
}
}}
