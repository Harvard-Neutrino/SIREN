#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"

namespace {
using namespace siren::injection;
class IntervalChannel : public PhaseSpaceChannel {
public:
    bool cubic;
    explicit IntervalChannel(bool c): cubic(c) {}
    void Sample(std::shared_ptr<siren::utilities::SIREN_random> random,
                std::shared_ptr<siren::detector::DetectorModel const>,
                siren::dataclasses::InteractionRecord & r) const override {
        double u=random->Uniform(0,1);
        r.interaction_parameters["x"] = cubic ? std::cbrt(u) : u;
    }
    double Density(std::shared_ptr<siren::detector::DetectorModel const>,
                   siren::dataclasses::InteractionRecord const & r) const override {
        double x=r.interaction_parameters.at("x");
        return cubic ? 3*x*x : 1.;
    }
    std::string Name() const override { return cubic ? "cubic" : "uniform"; }
    PhaseSpaceTopology Topology() const override { return PhaseSpaceTopology::Unspecified; }
    PhaseSpaceMeasure Measure() const override { return PhaseSpaceMeasure::Unspecified(); }
};

MultiChannelPhaseSpace MakeMixture() {
    return MultiChannelPhaseSpace({std::make_shared<IntervalChannel>(true),
                                  std::make_shared<IntervalChannel>(false)}, {.5,.5});
}

// Independent midpoint integration of variance for f(x)=3x^2, integral=1.
double Variance(std::vector<double> const & a) {
    double second=0;
    for(int i=0;i<100000;++i) {
        double x=(i+.5)/100000., f=3*x*x, g=a[0]*f+a[1];
        second+=f*f/g/100000.;
    }
    return second-1;
}

TEST(PhaseSpaceTuning, IndependentIntegralAndVariance) {
    auto mc=MakeMixture();
    auto rng=std::make_shared<siren::utilities::SIREN_random>(191);
    siren::dataclasses::InteractionRecord r;
    for(int round=0;round<20;++round) {
        for(int i=0;i<3000;++i) {
            mc.Sample(rng,nullptr,r);
            double x=r.interaction_parameters["x"];
            double w=3*x*x/mc.Density(nullptr,r);
            mc.Accumulate(nullptr,r,w);
        }
        mc.UpdateWeights("alpha_sqrt_W",1.,.005);
    }
    ASSERT_GE(mc.Weights()[1],.005);
    EXPECT_LT(Variance(mc.Weights()),.12*Variance({.5,.5}));
    // Production uses another stream and the frozen final proposal.
    rng=std::make_shared<siren::utilities::SIREN_random>(911);
    double sum=0, square=0;
    for(int i=0;i<20000;++i) {
        mc.Sample(rng,nullptr,r);
        double x=r.interaction_parameters["x"], w=3*x*x/mc.Density(nullptr,r);
        sum+=w; square+=w*w;
    }
    double mean=sum/20000, se=std::sqrt((square/20000-mean*mean)/19999);
    EXPECT_NEAR(mean,1.,6*se);
}

TEST(PhaseSpaceTuning, CommonMeasureAndExtremeWeights) {
    auto mc=MakeMixture();
    siren::dataclasses::InteractionRecord r;
    r.interaction_parameters["x"]=.5;
    mc.Accumulate(nullptr,r,1e300);
    mc.UpdateWeights("sqrt_W",1.,.01);
    double ratio=std::sqrt(.75);
    EXPECT_NEAR(mc.Weights()[0],.01+.98*ratio/(1+ratio),1e-13);
    EXPECT_THROW(mc.UpdateWeights("sqrt_W",1.,0.),std::invalid_argument);
    EXPECT_THROW(mc.Accumulate(nullptr,r,std::numeric_limits<double>::infinity()),std::runtime_error);
    EXPECT_EQ(mc.TuningCount(),0u);
}
}

namespace {
class AngularChannel : public PhaseSpaceChannel {
public:
    double slope;
    bool solid_angle;
    AngularChannel(double slope, bool solid): slope(slope), solid_angle(solid) {}
    void Sample(std::shared_ptr<siren::utilities::SIREN_random> random,
                std::shared_ptr<siren::detector::DetectorModel const>,
                siren::dataclasses::InteractionRecord & r) const override {
        double u=random->Uniform(0,1);
        r.interaction_parameters["cosine"] = slope == 0 ? 2*u-1
            : (-1+std::sqrt((1-slope)*(1-slope)+4*slope*u))/slope;
    }
    double Density(std::shared_ptr<siren::detector::DetectorModel const>,
                   siren::dataclasses::InteractionRecord const & r) const override {
        return .5*(1+slope*r.interaction_parameters.at("cosine"))
            / (solid_angle ? 2*std::acos(-1.) : 1.);
    }
    PhaseSpaceTopology Topology() const override { return PhaseSpaceTopology::Decay2Body; }
    PhaseSpaceMeasure Measure() const override {
        return solid_angle ? PhaseSpaceMeasure::SolidAngleRest() : PhaseSpaceMeasure::CosThetaRest();
    }
    std::string Name() const override { return "angular"; }
};
}

TEST(PhaseSpaceTuning, NestedDistinctMeasuresMatchIndependentQuadrature) {
    auto inner=std::make_shared<MultiChannelPhaseSpace>(
        std::vector<std::shared_ptr<PhaseSpaceChannel>>{
            std::make_shared<AngularChannel>(.9,false),std::make_shared<AngularChannel>(-.7,true)},
        std::vector<double>{.25,.75});
    MultiChannelPhaseSpace outer({std::make_shared<NestedMixtureChannel>(inner),
                                  std::make_shared<AngularChannel>(.4,false)}, {.4,.6});
    // Analytic cos(theta) densities: the outer proposal has slope .12.
    constexpr int n=100000;
    double wi[2]={0,0},wo[2]={0,0};
    siren::dataclasses::InteractionRecord r;
    for(int i=0;i<n;++i) {
        double x=-1+(i+.5)*2/n, a=.5*(1+.9*x), b=.5*(1-.7*x);
        double nested=.25*a+.75*b, c=.5*(1+.4*x), g=.4*nested+.6*c, f=.5*(1+.8*x);
        wi[0]+=2./n*f*f*a/(g*g);wi[1]+=2./n*f*f*b/(g*g);
        wo[0]+=2./n*f*f*nested/(g*g);wo[1]+=2./n*f*f*c/(g*g);
        double u=(i+.5)/n;
        x=(-1+std::sqrt(.88*.88+.48*u))/.12;
        g=.5*(1+.12*x);f=.5*(1+.8*x);
        r.interaction_parameters["cosine"]=x;
        outer.Accumulate(nullptr,r,f/g);
    }
    outer.UpdateWeights("alpha_sqrt_W",1.,.01);
    auto expected=[](double alpha,double a,double b) {
        double p=alpha*std::sqrt(a),q=(1-alpha)*std::sqrt(b);
        return .01+.98*p/(p+q);
    };
    EXPECT_NEAR(inner->Weights()[0],expected(.25,wi[0],wi[1]),1e-8);
    EXPECT_NEAR(outer.Weights()[0],expected(.4,wo[0],wo[1]),1e-8);
    EXPECT_GT(std::abs(inner->Weights()[0]-.25),.02);
}
