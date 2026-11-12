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
