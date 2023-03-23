
#include <cmath>
#include <tuple>
#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>

#include <gtest/gtest.h>

#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/distributions/primary/direction/PrimaryDirectionDistribution.h"
#include "LeptonInjector/distributions/primary/direction/Cone.h"

using namespace LI::math;
using namespace LI::distributions;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(Cone, Constructor) {
    Vector3D direction(0,0,1);
    double opening_angle = M_PI;
    Cone A(direction, opening_angle);
}

static const std::map<std::tuple<size_t, size_t>, double> delta_chi2_table = {
    {{1, 1}, 0.9999999999999993},
    {{1, 2}, 2.295748928898636},
    {{1, 3}, 3.5267403802617303},
    {{1, 4}, 4.719474460025883},
    {{1, 5}, 5.887595445915204},
    {{1, 6}, 7.038400923736641},
    {{1, 7}, 8.176236497856527},
    {{1, 8}, 9.30391276903717},
    {{1, 9}, 10.423363154355838},
    {{1, 10}, 11.535981713319316},
    {{1, 11}, 12.64281133339149},
    {{1, 12}, 13.744655587189282},
    {{1, 13}, 14.842148802786893},
    {{1, 14}, 15.935801892195538},
    {{1, 15}, 17.026033423371082},
    {{1, 16}, 18.11319133873574},
    {{1, 17}, 19.197568537049687},
    {{1, 18}, 20.279414307857422},
    {{1, 19}, 21.358942889620966},
    {{1, 20}, 22.436339987385537},
    {{1, 21}, 23.511767813648422},
    {{1, 22}, 24.585369041159232},
    {{1, 23}, 25.657269941145202},
    {{1, 24}, 26.72758290286756},
    {{1, 25}, 27.79640847713285},
    {{1, 26}, 28.863837049131934},
    {{1, 27}, 29.92995021950986},
    {{1, 28}, 30.99482195347944},
    {{1, 29}, 32.05851954383435},
    {{1, 30}, 33.121104423385034},
    {{1, 31}, 34.182632854602666},
    {{1, 32}, 35.24315651839821},
    {{1, 33}, 36.30272301948572},
    {{1, 34}, 37.361376322322016},
    {{1, 35}, 38.4191571289259},
    {{1, 36}, 39.47610320776851},
    {{1, 37}, 40.532249681257746},
    {{1, 38}, 41.587629278010624},
    {{1, 39}, 42.64227255504241},
    {{1, 40}, 43.696208094140864},
    {{1, 41}, 44.74946267599779},
    {{1, 42}, 45.80206143509929},
    {{1, 43}, 46.854027997910244},
    {{1, 44}, 47.90538460650226},
    {{1, 45}, 48.95615222945574},
    {{1, 46}, 50.00635066160044},
    {{1, 47}, 51.055998613936985},
    {{1, 48}, 52.10511379489487},
    {{1, 49}, 53.153712983925445},
    {{1, 50}, 54.20181209829501},
    {{2, 1}, 4.000000000000002},
    {{2, 2}, 6.180074306244173},
    {{2, 3}, 8.024881760266252},
    {{2, 4}, 9.715627154871333},
    {{2, 5}, 11.313855908361862},
    {{2, 6}, 12.848834791793395},
    {{2, 7}, 14.337110231671799},
    {{2, 8}, 15.789092974617745},
    {{2, 9}, 17.21182898078949},
    {{2, 10}, 18.610346565823498},
    {{2, 11}, 19.988381717650192},
    {{2, 12}, 21.348799569984315},
    {{2, 13}, 22.693854280452445},
    {{2, 14}, 24.025357063756637},
    {{2, 15}, 25.344789151124267},
    {{2, 16}, 26.653380234523553},
    {{2, 17}, 27.952164463248984},
    {{2, 18}, 29.24202137323326},
    {{2, 19}, 30.52370642400417},
    {{2, 20}, 31.797874195789504},
    {{2, 21}, 33.0650962935076},
    {{2, 22}, 34.325875362488986},
    {{2, 23}, 35.58065620044425},
    {{2, 24}, 36.82983466857519},
    {{2, 25}, 38.07376491213091},
    {{2, 26}, 39.31276526654019},
    {{2, 27}, 40.5471231301975},
    {{2, 28}, 41.77709901660786},
    {{2, 29}, 43.00292994871745},
    {{2, 30}, 44.22483232140679},
    {{2, 31}, 45.44300433056983},
    {{2, 32}, 46.65762804637739},
    {{2, 33}, 47.86887119242127},
    {{2, 34}, 49.076888680179},
    {{2, 35}, 50.28182393870856},
    {{2, 36}, 51.48381007201128},
    {{2, 37}, 52.682970870596726},
    {{2, 38}, 53.87942169908884},
    {{2, 39}, 55.07327027794646},
    {{2, 40}, 56.264617374339274},
    {{2, 41}, 57.45355741475761},
    {{2, 42}, 58.64017902992812},
    {{2, 43}, 59.82456554095781},
    {{2, 44}, 61.00679539427409},
    {{2, 45}, 62.18694255180141},
    {{2, 46}, 63.365076841879606},
    {{2, 47}, 64.54126427564665},
    {{2, 48}, 65.71556733295104},
    {{2, 49}, 66.88804522130543},
    {{2, 50}, 68.05875411092512},
    {{3, 1}, 8.999999999999986},
    {{3, 2}, 11.829158081900795},
    {{3, 3}, 14.156413609126675},
    {{3, 4}, 16.251340813956187},
    {{3, 5}, 18.205314008384093},
    {{3, 6}, 20.06208616571403},
    {{3, 7}, 21.846581673015194},
    {{3, 8}, 23.574591022671044},
    {{3, 9}, 25.25686586179292},
    {{3, 10}, 26.901119405801232},
    {{3, 11}, 28.513108748422717},
    {{3, 12}, 30.097266729568556},
    {{3, 13}, 31.657093249667593},
    {{3, 14}, 33.19540931314273},
    {{3, 15}, 34.714528445229774},
    {{3, 16}, 36.21637614028137},
    {{3, 17}, 37.70257539809229},
    {{3, 18}, 39.17450942666591},
    {{3, 19}, 40.6333685500316},
    {{3, 20}, 42.08018593012067},
    {{3, 21}, 43.515865201408936},
    {{3, 22}, 44.94120215061899},
    {{3, 23}, 46.356901939365294},
    {{3, 24}, 47.763592941568135},
    {{3, 25}, 49.16183797543894},
    {{3, 26}, 50.5521435059642},
    {{3, 27}, 51.934967249091855},
    {{3, 28}, 53.310724504519875},
    {{3, 29}, 54.679793467760035},
    {{3, 30}, 56.042519715732546},
    {{3, 31}, 57.399220017895495},
    {{3, 32}, 58.75018559292716},
    {{3, 33}, 60.095684906516844},
    {{3, 34}, 61.43596608694126},
    {{3, 35}, 62.7712590204},
    {{3, 36}, 64.10177717654342},
    {{3, 37}, 65.42771920549643},
    {{3, 38}, 66.7492703404089},
    {{3, 39}, 68.06660363372812},
    {{3, 40}, 69.37988105068075},
    {{3, 41}, 70.6892544396268},
    {{3, 42}, 71.99486639582527},
    {{3, 43}, 73.29685103258723},
    {{3, 44}, 74.59533467167566},
    {{3, 45}, 75.89043646305683},
    {{3, 46}, 77.18226894264781},
    {{3, 47}, 78.47093853547749},
    {{3, 48}, 79.7565460106578},
    {{3, 49}, 81.03918689368778},
    {{3, 50}, 82.31895184088492},
    {{4, 1}, 15.999999999999998},
    {{4, 2}, 19.333908611934685},
    {{4, 3}, 22.061320636960435},
    {{4, 4}, 24.502059115721934},
    {{4, 5}, 26.766257283708462},
    {{4, 6}, 28.907360065488362},
    {{4, 7}, 30.956146442791987},
    {{4, 8}, 32.932291600972206},
    {{4, 9}, 34.84929300935437},
    {{4, 10}, 36.71689548914631},
    {{4, 11}, 38.54241238232088},
    {{4, 12}, 40.331501650282},
    {{4, 13}, 42.08864918221147},
    {{4, 14}, 43.81748413327278},
    {{4, 15}, 45.520992665148356},
    {{4, 16}, 47.20166750104189},
    {{4, 17}, 48.86161542476853},
    {{4, 18}, 50.50263635578319},
    {{4, 19}, 52.1262826933321},
    {{4, 20}, 53.733904641936924},
    {{4, 21}, 55.32668537106272},
    {{4, 22}, 56.90566866824598},
    {{4, 23}, 58.4717809590671},
    {{4, 24}, 60.0258490380714},
    {{4, 25}, 61.56861449099265},
    {{4, 26}, 63.100745534044826},
    {{4, 27}, 64.62284681488875},
    {{4, 28}, 66.13546758902157},
    {{4, 29}, 67.63910858949752},
    {{4, 30}, 69.13422783680946},
    {{4, 31}, 70.62124558242385},
    {{4, 32}, 72.10054853901096},
    {{4, 33}, 73.57249351941898},
    {{4, 34}, 75.03741058248194},
    {{4, 35}, 76.49560576506124},
    {{4, 36}, 77.94736346502879},
    {{4, 37}, 79.39294852825779},
    {{4, 38}, 80.83260808340223},
    {{4, 39}, 82.26657316078415},
    {{4, 40}, 83.6950601256804},
    {{4, 41}, 85.11827195139763},
    {{4, 42}, 86.5363993535151},
    {{4, 43}, 87.94962180338156},
    {{4, 44}, 89.3581084362285},
    {{4, 45}, 90.76201886700397},
    {{4, 46}, 92.16150392514471},
    {{4, 47}, 93.55670631792795},
    {{4, 48}, 94.94776123071607},
    {{4, 49}, 96.3347968712858},
    {{4, 50}, 97.71793496448579},
    {{5, 1}, 25.00000000007513},
    {{5, 2}, 28.743702426935496},
    {{5, 3}, 31.812108348682955},
    {{5, 4}, 34.555046564970844},
    {{5, 5}, 37.0947940997286},
    {{5, 6}, 39.491406276472254},
    {{5, 7}, 41.77980967672749},
    {{5, 8}, 43.98251343814422},
    {{5, 9}, 46.115065924506006},
    {{5, 10}, 48.18875924131356},
    {{5, 11}, 50.212111633914695},
    {{5, 12}, 52.19174308533455},
    {{5, 13}, 54.13292330108691},
    {{5, 14}, 56.039930864816185},
    {{5, 15}, 57.916297685886235},
    {{5, 16}, 59.764980689340234},
    {{5, 17}, 61.58848565684374},
    {{5, 18}, 63.38895861696615},
    {{5, 19}, 65.16825463726117},
    {{5, 20}, 66.92799051237009},
    {{5, 21}, 68.66958574135634},
    {{5, 22}, 70.39429483491575},
    {{5, 23}, 72.10323310011631},
    {{5, 24}, 73.79739744736494},
    {{5, 25}, 75.47768334890382},
    {{5, 26}, 77.14489878672772},
    {{5, 27}, 78.79977581999604},
    {{5, 28}, 80.44298025157822},
    {{5, 29}, 82.07511976297646},
    {{5, 30}, 83.69675080483483},
    {{5, 31}, 85.30838446857548},
    {{5, 32}, 86.91049151784733},
    {{5, 33}, 88.50350672251678},
    {{5, 34}, 90.08783261008381},
    {{5, 35}, 91.6638427276561},
    {{5, 36}, 93.23188449048513},
    {{5, 37}, 94.79228167948162},
    {{5, 38}, 96.34533663927007},
    {{5, 39}, 97.89133221961062},
    {{5, 40}, 99.43053349595122},
    {{5, 41}, 100.96318929911892},
    {{5, 42}, 102.48953357944922},
    {{5, 43}, 104.00978662677949},
    {{5, 44}, 105.52415616452262},
    {{5, 45}, 107.03283833337696},
    {{5, 46}, 108.53601857800138},
    {{5, 47}, 110.03387244812349},
    {{5, 48}, 111.52656632397492},
    {{5, 49}, 113.01425807462581},
    {{5, 50}, 114.49709765666123}
};

class DistributionTest {
    size_t num_bins;
    size_t num_entries;
    std::vector<double> bin_edges;
    std::vector<double> bin_contents;
public:
    DistributionTest(double min_edge, double max_edge, size_t num_bins)
        : num_bins(num_bins), num_entries(0),
          bin_edges(num_bins + 1), bin_contents(num_bins, 0) {
        assert(num_bins > 0);
        double bin_range = max_edge - min_edge;
        double bin_delta = bin_range / num_bins;
        for(size_t j=0; j<num_bins + 1; ++j) {
            bin_edges[j] = bin_delta * j + min_edge;
        }
    }

    void AddValue(double value) {
        num_entries += 1;
        int bin_idx = std::distance(bin_edges.begin(), std::lower_bound(bin_edges.begin(), bin_edges.end(), value)) - 1;
        if(bin_idx < 0)
            return;
        if(bin_idx >= num_bins)
            return;
        bin_contents[bin_idx] += 1;
    }

    bool TestContents(size_t sigma, std::vector<double> const & expectation) {
        double chi2 = 0;
        for(size_t j=0; j<num_bins; ++j) {
            double contents = bin_contents[j];
            double expected_error = sqrt(expectation[j]);
            double error = std::abs(contents - expectation[j]);
            double term = (error / expected_error);
            chi2 += term*term;
        }
        double max_delta_chi2 = delta_chi2_table.at({sigma, num_bins});
        return chi2 <= max_delta_chi2;
    }

    bool TestFractionalContents(size_t sigma, std::vector<double> const & fractional_expectation) {
        double chi2 = 0;
        for(size_t j=0; j<num_bins; ++j) {
            double contents = bin_contents[j] / num_entries;
            double expected_error = sqrt(fractional_expectation[j] / num_entries);
            double error = std::abs(contents - fractional_expectation[j]);
            double term = (error / expected_error);
            chi2 += term*term;
        }
        double max_delta_chi2 = delta_chi2_table.at({sigma, num_bins});
        return chi2 <= max_delta_chi2;
    }
};

/*
TEST(Quaternion, Sample) {
    for(size_t i=0; i<100; ++i) {
        Vector3D start(0,0,1);
        Vector3D dir(RandomDouble() - 0.5, RandomDouble() - 0.5, RandomDouble() - 0.5);
        dir.normalize();
        LI::math::Vector3D r = cross_product(start, dir);
        // r.normalize();
        Quaternion rotation = LI::math::Quaternion(r);
        rotation.SetW(1.0 + dir.GetZ());
        rotation.normalize();
        //Quaternion rotation = LI::math::Quaternion(0,0,0,1);

        double angle = M_PI * RandomDouble();
        double c = RandomDouble() * (1.0 - cos(angle)) + cos(angle);
        double theta = acos(c);
        std::cout << angle << " " << c << " " << theta << std::endl;
        EXPECT_TRUE(theta < angle);
        double phi = RandomDouble() * 2.0 * M_PI;
        LI::math::Quaternion q;
        q.SetEulerAnglesZXZr(phi, theta, 0.0);
        Vector3D result = rotation.rotate(q.rotate(LI::math::Vector3D(0,0,1), false), false);

        //std::cout << "Start:" << std::endl;
        //std::cout << start << std::endl;
        //std::cout << std::endl;
        //std::cout << "Dir:" << std::endl;
        //std::cout << dir << std::endl;
        //std::cout << std::endl;
        //std::cout << "Rotation:" << std::endl;
        //std::cout << rotation << std::endl;
        //std::cout << std::endl;
        //std::cout << "phi: " << phi << std::endl;
        //std::cout << std::endl;
        //std::cout << "q:" << std::endl;
        //std::cout << q << std::endl;
        //std::cout << std::endl;
        //std::cout << "result:" << std::endl;
        //std::cout << result << std::endl;
        //std::cout << std::endl;
        //std::cout << "dot:" << std::endl;
        //std::cout << scalar_product(dir, result) << std::endl;
        //std::cout << std::endl;
        std::cout << "angle:" << std::endl;
        std::cout << theta << std::endl;
        std::cout << acos(scalar_product(dir, result)) << std::endl;
        std::cout << std::endl;
    }
}
*/

TEST(Cone, SampleBounds) {
    size_t N = 100;
    size_t M = 1000;
    std::shared_ptr<LI::utilities::LI_random> rand = std::make_shared<LI::utilities::LI_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        direction.normalize();

        double opening_angle = M_PI * RandomDouble();
        Cone A(direction, opening_angle);
        for(size_t j=0; j<M; ++j) {
            LI::dataclasses::InteractionRecord record;
            record.primary_momentum[0] = 1;
            record.primary_mass = 0;
            A.Sample(rand, nullptr, nullptr, record);
            Vector3D vec(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
            vec.normalize();
            double angle = acos(scalar_product(vec, direction));
            EXPECT_TRUE(angle <= opening_angle);
            EXPECT_TRUE(angle >= 0);
        }
    }
}

TEST(Cone, SampleDistribution) {
    size_t N = 100;
    size_t M = 10000;
    std::shared_ptr<LI::utilities::LI_random> rand = std::make_shared<LI::utilities::LI_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        direction.normalize();

        double opening_angle = M_PI * RandomDouble();
        Cone A(direction, opening_angle);
        size_t n_bins = (RandomDouble() * M / 500) + 1;
        double bin_max = 1.0;
        double bin_min = cos(opening_angle);
        DistributionTest test(bin_min, bin_max, n_bins);
        for(size_t j=0; j<M; ++j) {
            LI::dataclasses::InteractionRecord record;
            record.primary_momentum[0] = 1;
            record.primary_mass = 0;
            A.Sample(rand, nullptr, nullptr, record);
            Vector3D vec(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
            vec.normalize();
            double c = scalar_product(vec, direction);
            test.AddValue(c);
            double angle = acos(scalar_product(vec, direction));
            EXPECT_TRUE(angle <= opening_angle);
            EXPECT_TRUE(angle >= 0);
        }
        double expected_contents = double(M) / double(n_bins);
        std::vector<double> expect(n_bins, expected_contents);
        EXPECT_TRUE(test.TestContents(3, expect));
    }
}

TEST(Cone, GenerationProbability) {
    size_t N = 1000;
    size_t M = 10000;
    std::shared_ptr<LI::utilities::LI_random> rand = std::make_shared<LI::utilities::LI_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        while(true) {
            direction.normalize();
            double magnitude = direction.magnitude();
            if(std::isnan(direction.magnitude()) or magnitude == 0)
                direction = Vector3D(RandomDouble(), RandomDouble(), RandomDouble());
            else
                break;
        }

        double opening_angle = M_PI * RandomDouble();
        Cone A(direction, opening_angle);
        double expected_density = 1.0 / ((2.0 * M_PI) * (1.0 - cos(opening_angle)));
        for(size_t j=0; j<M; ++j) {
            double input_angle = RandomDouble() * opening_angle;
            Vector3D vec(RandomDouble(), RandomDouble(), RandomDouble());
            while(true) {
                vec.normalize();
                double magnitude = vec.magnitude();
                if(std::isnan(vec.magnitude()) or magnitude == 0)
                    vec = Vector3D(RandomDouble(), RandomDouble(), RandomDouble());
                else
                    break;
            }
            vec = cross_product(direction, vec).normalized();
            Quaternion q(vec * sin(input_angle));
            q.SetW(1.0 + cos(input_angle));
            q.normalize();
            vec = q.rotate(direction, false);

            EXPECT_TRUE(acos(scalar_product(vec, direction)) <= opening_angle);

            LI::dataclasses::InteractionRecord record;
            record.primary_momentum[1] = vec.GetX();
            record.primary_momentum[2] = vec.GetY();
            record.primary_momentum[3] = vec.GetZ();

            double c = scalar_product(direction.normalized(), vec.normalized());
            double output_angle;
            if(c > 1)
                output_angle = 0;
            else
                output_angle = acos(c);

            double resolution;
            if(input_angle == 0)
                resolution = 0;
            else
                resolution = 1.0 / sqrt(input_angle) * 1e-8;

            EXPECT_NEAR(input_angle, output_angle, resolution);

            double density = A.GenerationProbability(nullptr, nullptr, record);

            if(output_angle < opening_angle)
                EXPECT_NEAR(density, expected_density, expected_density * 1e-8);
            else
                EXPECT_TRUE(density == 0);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

