#include "SIREN/interactions/GammaDISFromSpline.h"

#include <math.h>   // before rk: geom3::rotationBetween calls an unqualified abs on doubles
#include <cmath>
#include <cctype>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include <Pythia8/Pythia.h>

#include "cnpy.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/Random.h"

namespace siren {
namespace interactions {

// lets Pythia draw from the SIREN generator
class GammaDISRndm : public Pythia8::RndmEngine {
public:
    std::shared_ptr<siren::utilities::SIREN_random> rng;
    double flat() override { return rng->Uniform(0.0, 1.0); }
};

namespace {
    constexpr double EARLY_SIGMA = 2.5;  // reject draws this many sigma too early or too narrow
    constexpr int    MAX_TRIES   = 100;
    constexpr double COV_RIDGE  = 1e-3;  // floor on the covariance eigenvalues
    constexpr int    ELECTRON_PDG = 11;
    // Start scale of the parton shower as a fraction of Q. 0.3 reproduces the hadron
    // multiplicity and leading-hadron energy fraction of full Pythia DIS events.
    constexpr double SHOWER_SCALE_FRACTION = 0.3;
    constexpr double MIN_SHOWER_SCALE  = 0.5;   // GeV
    constexpr double MIN_STRING_EXCESS = 0.3;   // GeV above the quark + diquark mass
    constexpr int    HADRONIZE_TRIES   = 20;

    // Same check as in DISFromSpline.cxx, where it is file-local.
    bool kinematicallyAllowed(double x, double y, double E, double M, double m) {
        if(x > 1) return false;
        if(x < ((m * m) / (2 * M * (E - m)))) return false;
        double d = 2 * (1 + (M * x) / (2 * E));
        double ad = 1 - m * m * ((1 / (2 * M * E * x)) + (1 / (2 * E * E)));
        double term = 1 - ((m * m) / (2 * M * E * x));
        double bd = std::sqrt(term * term - ((m * m) / (E * E)));
        return (ad - bd) <= d * y and d * y <= (ad + bd);
    }

    int charged_lepton_pdg(siren::dataclasses::ParticleType nu) {
        using PT = siren::dataclasses::ParticleType;
        if(nu == PT::NuE)      return 11;
        if(nu == PT::NuEBar)   return -11;
        if(nu == PT::NuMu)     return 13;
        if(nu == PT::NuMuBar)  return -13;
        if(nu == PT::NuTau)    return 15;
        if(nu == PT::NuTauBar) return -15;
        return 0;
    }

    // CC: the charged lepton. NC: the neutrino itself.
    int outgoing_lepton_pdg(siren::dataclasses::ParticleType nu, int interaction_type) {
        if(interaction_type == 1) return charged_lepton_pdg(nu);
        if(interaction_type == 2) return static_cast<int>(nu);
        throw std::runtime_error("GammaDISFromSpline: only CC (1) and NC (2) splines are supported, got interaction type "
                                 + std::to_string(interaction_type));
    }

    // Box-Muller
    double draw_normal(std::shared_ptr<siren::utilities::SIREN_random> r) {
        double u1 = std::max(r->Uniform(0.0, 1.0), 1e-300);
        double u2 = r->Uniform(0.0, 1.0);
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    }

    int pdg_to_model_pid(int pid, std::shared_ptr<siren::utilities::SIREN_random> r) {
        // K0/K0bar -> KS or KL with equal probability
        if(std::abs(pid) == 311) return (r->Uniform(0.0,1.0) < 0.5) ? 310 : 130;
        return pid;
    }

    // The shower model is parametrized in kinetic energy and the hadrons carry
    // total energy, so the model is evaluated at E - m.
    double rest_mass_GeV(int pid) {
        switch(std::abs(pid)) {
            case 11:   return 0.51099895e-3;   // e
            case 22:   return 0.0;             // photon
            case 111:  return 0.1349768;       // pi0
            case 211:  return 0.13957039;      // pi+-
            case 130: case 310: case 311:
                       return 0.497611;        // K0 / KS / KL
            case 321:  return 0.493677;        // K+-
            case 2212: return 0.93827208816;   // p, pbar
            case 2112: return 0.93956542052;   // n, nbar
            default:   return 0.0;
        }
    }

    // Eigen-decomposition of a small symmetric matrix (cyclic Jacobi).
    // a is n x n row-major; the eigenvectors are the columns of V.
    void symmetric_eigen(std::vector<double> a, int n,
                         std::vector<double> & w, std::vector<double> & V) {
        V.assign(n * n, 0.0);
        for(int i = 0; i < n; ++i) V[i*n+i] = 1.0;
        for(int sweep = 0; sweep < 100; ++sweep) {
            double off = 0.0;
            for(int p = 0; p < n; ++p)
                for(int q = p+1; q < n; ++q) off += a[p*n+q] * a[p*n+q];
            if(off < 1e-30) break;
            for(int p = 0; p < n; ++p) {
                for(int q = p+1; q < n; ++q) {
                    double apq = a[p*n+q];
                    if(std::abs(apq) < 1e-300) continue;
                    double th = (a[q*n+q] - a[p*n+p]) / (2.0 * apq);
                    double t = (th >= 0 ? 1.0 : -1.0) / (std::abs(th) + std::sqrt(th*th + 1.0));
                    double c = 1.0 / std::sqrt(t*t + 1.0), s = t * c;
                    for(int k = 0; k < n; ++k) {
                        double akp = a[k*n+p], akq = a[k*n+q];
                        a[k*n+p] = c*akp - s*akq;  a[k*n+q] = s*akp + c*akq;
                    }
                    for(int k = 0; k < n; ++k) {
                        double apk = a[p*n+k], aqk = a[q*n+k];
                        a[p*n+k] = c*apk - s*aqk;  a[q*n+k] = s*apk + c*aqk;
                    }
                    for(int k = 0; k < n; ++k) {
                        double vkp = V[k*n+p], vkq = V[k*n+q];
                        V[k*n+p] = c*vkp - s*vkq;  V[k*n+q] = s*vkp + c*vkq;
                    }
                }
            }
        }
        w.resize(n);
        for(int i = 0; i < n; ++i) w[i] = a[i*n+i];
    }
} // namespace

GammaDISFromSpline::GammaDISFromSpline() : DISFromSpline() {}

GammaDISFromSpline::~GammaDISFromSpline() = default;

GammaDISFromSpline::GammaDISFromSpline(
    std::string differential_filename, std::string total_filename,
    std::vector<siren::dataclasses::ParticleType> primary_types,
    std::vector<siren::dataclasses::ParticleType> target_types,
    std::string shower_model_path, std::string pythia_data_path,
    std::string units)
    : DISFromSpline(differential_filename, total_filename,
                    primary_types, target_types, units),
      shower_model_path_(shower_model_path),
      pythia_data_path_(pythia_data_path) {
    // Spline files without an INTERACTION key are read as CC by the base class.
    // The CSMS splines have no such key, so the file name decides CC or NC.
    std::string base = differential_filename.substr(differential_filename.find_last_of("/\\") + 1);
    for(char & ch : base) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    if(GetInteractionType() == 1 &&
       (base.find("_NC_") != std::string::npos || base.find("_NC.") != std::string::npos))
        interaction_type_override_ = 2;
    LoadShowerModel(shower_model_path_);
    InitializeGammaSignatures();
}

int GammaDISFromSpline::InteractionType() const {
    return interaction_type_override_ ? interaction_type_override_ : GetInteractionType();
}

void GammaDISFromSpline::SetInteractionType(int interaction_type) {
    if(interaction_type != 1 && interaction_type != 2)
        throw std::runtime_error("GammaDISFromSpline::SetInteractionType: use 1 (CC) or 2 (NC), got "
                                 + std::to_string(interaction_type));
    interaction_type_override_ = interaction_type;
    InitializeGammaSignatures();
}

GammaDISFromSpline::GammaDISFromSpline(
    std::string differential_filename, std::string total_filename,
    int interaction, double target_mass, double minimum_Q2,
    std::vector<siren::dataclasses::ParticleType> primary_types,
    std::vector<siren::dataclasses::ParticleType> target_types,
    std::string shower_model_path, std::string pythia_data_path,
    std::string units)
    : DISFromSpline(differential_filename, total_filename, interaction,
                    target_mass, minimum_Q2, primary_types, target_types, units),
      shower_model_path_(shower_model_path),
      pythia_data_path_(pythia_data_path) {
    LoadShowerModel(shower_model_path_);
    InitializeGammaSignatures();
}

std::string GammaDISFromSpline::Name() const { return "GammaDISFromSpline"; }

bool GammaDISFromSpline::equal(CrossSection const & other) const {
    const GammaDISFromSpline* x = dynamic_cast<const GammaDISFromSpline*>(&other);
    if(!x) return false;
    return shower_model_path_ == x->shower_model_path_ &&
           pythia_data_path_ == x->pythia_data_path_ &&
           InteractionType() == x->InteractionType() &&
           DISFromSpline::equal(other);
}

void GammaDISFromSpline::LoadShowerModel(std::string const & npz_path) {
    cnpy::npz_t npz;
    try { npz = cnpy::npz_load(npz_path); }
    catch(std::exception const & e) {
        throw std::runtime_error("GammaDISFromSpline: failed to read '" + npz_path + "': " + e.what());
    }
    species_.clear();

    auto get_d = [&](std::string const & k) -> cnpy::NpyArray & {
        auto it = npz.find(k);
        if(it == npz.end()) throw std::runtime_error("shower model missing key: " + k);
        return it->second;
    };

    cnpy::NpyArray & pids_arr = get_d("pids");
    const int64_t* pids = pids_arr.data<int64_t>();
    Kmax_ = static_cast<int>(get_d("Kmax").data<int64_t>()[0]);
    size_t P = pids_arr.shape[0];

    for(size_t ip = 0; ip < P; ++ip) {
        int pid = static_cast<int>(pids[ip]);
        std::string p = std::to_string(pid);
        SpeciesModel s;
        s.pid = pid;

        auto & logE = get_d(p + "|logE");
        s.logE.assign(logE.data<double>(), logE.data<double>() + logE.shape[0]);
        s.logE_lo = s.logE.front();
        s.logE_hi = s.logE.back();

        auto & z = get_d(p + "|z");
        s.z.assign(z.data<double>(), z.data<double>() + z.shape[0]);

        auto & pm = get_d(p + "|pm");       // [L, Kmax]
        s.pm.assign(pm.data<double>(), pm.data<double>() + pm.num_vals);

        auto & lN = get_d(p + "|logN");
        s.logN.assign(lN.data<double>(), lN.data<double>() + lN.shape[0]);
        auto & lNs = get_d(p + "|logNsig");
        s.logNsig.assign(lNs.data<double>(), lNs.data<double>() + lNs.shape[0]);

        auto & ma = get_d(p + "|m_avail");
        const int64_t* mv = ma.data<int64_t>();
        for(size_t j = 0; j < ma.shape[0]; ++j) {
            int m = static_cast<int>(mv[j]);
            s.m_avail.push_back(m);
            auto & mean = get_d(p + "|m" + std::to_string(m) + "|mean"); // [L, 3m-1]
            auto & cov  = get_d(p + "|m" + std::to_string(m) + "|cov");  // [L, 3m-1, 3m-1]
            s.mean[m].assign(mean.data<double>(), mean.data<double>() + mean.num_vals);
            s.cov[m].assign(cov.data<double>(),  cov.data<double>()  + cov.num_vals);
        }
        species_[pid] = std::move(s);
    }
}

// linear interpolation of column comp of y (row stride ny_stride) at lq, clamped to the grid
double GammaDISFromSpline::interp1d_(std::vector<double> const & grid,
                                     double const * y, int ny_stride, int comp,
                                     double lq) const {
    size_t n = grid.size();
    if(n == 1) return y[comp];
    if(lq <= grid.front()) return y[0*ny_stride + comp];
    if(lq >= grid.back())  return y[(n-1)*ny_stride + comp];
    size_t hi = std::lower_bound(grid.begin(), grid.end(), lq) - grid.begin();
    size_t lo = hi - 1;
    double t = (lq - grid[lo]) / (grid[hi] - grid[lo]);
    return y[lo*ny_stride + comp] * (1 - t) + y[hi*ny_stride + comp] * t;
}

double GammaDISFromSpline::lq_clamped_(SpeciesModel const & s, double E) const {
    double lq = std::log10(E);
    return std::min(std::max(lq, s.logE_lo), s.logE_hi);
}

// z holds (log centroid c, log width s) per component: alpha = (c/s)^2, beta = c/s^2
void GammaDISFromSpline::from_z_(std::vector<double> const & z, int m,
                                 std::vector<double> & w,
                                 std::vector<double> & alpha,
                                 std::vector<double> & beta) {
    w.assign(m, 1.0); alpha.assign(m, 1.0); beta.assign(m, 1.0);
    if(m == 1) {
        double c = std::exp(z[0]), s = std::exp(z[1]);
        alpha[0] = (c/s)*(c/s);
        beta[0]  = c/(s*s);
        return;
    }
    // z = [alr(m-1), lc(m), ls(m)]
    double denom = 1.0; // exp(0) for the pinned last component
    std::vector<double> e(m);
    for(int i = 0; i < m - 1; ++i) { e[i] = std::exp(z[i]); denom += e[i]; }
    e[m-1] = 1.0;
    for(int i = 0; i < m; ++i) w[i] = e[i] / denom;
    for(int i = 0; i < m; ++i) {
        double c = std::exp(z[(m-1) + i]);
        double s = std::exp(z[(m-1) + m + i]);
        alpha[i] = (c/s)*(c/s);
        beta[i]  = c/(s*s);
    }
}

double GammaDISFromSpline::YieldMean(int pid, double E_GeV) const {
    auto it = species_.find(pid);
    if(it == species_.end()) return 0.0;
    SpeciesModel const & s = it->second;
    double lq = lq_clamped_(s, E_GeV);
    return std::exp(interp1d_(s.logE, s.logN.data(), 1, 0, lq));
}

std::vector<GammaDISFromSpline::Component> GammaDISFromSpline::SampleShower(
    int pid, double E_GeV,
    std::shared_ptr<siren::utilities::SIREN_random> random) const {

    std::vector<Component> out;
    auto it = species_.find(pid);
    if(it == species_.end()) return out;
    SpeciesModel const & s = it->second;
    double lq = lq_clamped_(s, E_GeV);

    // number of components m, drawn from p(m)
    std::vector<double> p; p.reserve(s.m_avail.size());
    double psum = 0.0;
    for(int m : s.m_avail) {
        double pm = interp1d_(s.logE, s.pm.data(), Kmax_, m - 1, lq);
        pm = std::max(pm, 0.0);
        p.push_back(pm); psum += pm;
    }
    if(psum <= 0.0) { for(auto & v : p) v = 1.0 / p.size(); }
    else            { for(auto & v : p) v /= psum; }
    double u = random->Uniform(0.0, 1.0);
    int m = s.m_avail.back(); double acc = 0.0;
    for(size_t i = 0; i < p.size(); ++i) { acc += p[i]; if(u <= acc) { m = s.m_avail[i]; break; } }

    int dim = 3 * m - 1;

    // mean and covariance at lq
    std::vector<double> mean(dim), cov(dim * dim);
    auto const & meanbuf = s.mean.at(m);   // [L, dim]
    auto const & covbuf  = s.cov.at(m);    // [L, dim*dim]
    for(int c = 0; c < dim; ++c)
        mean[c] = interp1d_(s.logE, meanbuf.data(), dim, c, lq);
    for(int c = 0; c < dim*dim; ++c)
        cov[c] = interp1d_(s.logE, covbuf.data(), dim*dim, c, lq);

    // symmetrize, then floor the eigenvalues at COV_RIDGE
    for(int a = 0; a < dim; ++a)
        for(int b = a+1; b < dim; ++b) {
            double v = 0.5 * (cov[a*dim+b] + cov[b*dim+a]);
            cov[a*dim+b] = cov[b*dim+a] = v;
        }
    std::vector<double> evals, evecs;
    symmetric_eigen(cov, dim, evals, evecs);
    for(auto & e : evals) e = std::max(e, COV_RIDGE);

    // z = mean + V sqrt(lambda) g with g ~ N(0,1). A draw is repeated until every
    // component has alpha > 1 (the profile vanishes at the vertex) and is not more
    // than EARLY_SIGMA below the mean in log centroid or log width. The late and
    // wide side is not restricted.
    std::vector<double> sd(dim);
    for(int a = 0; a < dim; ++a) {
        double var = 0.0;
        for(int k = 0; k < dim; ++k) var += evecs[a*dim+k] * evecs[a*dim+k] * evals[k];
        sd[a] = std::sqrt(std::max(var, 0.0));
    }
    int shape_lo = (m == 1) ? 0 : (m - 1);               // first log-centroid entry
    std::vector<double> g(dim), z(dim), w, alpha, beta;
    bool ok = false;
    for(int tries = 0; tries < MAX_TRIES && !ok; ++tries) {
        for(int a = 0; a < dim; ++a) g[a] = draw_normal(random);
        for(int a = 0; a < dim; ++a) {
            double v = mean[a];
            for(int k = 0; k < dim; ++k)
                v += evecs[a*dim+k] * std::sqrt(evals[k]) * g[k];
            z[a] = v;
        }
        ok = true;
        for(int a = shape_lo; a < dim && ok; ++a)        // all lc and ls entries
            if(z[a] < mean[a] - EARLY_SIGMA * sd[a]) ok = false;
        if(!ok) continue;
        from_z_(z, m, w, alpha, beta);
        for(int i = 0; i < m && ok; ++i)
            if(!(alpha[i] > 1.0)) ok = false;
    }
    if(!ok) { z = mean; from_z_(z, m, w, alpha, beta); }

    // lognormal yield, linear in E outside the trained range
    double N = std::exp(interp1d_(s.logE, s.logN.data(), 1, 0, lq));
    double sig = std::max(interp1d_(s.logE, s.logNsig.data(), 1, 0, lq), 0.0);
    if(sig > 0.0) N *= std::exp(sig * draw_normal(random));
    N *= E_GeV / std::pow(10.0, lq);

    out.reserve(m);
    for(int i = 0; i < m; ++i)
        out.push_back({alpha[i], beta[i], N * w[i]});
    return out;
}

void GammaDISFromSpline::InitializeGammaSignatures() {
    gamma_signatures_.clear();
    gamma_signatures_by_parents_.clear();
    // top_k hadrons + 1 pi0 remainder + 1 electron, each up to m=2 components
    kMaxSub_ = (top_k_ + 2) * 2;

    std::vector<siren::dataclasses::ParticleType> primaries = GetPossiblePrimaries();
    std::vector<siren::dataclasses::ParticleType> targets   = GetPossibleTargets();
    for(siren::dataclasses::ParticleType primary : primaries) {
        int lep = outgoing_lepton_pdg(primary, InteractionType());
        bool shower_lepton = (std::abs(lep) == 11);   // only e+- are showered
        dataclasses::InteractionSignature sig;
        sig.primary_type = primary;
        if(!shower_lepton)
            sig.secondary_types.push_back(static_cast<siren::dataclasses::ParticleType>(lep));  // slot 0
        for(int i = 0; i < kMaxSub_; ++i)
            sig.secondary_types.push_back(siren::dataclasses::ParticleType::GammaSubshower);
        for(siren::dataclasses::ParticleType target : targets) {
            sig.target_type = target;
            gamma_signatures_.push_back(sig);
            gamma_signatures_by_parents_[{primary, target}].push_back(sig);
        }
    }
}

std::vector<siren::dataclasses::InteractionSignature>
GammaDISFromSpline::GetPossibleSignatures() const {
    return gamma_signatures_;
}

std::vector<siren::dataclasses::InteractionSignature>
GammaDISFromSpline::GetPossibleSignaturesFromParents(
    siren::dataclasses::ParticleType primary_type,
    siren::dataclasses::ParticleType target_type) const {
    auto it = gamma_signatures_by_parents_.find({primary_type, target_type});
    if(it != gamma_signatures_by_parents_.end()) return it->second;
    return {};
}

void GammaDISFromSpline::InitializePythia() const {
    if(pythia_data_path_.empty())
        pythia_ = std::make_unique<Pythia8::Pythia>("../share/Pythia8/xmldoc", false);
    else
        pythia_ = std::make_unique<Pythia8::Pythia>(pythia_data_path_, false);

    // hadron level only: Hadronize fills the event record with the partons
    pythia_->readString("ProcessLevel:all = off");
    pythia_->readString("111:mayDecay = off");
    pythia_->readString("Print:quiet = on");
    pythia_->readString("Next:numberCount = 0");
    if(!pythia_->init())
        throw std::runtime_error("GammaDISFromSpline: Pythia initialization failed");

    // after init(), which uses Pythia's own generator
    rndm_ = std::make_shared<GammaDISRndm>();
    pythia_->rndm.rndmEnginePtr(rndm_);
}

std::vector<std::pair<int,std::array<double,4>>> GammaDISFromSpline::Hadronize(
    std::array<double,4> const & p4, double Q,
    siren::dataclasses::ParticleType primary,
    siren::dataclasses::ParticleType target,
    std::shared_ptr<siren::utilities::SIREN_random> random) const {
    using PT = siren::dataclasses::ParticleType;

    if(!pythia_) InitializePythia();
    rndm_->rng = random;

    double W = std::sqrt(std::max(0.0, p4[0]*p4[0] - p4[1]*p4[1] - p4[2]*p4[2] - p4[3]*p4[3]));

    // an unspecified nucleon is a proton 10 times out of 18, as in PythiaDISCrossSection
    bool proton;
    if(target == PT::PPlus)        proton = true;
    else if(target == PT::Neutron) proton = false;
    else                           proton = random->Uniform(0.0, 1.0) < 10.0 / 18.0;

    bool cc = (InteractionType() == 1);
    bool nu = static_cast<int>(primary) > 0;

    // particles in the rest frame of the hadronic system: (PDG code, px, py, pz, mass)
    std::vector<std::array<double,5>> rest;

    // Valence picture: the struck quark and the diquark that is left span one string.
    // CC nu: d -> u. CC nubar: u -> d. NC: a valence quark, flavor unchanged.
    bool strike_u = cc ? !nu : (random->Uniform(0.0, 1.0) < (proton ? 2.0 / 3.0 : 1.0 / 3.0));
    int quark = cc ? (nu ? 2 : 1) : (strike_u ? 2 : 1);
    int diquark;
    if(proton == strike_u) diquark = (random->Uniform(0.0, 1.0) < 0.75) ? 2101 : 2103;   // ud
    else                   diquark = proton ? 2203 : 1103;                                // uu, dd

    double m1 = pythia_->particleData.m0(quark);
    double m2 = pythia_->particleData.m0(diquark);
    if(W > m1 + m2 + MIN_STRING_EXCESS) {
        double pz = std::sqrt((W*W - (m1+m2)*(m1+m2)) * (W*W - (m1-m2)*(m1-m2))) / (2.0 * W);
        double scale = std::max(SHOWER_SCALE_FRACTION * Q, MIN_SHOWER_SCALE);
        Pythia8::Event & ev = pythia_->event;
        bool ok = false;
        for(int tries = 0; tries < HADRONIZE_TRIES && !ok; ++tries) {
            ev.reset();
            // the struck quark moves along +z, which becomes the direction of p4
            ev.append(quark,   23, 101, 0, 0.0, 0.0,  pz, std::sqrt(pz*pz + m1*m1), m1);
            ev.append(diquark, 23, 0, 101, 0.0, 0.0, -pz, std::sqrt(pz*pz + m2*m2), m2);
            ev[1].scale(scale);
            ev[2].scale(scale);
            pythia_->forceTimeShower(1, 2, scale);
            ok = pythia_->next();
        }
        if(ok) {
            for(int i = 0; i < ev.size(); ++i) {
                if(!ev[i].isFinal()) continue;
                rest.push_back({double(ev[i].id()), ev[i].px(), ev[i].py(), ev[i].pz(), std::max(ev[i].m(), 0.0)});
            }
        }
    }

    // Too light for a string: nucleon + pion, isotropic, with the charge of the system.
    if(rest.empty()) {
        int charge = (proton ? 1 : 0) + (cc ? (nu ? 1 : -1) : 0);
        int nucleon, pion;
        bool neutral_pion = random->Uniform(0.0, 1.0) < 2.0 / 3.0;
        if(charge == 2)       { nucleon = 2212; pion = 211; }
        else if(charge == -1) { nucleon = 2112; pion = -211; }
        else if(charge == 1)  { nucleon = neutral_pion ? 2212 : 2112; pion = neutral_pion ? 111 : 211; }
        else                  { nucleon = neutral_pion ? 2112 : 2212; pion = neutral_pion ? 111 : -211; }
        double mN = pythia_->particleData.m0(nucleon);
        double mpi = pythia_->particleData.m0(pion);
        if(W > mN + mpi) {
            double p = std::sqrt((W*W - (mN+mpi)*(mN+mpi)) * (W*W - (mN-mpi)*(mN-mpi))) / (2.0 * W);
            double ct = random->Uniform(-1.0, 1.0);
            double st = std::sqrt(1.0 - ct*ct);
            double phi = random->Uniform(0.0, 2.0 * M_PI);
            double px = p * st * std::cos(phi), py = p * st * std::sin(phi), pz = p * ct;
            rest.push_back({double(nucleon),  px,  py,  pz, mN});
            rest.push_back({double(pion),    -px, -py, -pz, mpi});
        } else {
            rest.push_back({double(proton ? 2212 : 2112), 0.0, 0.0, 0.0, W});
        }
    }

    rk::P4 P(geom3::Vector3(p4[1], p4[2], p4[3]), W);
    rk::Boost toLab = P.labBoost();
    geom3::Rotation3 z_to_p = geom3::rotationBetween(geom3::UnitVector3::zAxis(), P.momentum().direction());

    std::vector<std::pair<int,std::array<double,4>>> out;
    out.reserve(rest.size());
    for(auto const & r : rest) {
        rk::P4 h(geom3::Vector3(r[1], r[2], r[3]), r[4]);
        h.rotate(z_to_p);
        h.boost(toLab);
        std::array<double,4> q = {h.e(), h.px(), h.py(), h.pz()};
        out.emplace_back(static_cast<int>(r[0]), q);
    }
    return out;
}

std::vector<std::pair<int,std::array<double,4>>> GammaDISFromSpline::SelectParticles(
    std::vector<std::pair<int,std::array<double,4>>> particles, double E_had,
    std::shared_ptr<siren::utilities::SIREN_random> random) const {

    std::sort(particles.begin(), particles.end(),
              [](auto const & a, auto const & b){ return a.second[0] > b.second[0]; });

    std::vector<std::pair<int,std::array<double,4>>> fs;
    std::array<double,4> remainder = {0.0, 0.0, 0.0, 0.0};
    int kept = 0;
    for(auto const & pe : particles) {
        int pid = pdg_to_model_pid(pe.first, random);
        if(kept < top_k_ && HasSpecies(pid)) { fs.push_back({pid, pe.second}); ++kept; }
        else for(int k = 0; k < 4; ++k) remainder[k] += pe.second[k];
    }
    if(remainder[0] > 1e-9 * E_had && HasSpecies(111))
        fs.push_back({111, remainder});
    return fs;
}

void GammaDISFromSpline::SampleFinalState(
    dataclasses::CrossSectionDistributionRecord & record,
    std::shared_ptr<siren::utilities::SIREN_random> random) const {

    siren::dataclasses::ParticleType primary = record.signature.primary_type;
    int lep_pdg = outgoing_lepton_pdg(primary, InteractionType());
    double m = DISFromSpline::GetLeptonMass(static_cast<siren::dataclasses::ParticleType>(lep_pdg));
    bool shower_lepton = (std::abs(lep_pdg) == ELECTRON_PDG);

    rk::P4 p1(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), record.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy

    double primary_energy = p1.e();
    double E1_lab = p1.e();
    double E2_lab = p2.e();
    double m1 = record.primary_mass; 
    double m3 = m;

    // (x, y) by Metropolis-Hastings on the differential cross section
    double yMax = 1.0 - m / primary_energy;
    double yMin = GetMinimumQ2() / (2.0 * E1_lab * E2_lab);
    double xMin = GetMinimumQ2() / (2.0 * E1_lab * E2_lab * yMax);
    auto dxs = [&](double x, double y) {
        if(!kinematicallyAllowed(x, y, primary_energy, record.target_mass, m)) return 0.0;
        // x*y is the Jacobian of the proposal, which is uniform in (log x, log y)
        return x * y * DISFromSpline::DifferentialCrossSection(primary_energy, x, y, m);
    };
    // initial valid point
    double x = 0.0, y = 0.0, cs = 0.0;
    for(int tries = 0; tries < 10000 && cs <= 0.0; ++tries) {
        x = std::pow(10.0, random->Uniform(std::log10(xMin), 0.0));
        y = std::pow(10.0, random->Uniform(std::log10(yMin), std::log10(yMax)));
        cs = dxs(x, y);
    }
    const int burnin = 40;                 // same as DISFromSpline
    for(int j = 0; j <= burnin; ++j) {
        double tx = std::pow(10.0, random->Uniform(std::log10(xMin), 0.0));
        double ty = std::pow(10.0, random->Uniform(std::log10(yMin), std::log10(yMax)));
        double tcs = dxs(tx, ty);
        double odds = (cs > 0.0) ? (tcs / cs) : 1.0;
        if(cs == 0.0 || odds > 1.0 || random->Uniform(0.0, 1.0) < odds) {
            x = tx; y = ty; cs = tcs;
        }
    }

    double final_x = x, final_y = y;
    record.interaction_parameters.clear();
    record.interaction_parameters["energy"]    = E1_lab;
    record.interaction_parameters["bjorken_x"] = final_x;
    record.interaction_parameters["bjorken_y"] = final_y;

    // scattered lepton p3 and hadronic system p4, as in DISFromSpline
    double Q2 = 2.0 * E1_lab * E2_lab * final_x * final_y;
    double p1x_lab = std::sqrt(p1.px()*p1.px() + p1.py()*p1.py() + p1.pz()*p1.pz());
    double pqx_lab = (m1*m1 + m3*m3 + 2.0*p1x_lab*p1x_lab + Q2 + 2.0*E1_lab*E1_lab*(final_y-1.0))
                     / (2.0 * p1x_lab);
    double momq_lab = std::sqrt(m1*m1 + p1x_lab*p1x_lab + Q2 + E1_lab*E1_lab*(final_y*final_y-1.0));
    double pqy_lab = (pqx_lab > momq_lab) ? 0.0 : std::sqrt(momq_lab*momq_lab - pqx_lab*pqx_lab);
    double Eq_lab = E1_lab * final_y;

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::UnitVector3 p1_dir = p1.momentum().direction();
    geom3::Rotation3 x_to_p1 = geom3::rotationBetween(x_dir, p1_dir);
    double phi = random->Uniform(0.0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(p1_dir, phi);

    rk::P4 pq_lab(Eq_lab, geom3::Vector3(pqx_lab, pqy_lab, 0.0));
    pq_lab.rotate(x_to_p1);
    pq_lab.rotate(rand_rot);

    rk::P4 p3 = rk::P4((p1 - pq_lab).momentum(), m3);
    rk::P4 p4 = p2 + pq_lab;

    // hadronic system and lepton; DifferentialCrossSection reads the lepton back
    record.interaction_parameters["had_E"]  = p4.e();
    record.interaction_parameters["had_px"] = p4.px();
    record.interaction_parameters["had_py"] = p4.py();
    record.interaction_parameters["had_pz"] = p4.pz();
    record.interaction_parameters["lep_E"]  = p3.e();
    record.interaction_parameters["lep_px"] = p3.px();
    record.interaction_parameters["lep_py"] = p3.py();
    record.interaction_parameters["lep_pz"] = p3.pz();

    auto & secondaries = record.GetSecondaryParticleRecords();
    size_t slot = 0;

    // mu/tau (CC) or the neutrino (NC) stays a regular secondary in slot 0
    if(!shower_lepton && !secondaries.empty()) {
        secondaries[slot].SetFourMomentum({p3.e(), p3.px(), p3.py(), p3.pz()});
        secondaries[slot].SetMass(p3.m());
        secondaries[slot].SetHelicity(record.primary_helicity);
        ++slot;
    }

    // The 4-momentum of each hadron goes into interaction_parameters as hadron<i>_*.
    // The GammaSubshower secondaries carry [amplitude, alpha, beta, parent index]
    // in their 4-momentum slots.
    std::array<double,4> p4arr = {p4.e(), p4.px(), p4.py(), p4.pz()};
    std::vector<std::pair<int,std::array<double,4>>> particles = SelectParticles(
        Hadronize(p4arr, std::sqrt(Q2), primary, record.signature.target_type, random), p4.e(), random);

    int hidx = 0;
    for(auto & pe : particles) {
        int pid = pe.first;
        std::array<double,4> const & q = pe.second;
        std::string pre = "hadron" + std::to_string(hidx);
        record.interaction_parameters[pre + "_pid"] = static_cast<double>(pid);
        record.interaction_parameters[pre + "_E"]   = q[0];
        record.interaction_parameters[pre + "_px"]  = q[1];
        record.interaction_parameters[pre + "_py"]  = q[2];
        record.interaction_parameters[pre + "_pz"]  = q[3];
        int parent_index = hidx;
        ++hidx;

        // parent index i refers to hadron<i>_*, -1 to the electron
        double E_kin = q[0] - rest_mass_GeV(pid);
        std::vector<Component> comps;
        if(E_kin > 0.0) comps = SampleShower(pid, E_kin, random);
        for(auto & c : comps) {
            if(slot >= secondaries.size()) break;   // the signature has kMaxSub_ slots
            secondaries[slot].SetFourMomentum(
                {c.amplitude, c.alpha, c.beta, static_cast<double>(parent_index)});
            secondaries[slot].SetMass(0.0);
            secondaries[slot].SetHelicity(0.0);
            ++slot;
        }
    }
    record.interaction_parameters["had_count"] = static_cast<double>(hidx);

    // CC nu_e: the scattered electron showers as well
    if(shower_lepton && HasSpecies(ELECTRON_PDG)) {
        std::vector<Component> comps = SampleShower(ELECTRON_PDG, p3.e() - rest_mass_GeV(ELECTRON_PDG), random);
        for(auto & c : comps) {
            if(slot >= secondaries.size()) break;
            secondaries[slot].SetFourMomentum(
                {c.amplitude, c.alpha, c.beta, -1.0});
            secondaries[slot].SetMass(0.0);
            secondaries[slot].SetHelicity(0.0);
            ++slot;
        }
    }
    // Finalize() reads mass and 4-momentum from every secondary slot, so the unused
    // ones are set to zero. Readers skip slots with amplitude 0.
    for(; slot < secondaries.size(); ++slot) {
        secondaries[slot].SetFourMomentum({0.0, 0.0, 0.0, 0.0});
        secondaries[slot].SetMass(0.0);
        secondaries[slot].SetHelicity(0.0);
    }
}

// Used by the Weighter. Same x, y, Q^2 as DISFromSpline, with the lepton taken from
// the stored lep_* parameters because the secondaries are GammaSubshowers.
double GammaDISFromSpline::DifferentialCrossSection(
    dataclasses::InteractionRecord const & interaction) const {
    using PT = siren::dataclasses::ParticleType;
    auto const & ip = interaction.interaction_parameters;
    auto get = [&](std::string const & k) {
        auto it = ip.find(k);
        if(it == ip.end())
            throw std::runtime_error("GammaDISFromSpline::DifferentialCrossSection: record has no '" + k + "'");
        return it->second;
    };

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2],
                             interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), interaction.target_mass);
    double primary_energy = interaction.primary_momentum[0];

    int lep_pdg = outgoing_lepton_pdg(interaction.signature.primary_type, InteractionType());
    double lepton_mass = DISFromSpline::GetLeptonMass(static_cast<PT>(lep_pdg));
    rk::P4 p3(geom3::Vector3(get("lep_px"), get("lep_py"), get("lep_pz")), lepton_mass);

    rk::P4 q = p1 - p3;
    double Q2 = -q.dot(q);
    double y = 1.0 - p2.dot(p3) / p2.dot(p1);
    double x = Q2 / (2.0 * p2.dot(q));
    return DISFromSpline::DifferentialCrossSection(primary_energy, x, y, lepton_mass, Q2);
}

} // namespace interactions
} // namespace siren
