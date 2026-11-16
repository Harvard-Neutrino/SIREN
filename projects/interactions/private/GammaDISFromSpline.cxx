#include "SIREN/interactions/GammaDISFromSpline.h"

#include <cmath>
#include <cctype>
#include <array>
#include <tuple>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include "cnpy.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/Random.h"

namespace siren {
namespace interactions {

namespace {
    constexpr double EARLY_SIGMA = 2.5;  // reject draws this many sigma too early or too narrow
    constexpr int    MAX_TRIES   = 100;
    constexpr double COV_RIDGE  = 1e-3;  // floor on the covariance eigenvalues
    constexpr int    ELECTRON_PDG = 11;

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

GammaDISFromSpline::GammaDISFromSpline(
    std::string differential_filename, std::string total_filename,
    std::vector<siren::dataclasses::ParticleType> primary_types,
    std::vector<siren::dataclasses::ParticleType> target_types,
    std::string shower_model_path, std::string pythia_events_path,
    std::string units)
    : DISFromSpline(differential_filename, total_filename,
                    primary_types, target_types, units),
      shower_model_path_(shower_model_path),
      pythia_events_path_(pythia_events_path) {
    // Spline files without an INTERACTION key are read as CC by the base class.
    // The CSMS splines have no such key, so the file name decides CC or NC.
    std::string base = differential_filename.substr(differential_filename.find_last_of("/\\") + 1);
    for(char & ch : base) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    if(GetInteractionType() == 1 &&
       (base.find("_NC_") != std::string::npos || base.find("_NC.") != std::string::npos))
        interaction_type_override_ = 2;
    LoadShowerModel(shower_model_path_);
    LoadPythiaEvents(pythia_events_path_);
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
    std::string shower_model_path, std::string pythia_events_path,
    std::string units)
    : DISFromSpline(differential_filename, total_filename, interaction,
                    target_mass, minimum_Q2, primary_types, target_types, units),
      shower_model_path_(shower_model_path),
      pythia_events_path_(pythia_events_path) {
    LoadShowerModel(shower_model_path_);
    LoadPythiaEvents(pythia_events_path_);
    InitializeGammaSignatures();
}

std::string GammaDISFromSpline::Name() const { return "GammaDISFromSpline"; }

bool GammaDISFromSpline::equal(CrossSection const & other) const {
    const GammaDISFromSpline* x = dynamic_cast<const GammaDISFromSpline*>(&other);
    if(!x) return false;
    return shower_model_path_ == x->shower_model_path_ &&
           pythia_events_path_ == x->pythia_events_path_ &&
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

void GammaDISFromSpline::LoadPythiaEvents(std::string const & npz_path) {
    cnpy::npz_t npz;
    try { npz = cnpy::npz_load(npz_path); }
    catch(std::exception const & e) {
        throw std::runtime_error("GammaDISFromSpline: failed to read '" + npz_path + "': " + e.what());
    }
    auto get = [&](std::string const & k) -> cnpy::NpyArray & {
        auto it = npz.find(k);
        if(it == npz.end()) throw std::runtime_error("pythia store missing key: " + k);
        return it->second;
    };
    cnpy::NpyArray & EH = get("E_had");          // [Nev]
    cnpy::NpyArray & PD = get("top20_pids");     // [Nev, 20] int64
    cnpy::NpyArray & EN = get("top20_energies"); // [Nev, 20] double
    size_t Nev = EH.shape[0];
    size_t cols = PD.shape.size() > 1 ? PD.shape[1] : 20;

    const double*  eh = EH.data<double>();
    const int64_t* pd = PD.data<int64_t>();
    const double*  en = EN.data<double>();

    pythia_Ehad_.assign(eh, eh + Nev);
    pythia_pids_.resize(Nev);
    pythia_energies_.resize(Nev);
    for(size_t i = 0; i < Nev; ++i) {
        for(size_t j = 0; j < cols; ++j) {
            pythia_pids_[i].push_back(static_cast<int>(pd[i*cols + j]));
            pythia_energies_[i].push_back(en[i*cols + j]);
        }
    }

    // top20_momenta is [Nev, 20, 4] and p_had is [Nev, 4]. Without them the
    // breakdown falls back to scaling the energies.
    have_momenta_ = (npz.find("top20_momenta") != npz.end()) &&
                    (npz.find("p_had") != npz.end());
    if(have_momenta_) {
        cnpy::NpyArray & MOM = get("top20_momenta");   // [Nev, 20, 4]
        cnpy::NpyArray & PH  = get("p_had");           // [Nev, 4]
        const double* mom = MOM.data<double>();
        const double* ph  = PH.data<double>();
        pythia_momenta_.resize(Nev);
        pythia_phad_.resize(Nev);
        for(size_t i = 0; i < Nev; ++i) {
            pythia_momenta_[i].resize(cols);
            for(size_t j = 0; j < cols; ++j)
                for(int k = 0; k < 4; ++k)
                    pythia_momenta_[i][j][k] = mom[(i*cols + j)*4 + k];
            for(int k = 0; k < 4; ++k) pythia_phad_[i][k] = ph[i*4 + k];
        }
    }

    // quantile bin edges in E_had
    std::vector<double> sorted = pythia_Ehad_;
    std::sort(sorted.begin(), sorted.end());
    bin_edges_.resize(n_bins_ + 1);
    for(int b = 0; b <= n_bins_; ++b) {
        double q = double(b) / n_bins_;
        double pos = q * (sorted.size() - 1);
        size_t lo = static_cast<size_t>(std::floor(pos));
        size_t hi = std::min(lo + 1, sorted.size() - 1);
        double frac = pos - lo;
        bin_edges_[b] = sorted[lo] * (1 - frac) + sorted[hi] * frac;
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

int GammaDISFromSpline::bin_index_(double E_had) const {
    if(bin_edges_.empty()) return 0;
    int b = 0;
    while(b < n_bins_ && E_had >= bin_edges_[b+1]) ++b;
    return std::min(std::max(b, 0), n_bins_ - 1);
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

// Pick a library event near E_had and Lorentz-transform its hadrons so that their
// total 4-momentum equals p4, the hadronic system from the DIS kinematics.
std::vector<std::pair<int,std::array<double,4>>> GammaDISFromSpline::BoostedBreakdown(
    std::array<double,4> const & p4,
    std::shared_ptr<siren::utilities::SIREN_random> random) const {

    std::vector<std::pair<int,std::array<double,4>>> fs;
    if(pythia_Ehad_.empty()) return fs;
    double E_had = p4[0];

    // nearest non-empty E_had bin, then a uniform event in it
    int b = bin_index_(E_had);
    std::vector<size_t> pool;
    for(int radius = 0; radius < n_bins_ && pool.empty(); ++radius) {
        for(int bb : {b - radius, b + radius}) {
            if(bb < 0 || bb >= n_bins_) continue;
            double lo = bin_edges_[bb], hi = bin_edges_[bb+1];
            for(size_t i = 0; i < pythia_Ehad_.size(); ++i)
                if(pythia_Ehad_[i] >= lo && (pythia_Ehad_[i] < hi || bb == n_bins_-1))
                    pool.push_back(i);
            if(radius == 0) break;
        }
    }
    if(pool.empty()) { for(size_t i=0;i<pythia_Ehad_.size();++i) pool.push_back(i); }
    size_t ievt = pool[std::min<size_t>(pool.size()-1,
                      static_cast<size_t>(random->Uniform(0.0,1.0) * pool.size()))];

    // ranked list of (energy, pid, boosted 4-mom [E,px,py,pz])
    std::vector<std::tuple<double,int,std::array<double,4>>> ranked;
    // rk::P4(E, p) takes sqrt(E^2 - p^2), which round-off can turn into NaN for
    // massless particles, so P4 is built from the momentum and a clamped mass.
    auto clamped_mass = [](double E, double px, double py, double pz) {
        return std::sqrt(std::max(0.0, E*E - px*px - py*py - pz*pz));
    };
    // transformed 4-momentum of the library particles that are not stored
    std::array<double,4> leftover = {0.0, 0.0, 0.0, 0.0};
    if(have_momenta_ && !pythia_momenta_[ievt].empty()) {
        std::array<double,4> lib_sum = {0.0, 0.0, 0.0, 0.0};   // stored particles, library frame
        rk::P4 P_tgt(geom3::Vector3(p4[1], p4[2], p4[3]),
                     clamped_mass(p4[0], p4[1], p4[2], p4[3]));
        auto const & ph = pythia_phad_[ievt];
        rk::P4 P_lib(geom3::Vector3(ph[1], ph[2], ph[3]),
                     clamped_mass(ph[0], ph[1], ph[2], ph[3]));
        rk::Boost toCM  = P_lib.restBoost();     // library frame -> CM
        rk::Boost toLab = P_tgt.labBoost();      // CM -> lab
        double Mlib = P_lib.m(), Mtgt = P_tgt.m();
        double lam  = (Mlib > 1e-9) ? (Mtgt / Mlib) : 1.0;   // match invariant mass
        for(size_t j = 0; j < pythia_pids_[ievt].size(); ++j) {
            int pid = pythia_pids_[ievt][j];
            auto const & q = pythia_momenta_[ievt][j];
            if(pid == 0 || !(q[0] > 0.0)) continue;
            double mj = clamped_mass(q[0], q[1], q[2], q[3]);
            rk::P4 h(geom3::Vector3(q[1], q[2], q[3]), mj);
            h.boost(toCM);
            h = rk::P4(h.momentum() * lam, mj * lam);            // rescale in CM
            h.boost(toLab);
            std::array<double,4> hq = {h.e(), h.px(), h.py(), h.pz()};
            for(int k = 0; k < 4; ++k) lib_sum[k] += q[k];
            ranked.push_back({h.e(), pid, hq});
        }
        // p_had sums all particles but only the leading ones are stored. The rest,
        // L = p_had - Sum(stored), gets the same transformation, so the hadrons
        // still add up to p4.
        double L[4] = {ph[0]-lib_sum[0], ph[1]-lib_sum[1], ph[2]-lib_sum[2], ph[3]-lib_sum[3]};
        if(L[0] > 1e-12 * ph[0]) {
            double mL = clamped_mass(L[0], L[1], L[2], L[3]);
            rk::P4 hL(geom3::Vector3(L[1], L[2], L[3]), mL);
            hL.boost(toCM);
            hL = rk::P4(hL.momentum() * lam, mL * lam);
            hL.boost(toLab);
            leftover = {hL.e(), hL.px(), hL.py(), hL.pz()};
        }
    } else {
        // no momenta in the library: scale the energies
        double eh_lib = pythia_Ehad_[ievt];
        double scale = (eh_lib > 0.0) ? (E_had / eh_lib) : 1.0;
        for(size_t j = 0; j < pythia_pids_[ievt].size(); ++j) {
            int pid = pythia_pids_[ievt][j];
            double E = pythia_energies_[ievt][j] * scale;
            if(pid == 0 || !(E > 0.0)) continue;
            ranked.push_back({E, pid, {E, 0.0, 0.0, 0.0}});
        }
    }

    std::sort(ranked.begin(), ranked.end(),
              [](auto const & a, auto const & b){ return std::get<0>(a) > std::get<0>(b); });
    std::array<double,4> remainder = {0.0, 0.0, 0.0, 0.0};
    int kept = 0;
    for(auto & t : ranked) {
        int pid = pdg_to_model_pid(std::get<1>(t), random);
        std::array<double,4> const & q4 = std::get<2>(t);
        if(kept < top_k_ && HasSpecies(pid)) { fs.push_back({pid, q4}); ++kept; }
        else for(int k = 0; k < 4; ++k) remainder[k] += q4[k];
    }
    // the particles that are not stored go into the remainder as well
    for(int k = 0; k < 4; ++k) remainder[k] += leftover[k];
    if(remainder[0] > 1e-9 * E_had && HasSpecies(111))
        fs.push_back({111, remainder});   // pi0 remainder
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
    std::vector<std::pair<int,std::array<double,4>>> particles = BoostedBreakdown(p4arr, random);

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
