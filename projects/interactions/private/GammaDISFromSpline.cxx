#include "SIREN/interactions/GammaDISFromSpline.h"

#include <math.h>   // before rk: geom3::rotationBetween calls an unqualified abs on doubles
#include <cmath>
#include <cctype>
#include <fstream>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include <sys/stat.h>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include <Pythia8/Pythia.h>

#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Errors.h"

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
    constexpr int    PYTHIA_TRIES = 100;

    // one species as stored in the shower model file; 2-d arrays are flattened row-major
    struct SpeciesFile {
        int pid = 0;
        std::vector<double> logE, pm, logN, logNsig;
        std::vector<int> m;
        std::vector<std::vector<double>> mean, cov;   // one entry per m
        template<typename Archive>
        void serialize(Archive & archive) {
            archive(cereal::make_nvp("pid", pid), cereal::make_nvp("logE", logE),
                    cereal::make_nvp("pm", pm), cereal::make_nvp("logN", logN),
                    cereal::make_nvp("logNsig", logNsig), cereal::make_nvp("m", m),
                    cereal::make_nvp("mean", mean), cereal::make_nvp("cov", cov));
        }
    };

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
    std::string pdf_set, std::string units)
    : DISFromSpline(differential_filename, total_filename,
                    primary_types, target_types, units),
      shower_model_path_(shower_model_path),
      pythia_data_path_(pythia_data_path),
      pdf_set_(pdf_set) {
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
    std::string pdf_set, std::string units)
    : DISFromSpline(differential_filename, total_filename, interaction,
                    target_mass, minimum_Q2, primary_types, target_types, units),
      shower_model_path_(shower_model_path),
      pythia_data_path_(pythia_data_path),
      pdf_set_(pdf_set) {
    LoadShowerModel(shower_model_path_);
    InitializeGammaSignatures();
}

std::string GammaDISFromSpline::Name() const { return "GammaDISFromSpline"; }

bool GammaDISFromSpline::equal(CrossSection const & other) const {
    const GammaDISFromSpline* x = dynamic_cast<const GammaDISFromSpline*>(&other);
    if(!x) return false;
    return shower_model_path_ == x->shower_model_path_ &&
           pythia_data_path_ == x->pythia_data_path_ &&
           pdf_set_ == x->pdf_set_ &&
           InteractionType() == x->InteractionType() &&
           DISFromSpline::equal(other);
}

void GammaDISFromSpline::LoadShowerModel(std::string const & path) {
    std::ifstream in(path);
    if(!in)
        throw std::runtime_error("GammaDISFromSpline: cannot open shower model '" + path + "'");
    std::vector<SpeciesFile> file;
    try {
        cereal::JSONInputArchive archive(in);
        archive(cereal::make_nvp("Kmax", Kmax_), cereal::make_nvp("species", file));
    } catch(std::exception const & e) {
        throw std::runtime_error("GammaDISFromSpline: failed to read '" + path + "': " + e.what());
    }

    species_.clear();
    for(auto & f : file) {
        if(f.logE.empty() || f.mean.size() != f.m.size() || f.cov.size() != f.m.size())
            throw std::runtime_error("GammaDISFromSpline: malformed shower model entry for pid "
                                     + std::to_string(f.pid));
        SpeciesModel s;
        s.pid = f.pid;
        s.logE = std::move(f.logE);
        s.logE_lo = s.logE.front();
        s.logE_hi = s.logE.back();
        s.pm = std::move(f.pm);
        s.logN = std::move(f.logN);
        s.logNsig = std::move(f.logNsig);
        s.m_avail = f.m;
        for(size_t j = 0; j < f.m.size(); ++j) {
            s.mean[f.m[j]] = std::move(f.mean[j]);
            s.cov[f.m[j]] = std::move(f.cov[j]);
        }
        species_[s.pid] = std::move(s);
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

Pythia8::Pythia & GammaDISFromSpline::GetPythia(int beam, int target, double E_nu,
    std::shared_ptr<siren::utilities::SIREN_random> random) const {
    PythiaInstance & inst = pythia_[{beam, target}];
    if(!inst.pythia || inst.energy != E_nu) {
        const char* lhapdf_path = std::getenv("LHAPDF_DATA_PATH");
        if(!lhapdf_path)
            throw std::runtime_error("GammaDISFromSpline: LHAPDF_DATA_PATH is not set");
        std::string set_name = pdf_set_.substr(pdf_set_.find(':') + 1);
        std::string set_dir = std::string(lhapdf_path) + "/" + set_name;
        struct stat st;
        if(stat(set_dir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode))
            throw std::runtime_error("GammaDISFromSpline: PDF set '" + set_name
                                     + "' not found under LHAPDF_DATA_PATH=" + lhapdf_path);

        inst.pythia = std::make_unique<Pythia8::Pythia>(
            pythia_data_path_.empty() ? std::string("../share/Pythia8/xmldoc") : pythia_data_path_, false);
        Pythia8::Pythia & p = *inst.pythia;
        p.readString(InteractionType() == 1 ? "WeakBosonExchange:ff2ff(t:W) = on"
                                            : "WeakBosonExchange:ff2ff(t:gmZ) = on");
        p.readString("Beams:frameType = 2");
        p.readString("Beams:idA = " + std::to_string(beam));
        p.readString("Beams:idB = " + std::to_string(target));
        p.readString("Beams:eA = " + std::to_string(E_nu));
        p.readString("Beams:eB = 0.");
        p.readString("PDF:pSet = " + pdf_set_);
        p.readString("PDF:useHard = on");
        p.readString("PDF:pHardSet = " + pdf_set_);
        p.readString("PartonLevel:MPI = off");
        // keeps the scattered lepton's kinematics when the initial-state shower recoils
        p.readString("SpaceShower:dipoleRecoil = on");
        p.readString("PhaseSpace:mHatMin = 0.0");
        p.readString("PhaseSpace:Q2Min = " + std::to_string(GetMinimumQ2()));
        // Particles with c*tau below 10 mm decay. pi0 has a shower model and the tau is
        // returned as a track, so both are kept.
        p.readString("ParticleDecays:limitTau0 = on");
        p.readString("ParticleDecays:tau0Max = 10.0");
        p.readString("111:mayDecay = off");
        p.readString("15:mayDecay = off");
        p.readString("Print:quiet = on");
        p.readString("Next:numberCount = 0");
        if(!p.init())
            throw std::runtime_error("GammaDISFromSpline: Pythia initialization failed");

        // after init(), which uses Pythia's own generator
        inst.rndm = std::make_shared<GammaDISRndm>();
        p.rndm.rndmEnginePtr(inst.rndm);
        inst.energy = E_nu;
    }
    inst.rndm->rng = random;
    return *inst.pythia;
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
    bool shower_lepton = (std::abs(lep_pdg) == ELECTRON_PDG);

    double E_nu = record.primary_momentum[0];
    rk::P4 p1(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);

    // an unspecified nucleon is a proton 10 times out of 18, as in PythiaDISCrossSection
    int target_pdg;
    if(record.signature.target_type == siren::dataclasses::ParticleType::PPlus)        target_pdg = 2212;
    else if(record.signature.target_type == siren::dataclasses::ParticleType::Neutron) target_pdg = 2112;
    else target_pdg = (random->Uniform(0.0, 1.0) < 10.0 / 18.0) ? 2212 : 2112;

    Pythia8::Pythia & pythia = GetPythia(static_cast<int>(primary), target_pdg, E_nu, random);
    Pythia8::Event & ev = pythia.event;

    // the outgoing lepton is the most energetic final-state particle of its type
    int i_lep = -1;
    for(int tries = 0; tries < PYTHIA_TRIES && i_lep < 0; ++tries) {
        if(!pythia.next()) continue;
        double best = -1.0;
        for(int i = 0; i < ev.size(); ++i) {
            if(ev[i].isFinal() && ev[i].id() == lep_pdg && ev[i].e() > best) {
                best = ev[i].e();
                i_lep = i;
            }
        }
    }
    if(i_lep < 0)
        throw siren::utilities::InjectionFailure("GammaDISFromSpline: no Pythia event with the expected lepton after "
                                                 + std::to_string(PYTHIA_TRIES) + " attempts");

    // x and y from the lepton, in Pythia's frame (neutrino along +z), as in PythiaDISCrossSection
    Pythia8::Vec4 q = Pythia8::Vec4(0.0, 0.0, E_nu, E_nu) - ev[i_lep].p();
    double Q2 = -q.m2Calc();
    double y = 1.0 - ev[i_lep].e() / E_nu;
    double x = (y > 0.0) ? Q2 / (2.0 * record.target_mass * E_nu * y) : 0.0;
    record.interaction_parameters.clear();
    record.interaction_parameters["energy"]     = E_nu;
    record.interaction_parameters["bjorken_x"]  = x;
    record.interaction_parameters["bjorken_y"]  = y;
    record.interaction_parameters["target_pdg"] = static_cast<double>(target_pdg);

    // Pythia's frame -> lab: rotate +z onto the neutrino direction
    geom3::Rotation3 rot = geom3::rotationBetween(geom3::UnitVector3::zAxis(), p1.momentum().direction());
    auto to_lab = [&](Pythia8::Particle const & p) {
        geom3::Vector3 mom = rot * geom3::Vector3(p.px(), p.py(), p.pz());
        return std::array<double,4>{p.e(), mom.x(), mom.y(), mom.z()};
    };

    std::array<double,4> lep = to_lab(ev[i_lep]);
    rk::P4 p3(geom3::Vector3(lep[1], lep[2], lep[3]), std::max(ev[i_lep].m(), 0.0));

    // Everything else is the hadronic system. Neutrinos from decays are invisible and
    // are kept out of the particle list, but stay in had_* for the conservation check.
    std::vector<std::pair<int,std::array<double,4>>> hadrons;
    std::array<double,4> had = {0.0, 0.0, 0.0, 0.0};
    std::array<double,4> invisible = {0.0, 0.0, 0.0, 0.0};
    for(int i = 0; i < ev.size(); ++i) {
        if(!ev[i].isFinal() || i == i_lep) continue;
        std::array<double,4> h = to_lab(ev[i]);
        for(int k = 0; k < 4; ++k) had[k] += h[k];
        int id = std::abs(ev[i].id());
        if(id == 12 || id == 14 || id == 16) {
            for(int k = 0; k < 4; ++k) invisible[k] += h[k];
            continue;
        }
        hadrons.emplace_back(ev[i].id(), h);
    }

    // hadronic system and lepton; DifferentialCrossSection reads the lepton back
    record.interaction_parameters["had_E"]  = had[0];
    record.interaction_parameters["had_px"] = had[1];
    record.interaction_parameters["had_py"] = had[2];
    record.interaction_parameters["had_pz"] = had[3];
    record.interaction_parameters["lep_E"]  = p3.e();
    record.interaction_parameters["lep_px"] = p3.px();
    record.interaction_parameters["lep_py"] = p3.py();
    record.interaction_parameters["lep_pz"] = p3.pz();
    record.interaction_parameters["invisible_E"]  = invisible[0];
    record.interaction_parameters["invisible_px"] = invisible[1];
    record.interaction_parameters["invisible_py"] = invisible[2];
    record.interaction_parameters["invisible_pz"] = invisible[3];

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
    std::vector<std::pair<int,std::array<double,4>>> particles = SelectParticles(hadrons, had[0], random);

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

double GammaDISFromSpline::FinalStateProbability(dataclasses::InteractionRecord const &) const {
    return 1.0;
}

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
