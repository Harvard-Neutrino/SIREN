// Reference generator for the QuarkDISFromSpline slow-rescaling charm-CC splines.
//
// QuarkDISFromSpline samples charm-DIS final states from FITS splines of the
// slow-rescaling differential cross section d2sigma/dxi dy and the total
// sigma(E). The production analysis supplies its own (Dutta-Kim CT14) splines;
// this standalone tool produces physically-validated splines from any installed
// LHAPDF set so the slow-rescaling sampler and normalization can be exercised
// locally (see tests/python/test_quarkdis_slow_rescaling.py).
//
// Physics (leading-order slow rescaling, lab frame, stationary nucleon), exactly
// matching QuarkDISFromSpline's kinematics (slowRescalingQ2 / kinematicallyAllowed):
//   xi  = struck-quark momentum fraction (the PDF argument)
//   Q^2 = 2 M E y xi - m_c^2,   W^2 = M^2 + 2 M E y (1-xi) + m_c^2
//   d2sigma/dxi dy = (G_F^2 M E / pi) (M_W^2/(Q^2+M_W^2))^2
//                    [ |V_cd|^2 xd(xi,Q^2) + |V_cs|^2 xs(xi,Q^2)
//                    + (|V_cd|^2 xdbar + |V_cs|^2 xsbar)(1-y)^2 ]
// PDFs are evaluated at xi (NOT at the measured Bjorken x). The integrated charm
// fraction reproduces ~4% at 100 GeV rising to ~6% at 10 TeV (literature band,
// cross-checks the PythiaDISCrossSection charm fraction ~6.5%).
//
// Build/run:
//   g++ -std=c++17 $(lhapdf-config --cflags --libs) \
//       generate_charm_slowrescaling_splines.cpp -o gen_sr
//   LHAPDF_DATA_PATH=<dir> ./gen_sr            # writes xidiff.txt + xisigma.txt
//   python fit_charm_slowrescaling_splines.py  # -> the two FITS splines
#include <LHAPDF/LHAPDF.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

int main(int argc, char** argv) {
    std::string pdfname = (argc > 1) ? argv[1] : "CT18NLO";
    LHAPDF::PDF* pdf = LHAPDF::mkPDF(pdfname, 0);

    const double GF = 1.1663787e-5, M = 0.9314943, MW = 80.379;
    const double mc = 1.280, MD0 = 1.86484, mmu = 0.105658;
    const double Vcd2 = 0.221 * 0.221, Vcs2 = 0.975 * 0.975;
    const double GEV2_TO_CM2 = 0.389379e-27;   // (hbar c)^2

    auto dsig = [&](double E, double xi, double y) -> double {
        double Q2 = 2 * M * E * y * xi - mc * mc;       if (Q2 <= 1.0) return 0.0;
        double W2 = M * M + 2 * M * E * y * (1 - xi) + mc * mc;
        if (W2 <= (M + MD0) * (M + MD0)) return 0.0;
        if (y >= 1 - mmu / E) return 0.0;
        double P1 = E, pqx = (mmu * mmu + 2 * P1 * P1 + Q2 + 2 * E * E * (y - 1)) / (2 * P1);
        double mq2 = P1 * P1 + Q2 + E * E * (y * y - 1);
        if (mq2 - pqx * pqx < 0) return 0.0;
        if (xi <= 0 || xi >= 1) return 0.0;
        double xd = pdf->xfxQ2(1, xi, Q2),  xs = pdf->xfxQ2(3, xi, Q2);
        double xdb = pdf->xfxQ2(-1, xi, Q2), xsb = pdf->xfxQ2(-3, xi, Q2);
        double pref = GF * GF * M * E / M_PI * std::pow(MW * MW / (Q2 + MW * MW), 2);
        double v = pref * ((Vcd2 * xd + Vcs2 * xs) + (Vcd2 * xdb + Vcs2 * xsb) * (1 - y) * (1 - y));
        return v > 0 ? v * GEV2_TO_CM2 : 0.0;       // cm^2
    };

    int nE = 20, nxi = 26, ny = 26;
    std::vector<double> Ev(nE), lxi(nxi), ly(ny);
    for (int i = 0; i < nE; i++)  Ev[i]  = std::pow(10, 1.0 + 5.0 * i / (nE - 1));     // 10 GeV - 1 PeV
    for (int i = 0; i < nxi; i++) lxi[i] = -4.0 + (std::log10(0.95) + 4.0) * i / (nxi - 1);
    for (int i = 0; i < ny; i++)  ly[i]  = -3.0 + (std::log10(0.95) + 3.0) * i / (ny - 1);

    FILE* fd = std::fopen("xidiff.txt", "w");
    std::fprintf(fd, "%d %d %d\n", nE, nxi, ny);
    for (double v : Ev)  std::fprintf(fd, "%.10e\n", v);
    for (double v : lxi) std::fprintf(fd, "%.10e\n", v);
    for (double v : ly)  std::fprintf(fd, "%.10e\n", v);
    for (int i = 0; i < nE; i++)
        for (int j = 0; j < nxi; j++)
            for (int k = 0; k < ny; k++) {
                double v = dsig(Ev[i], std::pow(10, lxi[j]), std::pow(10, ly[k]));
                if (v > 0) std::fprintf(fd, "%.8e 1\n", std::log10(v));
                else       std::fprintf(fd, "-99 0\n");
            }
    std::fclose(fd);

    FILE* fs = std::fopen("xisigma.txt", "w");
    for (double E : Ev) {
        int Nx = 400, Ny = 400; double lo = -4, hi = std::log10(0.95), tot = 0;
        double ymax = 1 - mmu / E - 1e-3, dy = (ymax - 1e-3) / Ny;
        for (int a = 0; a < Nx; a++) {
            double l0 = lo + (hi - lo) * a / Nx, l1 = lo + (hi - lo) * (a + 1) / Nx;
            double xi = std::pow(10, 0.5 * (l0 + l1)), dxi = std::pow(10, l1) - std::pow(10, l0), col = 0;
            for (int b = 0; b < Ny; b++) { double y = 1e-3 + (b + 0.5) * dy; col += dsig(E, xi, y) * dy; }
            tot += col * dxi;
        }
        std::fprintf(fs, "%.10e %.10e\n", E, tot);
    }
    std::fclose(fs);
    std::printf("wrote xidiff.txt (%dx%dx%d) + xisigma.txt with PDF %s\n", nE, nxi, ny, pdfname.c_str());
    return 0;
}
