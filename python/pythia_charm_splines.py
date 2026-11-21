"""Generate Pythia-derived charm-DIS splines for ``PythiaDISCrossSection``.

``PythiaDISCrossSection`` samples the final state from Pythia and uses splines for
the cross-section *rate* (and, optionally, a differential density). These helpers
build those splines directly from Pythia so they correspond to exactly the events
the sampler produces:

  * :func:`generate_total_spline`        -> 1D ``sigma(E)`` FITS (always required).
  * :func:`generate_differential_spline` -> 3D ``d2sigma/dx dy`` FITS (optional;
    supplying it gives ``FinalStateProbability`` a real density so weights stay
    correct under reweighting; omitting it makes ``FinalStateProbability`` a
    constant that cancels in the unbiased weight).

Both run Pythia via ``PythiaDISCrossSection.GeneratePythiaCharmSamples`` and fit
with photospline. ``LHAPDF_DATA_PATH`` must be set in the environment, and SIREN
must be built with ``SIREN_WITH_PYTHIA8=ON``.

The differential is fit in the muon-reconstructed Bjorken ``x`` (the same variable
``DifferentialCrossSection`` reconstructs), so it closes against the sampler.
"""
import numpy as np


def _knots(c, order=2):
    c = np.asarray(c, float)
    d = c[1] - c[0]
    return np.concatenate([c[0] - d * np.arange(order, 0, -1), c,
                           c[-1] + d * np.arange(1, order + 1)])


def _run_pythia(interaction_type, primary_pdg, target_pdg, target_mass, pdf_set,
                pythia_data_path, minimum_Q2, energies, n_events):
    import siren.interactions
    sigma_mb, E, x, y = siren.interactions.PythiaDISCrossSection.GeneratePythiaCharmSamples(
        int(interaction_type), int(primary_pdg), int(target_pdg), float(target_mass),
        str(pdf_set), str(pythia_data_path), float(minimum_Q2),
        [float(e) for e in energies], int(n_events))
    return np.asarray(sigma_mb), np.asarray(E), np.asarray(x), np.asarray(y)


def generate_total_spline(out_path, *, interaction_type, primary_pdg, target_pdg,
                          target_mass, pdf_set, pythia_data_path, energies,
                          minimum_Q2=1.0, n_events=20000, mb_to_cm2=1e-27):
    """Fit and write the 1D total cross-section spline ``log10(sigma_cm2)`` vs ``log10(E)``."""
    import photospline
    sigma_mb, _, _, _ = _run_pythia(interaction_type, primary_pdg, target_pdg,
                                    target_mass, pdf_set, pythia_data_path,
                                    minimum_Q2, energies, n_events)
    logE = np.log10(np.asarray(energies, float))
    logsig = np.log10(sigma_mb * mb_to_cm2)
    z, w = photospline.ndsparse.from_data(logsig, np.ones_like(logsig))
    spline = photospline.glam_fit(z, w, [logE], [_knots(logE)], [2], [1e-3], [2])
    spline.write(out_path)
    return out_path


def generate_differential_spline(out_path, *, interaction_type, primary_pdg, target_pdg,
                                 target_mass, pdf_set, pythia_data_path, energies,
                                 minimum_Q2=1.0, n_events=50000,
                                 logx_range=(-4.0, -0.1), logy_range=(-2.3, -0.004),
                                 nbins=12, mb_to_cm2=1e-27):
    """Fit and write the 3D differential spline ``log10(d2sigma/dx dy)`` vs ``(log10 E, log10 x, log10 y)``.

    The per-bin density is ``sigma(E) * dN/dx dy / N`` (log-space Jacobian folded
    in) so that ``FinalStateProbability = dsigma/sigma`` reproduces the sampled
    ``(x, y)`` density.
    """
    import photospline
    sigma_mb, E, x, y = _run_pythia(interaction_type, primary_pdg, target_pdg,
                                    target_mass, pdf_set, pythia_data_path,
                                    minimum_Q2, energies, n_events)
    energies = np.asarray(energies, float)
    logE_grid = np.log10(energies)
    logx_edges = np.linspace(logx_range[0], logx_range[1], nbins + 1)
    logy_edges = np.linspace(logy_range[0], logy_range[1], nbins + 1)
    logx_c = 0.5 * (logx_edges[:-1] + logx_edges[1:])
    logy_c = 0.5 * (logy_edges[:-1] + logy_edges[1:])
    dlx = logx_edges[1] - logx_edges[0]
    dly = logy_edges[1] - logy_edges[0]
    ln10 = np.log(10.0)
    nE, nx, ny = len(energies), len(logx_c), len(logy_c)
    Z = np.full((nE, nx, ny), -99.0)
    W = np.zeros((nE, nx, ny))
    for ie, Eval in enumerate(energies):
        sigma = sigma_mb[ie] * mb_to_cm2
        m = (E == Eval)
        Ntot = int(m.sum())
        if Ntot == 0:
            continue
        H, _, _ = np.histogram2d(np.log10(x[m]), np.log10(y[m]),
                                 bins=[logx_edges, logy_edges])
        for ix in range(nx):
            xc = 10.0 ** logx_c[ix]
            for iy in range(ny):
                c = H[ix, iy]
                if c <= 0:
                    continue
                yc = 10.0 ** logy_c[iy]
                dxdy = (xc * ln10 * dlx) * (yc * ln10 * dly)
                Z[ie, ix, iy] = np.log10(sigma * c / (Ntot * dxdy))
                W[ie, ix, iy] = c
    z, w = photospline.ndsparse.from_data(Z, W)
    spline = photospline.glam_fit(
        z, w, [logE_grid, logx_c, logy_c],
        [_knots(logE_grid), _knots(logx_c), _knots(logy_c)],
        [2, 2, 2], [5e-2, 5e-2, 5e-2], [2, 2, 2])
    spline.write(out_path)
    return out_path
