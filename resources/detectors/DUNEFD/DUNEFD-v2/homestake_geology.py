#!/usr/bin/env python3
"""
Homestake / SURF 4850L geology: schematic rock units with density + isotopic
composition, for the DUNE FD 3D overburden model.

Provides, per rock unit, a SIREN material (PDG isotope -> mass fraction) built by
expanding a representative major-oxide chemistry into natural isotopes, with the
measured U/Th/K assays folded in. Writes a SIREN materials.dat.

SCHEMATIC: the major-oxide chemistries are representative literature values for
the rock TYPES (mafic amphibolite, pelitic/carbonate phyllite, felsic rhyolite);
the U/Th/K are the measured SURF assays. Bulk densities are representative.
Stratigraphy/structure at SURF is multiply-deformed (Lead anticlinorium), so the
unit GEOMETRY (handled separately) is schematic too.

Sources: USGS Geolex (Poorman/Homestake/Yates/Ellison); SURF radioassays
(amphibolite 0.22/0.33 ppm U/Th, 0.96% K; Poorman 2.58/10.48 ppm, 2.12% K;
rhyolite 8.75/10.86 ppm, 4.17% K); DUNE DUSEL_Rock (duneggd) for the average.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# Standard atomic weights (for oxide -> element mass conversion)
# ---------------------------------------------------------------------------
AW = dict(H=1.008, C=12.011, N=14.007, O=15.999, Na=22.990, Mg=24.305,
          Al=26.982, Si=28.085, P=30.974, S=32.06, K=39.098, Ca=40.078,
          Ti=47.867, Mn=54.938, Fe=55.845, U=238.029, Th=232.038, Ar=39.948)
Z = dict(H=1, C=6, N=7, O=8, Na=11, Mg=12, Al=13, Si=14, P=15, S=16, K=19,
         Ca=20, Ti=22, Mn=25, Fe=26, U=92, Th=90, Ar=18)

# Natural isotopes: symbol -> [(A, atomic_mass_amu, mole_abundance_fraction), ...]
ISO = {
    "H":  [(1, 1.007825, 0.999885), (2, 2.014102, 0.000115)],
    "C":  [(12, 12.000000, 0.9893), (13, 13.003355, 0.0107)],
    "N":  [(14, 14.003074, 0.99636), (15, 15.000109, 0.00364)],
    "O":  [(16, 15.994915, 0.99757), (17, 16.999132, 0.00038), (18, 17.999160, 0.00205)],
    "Na": [(23, 22.989769, 1.0)],
    "Mg": [(24, 23.985042, 0.7899), (25, 24.985837, 0.1000), (26, 25.982593, 0.1101)],
    "Al": [(27, 26.981539, 1.0)],
    "Si": [(28, 27.976927, 0.92223), (29, 28.976495, 0.04685), (30, 29.973770, 0.03092)],
    "P":  [(31, 30.973762, 1.0)],
    "S":  [(32, 31.972071, 0.9499), (33, 32.971459, 0.0075), (34, 33.967867, 0.0425), (36, 35.967081, 0.0001)],
    "K":  [(39, 38.963707, 0.932581), (40, 39.963998, 0.000117), (41, 40.961825, 0.067302)],
    "Ca": [(40, 39.962591, 0.96941), (42, 41.958618, 0.00647), (43, 42.958766, 0.00135),
           (44, 43.955481, 0.02086), (46, 45.953693, 0.00004), (48, 47.952534, 0.00187)],
    "Ti": [(46, 45.952628, 0.0825), (47, 46.951759, 0.0744), (48, 47.947942, 0.7372),
           (49, 48.947866, 0.0541), (50, 49.944787, 0.0518)],
    "Mn": [(55, 54.938044, 1.0)],
    "Fe": [(54, 53.939609, 0.05845), (56, 55.934936, 0.91754), (57, 56.935393, 0.02119), (58, 57.933274, 0.00282)],
    "U":  [(234, 234.040952, 0.000054), (235, 235.043930, 0.007204), (238, 238.050788, 0.992742)],
    "Th": [(232, 232.038056, 1.0)],
    "Ar": [(36, 35.967545, 0.003336), (38, 37.962732, 0.000629), (40, 39.962383, 0.996035)],
}

# Oxide -> {element: count}
OXIDE = dict(SiO2={"Si":1,"O":2}, Al2O3={"Al":2,"O":3}, FeO={"Fe":1,"O":1},
            Fe2O3={"Fe":2,"O":3}, MgO={"Mg":1,"O":1}, CaO={"Ca":1,"O":1},
            Na2O={"Na":2,"O":1}, K2O={"K":2,"O":1}, TiO2={"Ti":1,"O":2},
            MnO={"Mn":1,"O":1}, P2O5={"P":2,"O":5}, CO2={"C":1,"O":2}, H2O={"H":2,"O":1})


def pdg(zz: int, aa: int) -> int:
    return 1000000000 + zz * 10000 + aa * 10


def _k2o_from_pct_K(pct_K: float) -> float:
    """wt% K2O that contains pct_K wt% elemental K."""
    return pct_K * (2 * AW["K"] + AW["O"]) / (2 * AW["K"])


# ---------------------------------------------------------------------------
# Rock units. Composition in wt% of components (oxides or free elements).
# trace: ppm by mass for U/Th. K given as elemental wt% -> converted to K2O.
# ---------------------------------------------------------------------------
def _poorman():
    K2O = _k2o_from_pct_K(2.12)
    comp = dict(SiO2=58.0, Al2O3=17.0, FeO=7.0, MgO=3.0, CaO=3.5, Na2O=1.5,
                K2O=K2O, TiO2=0.8, MnO=0.1, P2O5=0.2, CO2=4.5, H2O=1.0, C=0.3)
    return dict(name="Poorman_phyllite", density=2.80, comp=comp,
                trace_ppm=dict(U=2.58, Th=10.48),
                desc="Poorman Fm: carbonate/sericite/graphitic phyllite (bulk overburden)")


def _yates():
    K2O = _k2o_from_pct_K(0.96)
    comp = dict(SiO2=49.0, Al2O3=14.0, FeO=12.0, MgO=7.0, CaO=10.0, Na2O=2.5,
                K2O=K2O, TiO2=1.5, MnO=0.2, P2O5=0.2)
    return dict(name="Yates_amphibolite", density=2.95, comp=comp,
                trace_ppm=dict(U=0.22, Th=0.33),
                desc="Yates amphibolite: mafic metavolcanic (dense, radio-clean, hosts 4850L)")


def _rhyolite():
    K2O = _k2o_from_pct_K(4.17)
    comp = dict(SiO2=72.0, Al2O3=13.0, FeO=2.0, MgO=0.4, CaO=1.2, Na2O=3.2,
                K2O=K2O, TiO2=0.3)
    return dict(name="Rhyolite_dike", density=2.62, comp=comp,
                trace_ppm=dict(U=8.75, Th=10.86),
                desc="Tertiary rhyolite dike (felsic, high U/Th/K)")


def _dusel_avg():
    # DUNE-official average rock (duneggd), already element/oxide mix
    comp = dict(SiO2=52.67, FeO=11.74, Al2O3=10.25, MgO=4.73, CO2=4.22,
                CaO=3.82, C=2.40, S=1.86, Na2O=0.53, P2O5=0.07, O=7.71)
    return dict(name="DUSEL_Rock", density=2.82, comp=comp, trace_ppm={},
                desc="DUNE-official average cavern rock (fallback / generic)")


def _air():
    # dry air element mass fractions
    comp = dict(N=75.52, O=23.14, Ar=1.288, C=0.012)
    return dict(name="AIR", density=1.205e-3, comp=comp, trace_ppm={},
                desc="dry air")


def _lar():
    return dict(name="LAr", density=1.396, comp=dict(Ar=100.0), trace_ppm={},
                desc="liquid argon (detector placeholder)")


# --- PREM background materials (for the Earth shells the local model embeds in) ---
def _crust():   # generic upper/inner crust
    return dict(name="ROCK", density=2.65, comp=dict(SiO2=66.0, Al2O3=15.0, FeO=5.0,
                CaO=4.0, Na2O=3.5, K2O=3.0, MgO=2.5), trace_ppm={}, desc="generic crust (PREM)")


def _mantle():
    return dict(name="MANTLE", density=3.3, comp=dict(SiO2=45.0, MgO=38.0, FeO=8.0,
                Al2O3=4.5, CaO=3.5), trace_ppm={}, desc="pyrolite mantle (PREM)")


def _innercore():
    return dict(name="INNERCORE", density=13.0, comp=dict(Fe=100.0), trace_ppm={}, desc="Fe inner core")


def _outercore():
    return dict(name="OUTERCORE", density=11.0, comp=dict(Fe=100.0), trace_ppm={}, desc="Fe outer core")


UNITS = [_air(), _lar(), _dusel_avg(), _poorman(), _yates(), _rhyolite(),
         _crust(), _mantle(), _innercore(), _outercore()]


# ---------------------------------------------------------------------------
# Composition -> element mass fractions -> isotope (PDG) mass fractions
# ---------------------------------------------------------------------------
def element_mass_fractions(unit: dict) -> dict:
    em: dict[str, float] = {}
    for comp, wt in unit["comp"].items():
        if comp in OXIDE:
            fw = sum(n * AW[el] for el, n in OXIDE[comp].items())
            for el, n in OXIDE[comp].items():
                em[el] = em.get(el, 0.0) + wt * n * AW[el] / fw
        else:  # free element symbol
            em[comp] = em.get(comp, 0.0) + wt
    for el, ppm in unit.get("trace_ppm", {}).items():
        em[el] = em.get(el, 0.0) + ppm * 1e-4  # ppm -> wt%
    tot = sum(em.values())
    return {el: w / tot for el, w in em.items()}


def isotope_mass_fractions(unit: dict) -> dict:
    """PDG code -> mass fraction (sums to 1)."""
    em = element_mass_fractions(unit)
    out: dict[int, float] = {}
    for el, w in em.items():
        isos = ISO[el]
        denom = sum(ab * m for _, m, ab in isos)
        for A, m, ab in isos:
            out[pdg(Z[el], A)] = out.get(pdg(Z[el], A), 0.0) + w * ab * m / denom
    tot = sum(out.values())
    return {k: v / tot for k, v in out.items()}


def zoverA(unit: dict) -> tuple[float, float, float]:
    """Return (<Z>_mass, <A>_mass, Z/A) for sanity checks."""
    em = element_mass_fractions(unit)
    ZA = sum(w * Z[el] / AW[el] for el, w in em.items())          # electron fraction
    meanZ = sum(w * Z[el] for el, w in em.items())
    meanA = sum(w * AW[el] for el, w in em.items())
    return meanZ, meanA, ZA


def write_materials_dat(path: str) -> None:
    lines = ["# DUNE FD Homestake 3D overburden materials",
             "# PDG isotope -> mass fraction. See homestake_geology.py for provenance.",
             ""]
    for u in UNITS:
        iso = isotope_mass_fractions(u)
        iso = {k: v for k, v in iso.items() if v > 1e-12}
        lines.append(f"{u['name']} {len(iso)}")
        for code in sorted(iso):
            lines.append(f"{code} {iso[code]:.10f}")
        lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    import os
    out = os.path.join(os.path.dirname(__file__), "DUNEFD_homestake_materials.dat")
    write_materials_dat(out)
    print("wrote", out)
    print(f"\n{'unit':22s} {'rho':>5s} {'<Z>':>6s} {'<A>':>6s} {'Z/A':>6s}  niso")
    for u in UNITS:
        mz, ma, za = zoverA(u)
        n = len([v for v in isotope_mass_fractions(u).values() if v > 1e-12])
        print(f"{u['name']:22s} {u['density']:5.2f} {mz:6.2f} {ma:6.2f} {za:6.4f}  {n}")
