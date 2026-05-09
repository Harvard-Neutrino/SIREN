#!/usr/bin/env python3
"""
Smoke test for slow-rescaling (xi, y) sampling in QuarkDISFromSpline.
Exercises 100 events and verifies kinematic bounds for each.
"""
import sys
import traceback

# ---------------------------------------------------------------------------
# Constants (mirror C++ values exactly)
# ---------------------------------------------------------------------------
M_C    = 1.27       # Constants::CharmMass
M_D0   = 1.86484    # Constants::D0Mass
M_N    = (0.938272 + 0.939565) / 2   # isoscalar nucleon mass
M_MU   = 0.105658   # muon mass
Q2MIN  = 1.0        # GeV^2
N_EVENTS = 100
E_NU   = 100.0      # neutrino energy in GeV

# Spline files
SPLINE_DIR = (
    "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pzhelnin/"
    "DiMuons/Simulation/Resources/Splines/Maboi_M_Muon_SR"
)
DIFF_FILE  = SPLINE_DIR + "/dsdxidy_nu-N-cc-charm-CT14nlo_central.fits"
TOTAL_FILE = SPLINE_DIR + "/sigma_nu-N-cc-charm-CT14nlo_central.fits"

# ---------------------------------------------------------------------------
# Import SIREN
# ---------------------------------------------------------------------------
try:
    import siren
    import siren.interactions
    import siren.dataclasses
    import siren.utilities
except Exception as exc:
    print(f"IMPORT ERROR: {exc}", file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

PT = siren.dataclasses.Particle.ParticleType

# ---------------------------------------------------------------------------
# Construct cross-section object
# ---------------------------------------------------------------------------
try:
    xs = siren.interactions.QuarkDISFromSpline(
        DIFF_FILE,
        TOTAL_FILE,
        int(1),                                   # interaction_type: CC
        M_N,                                      # isoscalar mass
        int(1),                                   # minimum Q2
        [PT.NuMu],                                # primary types
        [PT.O16Nucleus],                          # target types
        "m",                                      # units
    )
except Exception as exc:
    print(f"CONSTRUCTOR ERROR: {exc}", file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Get a valid signature
# ---------------------------------------------------------------------------
try:
    sigs = list(xs.GetPossibleSignatures())
    assert len(sigs) > 0, "QuarkDISFromSpline returned no signatures"
    sig = sigs[0]
except Exception as exc:
    print(f"SIGNATURE ERROR: {exc}", file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# RNG
# ---------------------------------------------------------------------------
rng = siren.utilities.SIREN_random(1234)

# ---------------------------------------------------------------------------
# Expected kinematic bounds (mirror C++ Step 5)
# ---------------------------------------------------------------------------
y_max   = 1.0 - M_MU / E_NU
W2_thr  = (M_N + M_D0) ** 2
y_min   = (W2_thr - M_N**2 + Q2MIN) / (2.0 * M_N * E_NU)
xi_min  = (M_C**2 + Q2MIN) / (2.0 * M_N * E_NU * y_max)

print(f"Expected bounds:")
print(f"  y_min  = {y_min:.6f}")
print(f"  y_max  = {y_max:.6f}")
print(f"  xi_min = {xi_min:.6f}")
print(f"  W2_thr = {W2_thr:.6f}")
print(f"  Q2MIN  = {Q2MIN:.6f}")
print()

# ---------------------------------------------------------------------------
# Run N_EVENTS events
# ---------------------------------------------------------------------------
failures = []

for event_idx in range(N_EVENTS):
    try:
        ir = siren.dataclasses.InteractionRecord()
        ir.signature       = sig
        ir.primary_momentum = [E_NU, 0.0, 0.0, E_NU]   # massless nu along +z
        ir.primary_mass    = 0.0
        ir.target_mass     = M_N

        cdr = siren.dataclasses.CrossSectionDistributionRecord(ir)
        xs.SampleFinalState(cdr, rng)

        params = dict(cdr.interaction_parameters)
        xi = params["bjorken_xi"]
        y  = params["bjorken_y"]
        x  = params["bjorken_x"]

        # Derived quantities
        Q2 = 2.0 * M_N * E_NU * y * xi - M_C**2
        W2 = M_N**2 + 2.0 * M_N * E_NU * y * (1.0 - xi) + M_C**2

        # Check all assertions
        event_failures = []
        if not (xi_min <= xi <= 1.0):
            event_failures.append(
                f"xi={xi:.6g} not in [{xi_min:.6g}, 1.0]"
            )
        if not (y_min <= y <= y_max):
            event_failures.append(
                f"y={y:.6g} not in [{y_min:.6g}, {y_max:.6g}]"
            )
        if not (x > 0.0):
            event_failures.append(f"x={x:.6g} not > 0")
        if not (Q2 >= Q2MIN):
            event_failures.append(
                f"Q2={Q2:.6g} < Q2MIN={Q2MIN:.6g}"
            )
        if not (W2 > W2_thr):
            event_failures.append(
                f"W2={W2:.6g} not > W2_thr={W2_thr:.6g}"
            )

        if event_failures:
            msg = (
                f"Event {event_idx}: xi={xi:.6g} y={y:.6g} "
                f"x={x:.6g} Q2={Q2:.6g} W2={W2:.6g} | "
                + "; ".join(event_failures)
            )
            failures.append(msg)
            print(f"FAIL: {msg}")

    except Exception as exc:
        msg = f"Event {event_idx}: unexpected exception: {exc}"
        failures.append(msg)
        print(f"FAIL: {msg}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)

    if len(failures) >= 10:
        print(f"Stopping early after {len(failures)} failures.")
        break

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
if failures:
    n_fail = len(failures)
    n_ok   = N_EVENTS - n_fail
    print(f"\nFAIL {n_ok}/{N_EVENTS} ({n_fail} events failed bounds checks)")
    sys.exit(1)
else:
    print(f"OK {N_EVENTS}/{N_EVENTS}")
    sys.exit(0)
