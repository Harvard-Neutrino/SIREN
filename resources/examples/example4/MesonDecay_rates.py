"""
Meson three-body decay rate calculation.

Computes and plots the branching ratio enhancement for
pi/K -> l nu phi relative to the SM two-body decay,
as a function of mediator mass, for scalar and pseudoscalar
couplings.

Reference: Dutta et al., PRL 129, 111803 (2022), Fig. 4
           Carlson & Rislow, Phys.Rev.D 86, 035013 (2012)
"""

from pathlib import Path

import numpy as np
import siren

_resources_dir = Path(__file__).resolve().parents[2]
_meson_module = siren._util.load_module(
    "siren.resources.processes.DarkNewsTables.MesonProduction",
    str(_resources_dir / "processes" / "DarkNewsTables" / "MesonProduction.py"),
)
MesonThreeBodyDecay = _meson_module.MesonThreeBodyDecay

# Physical constants
M_PI = 0.13957    # GeV
M_K  = 0.49368    # GeV
M_MU = 0.10566    # GeV
M_E  = 0.000511   # GeV
GF   = 1.16638e-5 # GeV^-2
FPI  = 0.1307     # GeV
FK   = 0.1598     # GeV
VUD  = 0.9737
VUS  = 0.2245

g_mu = 1.0  # reference coupling

def sm_width_2body(m_M, m_l, f_M, V_Mq):
    r = (m_l / m_M)**2
    return GF**2 * f_M**2 * V_Mq**2 * m_M * m_l**2 * (1.0 - r)**2 / (8.0 * np.pi)


# Scan mediator mass
m_phi_vals = np.logspace(-3, np.log10(0.35), 100)  # 1 MeV to 350 MeV

configs = [
    ("pi, scalar",       M_PI, M_MU, FPI, VUD, "scalar"),
    ("pi, pseudoscalar", M_PI, M_MU, FPI, VUD, "pseudoscalar"),
    ("K, scalar",        M_K,  M_MU, FK,  VUS, "scalar"),
    ("K, pseudoscalar",  M_K,  M_MU, FK,  VUS, "pseudoscalar"),
]

print(f"{'m_phi [MeV]':>12s}", end="")
for label, *_ in configs:
    print(f"  {label:>20s}", end="")
print()

for m_phi in m_phi_vals:
    line = f"{m_phi*1e3:12.3f}"
    for label, m_M, m_l, f_M, V_Mq, mtype in configs:
        if m_phi >= m_M - m_l:
            line += f"  {'---':>20s}"
            continue
        try:
            decay = MesonThreeBodyDecay(m_M, m_l, m_phi, g_mu, mtype)
            w3 = decay.total_width()
            w2 = sm_width_2body(m_M, m_l, f_M, V_Mq)
            ratio = w3 / w2 if w2 > 0 else 0.0
            line += f"  {ratio:20.6e}"
        except ValueError:
            line += f"  {'forbidden':>20s}"
    print(line)
