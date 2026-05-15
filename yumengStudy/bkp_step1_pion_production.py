# BdNMC way to do Step 1: turn the CCM-200 p+W -> pi0 histogram into a pion list SIREN can read.
#
# input:  ../../DM_calculations/Experiments/CCM200/CCM_720pi0.txt
#         (T_pi [MeV], theta [rad], counts per 1e5 POT)  -- produced upstream by Geant4/BdNMC.
#         or: do everything in g4.
# output: pi0_list_CCM200.csv with columns E,m,x0,y0,z0,px,py,pz (GeV, m, pion decay vertex)
#         then step2(to dark photon) siren: to be input in siren.distributions.PrimaryExternalDistribution
#         (because my rare decay case & pick up one random pion from the list, so add pion decay vertex module
#         (pk currently module only starts where pion begins))

from __future__ import annotations
import csv
import numpy as np


# CCM-v2 world frame (m). See resources/detectors/CCM/CCM-v2/densities.dat.
LOWER_W_TARGET = np.array([0.0,  0.0, -0.241])
DETECTOR_ORIG  = np.array([23.0, 0.0, -0.65])

# proton beam direction at the target. CCM_720pi0.txt's theta is measured from
# this axis. detector sits at +x so +x as forward.
BEAM_AXIS = np.array([1.0, 0.0, 0.0])

M_PI0_GEV = 0.1349768   # PDG

DEFAULT_PI0_HISTOGRAM = (
    "/work/submit/yumeng1/workspaces/CCM/sources/DM_calculations/"
    "Experiments/CCM200/CCM_720pi0.txt"
)

# bin widths read off the file
DELTA_T_MEV = 1.0
DELTA_THETA = 0.01


def load_histogram(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # parse comma-separated (T_MeV, theta_rad, counts), skip blanks and '#' lines
    T, th, w = [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = [p.strip() for p in s.split(",")]
            if len(parts) != 3:
                continue
            T.append(float(parts[0]))
            th.append(float(parts[1]))
            w.append(float(parts[2]))
    return np.asarray(T), np.asarray(th), np.asarray(w)


def build(filename: str = "pi0_list_CCM200.csv",
          n: int = 100_000,
          seed: int = 0,
          histogram_path: str = DEFAULT_PI0_HISTOGRAM) -> None:
    """Resample the (T, theta) histogram into a flat pi0 list."""

    rng = np.random.default_rng(seed)

    # 1) load the upstream histogram
    T_bin, th_bin, counts = load_histogram(histogram_path)
    print(f"loaded {len(T_bin)} (T,theta) bins from {histogram_path}")
    print(f"  sum(counts) = {counts.sum():.0f}  per 1e5 POT")

    # 2) make a discrete pdf over bins, weighted by counts
    probs = counts / counts.sum()

    # 3) draw n bin indices with replacement
    idx = rng.choice(len(T_bin), size=n, p=probs)

    # 4) jitter inside each bin so events don't pile up on the bin centres
    T_MeV = T_bin[idx]  + rng.uniform(-0.5, 0.5, size=n) * DELTA_T_MEV
    T_MeV = np.clip(T_MeV, 0.0, None)
    theta = th_bin[idx] + rng.uniform(-0.5, 0.5, size=n) * DELTA_THETA
    theta = np.clip(theta, 0.0, np.pi)

    # 5) the histogram is integrated over phi (axisymmetric around the beam),
    #    so draw phi uniformly
    phi = rng.uniform(0.0, 2 * np.pi, size=n)

    # 6) (T, theta, phi) -> 4-momentum
    E_GeV = T_MeV * 1e-3 + M_PI0_GEV
    p_mag = np.sqrt(np.maximum(E_GeV**2 - M_PI0_GEV**2, 0.0))

    # rotate from "beam frame" into world: beam axis is z0, then build any
    # two perpendicular unit vectors x0, y0
    z0 = BEAM_AXIS / np.linalg.norm(BEAM_AXIS)
    helper = np.array([0.0, 1.0, 0.0]) if abs(z0[1]) < 0.9 else np.array([0.0, 0.0, 1.0])
    x0 = np.cross(z0, helper); x0 /= np.linalg.norm(x0)
    y0 = np.cross(z0, x0)

    sn, cs = np.sin(theta), np.cos(theta)
    dirs = (cs[:, None] * z0
            + sn[:, None] * np.cos(phi)[:, None] * x0
            + sn[:, None] * np.sin(phi)[:, None] * y0)
    p_vec = p_mag[:, None] * dirs

    # 7) vertex: pick a point uniformly inside the lower W target cylinder
    #    (Mark-III: R = 5 cm, L = 30 cm, centred at LOWER_W_TARGET)
    R, L = 0.05, 0.298
    r   = R * np.sqrt(rng.uniform(0, 1, size=n))   # sqrt for uniform-in-area
    ang = rng.uniform(0, 2 * np.pi, size=n)
    dz  = rng.uniform(-L / 2, L / 2, size=n)
    pos = LOWER_W_TARGET[None, :] + np.column_stack(
        [r * np.cos(ang), r * np.sin(ang), dz]
    )

    # 8) dump to csv in the column order PrimaryExternalDistribution expects
    with open(filename, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E", "m", "x0", "y0", "z0", "px", "py", "pz"])
        for i in range(n):
            w.writerow([f"{E_GeV[i]:.8f}", f"{M_PI0_GEV:.8f}",
                        f"{pos[i, 0]:.6f}", f"{pos[i, 1]:.6f}", f"{pos[i, 2]:.6f}",
                        f"{p_vec[i, 0]:.8f}", f"{p_vec[i, 1]:.8f}", f"{p_vec[i, 2]:.8f}"])

    print(f"wrote {n} pi0 records -> {filename}")
    print(f"  beam axis = {BEAM_AXIS}, target = {LOWER_W_TARGET}, detector = {DETECTOR_ORIG}")


if __name__ == "__main__":
    build()
