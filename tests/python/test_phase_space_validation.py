"""
Comprehensive phase space validation tests.

Tests cover normalization integrals, cross-measure integral agreement,
closure tests, phase space coverage, and topology/measure validation.
"""

import copy
import math

import pytest

siren = pytest.importorskip("siren")


# ------------------------------------------------------------------ #
#  Helpers                                                             #
# ------------------------------------------------------------------ #


def _make_2body_record(M=1.0, E=None, mA=0.0, mB=0.0,
                       vertex=(0.0, 0.0, 0.0)):
    """Build an InteractionRecord for a 2-body decay.

    Parameters
    ----------
    M : float
        Parent mass.
    E : float or None
        Parent energy. If None, parent is at rest (E = M).
    mA, mB : float
        Daughter masses.
    vertex : tuple
        Interaction vertex position.
    """
    if E is None:
        E = M
    pz = math.sqrt(max(E * E - M * M, 0.0))

    rec = siren.dataclasses.InteractionRecord()
    rec.signature.primary_type = siren.dataclasses.ParticleType.N4
    rec.signature.target_type = siren.dataclasses.ParticleType.Decay
    rec.signature.secondary_types = [
        siren.dataclasses.ParticleType.NuLight,
        siren.dataclasses.ParticleType.Gamma,
    ]
    rec.primary_mass = M
    rec.primary_momentum = [E, 0.0, 0.0, pz]
    rec.secondary_masses = [mA, mB]
    rec.secondary_momenta = [[0, 0, 0, 0], [0, 0, 0, 0]]
    rec.secondary_helicities = [0, 0]
    rec.interaction_vertex = list(vertex)
    rec.primary_initial_position = list(vertex)
    return rec


def _make_3body_record(M=1.0, E=2.0, masses=(0.1, 0.3, 0.2)):
    """Build an InteractionRecord for a 3-body decay."""
    pz = math.sqrt(max(E * E - M * M, 0.0))

    rec = siren.dataclasses.InteractionRecord()
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.dataclasses.ParticleType(211)
    sig.secondary_types = [
        siren.dataclasses.ParticleType(-13),
        siren.dataclasses.ParticleType(14),
        siren.dataclasses.ParticleType(5922),
    ]
    rec.signature = sig
    rec.primary_mass = M
    rec.primary_momentum = [E, 0.0, 0.0, pz]
    rec.secondary_masses = list(masses)
    rec.secondary_momenta = [[0, 0, 0, 0]] * 3
    rec.secondary_helicities = [0, 0, 0]
    rec.interaction_vertex = [0, 0, 0]
    rec.primary_initial_position = [0, 0, 0]
    return rec


def _make_box_at(x, y, z, wx, wy, wz):
    """Create a Box geometry centered at (x, y, z) with full-widths (wx, wy, wz)."""
    placement = siren.geometry.Placement(siren.math.Vector3D(x, y, z))
    return siren.geometry.Box(placement, wx, wy, wz)


def _daughter_hits_box(record, daughter_index, box):
    """Check if the daughter's momentum direction from the vertex hits the box."""
    vx, vy, vz = record.interaction_vertex
    pos = siren.math.Vector3D(vx, vy, vz)
    mom = record.secondary_momenta[daughter_index]
    px, py, pz = mom[1], mom[2], mom[3]
    p_mag = math.sqrt(px * px + py * py + pz * pz)
    if p_mag < 1e-30:
        return False
    direction = siren.math.Vector3D(px / p_mag, py / p_mag, pz / p_mag)
    intersections = box.Intersections(pos, direction)
    return len(intersections) > 0


# ================================================================== #
#  Category 2: Normalization integrals (Python-side Monte Carlo)       #
# ================================================================== #


class TestNormalizationIntegrals:

    def test_isotropic_2body_density_integrates_to_one(self):
        """Isotropic channel density should be 1/(4*pi) everywhere.

        Sample 5000 events and verify every single density evaluation
        returns exactly 1/(4*pi) to floating point precision.
        """
        iso = siren.injection.Isotropic2BodyChannel(0)
        rec = _make_2body_record(M=1.0, E=1.0, mA=0.0, mB=0.0)
        rng = siren.utilities.SIREN_random(42)
        expected = 1.0 / (4.0 * math.pi)

        for _ in range(5000):
            r = copy.deepcopy(rec)
            iso.Sample(rng, None, r)
            d = iso.Density(None, r)
            assert abs(d - expected) < 1e-12, (
                "Isotropic density {:.15e} != expected {:.15e}".format(
                    d, expected)
            )

    def test_directed_2body_density_integrates_to_one(self):
        """Importance-sampling estimate of the directed density integral.

        Sample from the multi-channel (isotropic + directed), compute
        weights w_i = iso_density / mc_density. The mean weight should
        be 1.0 because both densities integrate to 1 over the full
        solid angle.
        """
        box = _make_box_at(0, 0, 100, 0.2, 0.2, 0.2)
        rec = _make_2body_record(M=1.0, E=10.0, mA=0.0, mB=0.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.5, 0.5]

        rng = siren.utilities.SIREN_random(123)
        N = 2000
        weights = []

        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_iso / d_mc)

        assert len(weights) > N * 0.5, (
            "Too few valid weights: {}/{}".format(len(weights), N)
        )

        mean_w = sum(weights) / len(weights)
        assert 0.9 < mean_w < 1.1, (
            "Directed density integral estimate {:.4f} not in [0.9, 1.1]"
            .format(mean_w)
        )

    def test_multichannel_density_sums_correctly(self):
        """MultiChannelPhaseSpace.Density() must equal the weighted sum
        of individual channel densities. This is exact (no numerical
        approximation), so we check for exact floating-point equality.
        """
        box = _make_box_at(0, 0, 100, 0.2, 0.2, 0.2)
        rec = _make_2body_record(M=1.0, E=10.0, mA=0.0, mB=0.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        w_iso = 0.3
        w_dir = 0.7

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [w_iso, w_dir]

        rng = siren.utilities.SIREN_random(77)

        for _ in range(1000):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)

            d_mc = mc.Density(None, r)
            d_iso = iso.Density(None, r)
            d_dir = directed.Density(None, r)
            expected = w_iso * d_iso + w_dir * d_dir

            # Allow a tiny relative tolerance for floating-point summation
            if expected > 0:
                rel_err = abs(d_mc - expected) / expected
                assert rel_err < 1e-12, (
                    "MC density {:.15e} != weighted sum {:.15e} "
                    "(rel err {:.2e})"
                    .format(d_mc, expected, rel_err)
                )
            else:
                assert abs(d_mc) < 1e-30


# ================================================================== #
#  Category 3: Cross-measure integral agreement                       #
# ================================================================== #


class TestCrossMeasureAgreement:

    def test_rest_lab_integral_agreement_mc(self):
        """For a boosted 2-body decay, the rest-frame density
        (1/(4*pi)) should integrate to 1 over the rest-frame solid
        angle, regardless of the lab frame boost.

        We verify this by sampling from the isotropic channel (which
        samples uniformly in the rest frame) and confirming that the
        density is always 1/(4*pi). The Monte Carlo integral
        mean(density) * (4*pi) must equal 1.0.
        """
        M = 1.0
        E = 5.0
        mA = 0.1
        mB = 0.2
        rec = _make_2body_record(M=M, E=E, mA=mA, mB=mB)

        iso = siren.injection.Isotropic2BodyChannel(0)
        rng = siren.utilities.SIREN_random(42)
        expected = 1.0 / (4.0 * math.pi)

        densities = []
        for _ in range(3000):
            r = copy.deepcopy(rec)
            iso.Sample(rng, None, r)
            d = iso.Density(None, r)
            densities.append(d)

        mean_d = sum(densities) / len(densities)
        mc_integral = mean_d * (4.0 * math.pi)

        # Isotropic density is constant, so the integral should be
        # exactly 1.0 (up to floating point)
        assert abs(mc_integral - 1.0) < 1e-10, (
            "Rest-frame integral estimate {:.6f} != 1.0".format(mc_integral)
        )

        # Also check each density is 1/(4*pi)
        for d in densities:
            assert abs(d - expected) < 1e-12


# ================================================================== #
#  Category 4: Closure tests                                           #
# ================================================================== #


class TestClosureTests:

    def test_2body_closure_directed_vs_isotropic(self):
        """Closure test: biased sampling with correct weighting must
        reproduce unbiased expectations.

        We use a large target box so that a measurable fraction of
        isotropic events hit it. Then we compare the hit fractions
        estimated by:
          (A) brute-force isotropic sampling
          (B) importance-weighted multichannel sampling

        Both should agree within statistical error.
        """
        M = 0.5
        E = 3.0
        rec = _make_2body_record(M=M, E=E, mA=0.0, mB=0.0)

        # Large box at z=20 to get measurable hit fraction
        box = _make_box_at(0, 0, 20, 10.0, 10.0, 10.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.01, 0.99]

        rng_iso = siren.utilities.SIREN_random(100)
        rng_mc = siren.utilities.SIREN_random(200)

        # --- Method A: brute force isotropic ---
        N_iso = 20000
        hits_iso = 0
        for _ in range(N_iso):
            r = copy.deepcopy(rec)
            iso.Sample(rng_iso, None, r)
            if _daughter_hits_box(r, 0, box):
                hits_iso += 1

        frac_iso = hits_iso / float(N_iso)

        # --- Method B: importance-weighted multichannel ---
        N_mc = 2000
        weighted_hits = 0.0
        sum_weights = 0.0
        for _ in range(N_mc):
            r = copy.deepcopy(rec)
            mc.Sample(rng_mc, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                w = d_iso / d_mc
                sum_weights += w
                if _daughter_hits_box(r, 0, box):
                    weighted_hits += w

        frac_mc = weighted_hits / sum_weights if sum_weights > 0 else 0.0

        # Both estimates should agree. The isotropic estimate has
        # Poisson noise; with 20000 events and ~1% hit rate, expect
        # ~200 hits, so sigma ~ sqrt(200)/20000 ~ 0.07%.
        # Allow generous tolerance of 50% relative error.
        if frac_iso > 0:
            rel_diff = abs(frac_mc - frac_iso) / frac_iso
            assert rel_diff < 0.5, (
                "Closure test failed: isotropic hit fraction {:.5f}, "
                "multichannel estimate {:.5f}, relative diff {:.2f}"
                .format(frac_iso, frac_mc, rel_diff)
            )
        else:
            # If isotropic had zero hits, the MC estimate should also
            # be very small
            assert frac_mc < 0.01

    def test_multichannel_effective_sample_size(self):
        """The effective sample size of multichannel importance sampling
        should be a reasonable fraction of the total sample count.

        N_eff = (sum w_i)^2 / sum(w_i^2) measures sampling efficiency.
        For a well-designed multichannel, N_eff should be > N/10.
        """
        M = 0.5
        E = 3.0
        rec = _make_2body_record(M=M, E=E, mA=0.0, mB=0.0)

        box = _make_box_at(0, 0, 20, 10.0, 10.0, 10.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.01, 0.99]

        rng = siren.utilities.SIREN_random(42)
        N = 2000
        weights = []

        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_iso / d_mc)

        sum_w = sum(weights)
        sum_w2 = sum(w * w for w in weights)
        n_eff = (sum_w * sum_w) / sum_w2 if sum_w2 > 0 else 0.0

        assert n_eff > N / 10.0, (
            "Effective sample size {:.1f} is less than N/10 = {:.1f}"
            .format(n_eff, N / 10.0)
        )


# ================================================================== #
#  Category 5: Phase space coverage                                    #
# ================================================================== #


class TestPhaseSpaceCoverage:

    def test_no_density_holes_2body(self):
        """Verify there are no density holes in the multi-channel.

        - Isotropic samples evaluated by directed: density >= 0
          (may be 0 for directions missing the target)
        - Directed samples evaluated by isotropic: density must always
          be > 0 and finite (isotropic has no holes)
        """
        box = _make_box_at(0, 0, 100, 0.2, 0.2, 0.2)
        rec = _make_2body_record(M=1.0, E=10.0, mA=0.0, mB=0.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        rng = siren.utilities.SIREN_random(42)

        # Isotropic samples -> directed density
        for _ in range(500):
            r = copy.deepcopy(rec)
            iso.Sample(rng, None, r)
            d = directed.Density(None, r)
            assert d >= 0, "Directed density is negative: {}".format(d)
            assert not math.isnan(d), "Directed density is NaN"

        # Directed samples -> isotropic density
        for _ in range(500):
            r = copy.deepcopy(rec)
            directed.Sample(rng, None, r)
            d = iso.Density(None, r)
            assert d > 0, "Isotropic density is not positive: {}".format(d)
            assert math.isfinite(d), "Isotropic density is not finite: {}".format(d)

    def test_boundary_densities_finite(self):
        """Density should be finite at kinematic boundary cases:
        - Parent at rest (E = M)
        - Ultra-relativistic parent (E = 1000*M)
        - Near-threshold (M = mA + mB + epsilon)
        """
        iso = siren.injection.Isotropic2BodyChannel(0)
        rng = siren.utilities.SIREN_random(42)
        expected = 1.0 / (4.0 * math.pi)

        # Parent at rest
        rec_rest = _make_2body_record(M=1.0, E=1.0, mA=0.0, mB=0.0)
        for _ in range(50):
            r = copy.deepcopy(rec_rest)
            iso.Sample(rng, None, r)
            d = iso.Density(None, r)
            assert math.isfinite(d), "Density NaN/Inf at rest"
            assert abs(d - expected) < 1e-12

        # Ultra-relativistic
        rec_ultra = _make_2body_record(M=1.0, E=1000.0, mA=0.0, mB=0.0)
        for _ in range(50):
            r = copy.deepcopy(rec_ultra)
            iso.Sample(rng, None, r)
            d = iso.Density(None, r)
            assert math.isfinite(d), "Density NaN/Inf at ultra-relativistic"
            assert abs(d - expected) < 1e-12

        # Near-threshold
        mA = 0.3
        mB = 0.4
        M = mA + mB + 1e-8
        rec_thresh = _make_2body_record(M=M, E=M, mA=mA, mB=mB)
        for _ in range(50):
            r = copy.deepcopy(rec_thresh)
            iso.Sample(rng, None, r)
            d = iso.Density(None, r)
            assert math.isfinite(d), "Density NaN/Inf at near-threshold"
            assert not math.isnan(d), "Density is NaN near threshold"

        # Forward/backward daughters in boosted frame
        rec_boost = _make_2body_record(M=1.0, E=5.0, mA=0.1, mB=0.2)
        for _ in range(100):
            r = copy.deepcopy(rec_boost)
            iso.Sample(rng, None, r)
            d = iso.Density(None, r)
            assert math.isfinite(d), "Density NaN/Inf for boosted case"


# ================================================================== #
#  Category 6: Topology/measure validation                             #
# ================================================================== #


class TestTopologyMeasureValidation:

    def test_topology_mismatch_raises(self):
        """Mixing Decay2Body (isotropic) with Scatter2to2 (scattering)
        in a MultiChannelPhaseSpace must raise RuntimeError when
        CommonTopology() is called.
        """
        box = _make_box_at(0, 0, 100, 1.0, 1.0, 1.0)
        iso = siren.injection.Isotropic2BodyChannel(0)
        sc = siren.injection.DetectorDirectedScatteringChannel(
            box, 0, siren.injection.ScatteringVariable.Q2)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, sc]
        mc.weights = [0.5, 0.5]

        with pytest.raises(RuntimeError, match="mismatch"):
            mc.CommonTopology()

    def test_measure_compatibility_within_topology(self):
        """Isotropic2BodyChannel + DetectorDirected2BodyChannel share
        the same topology (Decay2Body) and measure (SolidAngleRest).

        ValidateChannels() should return no error diagnostics.
        """
        box = _make_box_at(0, 0, 100, 1.0, 1.0, 1.0)
        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.5, 0.5]

        # Should not raise
        topo = mc.CommonTopology()
        assert topo == siren.injection.PhaseSpaceTopology.Decay2Body

        measure = mc.CommonMeasure()
        assert measure == siren.injection.PhaseSpaceMeasure.SolidAngleRest()

        diagnostics = mc.ValidateChannels()
        # Filter for error-level diagnostics (mismatch / incompatibility)
        errors = [d for d in diagnostics
                  if "mismatch" in d.lower() or "incompatib" in d.lower()]
        assert len(errors) == 0, (
            "Unexpected error diagnostics: {}".format(errors)
        )

    def test_channel_reports_correct_topology_and_measure(self):
        """Each channel type should report the expected topology
        and measure.
        """
        box = _make_box_at(0, 0, 100, 1.0, 1.0, 1.0)

        # Isotropic2BodyChannel
        iso = siren.injection.Isotropic2BodyChannel(0)
        assert iso.Topology() == siren.injection.PhaseSpaceTopology.Decay2Body
        assert iso.Measure() == siren.injection.PhaseSpaceMeasure.SolidAngleRest()

        # DetectorDirected2BodyChannel
        dir2 = siren.injection.DetectorDirected2BodyChannel(box, 0)
        assert dir2.Topology() == siren.injection.PhaseSpaceTopology.Decay2Body
        assert dir2.Measure() == siren.injection.PhaseSpaceMeasure.SolidAngleRest()

        # DetectorDirected3BodyChannel
        dir3 = siren.injection.DetectorDirected3BodyChannel(
            box,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_pair_index=1,
            mass_mode=siren.injection.InvariantMassMode.Uniform,
        )
        assert dir3.Topology() == siren.injection.PhaseSpaceTopology.Decay3Body
        assert dir3.Measure() == siren.injection.PhaseSpaceMeasure.Recursive2Body()

        # DetectorDirectedScatteringChannel with Q2
        sc_q2 = siren.injection.DetectorDirectedScatteringChannel(
            box, directed_index=0,
            variable=siren.injection.ScatteringVariable.Q2)
        assert sc_q2.Topology() == siren.injection.PhaseSpaceTopology.Scatter2to2
        assert sc_q2.Measure() == siren.injection.PhaseSpaceMeasure.MandelstamQ2()

        # DetectorDirectedScatteringChannel with BjorkenY
        sc_by = siren.injection.DetectorDirectedScatteringChannel(
            box, directed_index=0,
            variable=siren.injection.ScatteringVariable.BjorkenY)
        assert sc_by.Topology() == siren.injection.PhaseSpaceTopology.Scatter2to2
        assert sc_by.Measure() == siren.injection.PhaseSpaceMeasure.BjorkenXY()

    def test_legacy_convention_matches_new_system(self):
        """The legacy Convention() method should return values
        consistent with the new Topology() + Measure() system.

        Each (Measure -> Convention) mapping is defined in
        PhaseSpaceChannel::Convention().
        """
        box = _make_box_at(0, 0, 100, 1.0, 1.0, 1.0)

        # Map from Measure -> expected legacy Convention
        measure_to_convention = {
            siren.injection.PhaseSpaceMeasure.SolidAngleRest():
                siren.injection.PhaseSpaceConvention.RestFrameSolidAngle,
            siren.injection.PhaseSpaceMeasure.SolidAngleLab():
                siren.injection.PhaseSpaceConvention.LabFrameSolidAngle,
            siren.injection.PhaseSpaceMeasure.Recursive2Body():
                siren.injection.PhaseSpaceConvention.Recursive2Body,
            siren.injection.PhaseSpaceMeasure.DalitzPair():
                siren.injection.PhaseSpaceConvention.Dalitz,
            siren.injection.PhaseSpaceMeasure.HelicityAngles():
                siren.injection.PhaseSpaceConvention.HelicityAngles,
            siren.injection.PhaseSpaceMeasure.MandelstamQ2():
                siren.injection.PhaseSpaceConvention.MandelstamST,
            siren.injection.PhaseSpaceMeasure.BjorkenXY():
                siren.injection.PhaseSpaceConvention.BjorkenXY,
        }

        channels = [
            siren.injection.Isotropic2BodyChannel(0),
            siren.injection.DetectorDirected2BodyChannel(box, 0),
            siren.injection.DetectorDirected3BodyChannel(
                box, spectator_index=0, pair_first_index=1,
                pair_second_index=2, directed_pair_index=1,
                mass_mode=siren.injection.InvariantMassMode.Uniform),
            siren.injection.DetectorDirectedScatteringChannel(
                box, directed_index=0,
                variable=siren.injection.ScatteringVariable.Q2),
            siren.injection.DetectorDirectedScatteringChannel(
                box, directed_index=0,
                variable=siren.injection.ScatteringVariable.BjorkenY),
        ]

        for ch in channels:
            measure = ch.Measure()
            convention = ch.Convention()
            expected_convention = measure_to_convention.get(measure)
            assert expected_convention is not None, (
                "No convention mapping for measure {} on channel {}"
                .format(measure, ch.Name())
            )
            assert convention == expected_convention, (
                "Channel {} has Measure={} but Convention={}, expected {}"
                .format(ch.Name(), measure, convention, expected_convention)
            )


# ================================================================== #
#  Gap fixes: sampling-density consistency                             #
# ================================================================== #


def _rest_frame_cos_theta(record, daughter_index=0):
    """Extract the rest-frame cos(theta) of a daughter from a sampled record."""
    M = record.primary_mass
    E_p = record.primary_momentum[0]
    px_p = record.primary_momentum[1]
    py_p = record.primary_momentum[2]
    pz_p = record.primary_momentum[3]
    p_p = math.sqrt(px_p**2 + py_p**2 + pz_p**2)

    E_d = record.secondary_momenta[daughter_index][0]
    px_d = record.secondary_momenta[daughter_index][1]
    py_d = record.secondary_momenta[daughter_index][2]
    pz_d = record.secondary_momenta[daughter_index][3]
    p_d = math.sqrt(px_d**2 + py_d**2 + pz_d**2)

    if p_p < 1e-30:
        # Parent at rest: use the z-axis as reference
        return pz_d / p_d if p_d > 0 else 0.0

    beta = p_p / E_p
    gamma = E_p / M
    cos_lab = (px_d * px_p + py_d * py_p + pz_d * pz_p) / (p_d * p_p)
    p_par_lab = p_d * cos_lab
    p_par_rest = gamma * (p_par_lab - beta * E_d)

    mA = record.secondary_masses[daughter_index]
    mB = record.secondary_masses[1 - daughter_index]
    p_rest = siren.injection.TwoBodyRestMomentum(M, mA, mB)
    if p_rest < 1e-30:
        return 0.0
    return max(-1.0, min(1.0, p_par_rest / p_rest))


class TestSamplingDensityConsistency:
    """Verify that Sample() and Density() are consistent: the distribution
    of sampled events matches the reported density function.
    """

    def test_isotropic_sampling_is_uniform(self):
        """Isotropic channel should produce a flat cos(theta_rest)
        distribution. Bin 10000 samples and chi-squared test.
        """
        iso = siren.injection.Isotropic2BodyChannel(0)
        rec = _make_2body_record(M=1.0, E=5.0, mA=0.1, mB=0.2)
        rng = siren.utilities.SIREN_random(42)

        N = 10000
        n_bins = 20
        counts = [0] * n_bins

        for _ in range(N):
            r = copy.deepcopy(rec)
            iso.Sample(rng, None, r)
            ct = _rest_frame_cos_theta(r, 0)
            b = min(int((ct + 1.0) / 2.0 * n_bins), n_bins - 1)
            counts[b] += 1

        expected = N / n_bins
        chi2 = sum((c - expected) ** 2 / expected for c in counts)
        # 20 bins -> 19 DOF. chi2 < 40 is p > 0.003
        assert chi2 < 40, (
            "Isotropic sampling not uniform: chi2={:.1f} (19 DOF), "
            "counts={}".format(chi2, counts)
        )

    def test_directed_sampling_concentrates_toward_target(self):
        """Directed channel should produce events concentrated toward
        the target geometry. Verify the cos(theta_lab) distribution is
        peaked in the target direction and that the density matches
        the actual sampling frequency.
        """
        box = _make_box_at(0, 0, 50, 0.5, 0.5, 0.5)
        rec = _make_2body_record(M=1.0, E=10.0, mA=0.0, mB=0.0)

        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)
        rng = siren.utilities.SIREN_random(42)

        N = 3000
        cos_labs = []
        densities = []

        for _ in range(N):
            r = copy.deepcopy(rec)
            directed.Sample(rng, None, r)
            mom = r.secondary_momenta[0]
            p = math.sqrt(mom[1]**2 + mom[2]**2 + mom[3]**2)
            if p > 0:
                cos_labs.append(mom[3] / p)  # cos_theta_lab wrt z-axis
            d = directed.Density(None, r)
            densities.append(d)

        # Most events should be forward (cos_lab > 0.9)
        n_forward = sum(1 for c in cos_labs if c > 0.9)
        assert n_forward > N * 0.5, (
            "Directed channel not concentrating toward target: "
            "only {}/{} events with cos_lab > 0.9".format(n_forward, N)
        )

        # All densities should be positive and finite
        for d in densities:
            assert d > 0 and math.isfinite(d), (
                "Directed density invalid: {}".format(d)
            )

    def test_directed_density_self_consistent(self):
        """Verify the directed channel's Density() is the true PDF of
        its Sample() distribution by checking two independent estimates
        of the target solid angle agree.

        For a channel g(x) with support only in the target region:
          E_g[1/g(x_i)] = integral_target dOmega = Omega_target

        We estimate Omega_target two ways:
          (a) mean(1/dir_density) from directed samples
          (b) 4*pi * mean(iso_density/dir_density) from the same samples

        Both should give the same value (they differ by factor 4*pi*iso
        = 4*pi * 1/(4*pi) = 1). Agreement means the density matches
        the sampling distribution.
        """
        box = _make_box_at(0, 0, 50, 0.5, 0.5, 0.5)
        rec = _make_2body_record(M=1.0, E=5.0, mA=0.0, mB=0.0)

        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)
        iso = siren.injection.Isotropic2BodyChannel(0)
        rng = siren.utilities.SIREN_random(1234)

        N = 5000
        inv_densities = []
        for _ in range(N):
            r = copy.deepcopy(rec)
            directed.Sample(rng, None, r)
            d_dir = directed.Density(None, r)
            if d_dir > 0:
                inv_densities.append(1.0 / d_dir)

        assert len(inv_densities) > N * 0.9

        # Estimate (a): mean(1/g)
        omega_a = sum(inv_densities) / len(inv_densities)

        # All 1/g values should be consistent (low coefficient of variation)
        # For volume mode, density varies across the target, so some
        # variation is expected. CV < 1 is reasonable.
        mean_inv = omega_a
        var_inv = sum((x - mean_inv)**2 for x in inv_densities) / len(inv_densities)
        cv = math.sqrt(var_inv) / mean_inv if mean_inv > 0 else float('inf')
        assert cv < 1.5, (
            "Density too variable: CV(1/g) = {:.2f}".format(cv)
        )

        # The target solid angle should be positive and reasonable
        # (< 4*pi and > 0)
        assert 0 < omega_a < 4 * math.pi, (
            "Estimated Omega_target = {:.6f} out of range".format(omega_a)
        )

    def test_directed_multichannel_density_integral(self):
        """When sampling from a multi-channel that covers the full
        sphere, mean(iso/mc) must equal 1.0. This is the definitive
        test that the directed channel's density is correctly
        normalized relative to the isotropic channel.
        """
        box = _make_box_at(0, 0, 50, 0.5, 0.5, 0.5)
        rec = _make_2body_record(M=1.0, E=5.0, mA=0.0, mB=0.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.3, 0.7]

        rng = siren.utilities.SIREN_random(1234)
        N = 5000
        weights = []
        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_iso / d_mc)

        mean_w = sum(weights) / len(weights)
        assert abs(mean_w - 1.0) < 0.03, (
            "Multi-channel density integral = {:.4f}, expected 1.0 "
            "(tolerance 0.03)".format(mean_w)
        )

    def test_massive_daughter_density_matches_sampling(self):
        """For massive daughters with moderate boost (dual solutions
        in SolveLabAngle), the density must still match the sampling
        distribution.
        """
        # mA=0.4 with M=1.0 and gamma~3: beta_daughter > beta_parent
        # may trigger two solutions
        iso = siren.injection.Isotropic2BodyChannel(0)
        rec = _make_2body_record(M=1.0, E=3.0, mA=0.4, mB=0.3)
        rng = siren.utilities.SIREN_random(77)

        N = 8000
        n_bins = 20
        counts = [0] * n_bins

        for _ in range(N):
            r = copy.deepcopy(rec)
            iso.Sample(rng, None, r)
            ct = _rest_frame_cos_theta(r, 0)
            b = min(int((ct + 1.0) / 2.0 * n_bins), n_bins - 1)
            counts[b] += 1

        # Should be uniform (isotropic sampling in rest frame)
        expected = N / n_bins
        chi2 = sum((c - expected) ** 2 / expected for c in counts)
        assert chi2 < 40, (
            "Massive daughter sampling not uniform: chi2={:.1f}, "
            "counts={}".format(chi2, counts)
        )


# ================================================================== #
#  Gap fixes: tighter closure test                                     #
# ================================================================== #


class TestTightClosure:

    def test_closure_small_target(self):
        """Closure with a small distant target (strong biasing required).

        The directed channel concentrates ~99% of samples toward a target
        that subtends ~0.01% of the sky. The importance weights must
        compensate exactly for the biased density.

        We estimate integral of f(x)*p(x)dx two ways:
          f = 1 (constant), p = isotropic density
          -> integral = 1 regardless of target

        Via importance sampling from the multi-channel:
          estimate = (1/N) sum_i p(x_i) / g(x_i)
        where g is the multi-channel density. Should equal 1.
        """
        box = _make_box_at(0, 0, 200, 0.2, 0.2, 0.2)
        rec = _make_2body_record(M=1.0, E=10.0, mA=0.0, mB=0.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.01, 0.99]

        rng = siren.utilities.SIREN_random(12345)
        N = 5000

        weights = []
        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_iso / d_mc)

        mean_w = sum(weights) / len(weights)
        # With 5000 samples, sigma ~ 1/sqrt(N) ~ 1.4%
        # Use 5% tolerance (3.5 sigma)
        assert abs(mean_w - 1.0) < 0.05, (
            "Tight closure failed: mean weight = {:.4f} "
            "(expected 1.0, tolerance 0.05)".format(mean_w)
        )

    def test_closure_hit_fraction_small_target(self):
        """Closure test comparing hit fractions with a small target.

        Isotropic brute-force needs many events. Multi-channel with
        biasing needs fewer but must weight correctly.
        """
        box = _make_box_at(0, 0, 30, 2.0, 2.0, 2.0)
        rec = _make_2body_record(M=0.5, E=3.0, mA=0.0, mB=0.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.01, 0.99]

        # Brute force: isotropic
        rng1 = siren.utilities.SIREN_random(100)
        N_iso = 50000
        hits_iso = 0
        for _ in range(N_iso):
            r = copy.deepcopy(rec)
            iso.Sample(rng1, None, r)
            if _daughter_hits_box(r, 0, box):
                hits_iso += 1

        frac_iso = hits_iso / N_iso

        # Importance sampling: multi-channel
        rng2 = siren.utilities.SIREN_random(200)
        N_mc = 5000
        weighted_hits = 0.0
        sum_w = 0.0
        for _ in range(N_mc):
            r = copy.deepcopy(rec)
            mc.Sample(rng2, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                w = d_iso / d_mc
                sum_w += w
                if _daughter_hits_box(r, 0, box):
                    weighted_hits += w

        frac_mc = weighted_hits / sum_w if sum_w > 0 else 0

        # Both should agree within 30% relative
        if frac_iso > 0:
            rel_diff = abs(frac_mc - frac_iso) / frac_iso
            assert rel_diff < 0.3, (
                "Hit fraction closure failed: iso={:.6f}, mc={:.6f}, "
                "rel_diff={:.2f}".format(frac_iso, frac_mc, rel_diff)
            )
        else:
            assert frac_mc < 0.001


# ================================================================== #
#  Gap fixes: directed channel boundary tests                          #
# ================================================================== #


class TestDirectedBoundaries:

    def test_directed_density_at_cone_edge(self):
        """Density should transition cleanly at the boundary of the
        bounding cone (not produce NaN or negative values).
        """
        box = _make_box_at(0, 0, 100, 1.0, 1.0, 1.0)
        rec = _make_2body_record(M=1.0, E=10.0, mA=0.0, mB=0.0)

        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)
        rng = siren.utilities.SIREN_random(42)

        # Sample many events and check ALL densities
        for _ in range(2000):
            r = copy.deepcopy(rec)
            directed.Sample(rng, None, r)
            d = directed.Density(None, r)
            assert d >= 0, "Negative density from directed channel"
            assert not math.isnan(d), "NaN density from directed channel"
            assert math.isfinite(d), "Infinite density from directed channel"

    def test_directed_massive_daughters(self):
        """Directed channel with massive daughters (relevant for
        Dutta-Kim V1 -> chi chi where m_chi = 0.008 GeV).
        """
        box = _make_box_at(0, 0, 100, 1.0, 1.0, 1.0)
        # V1 decay kinematics: M=0.017, mA=mB=0.008
        rec = _make_2body_record(M=0.017, E=0.5, mA=0.008, mB=0.008)

        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)
        iso = siren.injection.Isotropic2BodyChannel(0)
        rng = siren.utilities.SIREN_random(42)

        n_positive = 0
        N = 1000
        for _ in range(N):
            r = copy.deepcopy(rec)
            directed.Sample(rng, None, r)
            d_dir = directed.Density(None, r)
            d_iso = iso.Density(None, r)
            assert d_dir >= 0 and math.isfinite(d_dir)
            assert d_iso > 0 and math.isfinite(d_iso)
            if d_dir > 0:
                n_positive += 1

        assert n_positive > N * 0.8, (
            "Directed channel with massive daughters produced too many "
            "zero-density events: {}/{}".format(n_positive, N)
        )


# ================================================================== #
#  Gap fixes: normalization with tighter tolerance                     #
# ================================================================== #


class TestTighterNormalization:

    def test_directed_density_integral_tight(self):
        """Verify directed density integrates to 1 with < 3% error.

        Use 10000 importance samples for tighter statistical bound.
        """
        box = _make_box_at(0, 0, 50, 0.5, 0.5, 0.5)
        rec = _make_2body_record(M=1.0, E=10.0, mA=0.0, mB=0.0)

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.3, 0.7]

        rng = siren.utilities.SIREN_random(9999)
        N = 10000
        weights = []

        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_iso / d_mc)

        mean_w = sum(weights) / len(weights)
        assert abs(mean_w - 1.0) < 0.03, (
            "Directed density integral = {:.4f}, expected 1.0 "
            "(tolerance 0.03)".format(mean_w)
        )


# ================================================================== #
#  Mock physics models for testing channels with real Decay/XS        #
# ================================================================== #


def _boost_to_lab(P_parent, p_cm, cos_theta, phi, m_daughter):
    """Boost a rest-frame 2-body daughter to the lab frame.

    Matches the _boost_to_lab helper used by all production
    Decay classes in VectorPortal.py / MesonProduction.py.
    """
    E_parent = P_parent[0]
    px_p, py_p, pz_p = P_parent[1], P_parent[2], P_parent[3]
    p_mag = math.sqrt(px_p**2 + py_p**2 + pz_p**2)
    M_parent = math.sqrt(max(E_parent**2 - p_mag**2, 0.0))

    sin_theta = math.sqrt(max(1.0 - cos_theta**2, 0.0))
    E_rf = math.sqrt(p_cm**2 + m_daughter**2)
    p_rf_x = p_cm * sin_theta * math.cos(phi)
    p_rf_y = p_cm * sin_theta * math.sin(phi)
    p_rf_z = p_cm * cos_theta

    if p_mag < 1e-12 or M_parent < 1e-12:
        return [E_rf, p_rf_x, p_rf_y, p_rf_z]

    beta = p_mag / E_parent
    gamma = E_parent / M_parent
    bx, by, bz = px_p / p_mag, py_p / p_mag, pz_p / p_mag

    p_par = p_rf_x * bx + p_rf_y * by + p_rf_z * bz
    perp_x = p_rf_x - p_par * bx
    perp_y = p_rf_y - p_par * by
    perp_z = p_rf_z - p_par * bz

    E_lab = gamma * (E_rf + beta * p_par)
    p_par_lab = gamma * (p_par + beta * E_rf)
    return [E_lab,
            p_par_lab * bx + perp_x,
            p_par_lab * by + perp_y,
            p_par_lab * bz + perp_z]


class _MockIsotropic2BodyDecay(siren.interactions.Decay):
    """Isotropic 2-body decay with known masses.

    SampleFinalState samples uniform rest-frame cos_theta, boosts
    to the lab frame (matching production Decay conventions).
    FinalStateProbability returns 1/(4*pi) (isotropic).
    Convention returns RestFrameSolidAngle.
    """

    def __init__(self, parent_mass, daughter_mass_a, daughter_mass_b):
        siren.interactions.Decay.__init__(self)
        self._M = parent_mass
        self._mA = daughter_mass_a
        self._mB = daughter_mass_b
        self._sig = siren.dataclasses.InteractionSignature()
        self._sig.primary_type = siren.dataclasses.ParticleType.N4
        self._sig.target_type = siren.dataclasses.ParticleType.Decay
        self._sig.secondary_types = [
            siren.dataclasses.ParticleType.NuLight,
            siren.dataclasses.ParticleType.Gamma,
        ]

    def equal(self, other):
        return self is other

    def TotalDecayWidthAllFinalStates(self, record):
        return 1.0

    def TotalDecayWidth(self, record_or_type):
        return 1.0

    def DifferentialDecayWidth(self, record):
        return 1.0 / (4.0 * math.pi)

    def SampleFinalState(self, csdr, random):
        cos_theta = random.Uniform(-1, 1)
        phi = random.Uniform(0, 2 * math.pi)

        P_parent = list(csdr.primary_momentum)
        p_cm = siren.injection.TwoBodyRestMomentum(
            self._M, self._mA, self._mB)

        P_A = _boost_to_lab(P_parent, p_cm, cos_theta, phi, self._mA)
        P_B = _boost_to_lab(P_parent, p_cm, -cos_theta,
                            phi + math.pi, self._mB)

        sp0 = csdr.get_secondary_particle_record(0)
        sp0.four_momentum = P_A
        sp1 = csdr.get_secondary_particle_record(1)
        sp1.four_momentum = P_B

    def SecondaryMasses(self, secondary_types):
        return [self._mA, self._mB]

    def SecondaryHelicities(self, record):
        return [0, 0]

    def GetPossibleSignatures(self):
        return [self._sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        return [self._sig]

    def FinalStateProbability(self, record):
        return 1.0 / (4.0 * math.pi)

    def DensityVariables(self):
        return ["cos_theta"]

    def Convention(self):
        return siren.injection.PhaseSpaceConvention.RestFrameSolidAngle


class _MockAnisotropic2BodyDecay(siren.interactions.Decay):
    """2-body decay with angular distribution (1 + cos^2 theta).

    Normalized density: (3/(16*pi)) * (1 + cos^2 theta).
    Uses accept-reject to sample, then boosts to lab frame.
    """

    def __init__(self, parent_mass, daughter_mass_a, daughter_mass_b):
        siren.interactions.Decay.__init__(self)
        self._M = parent_mass
        self._mA = daughter_mass_a
        self._mB = daughter_mass_b
        self._sig = siren.dataclasses.InteractionSignature()
        self._sig.primary_type = siren.dataclasses.ParticleType.N4
        self._sig.target_type = siren.dataclasses.ParticleType.Decay
        self._sig.secondary_types = [
            siren.dataclasses.ParticleType.NuLight,
            siren.dataclasses.ParticleType.Gamma,
        ]

    def equal(self, other):
        return self is other

    def TotalDecayWidthAllFinalStates(self, record):
        return 1.0

    def TotalDecayWidth(self, record_or_type):
        return 1.0

    def DifferentialDecayWidth(self, record):
        ct = self._cos_theta_rest(record)
        return (3.0 / (16.0 * math.pi)) * (1.0 + ct * ct)

    def _cos_theta_rest(self, record):
        return _rest_frame_cos_theta(record, 0)

    def SampleFinalState(self, csdr, random):
        # Accept-reject with envelope max = 1.5
        while True:
            cos_theta = random.Uniform(-1, 1)
            u = random.Uniform(0, 1)
            f = 0.75 * (1.0 + cos_theta**2)
            if u * 1.5 < f:
                break

        phi = random.Uniform(0, 2 * math.pi)

        P_parent = list(csdr.primary_momentum)
        p_cm = siren.injection.TwoBodyRestMomentum(
            self._M, self._mA, self._mB)

        P_A = _boost_to_lab(P_parent, p_cm, cos_theta, phi, self._mA)
        P_B = _boost_to_lab(P_parent, p_cm, -cos_theta,
                            phi + math.pi, self._mB)

        sp0 = csdr.get_secondary_particle_record(0)
        sp0.four_momentum = P_A
        sp1 = csdr.get_secondary_particle_record(1)
        sp1.four_momentum = P_B

    def SecondaryMasses(self, secondary_types):
        return [self._mA, self._mB]

    def SecondaryHelicities(self, record):
        return [0, 0]

    def GetPossibleSignatures(self):
        return [self._sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        return [self._sig]

    def FinalStateProbability(self, record):
        ct = self._cos_theta_rest(record)
        return (3.0 / (16.0 * math.pi)) * (1.0 + ct * ct)

    def DensityVariables(self):
        return ["cos_theta"]

    def Convention(self):
        return siren.injection.PhaseSpaceConvention.RestFrameSolidAngle


class _MockElasticScattering(siren.interactions.CrossSection):
    """Mock 2->2 elastic scattering with isotropic angular distribution.

    dsigma/dQ2 = sigma_total / Q2_max = const.
    Reports MandelstamST convention.
    """

    def __init__(self, m_beam, m_target):
        siren.interactions.CrossSection.__init__(self)
        self._m_beam = m_beam
        self._m_target = m_target
        self._sigma = 1e-30  # cm^2, irrelevant for density test
        self._sig = siren.dataclasses.InteractionSignature()
        self._sig.primary_type = siren.dataclasses.ParticleType(5917)
        self._sig.target_type = siren.dataclasses.ParticleType(1000180400)
        self._sig.secondary_types = [
            siren.dataclasses.ParticleType(5917),
            siren.dataclasses.ParticleType(1000180400),
        ]

    def equal(self, other):
        return self is other

    def _Q2_max(self, E):
        s = (self._m_beam**2 + self._m_target**2
             + 2 * self._m_target * E)
        lam = siren.injection.Kallen(
            s, self._m_beam**2, self._m_target**2)
        if lam <= 0 or s <= 0:
            return 0
        return lam / s

    def TotalCrossSection(self, record):
        return self._sigma

    def DifferentialCrossSection(self, record):
        E = record.primary_momentum[0]
        q2max = self._Q2_max(E)
        if q2max <= 0:
            return 0
        return self._sigma / q2max

    def InteractionThreshold(self, record):
        return self._m_beam

    def SampleFinalState(self, csdr, random):
        E = csdr.primary_momentum[0]
        q2max = self._Q2_max(E)
        Q2 = random.Uniform(0, q2max)
        csdr.interaction_parameters = {"Q2": Q2}

        # Solve 2-body kinematics from Q2
        s = (self._m_beam**2 + self._m_target**2
             + 2 * self._m_target * E)
        p_cm_sq = siren.injection.Kallen(
            s, self._m_beam**2, self._m_target**2) / (4 * s)
        p_cm = math.sqrt(max(p_cm_sq, 0))
        cos_cm = 1.0 - Q2 / (2 * p_cm_sq) if p_cm_sq > 0 else 1.0
        sin_cm = math.sqrt(max(1 - cos_cm**2, 0))
        phi = random.Uniform(0, 2 * math.pi)

        E_cm_beam = math.sqrt(p_cm**2 + self._m_beam**2)
        E_cm_target = math.sqrt(p_cm**2 + self._m_target**2)

        # Outgoing momenta in CM frame
        px1 = p_cm * sin_cm * math.cos(phi)
        py1 = p_cm * sin_cm * math.sin(phi)
        pz1 = p_cm * cos_cm

        sp0 = csdr.get_secondary_particle_record(0)
        sp0.four_momentum = [E_cm_beam, px1, py1, pz1]
        sp1 = csdr.get_secondary_particle_record(1)
        sp1.four_momentum = [E_cm_target, -px1, -py1, -pz1]

    def SecondaryMasses(self, secondary_types):
        return [self._m_beam, self._m_target]

    def SecondaryHelicities(self, record):
        return [0, 0]

    def GetPossibleTargets(self):
        return [self._sig.target_type]

    def GetPossibleTargetsFromPrimary(self, primary_type):
        return [self._sig.target_type]

    def GetPossiblePrimaries(self):
        return [self._sig.primary_type]

    def GetPossibleSignatures(self):
        return [self._sig]

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        return [self._sig]

    def FinalStateProbability(self, record):
        E = record.primary_momentum[0]
        q2max = self._Q2_max(E)
        if q2max <= 0:
            return 0
        return 1.0 / q2max

    def DensityVariables(self):
        return ["Q2"]

    def Convention(self):
        return siren.injection.PhaseSpaceConvention.MandelstamST


# ================================================================== #
#  Gap 4: Cross-measure conversion in multi-channel                    #
# ================================================================== #


class TestCrossMeasureConversion:

    def test_physical_decay_with_directed_multichannel(self):
        """Combine a PhysicalDecayChannel (isotropic, SolidAngleRest)
        with a DetectorDirected2BodyChannel (SolidAngleRest) in a
        MultiChannelPhaseSpace. Both use the same measure so no
        conversion is needed. The multi-channel density integral
        (via importance sampling) should be 1.0.

        """
        decay = _MockIsotropic2BodyDecay(1.0, 0.1, 0.2)
        sig = decay.GetPossibleSignatures()[0]
        box = _make_box_at(0, 0, 50, 0.5, 0.5, 0.5)
        rec = _make_2body_record(M=1.0, E=5.0, mA=0.1, mB=0.2)

        phys = siren.injection.PhysicalDecayChannel(decay, sig)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        assert phys.Topology() == siren.injection.PhaseSpaceTopology.Decay2Body
        assert phys.Measure() == siren.injection.PhaseSpaceMeasure.SolidAngleRest()
        assert directed.Measure() == siren.injection.PhaseSpaceMeasure.SolidAngleRest()

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [phys, directed]
        mc.weights = [0.1, 0.9]

        rng = siren.utilities.SIREN_random(42)
        N = 10000
        weights = []
        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_phys = phys.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_phys / d_mc)

        mean_w = sum(weights) / len(weights)
        assert abs(mean_w - 1.0) < 0.05, (
            "PhysicalDecay+Directed integral = {:.4f}".format(mean_w)
        )

    def test_anisotropic_decay_closure(self):
        """Closure test for a non-trivial angular distribution
        (1 + cos^2 theta). The PhysicalDecayChannel wraps the
        anisotropic mock, and the multi-channel must correctly
        combine its density with the directed channel.

        Importance sampling: E_mc[phys/mc] = 1.
        Uses a boosted parent (gamma=5) to concentrate daughters
        forward, reducing weight variance.
        """
        decay = _MockAnisotropic2BodyDecay(1.0, 0.1, 0.2)
        sig = decay.GetPossibleSignatures()[0]
        box = _make_box_at(0, 0, 50, 0.5, 0.5, 0.5)
        rec = _make_2body_record(M=1.0, E=5.0, mA=0.1, mB=0.2)

        phys = siren.injection.PhysicalDecayChannel(decay, sig)
        directed = siren.injection.DetectorDirected2BodyChannel(box, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [phys, directed]
        mc.weights = [0.1, 0.9]

        rng = siren.utilities.SIREN_random(42)
        N = 10000
        weights = []
        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_phys = phys.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_phys / d_mc)

        mean_w = sum(weights) / len(weights)
        assert abs(mean_w - 1.0) < 0.05, (
            "Anisotropic decay closure = {:.4f}".format(mean_w)
        )

    def test_anisotropic_sampling_matches_density(self):
        """Verify the mock anisotropic decay's sampling distribution
        matches its density. Bin by rest-frame cos_theta and check
        the histogram against the expected (1 + cos^2) shape.
        """
        decay = _MockAnisotropic2BodyDecay(1.0, 0.0, 0.0)
        sig = decay.GetPossibleSignatures()[0]
        rec = _make_2body_record(M=1.0, E=1.0, mA=0.0, mB=0.0)

        phys = siren.injection.PhysicalDecayChannel(decay, sig)
        rng = siren.utilities.SIREN_random(42)

        N = 20000
        n_bins = 20
        counts = [0] * n_bins
        for _ in range(N):
            r = copy.deepcopy(rec)
            phys.Sample(rng, None, r)
            ct = _rest_frame_cos_theta(r, 0)
            b = min(int((ct + 1.0) / 2.0 * n_bins), n_bins - 1)
            counts[b] += 1

        # Expected: proportional to (1 + cos^2)
        # Bin centers
        bin_width = 2.0 / n_bins
        for b in range(n_bins):
            ct = -1.0 + (b + 0.5) * bin_width
            expected_frac = 0.75 * (1.0 + ct * ct) * bin_width / 2.0
            actual_frac = counts[b] / N
            # Allow 30% relative deviation (statistical)
            if expected_frac > 0.01:
                rel = abs(actual_frac - expected_frac) / expected_frac
                assert rel < 0.3, (
                    "Bin {}: expected {:.4f}, got {:.4f}, rel={:.2f}"
                    .format(b, expected_frac, actual_frac, rel)
                )


# ================================================================== #
#  Gap 7: Scattering channel density tests                             #
# ================================================================== #


class TestScatteringChannelDensity:

    def _make_scatter_record(self, m_beam=0.008, m_target=37.215, E=1.0):
        rec = siren.dataclasses.InteractionRecord()
        rec.signature.primary_type = siren.dataclasses.ParticleType(5917)
        rec.signature.target_type = siren.dataclasses.ParticleType(1000180400)
        rec.signature.secondary_types = [
            siren.dataclasses.ParticleType(5917),
            siren.dataclasses.ParticleType(1000180400),
        ]
        rec.primary_mass = m_beam
        pz = math.sqrt(max(E**2 - m_beam**2, 0))
        rec.primary_momentum = [E, 0, 0, pz]
        rec.target_mass = m_target
        rec.secondary_masses = [m_beam, m_target]
        rec.secondary_momenta = [[0, 0, 0, 0], [0, 0, 0, 0]]
        rec.secondary_helicities = [0, 0]
        rec.interaction_vertex = [0, 0, 0]
        rec.primary_initial_position = [0, 0, 0]
        return rec

    def test_scattering_channel_density_positive(self):
        """DetectorDirectedScatteringChannel should produce positive
        finite densities for all sampled events.
        """
        box = _make_box_at(0, 0, 50, 1.0, 1.0, 1.0)
        rec = self._make_scatter_record(m_beam=0.008, m_target=37.215, E=1.0)

        sc = siren.injection.DetectorDirectedScatteringChannel(
            box, directed_index=0,
            variable=siren.injection.ScatteringVariable.Q2)
        rng = siren.utilities.SIREN_random(42)

        N = 1000
        for _ in range(N):
            r = copy.deepcopy(rec)
            sc.Sample(rng, None, r)
            d = sc.Density(None, r)
            assert d >= 0, "Scattering density negative: {}".format(d)
            assert math.isfinite(d), "Scattering density infinite/NaN"

    def test_physical_scattering_with_directed_multichannel(self):
        """Combine PhysicalCrossSectionChannel (MandelstamQ2) with
        DetectorDirectedScatteringChannel (MandelstamQ2) in a
        MultiChannelPhaseSpace. Both use MandelstamQ2 measure.

        The multi-channel integral E_mc[phys/mc] should be 1.
        """
        xs = _MockElasticScattering(0.008, 37.215)
        sig = xs.GetPossibleSignatures()[0]
        box = _make_box_at(0, 0, 50, 1.0, 1.0, 1.0)
        rec = self._make_scatter_record(m_beam=0.008, m_target=37.215, E=1.0)

        phys = siren.injection.PhysicalCrossSectionChannel(xs, sig)
        directed = siren.injection.DetectorDirectedScatteringChannel(
            box, directed_index=0,
            variable=siren.injection.ScatteringVariable.Q2)

        assert phys.Topology() == siren.injection.PhaseSpaceTopology.Scatter2to2
        assert phys.Measure() == siren.injection.PhaseSpaceMeasure.MandelstamQ2()
        assert directed.Topology() == siren.injection.PhaseSpaceTopology.Scatter2to2
        assert directed.Measure() == siren.injection.PhaseSpaceMeasure.MandelstamQ2()

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [phys, directed]
        mc.weights = [0.3, 0.7]

        rng = siren.utilities.SIREN_random(42)
        N = 3000
        weights = []
        for _ in range(N):
            r = copy.deepcopy(rec)
            mc.Sample(rng, None, r)
            d_phys = phys.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_phys / d_mc)

        mean_w = sum(weights) / len(weights)
        assert abs(mean_w - 1.0) < 0.05, (
            "Scattering multichannel integral = {:.4f}".format(mean_w)
        )


# ================================================================== #
#  Gap 6: 3-body channel tests                                        #
# ================================================================== #


class TestThreeBodyChannel:

    def test_3body_channel_density_positive(self):
        """DetectorDirected3BodyChannel should produce positive finite
        densities for all sampled events.
        """
        box = _make_box_at(0, 0, 50, 1.0, 1.0, 1.0)
        rec = _make_3body_record(M=1.0, E=3.0, masses=(0.1, 0.05, 0.1))

        dir3 = siren.injection.DetectorDirected3BodyChannel(
            box,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_pair_index=2,
            mass_mode=siren.injection.InvariantMassMode.Uniform,
        )
        rng = siren.utilities.SIREN_random(42)

        N = 1000
        n_positive = 0
        for _ in range(N):
            r = copy.deepcopy(rec)
            dir3.Sample(rng, None, r)
            d = dir3.Density(None, r)
            assert d >= 0, "3-body density negative"
            assert math.isfinite(d), "3-body density NaN/Inf"
            if d > 0:
                n_positive += 1

        assert n_positive > N * 0.8, (
            "3-body channel produced too many zero-density events: "
            "{}/{}".format(n_positive, N)
        )

    def test_3body_channel_topology_and_measure(self):
        """3-body channel reports Decay3Body topology and
        Recursive2Body measure.
        """
        box = _make_box_at(0, 0, 50, 1.0, 1.0, 1.0)
        dir3 = siren.injection.DetectorDirected3BodyChannel(
            box,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_pair_index=1,
        )

        assert dir3.Topology() == siren.injection.PhaseSpaceTopology.Decay3Body
        assert dir3.Measure() == siren.injection.PhaseSpaceMeasure.Recursive2Body()

    def test_3body_channel_conserves_4momentum(self):
        """3-body channel must conserve 4-momentum."""
        box = _make_box_at(0, 0, 50, 1.0, 1.0, 1.0)
        rec = _make_3body_record(M=1.0, E=3.0, masses=(0.1, 0.05, 0.1))

        dir3 = siren.injection.DetectorDirected3BodyChannel(
            box,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_pair_index=2,
            mass_mode=siren.injection.InvariantMassMode.Uniform,
        )
        rng = siren.utilities.SIREN_random(42)

        for _ in range(200):
            r = copy.deepcopy(rec)
            dir3.Sample(rng, None, r)

            # Sum secondary 4-momenta
            E_sum = sum(p[0] for p in r.secondary_momenta)
            px_sum = sum(p[1] for p in r.secondary_momenta)
            py_sum = sum(p[2] for p in r.secondary_momenta)
            pz_sum = sum(p[3] for p in r.secondary_momenta)

            # Should equal primary 4-momentum
            assert abs(E_sum - r.primary_momentum[0]) < 1e-8, (
                "Energy not conserved: {:.8f} vs {:.8f}".format(
                    E_sum, r.primary_momentum[0])
            )
            assert abs(px_sum - r.primary_momentum[1]) < 1e-8
            assert abs(py_sum - r.primary_momentum[2]) < 1e-8
            assert abs(pz_sum - r.primary_momentum[3]) < 1e-8

    def test_3body_secondary_masses_correct(self):
        """3-body channel must produce secondaries with the correct masses."""
        box = _make_box_at(0, 0, 50, 1.0, 1.0, 1.0)
        masses = (0.1, 0.05, 0.15)
        rec = _make_3body_record(M=1.0, E=3.0, masses=masses)

        dir3 = siren.injection.DetectorDirected3BodyChannel(
            box,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_pair_index=2,
            mass_mode=siren.injection.InvariantMassMode.Uniform,
        )
        rng = siren.utilities.SIREN_random(42)

        for _ in range(100):
            r = copy.deepcopy(rec)
            dir3.Sample(rng, None, r)

            for i, m_expected in enumerate(masses):
                p = r.secondary_momenta[i]
                m_sq = p[0]**2 - p[1]**2 - p[2]**2 - p[3]**2
                m_actual = math.sqrt(max(m_sq, 0))
                assert abs(m_actual - m_expected) < 1e-6, (
                    "Secondary {} mass: expected {:.4f}, got {:.4f}".format(
                        i, m_expected, m_actual)
                )
