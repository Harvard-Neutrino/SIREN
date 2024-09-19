import os
import numpy as np
import pickle

from siren import _util

base_path = os.path.dirname(os.path.abspath(__file__))
logger_file = os.path.join(base_path, "logger.py")
_util.load_module("logger", logger_file)

# SIREN methods
from siren.interactions import DarkNewsDecay
from siren import dataclasses
from siren.dataclasses import Particle

# DarkNews methods
import DarkNews
from DarkNews.processes import FermionDileptonDecay, FermionSinglePhotonDecay
from DarkNews import processes as proc
from DarkNews import Cfourvec as Cfv
from DarkNews import phase_space

def get_decay_momenta_from_vegas_samples(vsamples, MC_case, decay_case, PN_LAB):
    """
    Construct the four momenta of all final state particles in the decay process from the
    vegas weights.

    Args:
            vsamples (np.ndarray): integration samples obtained from vegas
                            as hypercube coordinates. Always in the interval [0,1].

            MC_case (DarkNews.process.dec_case): the decay class of DarkNews

            PN_LAB (np.ndarray): four-momentum of the upscattered N in the lab frame: [E, pX, pY, pZ]

    Returns:
            dict: each key corresponds to a set of four momenta for a given particle involved,
                    so the values are 2D np.ndarrays with each row a different event and each column a different
                    four momentum component. Contains also the weights.
    """

    four_momenta = {}

    # N boost parameters
    boost_scattered_N = {
        "EP_LAB": PN_LAB.T[0],
        "costP_LAB": Cfv.get_cosTheta(PN_LAB),
        "phiP_LAB": np.arctan2(PN_LAB.T[2], PN_LAB.T[1]),
    }

    #######################
    # DECAY PROCESSES

    if type(decay_case) == proc.FermionDileptonDecay:

        mh = decay_case.m_parent
        mf = decay_case.m_daughter
        mm = decay_case.mm
        mp = decay_case.mm

        if decay_case.vector_on_shell or decay_case.scalar_on_shell:

            if decay_case.vector_on_shell and decay_case.scalar_off_shell:
                m_mediator = decay_case.mzprime
            elif decay_case.vector_off_shell and decay_case.scalar_on_shell:
                m_mediator = decay_case.mhprime
            else:
                raise NotImplementedError("Both mediators on-shell is not yet implemented.")

            ########################
            ### HNL decay
            N_decay_samples = {"unit_cost": np.array(vsamples[0])}
            # Ni (k1) --> Nj (k2)  Z' (k3)
            masses_decay = {
                "m1": mh,  # Ni
                "m2": mf,  # Nj
                "m3": m_mediator,  # Z'
            }
            # Phnl, Phnl_daughter, Pz'
            P1LAB_decay, P2LAB_decay, P3LAB_decay = phase_space.two_body_decay(N_decay_samples, boost=boost_scattered_N, **masses_decay, rng=MC_case.rng)

            # Z' boost parameters
            boost_Z = {
                "EP_LAB": P3LAB_decay.T[0],
                "costP_LAB": Cfv.get_cosTheta(P3LAB_decay),
                "phiP_LAB": np.arctan2(P3LAB_decay.T[2], P3LAB_decay.T[1]),
            }

            ########################
            ### Z' decay
            Z_decay_samples = {}  # all uniform
            # Z'(k1) --> ell- (k2)  ell+ (k3)
            masses_decay = {
                "m1": m_mediator,  # Ni
                "m2": mp,  # \ell+
                "m3": mm,  # \ell-
            }
            # PZ', pe-, pe+
            P1LAB_decayZ, P2LAB_decayZ, P3LAB_decayZ = phase_space.two_body_decay(Z_decay_samples, boost=boost_Z, **masses_decay, rng=MC_case.rng)

            four_momenta["P_decay_N_parent"] = P1LAB_decay
            four_momenta["P_decay_N_daughter"] = P2LAB_decay
            four_momenta["P_decay_ell_minus"] = P2LAB_decayZ
            four_momenta["P_decay_ell_plus"] = P3LAB_decayZ

        elif decay_case.vector_off_shell and decay_case.scalar_off_shell:

            ########################
            # HNL decay
            N_decay_samples = {
                "unit_t": vsamples[0],
                "unit_u": vsamples[1],
                "unit_c3": vsamples[2],
                "unit_phi34": vsamples[3],
            }

            # Ni (k1) --> ell-(k2)  ell+(k3)  Nj(k4)
            masses_decay = {
                "m1": mh,  # Ni
                "m2": mm,  # ell-
                "m3": mp,  # ell+
                "m4": mf,
            }  # Nj
            # Phnl, pe-, pe+, pnu
            (
                P1LAB_decay,
                P2LAB_decay,
                P3LAB_decay,
                P4LAB_decay,
            ) = phase_space.three_body_decay(N_decay_samples, boost=boost_scattered_N, **masses_decay, rng=MC_case.rng)

            four_momenta["P_decay_N_parent"] = P1LAB_decay
            four_momenta["P_decay_ell_minus"] = P2LAB_decay
            four_momenta["P_decay_ell_plus"] = P3LAB_decay
            four_momenta["P_decay_N_daughter"] = P4LAB_decay

    elif type(decay_case) == proc.FermionSinglePhotonDecay:

        mh = decay_case.m_parent
        mf = decay_case.m_daughter

        ########################
        ### HNL decay
        N_decay_samples = {"unit_cost": np.array(vsamples[0])}
        # Ni (k1) --> Nj (k2)  gamma (k3)
        masses_decay = {
            "m1": mh,  # Ni
            "m2": mf,  # Nj
            "m3": 0.0,  # gamma
        }
        # Phnl, Phnl', Pgamma
        P1LAB_decay, P2LAB_decay, P3LAB_decay = phase_space.two_body_decay(N_decay_samples, boost=boost_scattered_N, **masses_decay, rng=MC_case.rng)

        four_momenta["P_decay_N_parent"] = P1LAB_decay
        four_momenta["P_decay_N_daughter"] = P2LAB_decay
        four_momenta["P_decay_photon"] = P3LAB_decay

    return four_momenta


class _FakeMCInterface:
    def __init__(self, random):
        self.random = random
        self.rng_func = np.frompyfunc(lambda x: self.random.Uniform(0, 1), 1, 1)
        self.rng = lambda x: self.rng_func(np.empty(x)).astype(float)


# A class representing a single decay_case DarkNews class
# Only handles methods concerning the decay part
class PyDarkNewsDecay(DarkNewsDecay):
    def __init__(self, dec_case):
        DarkNewsDecay.__init__(self)  # C++ constructor
        self.dec_case = dec_case

        # Some variables for storing the decay phase space integrator
        self.decay_integrator = None
        self.decay_norm = None
        self.PS_samples = None
        self.PS_weights = None
        self.PS_weights_CDF = None
        self.total_width = None

    def load_from_table(self, table_dir):
        # Make the table directory where will we store cross section integrators
        if not os.path.exists(table_dir):
            try:
                os.makedirs(table_dir, exist_ok=False)
            except OSError as error:
                raise RuntimeError("Directory '%s' cannot be created" % table_dir)

        # Try to find the decay integrator
        decay_file = os.path.join(table_dir, "decay.pkl")
        if os.path.isfile(decay_file):
            with open(decay_file, "rb") as f:
                self.decay_norm, self.decay_integrator = pickle.load(f)

    def save_to_table(self, table_dir):
        with open(os.path.join(table_dir, "decay.pkl")) as f:
            pickle.dump(f, {
                "decay_integrator": self.decay_integrator,
                "decay_norm": self.decay_norm
            })

    # serialization method
    def get_representation(self):
        return {"decay_integrator": self.decay_integrator,
                "decay_norm": self.decay_norm,
                "dec_case": self.dec_case,
                "PS_samples": self.PS_samples,
                "PS_weights": self.PS_weights,
                "PS_weights_CDF": self.PS_weights_CDF,
                "total_width": self.total_width,
               }

    def SetIntegratorAndNorm(self, decay_norm, decay_integrator):
        self.decay_norm = decay_norm
        self.decay_integrator = decay_integrator

    def GetPossibleSignatures(self):
        signature = dataclasses.InteractionSignature()
        signature.primary_type = Particle.ParticleType(self.dec_case.nu_parent.pdgid)
        signature.target_type = Particle.ParticleType.Decay
        secondary_types = []
        secondary_types.append(Particle.ParticleType(self.dec_case.nu_daughter.pdgid))
        for secondary in self.dec_case.secondaries:
            secondary_types.append(Particle.ParticleType(secondary.pdgid))
        signature.secondary_types = secondary_types
        return [signature]

    def GetPossibleSignaturesFromParent(self, primary_type):
        if Particle.ParticleType(self.dec_case.nu_parent.pdgid) == primary_type:
            signature = dataclasses.InteractionSignature()
            signature.primary_type = Particle.ParticleType(
                self.dec_case.nu_parent.pdgid
            )
            signature.target_type = Particle.ParticleType.Decay
            secondary_types = []
            secondary_types.append(
                Particle.ParticleType(self.dec_case.nu_daughter.pdgid)
            )
            for secondary in self.dec_case.secondaries:
                secondary_types.append(Particle.ParticleType(secondary.pdgid))
            signature.secondary_types = secondary_types
            return [signature]
        return []

    def DifferentialDecayWidth(self, record):
        # Momentum variables of HNL necessary for calculating decay phase space
        PN = np.array(record.primary_momentum)

        if isinstance(self.dec_case, FermionSinglePhotonDecay):
            gamma_idx = 0
            for secondary in record.signature.secondary_types:
                if secondary == dataclasses.Particle.ParticleType.Gamma:
                    break
                gamma_idx += 1
            if gamma_idx >= len(record.signature.secondary_types):
                print("No gamma found in the list of secondaries!")
                exit(0)

            Pgamma = np.array(record.secondary_momenta[gamma_idx])
            momenta = np.expand_dims(PN, 0), np.expand_dims(Pgamma, 0)

        elif isinstance(self.dec_case, FermionDileptonDecay):
            lepminus_idx = -1
            lepplus_idx = -1
            nu_idx = -1
            for idx, secondary in enumerate(record.signature.secondary_types):
                if secondary in [
                    dataclasses.Particle.ParticleType.EMinus,
                    dataclasses.Particle.ParticleType.MuMinus,
                    dataclasses.Particle.ParticleType.TauMinus,
                ]:
                    lepminus_idx = idx
                elif secondary in [
                    dataclasses.Particle.ParticleType.EPlus,
                    dataclasses.Particle.ParticleType.MuPlus,
                    dataclasses.Particle.ParticleType.TauPlus,
                ]:
                    lepplus_idx = idx
                else:
                    nu_idx = idx
            if -1 in [lepminus_idx, lepplus_idx, nu_idx]:
                print("Couldn't find two leptons and a neutrino in the final state!")
                exit(0)
            Pnu = np.array(record.secondary_momenta[nu_idx])
            Plepminus = np.array(record.secondary_momenta[lepminus_idx])
            Plepplus = np.array(record.secondary_momenta[lepplus_idx])
            momenta = (
                np.expand_dims(PN, 0),
                np.expand_dims(Plepminus, 0),
                np.expand_dims(Plepplus, 0),
                np.expand_dims(Pnu, 0),
            )
        else:
            print("%s is not a valid decay class type!" % type(self.dec_case))
            exit(0)
        return self.dec_case.differential_width(momenta)

    def TotalDecayWidth(self, arg1):
        if isinstance(arg1, dataclasses.InteractionRecord):
            primary = arg1.signature.primary_type
        elif isinstance(arg1, dataclasses.Particle.ParticleType):
            primary = arg1
        else:
            print("Incorrect function call to TotalDecayWidth!")
            exit(0)
        if int(primary) != self.dec_case.nu_parent:
            return 0
        if self.total_width is None:
            # Need to set the total width
            if isinstance(self.dec_case, FermionDileptonDecay) and (
                self.dec_case.vector_off_shell and self.dec_case.scalar_off_shell
            ):
                # total width calculation requires evaluating an integral
                if self.decay_integrator is None or self.decay_norm is None:
                    # We need to initialize a new VEGAS integrator in DarkNews
                    self.total_width, dec_norm, dec_integrator = self.dec_case.total_width(
                        return_norm=True, return_dec=True
                    )
                    self.SetIntegratorAndNorm(dec_norm, dec_integrator)
                else:
                    self.total_width = (
                        self.decay_integrator["diff_decay_rate_0"].mean
                        * self.decay_norm["diff_decay_rate_0"]
                    )
            else:
                self.total_width = self.dec_case.total_width()
        return self.total_width

    def TotalDecayWidthForFinalState(self, record):
        sig = self.GetPossibleSignatures()[0]
        if (
            (record.signature.primary_type != sig.primary_type)
            or (record.signature.target_type != sig.target_type)
            or (len(record.signature.secondary_types) != len(sig.secondary_types))
            or (
                np.any(
                    [
                        record.signature.secondary_types[i] != sig.secondary_types[i]
                        for i in range(len(sig.secondary_types))
                    ]
                )
            )
        ):
            return 0
        ret = self.dec_case.total_width()
        return ret

    def DensityVariables(self):
        if isinstance(self.dec_case, FermionSinglePhotonDecay):
            return "cost"
        elif isinstance(self.dec_case, FermionDileptonDecay):
            if self.dec_case.vector_on_shell and self.dec_case.scalar_on_shell:
                print("Can't have both the scalar and vector on shell")
                exit(0)
            elif (self.dec_case.vector_on_shell and self.dec_case.scalar_off_shell) or (
                self.dec_case.vector_off_shell and self.dec_case.scalar_on_shell
            ):
                return "cost"
            elif self.dec_case.vector_off_shell and self.dec_case.scalar_off_shell:
                return "t,u,c3,phi34"
        else:
            print("%s is not a valid decay class type!" % type(self.dec_case))
            exit(0)
        return ""

    def GetPSSample(self, random):
        # Make the PS weight CDF if that hasn't been done
        if self.PS_weights_CDF is None:
            self.PS_weights_CDF = np.cumsum(self.PS_weights)

        # Random number to determine
        x = random.Uniform(0, self.PS_weights_CDF[-1])

        # find first instance of a CDF entry greater than x
        PSidx = np.argmax(x - self.PS_weights_CDF <= 0)
        return self.PS_samples[:, PSidx]

    def GetPSSample(self, random):
        # Make the PS weight CDF if that hasn't been done
        if self.PS_weights_CDF is None:
            self.PS_weights_CDF = np.cumsum(self.PS_weights)

        # Check that the CDF makes sense
        total_weight = self.PS_weights_CDF[-1]
        if total_weight == 0:
            raise ValueError("Total weight is zero, cannot sample")

        # Random number to determine
        x = random.Uniform(0, total_weight)

        # find first instance of a CDF entry greater than x
        PSidx = np.argmax(x - self.PS_weights_CDF <= 0)
        return self.PS_samples[:, PSidx]

    def SampleRecordFromDarkNews(self, record, random):
        # First, make sure we have PS samples and weights
        if self.PS_samples is None or self.PS_weights is None:
            # We need to generate new PS samples
            if self.decay_integrator is None or self.decay_norm is None:
                # We need to initialize a new VEGAS integrator in DarkNews
                (self.PS_samples, PS_weights_dict), dec_norm, dec_integrator = self.dec_case.SamplePS(
                    return_norm=True, return_dec=True
                )
                self.PS_weights = PS_weights_dict["diff_decay_rate_0"]
                self.SetIntegratorAndNorm(dec_norm, dec_integrator)
            else:
                # We already have an integrator, we just need new PS samples
                self.PS_samples, PS_weights_dict = self.dec_case.SamplePS(
                    existing_integrator=self.decay_integrator
                )
                self.PS_weights = PS_weights_dict["diff_decay_rate_0"]

        # Now we must sample an PS point on the hypercube
        PS = self.GetPSSample(random)

        # Find the four-momenta associated with this point
        # Expand dims required to call DarkNews function on signle sample
        four_momenta = get_decay_momenta_from_vegas_samples(
            np.expand_dims(PS, 0),
            _FakeMCInterface(random),
            self.dec_case,
            np.expand_dims(np.array(record.primary_momentum), 0),
        )

        secondaries = record.get_secondary_particle_records()

        if isinstance(self.dec_case, FermionSinglePhotonDecay):
            gamma_idx = 0
            for secondary in record.signature.secondary_types:
                if secondary == dataclasses.Particle.ParticleType.Gamma:
                    break
                gamma_idx += 1
            if gamma_idx >= len(record.signature.secondary_types):
                print("No gamma found in the list of secondaries!")
                exit(0)
            nu_idx = 1 - gamma_idx
            secondaries[gamma_idx].four_momentum = np.squeeze(four_momenta["P_decay_photon"])
            secondaries[gamma_idx].mass = 0
            secondaries[nu_idx].four_momentum = np.squeeze(four_momenta["P_decay_N_daughter"])
            secondaries[nu_idx].mass = 0

        elif isinstance(self.dec_case, FermionDileptonDecay):
            lepminus_idx = -1
            lepplus_idx = -1
            nu_idx = -1
            for idx, secondary in enumerate(record.signature.secondary_types):
                if secondary in [
                    dataclasses.Particle.ParticleType.EMinus,
                    dataclasses.Particle.ParticleType.MuMinus,
                    dataclasses.Particle.ParticleType.TauMinus,
                ]:
                    lepminus_idx = idx
                elif secondary in [
                    dataclasses.Particle.ParticleType.EPlus,
                    dataclasses.Particle.ParticleType.MuPlus,
                    dataclasses.Particle.ParticleType.TauPlus,
                ]:
                    lepplus_idx = idx
                else:
                    nu_idx = idx
            if -1 in [lepminus_idx, lepplus_idx, nu_idx]:
                print([lepminus_idx, lepplus_idx, nu_idx])
                print(record.signature.secondary_types)
                print("Couldn't find two leptons and a neutrino in the final state!")
                exit(0)
            secondaries[lepminus_idx].four_momentum = (
                np.squeeze(four_momenta["P_decay_ell_minus"])
            )
            secondaries[lepplus_idx].four_momentum = (
                np.squeeze(four_momenta["P_decay_ell_plus"])
            )
            secondaries[nu_idx].four_momentum = (
                np.squeeze(four_momenta["P_decay_N_daughter"])
            )
        return record

