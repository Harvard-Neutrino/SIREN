import os
import numpy as np
import functools

base_path = os.path.dirname(os.path.abspath(__file__))
loader_file = os.path.join(base_path, "loader.py")
siren._util.load_module("loader", loader_file)

# SIREN methods
from siren.interactions import DarkNewsDecay
from siren import dataclasses
from siren.dataclasses import Particle

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
        if table_dir is None:
            print(
                "No table_dir specified; will sample from new VEGAS integrator for each decay"
            )
            print("WARNING: this will siginficantly slow down event generation")
            return

        # Make the table directory where will we store cross section integrators
        table_dir_exists = False
        if os.path.exists(table_dir):
            # print("Directory '%s' already exists"%table_dir)
            table_dir_exists = True
        else:
            try:
                os.makedirs(table_dir, exist_ok=False)
                print("Directory '%s' created successfully" % table_dir)
            except OSError as error:
                print("Directory '%s' cannot be created" % table_dir)
                exit(0)

        # Try to find the decay integrator
        int_file = os.path.join(table_dir, "decay_integrator.pkl")
        if os.path.isfile(int_file):
            with open(int_file, "rb") as ifile:
                self.decay_integrator = pickle.load(ifile)
        # Try to find the normalization information
        norm_file = os.path.join(table_dir, "decay_norm.json")
        if os.path.isfile(norm_file):
            with open(
                norm_file,
            ) as nfile:
                self.decay_norm = json.load(nfile)


    # serialization method
    def get_representation(self):
        return {"decay_integrator":self.decay_integrator,
                "decay_norm":self.decay_norm,
                "dec_case":self.dec_case,
                "PS_samples":self.PS_samples,
                "PS_weights":self.PS_weights,
                "PS_weights_CDF":self.PS_weights_CDF,
                "total_width":self.total_width,
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

        if type(self.dec_case) == FermionSinglePhotonDecay:
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

        elif type(self.dec_case) == FermionDileptonDecay:
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
        if type(arg1) == dataclasses.InteractionRecord:
            primary = arg1.signature.primary_type
        elif type(arg1) == dataclasses.Particle.ParticleType:
            primary = arg1
        else:
            print("Incorrect function call to TotalDecayWidth!")
            exit(0)
        if int(primary) != self.dec_case.nu_parent:
            return 0
        if self.total_width is None:
            # Need to set the total width
            if type(self.dec_case) == FermionDileptonDecay and (
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
        if type(self.dec_case) == FermionSinglePhotonDecay:
            return "cost"
        elif type(self.dec_case) == FermionDileptonDecay:
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
            self.dec_case,
            np.expand_dims(np.array(record.primary_momentum), 0),
        )

        secondaries = record.GetSecondaryParticleRecords()

        if type(self.dec_case) == FermionSinglePhotonDecay:
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

        elif type(self.dec_case) == FermionDileptonDecay:
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

