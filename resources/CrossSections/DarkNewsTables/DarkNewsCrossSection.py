import os
import numpy as np
import functools
from scipy.interpolate import LinearNDInterpolator, PchipInterpolator

base_path = os.path.dirname(os.path.abspath(__file__))
loader_file = os.path.join(base_path, "loader.py")
siren._util.load_module("loader", loader_file)

# SIREN methods
from siren.interactions import DarkNewsCrossSection
from siren import dataclasses
from siren.dataclasses import Particle

# DarkNews methods
from DarkNews import phase_space

# A class representing a single ups_case DarkNews class
# Only handles methods concerning the upscattering part
class PyDarkNewsCrossSection(DarkNewsCrossSection):
    def __init__(
        self,
        ups_case,  # DarkNews UpscatteringProcess instance
        tolerance=1e-6,  # supposed to represent machine epsilon
        interp_tolerance=5e-2,  # relative interpolation tolerance
        always_interpolate=True,  # bool whether to always interpolate the total/differential cross section
    ):
        DarkNewsCrossSection.__init__(self)  # C++ constructor

        self.ups_case = ups_case
        self.tolerance = tolerance
        self.interp_tolerance = interp_tolerance
        self.always_interpolate = always_interpolate

        # 2D table in E, sigma
        self.total_cross_section_table = np.empty((0, 2), dtype=float)
        # 3D table in E, z, dsigma/dQ2 where z = (Q2 - Q2min) / (Q2max - Q2min)
        self.differential_cross_section_table = np.empty((0, 3), dtype=float)

    def load_from_table(self, table_dir):
        # Make the table directory where will we store cross section tables
        if not os.path.exists(table_dir):
            try:
                os.makedirs(table_dir, exist_ok=False)
            except OSError as error:
                raise RuntimeError("Directory '%s' cannot be created" % table_dir)

        # Look in table dir and check whether total/differential xsec tables exist
        total_xsec_file = os.path.join(table_dir, "total_cross_sections.npy")
        if os.path.exists(total_xsec_file):
            self.total_cross_section_table = np.load(total_xsec_file)
        diff_xsec_file = os.path.join(
            table_dir, "differential_cross_sections.npy"
        )
        if os.path.exists(diff_xsec_file):
            self.differential_cross_section_table = np.load(diff_xsec_file)

        self.configure()

    def save_to_table(self, table_dir, total=True, diff=True):
        if total:
            self._redefine_interpolation_objects(total=True)
            with open(
                os.path.join(table_dir, "total_cross_sections.npy"), "wb"
            ) as f:
                np.save(f, self.total_cross_section_table)
        if diff:
            self._redefine_interpolation_objects(diff=True)
            with open(
                os.path.join(table_dir, "differential_cross_sections.npy"), "wb"
            ) as f:
                np.save(f, self.differential_cross_section_table)

    # serialization method
    def get_representation(self):
        return {
            "total_cross_section_table": self.total_cross_section_table,
            "differential_cross_section_table": self.differential_cross_section_table,
            "ups_case": self.ups_case,
            "tolerance": self.tolerance,
            "interp_tolerance": self.interp_tolerance,
            "always_interpolate": self.always_interpolate,
            "is_configured": False,
        }

    # Configure function to set up member variables
    # assumes we have defined the following:
    #   ups_case, total_cross_section_table, differential_cross_section_table,
    #   tolerance, interp_tolerance, always_interpolate
    #   kwargs argument can be used to set any of these
    def configure(self, **kwargs):

        for k, v in kwargs.items():
            self.__setattr__(k, v)

        # Define the target particle
        # make sure protons are stored as H nuclei
        self.target_type = Particle.ParticleType(self.ups_case.nuclear_target.pdgid)
        if self.target_type == Particle.ParticleType.PPlus:
            self.target_type = Particle.ParticleType.HNucleus

        # Initialize interpolation objects
        self.total_cross_section_interpolator = None
        self.differential_cross_section_interpolator = None
        self._redefine_interpolation_objects(total=True, diff=True)
        self.is_configured = True

    # Sorts and redefines scipy interpolation objects
    def _redefine_interpolation_objects(self, total=False, diff=False):
        if total:
            if len(self.total_cross_section_table) <= 1:
                return
            idxs = np.argsort(self.total_cross_section_table[:, 0])
            self.total_cross_section_table = self.total_cross_section_table[idxs]
            self.total_cross_section_interpolator = PchipInterpolator(
                self.total_cross_section_table[:, 0],
                self.total_cross_section_table[:, 1],
            )
        if diff:
            if len(self.differential_cross_section_table) <= 1:
                return
            idxs = np.lexsort(
                (
                    self.differential_cross_section_table[:, 1],
                    self.differential_cross_section_table[:, 0],
                )
            )
            self.differential_cross_section_table = (
                self.differential_cross_section_table[idxs]
            )
            # If we only have two energy points, don't try to construct interpolator
            if len(np.unique(self.differential_cross_section_table[:, 0])) <= 2:
                return
            self.differential_cross_section_interpolator = LinearNDInterpolator(
                self.differential_cross_section_table[:, :2],
                self.differential_cross_section_table[:, 2],
                rescale=True,
            )

    # Check whether we have close-enough entries in the intrepolation tables
    def _interpolation_flags(self, inputs, mode):
        #
        # returns UseSinglePoint,Interpolate,closest_idx
        # UseSinglePoint: whether to use a single point in table
        # Interpolate: whether to interpolate bewteen different points
        # closest_idx: index of closest point in table (for UseSinglePoint)

        # Determine which table we are using
        if mode == "total":
            interp_table = self.total_cross_section_table
        elif mode == "differential":
            interp_table = self.differential_cross_section_table
        else:
            print("Invalid interpolation table mode %s" % mode)
            exit(0)

        # first check if we have saved table points already
        if len(interp_table) == 0:
            return False, False, -1

        # bools to keep track of whether to use a single point or interpolate
        UseSinglePoint = False
        Interpolate = True
        # order events by the relative difference
        rel_diff = np.abs((interp_table[:, :-1] - inputs) / inputs)
        rel_diff_length = np.sqrt(np.sum(rel_diff**2, axis=-1))
        closest_idx_abs = np.argmin(rel_diff_length, axis=-1)
        # First check whether we have a close-enough single point
        if np.all(np.abs(rel_diff[closest_idx_abs]) < self.tolerance):
            UseSinglePoint = True
        # Ensure we have enough points to interpolate
        if len(interp_table) < len(inputs) + 1:
            Interpolate = False
        # Require that we have at least len(inputs)+1 close points to interpolate
        else:
            close = np.all(rel_diff < self.interp_tolerance, axis=-1)
            if sum(close) < len(inputs) + 1:
                Interpolate = False
        return UseSinglePoint, Interpolate, closest_idx_abs

    # return entries in interpolation table if we have inputs
    def _query_interpolation_table(self, inputs, mode):
        #
        # returns:
        # 0 if we are not close enough to any points in the interpolation table
        # otherwise, returns the desired interpolated value

        # First make sure we are configured
        self._ensure_configured()

        # Determine which table we are using
        if mode == "total":
            interp_table = self.total_cross_section_table
            interpolator = self.total_cross_section_interpolator
        elif mode == "differential":
            interp_table = self.differential_cross_section_table
            interpolator = self.differential_cross_section_interpolator
        else:
            print("Invalid interpolation table mode %s" % mode)
            exit(0)

        if self.always_interpolate:
            # check if energy is within table range

            if interpolator is None or inputs[0] > interp_table[-1, 0]:
                print(
                    "Requested interpolation at %2.2f GeV. Either this is above the table boundary or the interpolator doesn't yet exist. Filling %s table"
                    % (inputs[0], mode)
                )
                n = self.FillInterpolationTables(
                    total=(mode == "total"),
                    diff=(mode == "differential"),
                    Emax=(1 + self.interp_tolerance) * inputs[0],
                )
                print("Added %d points" % n)
                if mode == "total":
                    interpolator = self.total_cross_section_interpolator
                elif mode == "differential":
                    interpolator = self.differential_cross_section_interpolator
            elif inputs[0] < interp_table[0, 0]:
                print(
                    "Requested interpolation at %2.2f GeV below table boundary. Requring calculation"
                    % inputs[0]
                )
                return 0
            val = max(0, interpolator(inputs))
            if val < 0:
                print(
                    "WARNING: negative interpolated value for %s-%s %s cross section at,"
                    % (
                        self.ups_case.nuclear_target.name,
                        self.ups_case.scattering_regime,
                        mode,
                    ),
                    inputs,
                )
            return val

        UseSinglePoint, Interpolate, closest_idx = self._interpolation_flags(
            inputs, mode
        )

        if UseSinglePoint:
            if closest_idx < 0:
                print(
                    "Trying to use a single table point, but no closest idx found. Exiting..."
                )
                exit(0)
            return interp_table[closest_idx, -1]
        elif Interpolate:
            return interpolator(inputs)
        else:
            return -1

    def FillTableAtEnergy(self, E, total=True, diff=True, factor=0.8):
        num_added_points = 0
        if total:
            xsec = self.ups_case.total_xsec(E)
            self.total_cross_section_table = np.append(
                self.total_cross_section_table, [[E, xsec]], axis=0
            )
            num_added_points += 1
        if diff:
            interaction = dataclasses.InteractionRecord()
            interaction.signature.primary_type = self.GetPossiblePrimaries()[
                0
            ]  # only one primary
            interaction.signature.target_type = self.GetPossibleTargets()[
                0
            ]  # only one target
            interaction.target_mass = self.ups_case.MA
            interaction.primary_momentum = [E, 0, 0, 0]
            zmin, zmax = self.tolerance, 1
            Q2min = self.Q2Min(interaction)
            Q2max = self.Q2Max(interaction)
            z = zmin
            while z < zmax:
                Q2 = Q2min + z * (Q2max - Q2min)
                dxsec = self.ups_case.diff_xsec_Q2(E, Q2).item()
                self.differential_cross_section_table = np.append(
                    self.differential_cross_section_table,
                    [[E, z, dxsec]],
                    axis=0,
                )
                num_added_points += 1
                z *= 1 + factor * self.interp_tolerance
        self._redefine_interpolation_objects(total=total, diff=diff)
        return num_added_points

    # Fills the total and differential cross section tables within interp_tolerance
    def FillInterpolationTables(self, total=True, diff=True, factor=0.8, Emax=None):
        increment_factor = 0.5 * factor * self.interp_tolerance
        Emin = (1.0 + self.tolerance) * self.ups_case.Ethreshold
        if Emax is None:
            if (
                len(self.total_cross_section_table)
                + len(self.differential_cross_section_table)
            ) <= 0:
                return 0
            Emax = max(
                np.max([0] + list(self.total_cross_section_table[:, 0])),
                np.max([0] + list(self.differential_cross_section_table[:, 0])),
            )
        num_added_points = 0
        E = Emin
        E_existing_total = np.unique(self.total_cross_section_table[:, 0])
        E_existing_diff = np.unique(self.differential_cross_section_table[:, 0])
        while E < Emax:
            # sample more coarsely past 1.5*threshold
            if E > 1.5 * self.ups_case.Ethreshold:
                increment_factor = factor * self.interp_tolerance
            n = self.FillTableAtEnergy(
                E,
                total=(total and (E not in E_existing_total)),
                diff=(diff and (E not in E_existing_diff)),
                factor=factor,
            )
            num_added_points += n
            E *= 1 + increment_factor
        self._redefine_interpolation_objects(total=total, diff=diff)
        return num_added_points

    def GetPossiblePrimaries(self):
        return [Particle.ParticleType(self.ups_case.nu_projectile.pdgid)]

    def _ensure_configured(self):
        if not self.is_configured:
            self.configure()

    def GetPossibleTargetsFromPrimary(self, primary_type):
        self._ensure_configured()
        if Particle.ParticleType(self.ups_case.nu_projectile.pdgid) == primary_type:
            return [self.target_type]
        return []

    def GetPossibleTargets(self):
        self._ensure_configured()
        return [self.target_type]

    def GetPossibleSignatures(self):
        self._ensure_configured()
        signature = dataclasses.InteractionSignature()
        signature.primary_type = Particle.ParticleType(
            self.ups_case.nu_projectile.pdgid
        )
        signature.target_type = self.target_type
        signature.secondary_types = []
        signature.secondary_types.append(
            Particle.ParticleType(self.ups_case.nu_upscattered.pdgid)
        )
        signature.secondary_types.append(self.target_type)
        return [signature]

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        if (
            Particle.ParticleType(self.ups_case.nu_projectile.pdgid) == primary_type
        ) and ((self.target_type == target_type)):
            signature = dataclasses.InteractionSignature()
            signature.primary_type = Particle.ParticleType(
                self.ups_case.nu_projectile.pdgid
            )
            signature.target_type = self.target_type
            secondary_types = []
            secondary_types.append(
                Particle.ParticleType(self.ups_case.nu_upscattered.pdgid)
            )
            secondary_types.append(
                Particle.ParticleType(self.ups_case.nuclear_target.pdgid)
            )
            signature.secondary_types = secondary_types
            return [signature]
        return []

    def DifferentialCrossSection(self, arg1, target=None, energy=None, Q2=None):
        if type(arg1) == dataclasses.InteractionRecord:
            interaction = arg1
            # Calculate Q2 assuming we are in the target rest frame
            m1sq = interaction.primary_momentum[0] ** 2 - np.sum(
                [p**2 for p in interaction.primary_momentum[1:]]
            )
            m3sq = interaction.secondary_momenta[0][0] ** 2 - np.sum(
                [p**2 for p in interaction.secondary_momenta[0][1:]]
            )
            p1p3 = interaction.primary_momentum[0] * interaction.secondary_momenta[0][
                0
            ] - np.sum(
                p1 * p3
                for p1, p3 in zip(
                    interaction.primary_momentum[1:],
                    interaction.secondary_momenta[0][1:],
                )
            )
            Q2 = -(m1sq + m3sq - 2 * p1p3)
            energy = interaction.primary_momentum[0]
        else:
            primary = arg1
            interaction = dataclasses.InteractionRecord()
            interaction.signature.primary_type = primary
            interaction.signature.target_type = target
            interaction.primary_momentum = [energy, 0, 0, 0]
            interaction.target_mass = self.ups_case.MA
        if interaction.signature.primary_type != Particle.ParticleType(
            self.ups_case.nu_projectile.pdgid
        ):
            return 0
        if interaction.primary_momentum[0] < self.InteractionThreshold(interaction):
            return 0
        Q2min = self.Q2Min(interaction)
        Q2max = self.Q2Max(interaction)
        if Q2 < Q2min or Q2 > Q2max:
            return 0
        z = (Q2 - Q2min) / (Q2max - Q2min)

        if self.always_interpolate:
            # Check if we can interpolate
            val = self._query_interpolation_table([energy, z], mode="differential")
            if val >= 0:
                # we have recovered the differential cross section from the interpolation table
                return val

        # If we have reached this block, we must compute the differential cross section using DarkNews
        dxsec = self.ups_case.diff_xsec_Q2(energy, Q2).item()
        return dxsec

    def TargetMass(self, target_type):
        target_mass = self.ups_case.MA
        return target_mass

    def SecondaryMasses(self, secondary_types):
        secondary_masses = []
        secondary_masses.append(self.ups_case.m_ups)
        secondary_masses.append(self.ups_case.MA)
        return secondary_masses

    def SecondaryHelicities(self, record):
        secondary_helicities = []
        secondary_helicities.append(
            self.ups_case.h_upscattered * record.primary_helicity
        )
        secondary_helicities.append(record.target_helicity)
        self.h_ups = self.ups_case.m_ups
        self.h_target = self.ups_case.MA
        return secondary_helicities

    def TotalCrossSection(self, arg1, energy=None, target=None):
        # Handle overloaded arguments
        if type(arg1) == dataclasses.InteractionRecord:
            primary = arg1.signature.primary_type
            energy = arg1.primary_momentum[0]
            target = arg1.signature.target_type
        elif energy is not None and target is not None:
            primary = arg1
        else:
            print("Incorrect function call to TotalCrossSection!")
            exit(0)
        if int(primary) != self.ups_case.nu_projectile:
            return 0
        interaction = dataclasses.InteractionRecord()
        interaction.signature.primary_type = primary
        interaction.signature.target_type = target
        interaction.primary_momentum[0] = energy
        if energy < self.InteractionThreshold(interaction):
            # print("Python: energy %2.2f < self.InteractionThreshold(interaction) %2.2f"%(energy,self.InteractionThreshold(interaction)))
            return 0

        # Check if we can interpolate
        val = self._query_interpolation_table([energy], mode="total")
        if val >= 0:
            # we have recovered the cross section from the interpolation table
            return val

        # If we have reached this block, we must compute the cross section using DarkNews
        xsec = self.ups_case.total_xsec(energy)
        self.total_cross_section_table = np.append(
            self.total_cross_section_table, [[energy, xsec]], axis=0
        )
        self._redefine_interpolation_objects(total=True)
        return xsec

    def InteractionThreshold(self, interaction):
        return self.ups_case.Ethreshold

    def Q2Min(self, interaction):
        return phase_space.upscattering_Q2min(
            interaction.primary_momentum[0],
            self.ups_case.m_ups,
            self.ups_case.MA,
        )

    def Q2Max(self, interaction):
        return phase_space.upscattering_Q2max(
            interaction.primary_momentum[0],
            self.ups_case.m_ups,
            self.ups_case.MA,
        )
