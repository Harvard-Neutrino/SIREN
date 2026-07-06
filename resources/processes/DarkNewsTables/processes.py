import os
from typing import Tuple, List, Any, Optional
import siren
import collections
import pickle

base_path = os.path.dirname(os.path.abspath(__file__))
logger_file = os.path.join(base_path, "logger.py")
siren._util.load_module("logger", logger_file)

def _load_local_module(module_name):
    return siren._util.load_module(
        f"siren.resources.processes.DarkNewsTables.{module_name}",
        os.path.join(base_path, f"{module_name}.py"),
    )


try:
    from siren.DNModelContainer import ModelContainer
    from DarkNews.nuclear_tools import NuclearTarget

    DarkNewsDecay = _load_local_module("DarkNewsDecay")
    DarkNewsCrossSection = _load_local_module("DarkNewsCrossSection")

    PyDarkNewsCrossSection = DarkNewsCrossSection.PyDarkNewsCrossSection
    PyDarkNewsDecay = DarkNewsDecay.PyDarkNewsDecay
    _DARKNEWS_AVAILABLE = True
except ImportError:
    ModelContainer = None

    class NuclearTarget:
        pass

    PyDarkNewsCrossSection = None
    PyDarkNewsDecay = None
    _DARKNEWS_AVAILABLE = False

xs_path = siren.utilities.get_processes_model_path(
    f"DarkNewsTables-v{siren.utilities.darknews_version()}", must_exist=False
)

def GetDetectorModelTargets(detector_model: siren.detector.DetectorModel) -> Tuple[List[siren.dataclasses.Particle.ParticleType], List[str]]:
    """
    Determines the targets that exist inside the detector model.

    Args:
        detector_model (siren.detector.DetectorModel): The detector model object.

    Returns:
        Tuple[List[siren.dataclasses.Particle.ParticleType], List[str]]: A tuple containing two lists:
            - List of target objects (ParticleType)
            - List of target strings
    """
    count = 0
    targets = []
    target_strs = []
    while detector_model.Materials.HasMaterial(count):
        for target in detector_model.Materials.GetMaterialTargets(count):
            if target not in targets:
                targets.append(target)
            if str(target).find("Nucleus") == -1:
                continue
            else:
                target_str = str(target)[
                    str(target).find("Type") + 5 : str(target).find("Nucleus")
                ]
                if target_str == "H":
                    target_str = "H1"
                if target_str not in target_strs:
                    target_strs.append(target_str)
        count += 1
    return targets, target_strs


def load_cross_section(
    model_container: ModelContainer,
    upscattering_key: Any,
    tolerance: float = 1e-6,
    interp_tolerance: float = 5e-2,
    always_interpolate: bool = True,
) -> PyDarkNewsCrossSection:
    """
    Loads a cross-section object based on the given parameters.

    Args:
        model_container (ModelContainer): The model container object.
        upscattering_key (Any): The key for the upscattering model.
        tolerance (float, optional): Tolerance for calculations. Defaults to 1e-6.
        interp_tolerance (float, optional): Interpolation tolerance. Defaults to 5e-2.
        always_interpolate (bool, optional): Whether to always interpolate. Defaults to True.

    Returns:
        PyDarkNewsCrossSection: The loaded cross-section object.

    Raises:
        KeyError: If the upscattering key is not present in model_container.ups_cases.
    """
    if upscattering_key not in model_container.ups_cases:
        raise KeyError(
            f'Upscattering key "{upscattering_key}" not present in model_container.ups_cases'
        )
    upscattering_model = model_container.ups_cases[upscattering_key]
    return PyDarkNewsCrossSection(
        upscattering_model,
        tolerance=tolerance,
        interp_tolerance=interp_tolerance,
        always_interpolate=always_interpolate,
    )


def load_cross_section_from_table(
    model_container,
    upscattering_key,
    table_dir,
    tolerance=1e-6,
    interp_tolerance=5e-2,
    always_interpolate=True,
):
    # subdir = "_".join(["CrossSection"] + [str(x) if type(x)!=NuclearTarget else str(x.name)  for x in upscattering_key])
    # table_subdir = os.path.join(table_dir, subdir)

    cross_section = load_cross_section(
        model_container,
        upscattering_key,
        tolerance=tolerance,
        interp_tolerance=interp_tolerance,
        always_interpolate=always_interpolate,
    )
    cross_section.load_from_table(table_dir)
    return cross_section


def load_cross_section_from_pickle(
    upscattering_key,
    table_dir,
    tolerance=1e-6,
    interp_tolerance=5e-2,
    always_interpolate=True,
):
    import pickle
    # subdir = "_".join(["CrossSection"] + [str(x) if type(x)!=NuclearTarget else str(x.name)  for x in upscattering_key])
    # table_subdir = os.path.join(table_dir, subdir)
    fname = os.path.join(table_dir, "xs_object.pkl")
    with open(fname, "rb") as f:
        xs_obj = pickle.load(f)
        xs_obj.configure(
            tolerance=tolerance,
            interp_tolerance=interp_tolerance,
            always_interpolate=always_interpolate,
        )
        return xs_obj


def attempt_to_load_cross_section(
    models: ModelContainer,
    ups_key: Any,
    table_dir: str,
    preferences: List[str],
    tolerance: float = 1e-6,
    interp_tolerance: float = 5e-2,
    always_interpolate: bool = True,
) -> PyDarkNewsCrossSection:
    """
    Attempts to load a cross-section object using different strategies based on preferences.

    Args:
        models (ModelContainer): The model container object.
        ups_key (Any): The key for the upscattering model.
        table_dir (str): Directory path for the tables.
        preferences (List[str]): List of loading preferences (e.g., ["table", "pickle", "normal"]).
        tolerance (float, optional): Tolerance for calculations. Defaults to 1e-6.
        interp_tolerance (float, optional): Interpolation tolerance. Defaults to 5e-2.
        always_interpolate (bool, optional): Whether to always interpolate. Defaults to True.

    Returns:
        PyDarkNewsCrossSection: The loaded cross-section object.

    Raises:
        ValueError: If preferences list is empty.
        RuntimeError: If unable to load the cross-section with any strategy.
    """
    if len(preferences) == 0:
        raise ValueError("preferences must have at least one entry")

    subdir = "_".join(["CrossSection"] + [str(x) if type(x)!=NuclearTarget else str(x.name)  for x in ups_key])
    loaded = False
    cross_section = None
    for p in preferences:
        if p == "table":
            table_subdir = os.path.join(table_dir, subdir)
            if os.path.isdir(table_subdir):
                try:
                    cross_section = load_cross_section_from_table(
                        models,
                        ups_key,
                        table_subdir,
                        tolerance=tolerance,
                        interp_tolerance=interp_tolerance,
                        always_interpolate=always_interpolate,
                    )
                    loaded = True
                except Exception as e:
                    print(
                        "Encountered exception while loading DN cross section from table"
                    )
                    raise e from None
                break
        elif p == "pickle":
            table_subdir = os.path.join(table_dir, subdir)
            if os.path.isdir(table_subdir):
                try:
                    cross_section = load_cross_section_from_pickle(
                        ups_key,
                        table_subdir,
                        tolerance=tolerance,
                        interp_tolerance=interp_tolerance,
                        always_interpolate=always_interpolate,
                    )
                    loaded = True
                except Exception as e:
                    print(
                        "Encountered exception while loading DN cross section from pickle"
                    )
                    raise e from None
                break
        elif p == "normal":
            try:
                cross_section = load_cross_section(
                    models,
                    ups_key,
                    tolerance=tolerance,
                    interp_tolerance=interp_tolerance,
                    always_interpolate=always_interpolate,
                )
                loaded = True
            except Exception as e:
                print("Encountered exception while loading DN cross section normally")
                raise e from None
            break

    if not loaded:
        raise RuntimeError("Not able to load DN cross section with any strategy")
    return cross_section


def load_cross_sections(
    models,
    table_dir=None,
    tolerance=1e-6,
    interp_tolerance=5e-2,
    always_interpolate=True,
    preferences=None,
):
    if preferences is None:
        preferences = ["table", "pickle", "normal"]

    if table_dir is None:
        table_dir = ""

    cross_sections = {}
    for ups_key, ups_case in models.ups_cases.items():
        cross_sections[ups_key] = (
            attempt_to_load_cross_section(models, ups_key,
                                          table_dir,
                                          preferences,
                                          tolerance,
                                          interp_tolerance,
                                          always_interpolate)
        )

    return cross_sections


def load_decay(
    model_container,
    decay_key,
):
    if decay_key not in model_container.dec_cases:
        raise KeyError(
            f'Decay key "{decay_key}" not present in model_container.dec_cases'
        )
    decay_model = model_container.dec_cases[decay_key]
    return PyDarkNewsDecay(
        decay_model,
    )


def load_decay_from_table(
    model_container,
    decay_key,
    table_dir,
):
    subdir = "_".join(["Decay"] + [str(x) for x in decay_key])
    table_subdir = os.path.join(table_dir, subdir)

    decay = load_decay(
        model_container,
        decay_key,
    )
    decay.load_from_table(table_subdir)
    return decay


def load_decay_from_pickle(
    decay_key,
    table_dir,
):
    import pickle
    subdir = "_".join(["Decay"] + [str(x) for x in decay_key])
    table_subdir = os.path.join(table_dir, subdir)
    fname = os.path.join(table_dir, "dec_object.pkl")
    with open(fname, "rb") as f:
        dec_obj = pickle.load(f)
        return dec_obj


def attempt_to_load_decay(
    models,
    decay_key,
    table_dir,
    preferences,
):
    if len(preferences) == 0:
        raise ValueError("preferences must have at least one entry")

    subdir = "_".join(["Decay"] + [str(x) for x in decay_key])
    loaded = False
    decay = None
    for p in preferences:
        if p == "table":
            table_subdir = os.path.join(table_dir, subdir)
            if os.path.isdir(table_subdir):
                try:
                    decay = load_decay_from_table(
                        models,
                        decay_key,
                        table_subdir,
                    )
                    loaded = True
                except Exception as e:
                    print("Encountered exception while loading DN decay from table")
                    raise e from None
                break
        elif p == "pickle":
            table_subdir = os.path.join(table_dir, subdir)
            if os.path.isdir(table_subdir):
                try:
                    decay = load_decay_from_pickle(
                        decay_key,
                        table_dir,
                    )
                    loaded = True
                except Exception as e:
                    print("Encountered exception while loading DN decay from pickle")
                    raise e from None
                break
        elif p == "normal":
            try:
                decay = load_decay(
                    models,
                    decay_key,
                )
                loaded = True
            except Exception as e:
                print("Encountered exception while loading DN decay normally")
                raise e from None
            break

    if not loaded:
        raise RuntimeError("Not able to load DN decay with any strategy")
    return decay


def load_decays(
    models,
    table_dir=None,
    preferences=None,
):
    if preferences is None:
        preferences = ["table", "pickle", "normal"]

    if table_dir is None:
        table_dir = ""

    decays = {}
    for decay_key, dec_case in models.dec_cases.items():
        decays[decay_key] = attempt_to_load_decay(models, decay_key, table_dir, preferences)

    return decays

def load_processes(
    primary_type: Optional[Any] = None,
    target_types: Optional[List[Any]] = None,
    fill_tables_at_start: bool = False,
    Emax: Optional[float] = None,
    nuclear_targets: Optional[List[str]] = None,
    detector_model: Optional[Any] = None,
    tolerance: float = 1e-6,
    interp_tolerance: float = 5e-2,
    always_interpolate: bool = True,
    table_name: Optional[str] = None,
    **model_kwargs,
    # m4: Optional[float] = None,
    # mu_tr_mu4: Optional[float] = None,
    # UD4: float = 0,
    # Umu4: float = 0,
    # epsilon: float = 0.0,
    # gD: float = 0.0,
    # decay_product: str = "photon",
    # noHC: bool = True,
    # HNLtype: str = "dirac",

) -> List[Any]:
    """
    Loads and returns a list of cross-section and decay objects based on the given parameters.

    Args:
        primary_type (Optional[Any]): The primary particle type.
        target_types (Optional[List[Any]]): List of target particle types.
        fill_tables_at_start (bool): Whether to fill interpolation tables at start.
        Emax (Optional[float]): Maximum energy for table filling.
        detector_model (Optional[Any]): Detector model object.
        tolerance (float): Tolerance for calculations.
        interp_tolerance (float): Interpolation tolerance.
        always_interpolate (bool): Whether to always interpolate.
        table_name: Optional[str] = None,
        **model_kwargs: dictionary of DarkNews model arguments

    Returns:
        List[Any]: A list of loaded cross-section and decay objects.

    Raises:
        ValueError: If neither nuclear_targets nor detector_model is provided.
    """

    if nuclear_targets is None and detector_model is None:
        raise ValueError(
            'Either "nuclear_targets" or "detector_model" must be provided'
        )
    if not _DARKNEWS_AVAILABLE:
        raise ImportError(
            "DarkNews is required for the generic DarkNewsTables loader. "
            "Use load_vector_portal_offshell or load_meson_scalar for the "
            "self-contained Dutta-Kim models."
        )

    if nuclear_targets is None:
        nuclear_targets = GetDetectorModelTargets(detector_model)[1]
    model_kwargs["nuclear_targets"] = list(nuclear_targets)
    if target_types: model_kwargs["nuclear_targets"]+=list(target_types)

    base_path = os.path.dirname(os.path.abspath(__file__))
    table_dir = os.path.join(base_path, table_name)

    models = ModelContainer(
        **model_kwargs
    )

    cross_sections = load_cross_sections(
        models,
        table_dir=table_dir,
        tolerance=tolerance,
        interp_tolerance=interp_tolerance,
        always_interpolate=always_interpolate,
    )

    decays = load_decays(
        models,
        table_dir=table_dir,
    )

    cross_sections = {k:xs for k,xs in cross_sections.items() if len([s for s in xs.GetPossibleSignatures() if s.primary_type == primary_type])>0}

    if fill_tables_at_start:
        if Emax is None:
            print(
                "WARNING: Cannot fill cross section tables without specifying a maximum energy"
            )
        else:
            for cross_section in cross_sections:
                cross_section.FillInterpolationTables(Emax=Emax)

    primary_processes = collections.defaultdict(list)
    primary_ups_keys = collections.defaultdict(list)
    # Loop over available cross sections and save those which match primary type
    for ups_key,cross_section in cross_sections.items():
        if primary_type == siren.dataclasses.Particle.ParticleType(
            cross_section.ups_case.nu_projectile.pdgid
        ):
            primary_processes[primary_type].append(cross_section)
            primary_ups_keys[primary_type].append(ups_key)

    secondary_processes = collections.defaultdict(list)
    secondary_dec_keys = collections.defaultdict(list)
    # Loop over available decays, group by parent type
    for dec_key,decay in decays.items():
        secondary_type = siren.dataclasses.Particle.ParticleType(
            decay.dec_case.nu_parent.pdgid
        )
        secondary_processes[secondary_type].append(decay)
        secondary_dec_keys[secondary_type].append(dec_key)

    return dict(primary_processes), dict(secondary_processes), dict(primary_ups_keys), dict(secondary_dec_keys)

def SaveDarkNewsProcesses(table_dir,
                          primary_processes,
                          primary_ups_keys,
                          secondary_processes,
                          secondary_dec_keys,
                          pickles=True):
    for primary in primary_processes.keys():
        for xs, ups_key in zip(primary_processes[primary], primary_ups_keys[primary]):
            subdir = "_".join(["CrossSection"] + [str(x) if type(x)!=NuclearTarget else str(x.name) for x in ups_key])
            table_subdir = os.path.join(table_dir, subdir)
            os.makedirs(table_subdir,exist_ok=True)
            print("Saving cross section table at %s" % table_subdir)
            xs.FillInterpolationTables()
            xs.save_to_table(table_subdir)
            # if pickles:
            #     with open(os.path.join(table_subdir, "xs_object.pkl"),"wb") as f:
            #         pickle.dump(xs,f)
    for secondary in secondary_processes.keys():
        for dec, dec_key in zip(secondary_processes[secondary],secondary_dec_keys[secondary]):
            subdir = "_".join(["Decay"] + [str(x) if type(x)!=NuclearTarget else str(x.name) for x in dec_key])
            table_subdir = os.path.join(table_dir, subdir)
            os.makedirs(table_subdir,exist_ok=True)
            print("Saving decay object at %s" % table_subdir)
            dec.save_to_table(table_subdir)
            if pickles:
                with open(os.path.join(table_subdir, "dec_object.pkl"),"wb") as f:
                    pickle.dump(dec,f)

def SaveDarkNewsProcesses(table_dir,
                          primary_processes,
                          primary_ups_keys,
                          secondary_processes,
                          secondary_dec_keys,
                          pickles=True):
    for primary in primary_processes.keys():
        for xs, ups_key in zip(primary_processes[primary], primary_ups_keys[primary]):
            subdir = "_".join(["CrossSection"] + [str(x) if type(x)!=NuclearTarget else str(x.name) for x in ups_key])
            table_subdir = os.path.join(table_dir, subdir)
            os.makedirs(table_subdir,exist_ok=True)
            print("Saving cross section table at %s" % table_subdir)
            xs.FillInterpolationTables()
            xs.save_to_table(table_subdir)
            # if pickles:
            #     with open(os.path.join(table_subdir, "xs_object.pkl"),"wb") as f:
            #         pickle.dump(xs,f)
    for secondary in secondary_processes.keys():
        for dec, dec_key in zip(secondary_processes[secondary],secondary_dec_keys[secondary]):
            subdir = "_".join(["Decay"] + [str(x) if type(x)!=NuclearTarget else str(x.name) for x in dec_key])
            table_subdir = os.path.join(table_dir, subdir)
            os.makedirs(table_subdir,exist_ok=True)
            print("Saving decay object at %s" % table_subdir)
            dec.save_to_table(table_subdir)
            if pickles:
                with open(os.path.join(table_subdir, "dec_object.pkl"),"wb") as f:
                    pickle.dump(dec,f)


# ===================================================================
#  Factory functions for Dutta-Kim models (no DarkNews dependency)
# ===================================================================

def load_vector_portal(
    m_chi,
    m_chi_prime,
    m_V1,
    m_V2,
    g_D,
    epsilon_1,
    epsilon_2,
    detector_model,
    *,
    pdgid_chi=5917,
    pdgid_chi_prime=5918,
    pdgid_V1=5922,
    flux_tag="pion_numu",
    m_meson=0.13957,
    m_lepton=0.10566,
    min_energy=None,
    max_energy=3.0,
    table_dir=None,
    physically_normalized=True,
):
    """
    Construct primary and secondary process dicts for the vector-portal
    dark matter model:  chi N -> chi' N,  chi' -> chi V1,  V1 -> e+ e-.

    Parameters
    ----------
    m_chi, m_chi_prime, m_V1, m_V2 : float
        Particle masses in GeV.  m_V2 is the heavy mediator for upscattering.
    g_D : float
        Dark gauge coupling.
    epsilon_1, epsilon_2 : float
        Kinetic mixing parameters for V1 and V2.
    detector_model : siren.detector.DetectorModel
        Detector model (used to determine target nuclei).
    flux_tag : str
        Flux table tag, e.g. "pion_numu", "kaon_numu".
    m_meson, m_lepton : float
        Parent meson and lepton masses for the production channel.

    Returns
    -------
    primary_processes : dict
        {chi_type: [PyDarkNewsCrossSection]}
    secondary_processes : dict
        {chi_prime_type: [ChiPrimeDecay], v1_type: [DarkPhotonDecay]}
    chi_flux : siren.distributions.TabulatedFluxDistribution
    """
    if not _DARKNEWS_AVAILABLE:
        raise ImportError(
            "load_vector_portal uses PyDarkNewsCrossSection and requires DarkNews. "
            "Use load_vector_portal_offshell for the self-contained Dutta-Kim model."
        )

    VectorPortal = _load_local_module("VectorPortal")
    (
        VectorPortalUpsCase,
        ChiPrimeDecay,
        DarkPhotonDecay,
        compute_chi_flux,
    ) = (
        VectorPortal.VectorPortalUpsCase,
        VectorPortal.ChiPrimeDecay,
        VectorPortal.DarkPhotonDecay,
        VectorPortal.compute_chi_flux,
    )

    _, target_strs = GetDetectorModelTargets(detector_model)

    chi_type = siren.dataclasses.Particle.ParticleType(pdgid_chi)
    chi_prime_type = siren.dataclasses.Particle.ParticleType(pdgid_chi_prime)
    v1_type = siren.dataclasses.Particle.ParticleType(pdgid_V1)

    target_db = {
        "Ar40": (1000180400, 37.215, 40, 18),
        "C12":  (1000060120, 11.178, 12, 6),
        "O16":  (1000080160, 14.899, 16, 8),
        "H1":   (1000010010,  0.938,  1, 1),
        "Fe56": (1000260560, 52.103, 56, 26),
        "Si28": (1000140280, 26.060, 28, 14),
        "Ca40": (1000200400, 37.225, 40, 20),
        "N14":  (1000070140, 13.044, 14, 7),
    }

    primary_processes = {chi_type: []}
    for tname in target_strs:
        if tname not in target_db:
            continue
        nuc_pdgid, nuc_mass, nuc_A, nuc_Z = target_db[tname]
        ups_case = VectorPortalUpsCase(
            m_chi=m_chi,
            m_chi_prime=m_chi_prime,
            m_V=m_V2,
            g_D=g_D,
            epsilon=epsilon_2,
            pdgid_chi=pdgid_chi,
            pdgid_chi_prime=pdgid_chi_prime,
            nuclear_pdgid=nuc_pdgid,
            nuclear_mass=nuc_mass,
            nuclear_name=tname,
            A=nuc_A,
            Z=nuc_Z,
        )
        xs = PyDarkNewsCrossSection(ups_case, always_interpolate=True)
        primary_processes[chi_type].append(xs)

    td = table_dir or "."
    chi_prime_decay = ChiPrimeDecay(
        m_chi, m_chi_prime, m_V1, g_D,
        pdgid_chi_prime=pdgid_chi_prime,
        pdgid_chi=pdgid_chi,
        pdgid_V1=pdgid_V1,
        table_dir=td,
    )
    v1_decay = DarkPhotonDecay(
        m_V1, epsilon_1,
        pdgid_V1=pdgid_V1,
        table_dir=td,
    )

    secondary_processes = {
        chi_prime_type: [chi_prime_decay],
        v1_type: [v1_decay],
    }

    if min_energy is None:
        min_energy = ups_case.Ethreshold if primary_processes[chi_type] else 0.001

    chi_flux = compute_chi_flux(
        m_meson=m_meson,
        m_lepton=m_lepton,
        m_V1=m_V1,
        m_chi=m_chi,
        m_chi_prime=m_chi_prime,
        g_D=g_D,
        epsilon=epsilon_1,
        flux_tag=flux_tag,
        min_energy=min_energy,
        max_energy=max_energy,
        physically_normalized=physically_normalized,
    )

    return primary_processes, secondary_processes, chi_flux


def load_vector_portal_offshell(
    m_chi,
    m_chi_prime,
    m_V1,
    m_V2,
    g_D,
    epsilon_1,
    epsilon_2,
    detector_model,
    *,
    pdgid_chi=5917,
    pdgid_V1=5922,
    nuclear_pdgid=1000180400,
    nuclear_mass=37.215,
    nuclear_name="Ar40",
    A=40,
    Z=18,
    table_dir=None,
):
    """Construct self-contained Dutta-Kim off-shell scattering processes."""
    VectorPortal = _load_local_module("VectorPortal")

    chi_type = siren.dataclasses.Particle.ParticleType(pdgid_chi)
    v1_type = siren.dataclasses.Particle.ParticleType(pdgid_V1)

    xs = VectorPortal.VectorPortalOffShellXS(
        m_chi=m_chi,
        m_chi_prime=m_chi_prime,
        m_V1=m_V1,
        m_V2=m_V2,
        g_D=g_D,
        epsilon=epsilon_2,
        pdgid_chi=pdgid_chi,
        pdgid_V1=pdgid_V1,
        nuclear_pdgid=nuclear_pdgid,
        nuclear_mass=nuclear_mass,
        nuclear_name=nuclear_name,
        A=A,
        Z=Z,
    )
    v1_decay = VectorPortal.DarkPhotonDecay(
        m_V1,
        epsilon_1,
        pdgid_V1=pdgid_V1,
        table_dir=table_dir or ".",
    )

    return {chi_type: [xs]}, {v1_type: [v1_decay]}


def load_vector_portal_onshell(
    m_chi,
    m_chi_prime,
    m_V1,
    m_V2,
    g_D,
    epsilon_1,
    epsilon_2,
    detector_model,
    *,
    pdgid_chi=5917,
    pdgid_chi_prime=5918,
    pdgid_V1_prod=5922,
    pdgid_V1_signal=5923,
    nuclear_pdgid=1000180400,
    nuclear_mass=37.215,
    nuclear_name="Ar40",
    A=40,
    Z=18,
    table_dir=None,
):
    """Construct topology-preserving Dutta-Kim on-shell processes."""
    VectorPortal = _load_local_module("VectorPortal")

    chi_type = siren.dataclasses.Particle.ParticleType(pdgid_chi)
    chi_prime_type = siren.dataclasses.Particle.ParticleType(pdgid_chi_prime)
    v1_prod_type = siren.dataclasses.Particle.ParticleType(pdgid_V1_prod)
    v1_signal_type = siren.dataclasses.Particle.ParticleType(pdgid_V1_signal)

    v1_to_chi = VectorPortal.DarkPhotonToChiDecay(
        m_V1,
        m_chi,
        g_D,
        pdgid_V1=pdgid_V1_prod,
        pdgid_chi=pdgid_chi,
        table_dir=table_dir or ".",
    )
    upscatter = VectorPortal.VectorPortalUpscatteringXS(
        m_chi=m_chi,
        m_chi_prime=m_chi_prime,
        m_V2=m_V2,
        g_D=g_D,
        epsilon=epsilon_2,
        pdgid_chi=pdgid_chi,
        pdgid_chi_prime=pdgid_chi_prime,
        nuclear_pdgid=nuclear_pdgid,
        nuclear_mass=nuclear_mass,
        nuclear_name=nuclear_name,
        A=A,
        Z=Z,
    )
    chi_prime_decay = VectorPortal.ChiPrimeDecay(
        m_chi,
        m_chi_prime,
        m_V1,
        g_D,
        pdgid_chi_prime=pdgid_chi_prime,
        pdgid_chi=pdgid_chi,
        pdgid_V1=pdgid_V1_signal,
        table_dir=table_dir or ".",
    )
    visible_decay = VectorPortal.DarkPhotonDecay(
        m_V1,
        epsilon_1,
        pdgid_V1=pdgid_V1_signal,
        table_dir=table_dir or ".",
    )

    primary_processes = {chi_type: [upscatter]}
    secondary_processes = {
        v1_prod_type: [v1_to_chi],
        chi_prime_type: [chi_prime_decay],
        v1_signal_type: [visible_decay],
    }
    return primary_processes, secondary_processes


def load_meson_scalar(
    m_meson,
    m_lepton,
    m_phi,
    g_mu,
    mediator_type="scalar",
    flux_tag="pion_numu",
    min_energy=0.0,
    max_energy=3.0,
    n_bins=50,
    physically_normalized=True,
):
    """
    Construct the scalar/pseudoscalar mediator flux from charged meson
    three-body decay.

    Returns
    -------
    phi_flux : siren.distributions.TabulatedFluxDistribution
    decay_info : MesonThreeBodyDecay
        The decay object (for inspecting widths, branching ratios, etc.)
    """
    MesonProduction = _load_local_module("MesonProduction")
    MesonThreeBodyDecay = MesonProduction.MesonThreeBodyDecay
    build_phi_flux = MesonProduction.build_phi_flux

    decay = MesonThreeBodyDecay(m_meson, m_lepton, m_phi, g_mu, mediator_type)

    phi_flux = build_phi_flux(
        m_meson=m_meson,
        m_lepton=m_lepton,
        m_phi=m_phi,
        g_mu=g_mu,
        mediator_type=mediator_type,
        flux_tag=flux_tag,
        min_energy=min_energy,
        max_energy=max_energy,
        n_bins=n_bins,
        physically_normalized=physically_normalized,
    )

    return phi_flux, decay
