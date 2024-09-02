import os
from typing import Tuple, List, Any, Optional
import siren
import collections

base_path = os.path.dirname(os.path.abspath(__file__))
logger_file = os.path.join(base_path, "logger.py")
siren._util.load_module("logger", logger_file)

from siren.DNModelContainer import ModelContainer

# Import PyDarkNewsDecay and PyDarkNewsCrossSection
decay_file = os.path.join(base_path, "DarkNewsDecay.py")
cross_section_file = os.path.join(base_path, "DarkNewsCrossSection.py")
DarkNewsDecay = siren._util.load_module("DarkNewsDecay", decay_file)
DarkNewsCrossSection = siren._util.load_module("DarkNewsCrossSection", cross_section_file)

PyDarkNewsCrossSection = DarkNewsCrossSection.PyDarkNewsCrossSection
PyDarkNewsDecay = DarkNewsDecay.PyDarkNewsDecay

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
    subdir = "_".join(["CrossSection"] + [str(x) for x in upscattering_key])
    table_subdir = os.path.join(table_dir, subdir)

    cross_section = load_cross_section(
        model_container,
        upscattering_key,
        tolerance=tolerance,
        interp_tolerance=interp_tolerance,
        always_interpolate=always_interpolate,
    )
    cross_section.load_from_table(table_subdir)
    return cross_section


def load_cross_section_from_pickle(
    upscattering_key,
    table_dir,
    tolerance=1e-6,
    interp_tolerance=5e-2,
    always_interpolate=True,
):
    import pickle
    subdir = "_".join(["CrossSection"] + [str(x) for x in upscattering_key])
    table_subdir = os.path.join(table_dir, subdir)
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

    subdir = "_".join(["CrossSection"] + [str(x) for x in ups_key])
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

    cross_sections = []
    for ups_key, ups_case in models.ups_cases.items():
        cross_sections.append(
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

    decays = []
    for decay_key, dec_case in models.dec_cases.items():
        decays.append(attempt_to_load_decay(models, decay_key, table_dir, preferences))

    return decays

# This class is a hacky workaround for an issue with the python reference counting of classes derived
# from a pybind11 trampoline class i.e. python cross-section classes and python decay classes that
# inherit from siren.interactions.CrossSection or siren.interactions.Decay. If these python classes
# are passed to the InteractionCollection constructor, but a python-side reference to them is not
# maintained, then their python side state/memory will be destroyed/deallocated. This class maintains
# a python-side reference to all PyDarkNewsCrossSection and PyDarkNewsDecay instances created by
# load_processes(...) to avoid this issue
class Holder:
    holders = []
    def __init__(self):
        Holder.holders.append(self)

def load_processes(
    primary_type: Optional[Any] = None,
    target_types: Optional[List[Any]] = None,
    fill_tables_at_start: bool = False,
    Emax: Optional[float] = None,
    m4: Optional[float] = None,
    mu_tr_mu4: Optional[float] = None,
    UD4: float = 0,
    Umu4: float = 0,
    epsilon: float = 0.0,
    gD: float = 0.0,
    decay_product: str = "photon",
    noHC: bool = True,
    HNLtype: str = "dirac",
    nuclear_targets: Optional[List[str]] = None,
    detector_model: Optional[Any] = None,
    tolerance: float = 1e-6,
    interp_tolerance: float = 5e-2,
    always_interpolate: bool = True,
) -> List[Any]:
    """
    Loads and returns a list of cross-section and decay objects based on the given parameters.

    Args:
        primary_type (Optional[Any]): The primary particle type.
        target_types (Optional[List[Any]]): List of target particle types.
        fill_tables_at_start (bool): Whether to fill interpolation tables at start.
        Emax (Optional[float]): Maximum energy for table filling.
        m4 (Optional[float]): Mass parameter.
        mu_tr_mu4 (Optional[float]): Transition magnetic moment parameter.
        UD4 (float): UD4 parameter.
        Umu4 (float): Umu4 parameter.
        epsilon (float): Epsilon parameter.
        gD (float): gD parameter.
        decay_product (str): Type of decay product.
        noHC (bool): noHC parameter.
        HNLtype (str): Type of HNL (e.g., "dirac").
        nuclear_targets (Optional[List[str]]): List of nuclear targets.
        detector_model (Optional[Any]): Detector model object.
        tolerance (float): Tolerance for calculations.
        interp_tolerance (float): Interpolation tolerance.
        always_interpolate (bool): Whether to always interpolate.

    Returns:
        List[Any]: A list of loaded cross-section and decay objects.

    Raises:
        ValueError: If neither nuclear_targets nor detector_model is provided.
    """

    if nuclear_targets is None and detector_model is None:
        raise ValueError(
            'Either "nuclear_targets" or "detector_model" must be provided'
        )

    if nuclear_targets is None:
        nuclear_targets = GetDetectorModelTargets(detector_model)[1]

    base_path = os.path.dirname(os.path.abspath(__file__))
    table_dir = os.path.join(base_path, "Dipole_M%2.2e_mu%2.2e" % (m4, mu_tr_mu4))

    models = ModelContainer(
        nuclear_targets=nuclear_targets,
        m4=m4,
        mu_tr_mu4=mu_tr_mu4,
        UD4=UD4,
        Umu4=Umu4,
        epsilon=epsilon,
        gD=gD,
        decay_product=decay_product,
        noHC=noHC,
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

    cross_sections = [xs for xs in cross_sections if len([s for s in xs.GetPossibleSignatures() if s.primary_type == primary_type])>0]

    if fill_tables_at_start:
        if Emax is None:
            print(
                "WARNING: Cannot fill cross section tables without specifying a maximum energy"
            )
        else:
            for cross_section in cross_sections:
                cross_section.FillInterpolationTables(Emax=Emax)

    primary_processes = collections.defaultdict(list)
    # Loop over available cross sections and save those which match primary type
    for cross_section in cross_sections:
        if primary_type == siren.dataclasses.Particle.ParticleType(
            cross_section.ups_case.nu_projectile.pdgid
        ):
            primary_processes[primary_type].append(cross_section)

    secondary_processes = collections.defaultdict(list)
    # Loop over available decays, group by parent type
    for decay in decays:
        secondary_type = siren.dataclasses.Particle.ParticleType(
            decay.dec_case.nu_parent.pdgid
        )
        secondary_processes[secondary_type].append(decay)


    holder = Holder()
    holder.primary_processes = primary_processes
    holder.secondary_processes = secondary_processes

    return dict(primary_processes), dict(secondary_processes)

