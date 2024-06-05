import os
import siren

base_path = os.path.dirname(os.path.abspath(__file__))
loader_file = os.path.join(base_path, "loader.py")
siren._util.load_module("loader", loader_file)

from DarkNews.ModelContainer import ModelContainer

xs_path = siren.utilities.get_cross_section_model_path(
    f"DarkNewsTables-v{siren.utilities.darknews_version()}", must_exist=False
)

def GetDetectorModelTargets(detector_model):
    """
    Determines the targets that exist inside the detector model
    :return: lists of targets and strings
    :rtype: (list<ParticleType>, list<str>)
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
    model_container,
    upscattering_key,
    tolerance=1e-6,
    interp_tolerance=5e-2,
    always_interpolate=True,
):
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
    models,
    ups_key,
    tabel_dir,
    preferences,
):
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
                    cross_section = append(
                        load_cross_section_from_table(
                            models,
                            ups_key,
                            table_subdir,
                            tolerance=tolerance,
                            interp_tolerance=interp_tolerance,
                            always_interpolate=always_interpolate,
                        )
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
                    cross_section = append(
                        load_cross_section_from_pickle(
                            ups_key,
                            table_subdir,
                            tolerance=tolerance,
                            interp_tolerance=interp_tolerance,
                            always_interpolate=always_interpolate,
                        )
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
                cross_sections = append(
                    load_cross_section(
                        models,
                        ups_key,
                        tolerance=tolerance,
                        interp_tolerance=interp_tolerance,
                        always_interpolate=always_interpolate,
                    )
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
    model_kwargs,
    table_dir=None,
    tolerance=1e-6,
    interp_tolerance=5e-2,
    always_interpolate=True,
    preferences=None,
):
    if preferences is None:
        preferences = ["table", "pickle", "normal"]

    models = ModelContainer(**model_kwargs)

    if table_dir is None:
        table_dir = ""

    cross_sections = []
    for ups_key, ups_case in models.ups_cases.items():
        cross_sections.append(
            attempt_to_load_cross_section(models, ups_key, table_dir, preferences)
        )

    return cross_sections


def load_processes(
    primary_type=None,
    target_types=None,
    fill_tables_at_start=False,
    Emax=None,
    m4=None,
    mu_tr_mu4=None,  # GeV^-1
    UD4=0,
    Umu4=0,
    epsilon=0.0,
    gD=0.0,
    decay_product="photon",
    noHC=True,
    HNLtype="dirac",
    nuclear_targets=None,
    detector_model=None,
    tolerance=1e-6,  # supposed to represent machine epsilon
    interp_tolerance=5e-2,  # relative interpolation tolerance
    always_interpolate=True,  # bool whether to always interpolate the total/differential cross section
):

    if nuclear_targets is None and detector_model is None:
        raise ValueError(
            'Either "nuclear_targets" or "detector_model" must be provided'
        )

    if nuclear_targets is None:
        nuclear_targets = GetDetectorModelTargets(detector_model)[1]

    base_path = os.path.dirname(os.path.abspath(__file__))
    table_dir = os.path.join(base_path, "Dipole_M%2.2e_mu%2.2e" % (m4, mu_tr_mu4))

    model_kwargs = {
        "m4": m4,
        "mu_tr_mu4": mu_tr_mu4,
        "UD4": UD4,
        "Umu4": Umu4,
        "epsilon": epsilon,
        "gD": gD,
        "decay_product": decay_product,
        "noHC": noHC,
    }

    cross_sections = load_cross_sections(
       model_kwargs,
       table_dir=None,
       tolerance=tolerance,
       interp_tolerance=interp_tolerance,
       always_interpolate=always_interpolate,
    )

    if fill_tables_at_start:
        if Emax is None:
            print(
                "WARNING: Cannot fill cross section tables without specifying a maximum energy"
            )
        else:
            for cross_section in cross_sections:
                cross_section.FillInterpolationTables(Emax=Emax)

    # Initialize primary InteractionCollection
    # Loop over available cross sections and save those which match primary type
    primary_cross_sections = []
    for cross_section in self.DN_processes.cross_sections:
        if primary_type == _dataclasses.Particle.ParticleType(
            cross_section.ups_case.nu_projectile.pdgid
        ):
            primary_cross_sections.append(cross_section)
    primary_interaction_collection = _interactions.InteractionCollection(
        primary_type, primary_cross_sections
    )

    # Initialize secondary processes and define secondary InteractionCollection objects
    secondary_decays = {}
    # Also keep track of the minimum decay width for defining the position distribution later
    self.DN_min_decay_width = np.inf
    # Loop over available decays, group by parent type
    for decay in self.DN_processes.decays:
        secondary_type = _dataclasses.Particle.ParticleType(
            decay.dec_case.nu_parent.pdgid
        )
        if secondary_type not in secondary_decays.keys():
            secondary_decays[secondary_type] = []
        secondary_decays[secondary_type].append(decay)
        total_decay_width = decay.TotalDecayWidth(secondary_type)
        if total_decay_width < self.DN_min_decay_width:
            self.DN_min_decay_width = total_decay_width
    # Now make the list of secondary cross section collections
    # Add new secondary injection and physical processes at the same time
    secondary_interaction_collections = []
    for secondary_type, decay_list in secondary_decays.items():
        # Define a sedcondary injection distribution
        secondary_injection_process = _injection.SecondaryInjectionProcess()
        secondary_physical_process = _injection.PhysicalProcess()
        secondary_injection_process.primary_type = secondary_type
        secondary_physical_process.primary_type = secondary_type

        # Add the secondary position distribution
        if self.fid_vol is not None:
            secondary_injection_process.AddSecondaryInjectionDistribution(
                _distributions.SecondaryBoundedVertexDistribution(self.fid_vol)
            )
        else:
            secondary_injection_process.AddSecondaryInjectionDistribution(
                _distributions.SecondaryPhysicalVertexDistribution()
            )

        self.secondary_injection_processes.append(secondary_injection_process)
        self.secondary_physical_processes.append(secondary_physical_process)

        secondary_interaction_collections.append(
            _interactions.InteractionCollection(secondary_type, decay_list)
        )

    self.SetInteractions(
        primary_interaction_collection, secondary_interaction_collections
    )


def GetFiducialVolume(self):
    """
    :return: identified fiducial volume for the experiment, None if not found
    """
    detector_model_file = _util.get_detector_model_path(self.experiment)
    with open(detector_model_file) as file:
        fiducial_line = None
        detector_line = None
        for line in file:
            data = line.split()
            if len(data) <= 0:
                continue
            elif data[0] == "fiducial":
                fiducial_line = line
            elif data[0] == "detector":
                detector_line = line
        if fiducial_line is None or detector_line is None:
            return None
        return _detector.DetectorModel.ParseFiducialVolume(fiducial_line, detector_line)
    return None
