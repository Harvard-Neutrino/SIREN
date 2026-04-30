import numpy as np
import math

xs_file = os.path.join(base_path, "differential_xs.py")
differential_xs = siren._util.load_module("siren.resources.processes.dark_scalar.differential_xs", xs_file)

neutrinos = [
    siren.dataclasses.Particle.ParticleType.NuE,
    siren.dataclasses.Particle.ParticleType.NuMu,
    siren.dataclasses.Particle.ParticleType.NuTau,
]
antineutrinos = [
    siren.dataclasses.Particle.ParticleType.NuEBar,
    siren.dataclasses.Particle.ParticleType.NuMuBar,
    siren.dataclasses.Particle.ParticleType.NuTauBar,
]

sm_masses = {
    11: 0.00051099892,  # Electron
    13: 0.105658369,  # Muon
    15: 1.77699,  # Tau
    12: 0,  # Neutrino
    14: 0,  # Neutrino
    16: 0,  # Neutrino
    22: 0,  # Photon
    23: 91.1876,  # Z boson
    24: 80.385,  # W boson
    111: 0.1349766,  # Neutral pion
    211: 0.13957018,  # Charged pion
    221: 0.54751,  # Eta meson
    311: 0.497648,  # Neutral kaon
    321: 0.493677,  # Charged kaon
    2212: 0.93827203,  # Proton
    2112: 0.93956536,  # Neutron
}

m_e = sm_masses[11]

decay_lengths = {
    111: 8.5e-17,
    211: 2.603e-8,
    321: 1.238e-8,
}

piplus_perpot = 0.1 #0.0520
kaplus_perpot = 0

# Necessary for three body decay calculation
Gf = 1.166e-5  #Fermi's constant in GeV**-2
fK = 0.11  # Kaon Decay constant in GeV
fpi = 0.093
Vus = 0.2253  #CKM Matrix element between u and s
Vud = 0.99
GamKPDG = 5.0513e-17  #pdg Kaon decay width in GeV
GamPiPDG = 2.52e-17

alphaem = 1/137
eqed = np.sqrt(4*math.pi*alphaem)

HBARC = 197 #MeVfm
m_e = 0.000511

def sm_masses(pid):
    #mass of common particles in GeV
    #leptons
    if pid in [11, -11]: return 0.000511
    elif pid in [13, -13]: return 0.105658
    elif pid in [15, -15]: return 1.777
    elif pid in [12, -12, 14, -14, 16, -16]: return 0

    #gauge bosons
    if pid in [22]: return 0
    elif pid in [23]: return 91.1876
    elif pid in [24]: return 80.385

    #mesons
    if pid in [111]: return 0.135  #pi0
    elif pid in [211, -211]: return 0.13957  #pi+-
    elif pid in [221]: return 0.547862  #eta
    elif pid in [311]: return 0.497611  #K0
    elif pid in [321, -321]: return 0.493677  #K+-

    #hadrons
    if pid in [2212, -2212]: return 0.938272
    elif pid in [2112]: return 0.939565

def mass_sm(pid):
    """
    Return the mass of common particles based on their PDG ID.

    Args:
        pid (int): PDG ID of the particle.

    Returns:
        float: Mass of the particle in GeV.
    """
    if abs(pid) in sm_masses:
        return sm_masses[abs(pid)]
    else:
        raise ValueError(f"Unsupported particle ID: {pid}")
        
def scaling(pid0):
    if pid0 in [321]: return kaplus_perpot
    elif pid0 in [211]: return piplus_perpot

def BFpdg(pid):
    # pions
    if pid in [211, -211]: return 2.5256e-17
    # kaons
    elif pid in [321, -321]: return 5.0513e-17

# c * lifetime in CM frame
def ctau(pid):
    if pid in [111]: return 3e8*8.5e-17
    elif pid in [211, -211]: return 3e8*2.603e-8
    elif pid in [321, -321]: return 3e8*1.238e-8
    else: return 0

# factor for three body decays
def Gfac(pid):
    if pid in [211, -211]: return Gf*Vud*fpi
    elif pid in [321, -321]: return Gf*Vus*fK

# Sample decay position for a given particle and its momentum
def decay_point2(pid, momentum):
    ct = ctau(pid)
    # Boost to lab frame
    dbarz = ct * momentum/mass_sm(pid)
    p = random.uniform(0,1)
    z = -dbarz*np.log(1-p)
    return z # meters

# Differential and total cross sections for scattering processes

## Shell model functions

a0 = 0.523 #fm
def r0_func(An):  #fm
    return 1.23*An**(1/3)

def R1val(An): #in fm
    c = r0_func(An) - 0.6
    a = 0.52
    s = 0.9
    return np.sqrt(c**2 + 7*math.pi**2*a**2/3 - 5*s**2)

def atomic_formfac(Z, Q2):
    a = 184.15*np.e**(-0.5)*Z**(-1/3)/m_e
    #print(a**2*Q2)
    return Z*a**2*Q2/(1 + a**2*Q2)


#Compute mandelstam t and prefactor values for center of mass differential cross section calculations
def get_CM_factors(E1, m_photon, m_nucleus, m_scalar_mediator, m_nucleus):
    # Compute the center-of-mass energy squared (s)
    s = m_photon**2 + m_nucleus**2 + 2 * self.E_photon * m_nucleus
    
    # Compute energies and momenta in center-of-mass frame
    E1star = (s + m_photon**2 - m_nucleus**2) / (2 * np.sqrt(s))
    E3star = (s + m_scalar_mediator**2 - m_nucleus**2) / (2 * np.sqrt(s))
    p1star = np.sqrt(E1star**2 - m_photon**2)
    p3star = np.sqrt(E3star**2 - m_scalar_mediator**2)
    
    # Compute Mandelstam t
    tval = m_photon**2 + m_scalar_mediator**2 - 2 * (E1star * E3star - p1star * p3star * self.cos_theta)
    
    # Compute prefactor (phase-space and flux factor)
    prefactor = 1 / (16 * math.pi * (s - (m_photon + m_nucleus)**2) * (s - (m_photon - m_nucleus)**2))

    return tval, prefactor

# j_1
def bes(x):
    # Some issue with the scipy builtin j_1
    return np.sin(x)/(x**2) - np.cos(x)/x

def helm_formfac(Z, An, Q2):
    R1 = R1val(An)
    s = 0.9
    Q2val = Q2/(0.197)**2
    x = np.sqrt(Q2val)*R1
    #print(r0**2*Q2val)
    return Z*np.exp(-s**2*Q2val/2)*(3*bes(x)/x)

# Final form factor combines atomic, electron cloud (-Z), and helm form factors
def FormFacSq(Z, An, Q2):  #1807.10973
    return (atomic_formfac(Z, Q2) - Z + helm_formfac(Z, An, Q2))**2 #

def Primakoff_Scal_matrixelement(gpg, s, t, m_scalar_mediator, m_nucleus, Z):
    # Nc, Qf = 1, -1
    # beta = m_scalar_mediator**2/(4*mloop**2)
    pref = eqed*gpg  # gphigamma_loop(gll, ml, mphi) + gpg #Yee*(alphaem*0.303/(math.pi*mloop))*Flepton2(beta, Nc, Qf)

    return (0.5)*FormFacSq(Z, m_nucleus, -t)*pref**2*((- m_scalar_mediator**4*(2*m_nucleus**2 + t) ) + 2*m_scalar_mediator**2*t*(m_nucleus**2 +
            s + t) + t*(-2*(m_nucleus**2 - s)**2 - 2*s*t - t**2))/(t**2)

    
class ScalarDMCrossSection(siren.interactions.CrossSection):
    def __init__(self, g_phi_e_e, scalar_mediator_mass, E_photon, Z, cos_theta, process="Primakoff"):
        self.g_phi_e_e = g_phi_e_e  # Coupling to electrons
        self.mediator_mass = scalar_mediator_mass  # Mass of scalar mediator
        self.E_photon = E_photon  # Incoming photon energy
        self.Z = Z  # Atomic number
        self.cos_theta = cos_theta  # Scattering angle cosine
        self.process = process  # Process type

    def DifferentialCrossSection(self, interaction):
        """Computes the differential cross section $\frac{d\sigma}{d\cos\theta} for the selected process in COM frame"""
        m_scalar_mediator = self.mediator_mass
        m_photon = mass_sm(22)
        m_nucleus = 37.26  #siren.dataclasses.Particle.ParticleType.Ar40Nucleus #mass of Ar40Nucleus in GeV/c^2
        E_photon = self.E_photon

        s = m_photon**2 + m_nucleus**2 + 2 * self.E_photon * m_nucleus
    
        # Compute energies and momenta in center-of-mass frame
        E1star = (s + m_photon**2 - m_nucleus**2) / (2 * np.sqrt(s))
        E3star = (s + m_scalar_mediator**2 - m_nucleus**2) / (2 * np.sqrt(s))
        p1star = np.sqrt(E1star**2 - m_photon**2)
        p3star = np.sqrt(E3star**2 - m_scalar_mediator**2)
        
        # Compute Mandelstam t
        t_mandelstam = m_photon**2 + m_scalar_mediator**2 - 2 * (E1star * E3star - p1star * p3star * self.cos_theta)
        
        # Compute prefactor (phase-space and flux factor)
        prefactor = 1 / (16 * math.pi * (s - (m_photon + m_nucleus)**2) * (s - (m_photon - m_nucleus)**2))

        # Compute the differential cross section based on the selected process
        dsigdt = Primakoff_Scal_matrixelement(self.g_gamgam, s, t_mandelstam, m_scalar_mediator, m_nucleus, self.Z) * prefactor

        # Compute Jacobian factor
        jacobian = 2 * p1star * p3star  # dt/dcosθ in GeV²

        # Convert to cm² using conversion factor (HBARC in GeV·cm)
        return dsigdt * jacobian * (HBARC * 1e-3 * 1e-13)**2

