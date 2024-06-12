import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from parse_output import analysis
filename = "output/FullSim.parquet"
plt.style.use('paper.mplstyle')


sim = analysis(filename)


c = 3 * 1e8 # m/s
m_D0 = 1.86962 # GeV
m_Dp = 1.86484
t_Dp = 1040 * 1e-15 # s
t_D0 = 410 * 1e-15
m_ice = 18.02 # g/mol
N = 6.02214 * 1e23 #mol^-1
rho = 0.917 # g/cm^3


def normalize(hist, xbins, ybins):
    normed_hist = np.zeros_like(hist)
    for i in range(len(xbins) - 1):
        tot = 0
        for j in range(len(ybins) - 1):
            tot += hist[i][j]
        for j in range(len(ybins) - 1):
            if tot != 0:
                normed_hist[i][j] = hist[i][j] / tot 
            else:
                normed_hist[i][j] = 0
    return normed_hist

def analytic_decay_length(E, t, m):
    return E * t / ((m/(c ** 2)) * c)

def xsec(E):
    return (np.exp(1.891 + 0.2095 * np.log10(E)) - 2.157 + 1.263 * np.log10(E)) * 1e-27 # convert to cm^2

def analytic_free_path(E):
    return (m_ice / (rho * N * xsec(E))) / 100 # convert to m

def plot_separation_distribution(analysis_):
    D0_energies, D0_separations, Dp_energies, Dp_separations = analysis_.separation_analysis()
    min_eng = 1e1
    max_eng = 1e9
    energy_bins = np.logspace(np.log10(min_eng), np.log10(max_eng), 20)

    energy_bins_centers = np.zeros((len(energy_bins) - 1,))
    for i in range(len(energy_bins_centers)):
        energy_bins_centers[i] = np.sqrt(energy_bins[i] * energy_bins[i + 1])
    D0_analytic_lengths = analytic_decay_length(energy_bins_centers, t_D0, m_D0)
    Dp_analytic_lengths = analytic_decay_length(energy_bins_centers, t_Dp, m_Dp)

    min_sep = 1e-3
    max_sep = 50000
    log_separation_bins = np.logspace(np.log10(min_sep), np.log10(max_sep), 20)

    X2, Y2 = np.meshgrid(energy_bins, log_separation_bins)
    log_hist_D0, _, _ = np.histogram2d(D0_energies, D0_separations, bins = (energy_bins, log_separation_bins))
    log_hist_Dp, _, _ = np.histogram2d(Dp_energies, Dp_separations, bins = (energy_bins, log_separation_bins))
    log_hist_D0 = normalize(log_hist_D0, energy_bins, log_separation_bins)
    log_hist_Dp = normalize(log_hist_Dp, energy_bins, log_separation_bins)

    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (11, 5))

    log_im1 = axes[0].pcolor(X2, Y2, log_hist_D0.T, cmap="plasma", alpha = 0.7, vmin=0, vmax=1)
    log_im2 = axes[1].pcolor(X2, Y2, log_hist_Dp.T, cmap="plasma", alpha = 0.7, vmin=0, vmax=1)

    # divider1 = make_axes_locatable(axes[0])
    # cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    divider2 = make_axes_locatable(axes[1])
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(log_im2, cax=cax2, orientation='vertical', alpha = 0.7)
    # fig.colorbar(log_im1, cax=cax1, orientation='vertical', alpha = 0.7)

    axes[0].set_title(r"$D^0$ Separation")
    axes[1].set_title(r"$D^+$ Separation")

    axes[0].set_xlabel(r"$E_{D^0}$ [GeV]")
    axes[1].set_xlabel(r"$E_{D^+}$ [GeV]")

    axes[0].set_ylabel("Separation Length [m]")

    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')

    axes[0].set_ylim(min_sep, max_sep)
    axes[1].set_ylim(min_sep, max_sep)
    axes[0].set_xlim(min_eng, max_eng)
    axes[1].set_xlim(min_eng, max_eng)

    # also plot the analytic lines
    axes[0].plot(energy_bins_centers, D0_analytic_lengths, color = '#FEF3E8', alpha = 0.7)
    axes[1].plot(energy_bins_centers, Dp_analytic_lengths, label = r"$d = \frac{E \tau}{mc}$", color = '#FEF3E8', alpha = 0.7)

    legend = axes[1].legend(loc = 'upper left')
    for text in legend.get_texts():
        text.set_color('#FEF3E8')

    fig.savefig("./plots/Separation_Length_Distribution", bbox_inches = 'tight')

# plot_separation_distribution(sim)

def plot_2d_energy_loss(analysis_):
    E_D0, E_Dp, n_D0, n_Dp = analysis_.energy_loss_analysis_2d()
    energy_bins = np.logspace(np.log10(min(min(E_D0), min(E_Dp))), np.log10(max(max(E_D0), max(E_Dp))), 20)
    num_bins = np.linspace(-0.01, 7.99, 9)
    X, Y = np.meshgrid(energy_bins, num_bins)
    hist_D0, _, _ = np.histogram2d(E_D0, n_D0, bins = (energy_bins, num_bins))
    hist_Dp, _, _ = np.histogram2d(E_Dp, n_Dp, bins = (energy_bins, num_bins))

    hist_D0 = normalize(hist_D0, energy_bins, num_bins)
    hist_Dp = normalize(hist_Dp, energy_bins, num_bins)

    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (11, 5))
    im1 = axes[0].pcolor(X, Y, hist_D0.T, cmap="plasma", alpha = 0.7)
    im2 = axes[1].pcolor(X, Y, hist_Dp.T, cmap="plasma", alpha = 0.7)
    divider2 = make_axes_locatable(axes[1])
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax2, orientation='vertical', alpha = 0.7)
    axes[0].axvline(x = 53 * 1e3, color = '#FEF3E8', alpha = 0.7, label = r"$d_{D^0} = l_{D^0}$")
    axes[1].axvline(x = 22 * 1e3, color = '#FEF3E8', alpha = 0.7, label = r"$d_{D^+} = l_{D^+}$")

    axes[0].set_title(r"$D^0-p$ Collision")
    axes[1].set_title(r"$D^+-p$ Collision")

    axes[0].set_xlabel(r"$E_{D^0}$ [GeV]")
    axes[1].set_xlabel(r"$E_{D^+}$ [GeV]")

    axes[0].set_ylim(0, 8)
    axes[1].set_ylim(0, 8)

    axes[0].set_ylabel(r"$n_{\textrm{Elastic Collision}}$")
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    legend0 = axes[0].legend()
    legend1 = axes[1].legend()
    for text in legend0.get_texts():
        text.set_color('#FEF3E8')
    for text in legend1.get_texts():
        text.set_color('#FEF3E8')


    fig.savefig("./plots/Energy_loss_2d_Distribution", bbox_inches = 'tight')

plot_2d_energy_loss(sim)
exit(0)

def plot_free_path_distribution():
    D0_E_list, D0_free_path_list, Dp_E_list, Dp_free_path_list = analysis_.free_path_analysis()
    energy_bins = np.logspace(1.5, 9, 20)
    distance_bins = np.logspace(-3, np.log10(5000), 20)
    X, Y = np.meshgrid(energy_bins, distance_bins)
    hist_D0, _, _ = np.histogram2d(D0_E_list, D0_free_path_list, bins = (energy_bins, distance_bins))
    hist_Dp, _, _ = np.histogram2d(Dp_E_list, Dp_free_path_list, bins = (energy_bins, distance_bins))

    hist_D0 = normalize(hist_D0, energy_bins, distance_bins)
    hist_Dp = normalize(hist_Dp, energy_bins, distance_bins)

    energy_bins_centers = np.zeros((len(energy_bins) - 1,))
    for i in range(len(energy_bins_centers)):
        energy_bins_centers[i] = np.sqrt(energy_bins[i] * energy_bins[i + 1])
    D0_analytic_lengths = analytic_free_path(energy_bins_centers)
    Dp_analytic_lengths = analytic_free_path(energy_bins_centers)


    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (11, 5))
    im1 = axes[0].pcolor(X, Y, hist_D0.T, cmap="plasma", alpha = 0.7, vmin=0, vmax=1)
    im2 = axes[1].pcolor(X, Y, hist_Dp.T, cmap="plasma", alpha = 0.7, vmin=0, vmax=1)
    divider2 = make_axes_locatable(axes[1])
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax2, orientation='vertical', alpha = 0.7)

    axes[0].set_title(r"$D^0-p$ Free Path")
    axes[1].set_title(r"$D^+-p$ Free Path")

    axes[0].set_xlabel(r"$E_{D^0}$ [GeV]")
    axes[1].set_xlabel(r"$E_{D^+}$ [GeV]")

    axes[0].plot(energy_bins_centers, D0_analytic_lengths, color = '#FEF3E8', alpha = 0.7)
    axes[1].plot(energy_bins_centers, Dp_analytic_lengths, label = r"$l = \frac{m_{\textrm{ice}}}{\rho N_A \sigma(E)}$", color = '#FEF3E8', alpha = 0.7)

    # also plot the decay lengths to explain low energy increase
    D0_decay_analytic_lengths = analytic_decay_length(energy_bins_centers, t_D0, m_D0)
    Dp_decay_analytic_lengths = analytic_decay_length(energy_bins_centers, t_Dp, m_Dp)

    axes[0].plot(energy_bins_centers, D0_decay_analytic_lengths, label = r"$d = \frac{E \tau}{mc}$", color = '#A597B6', alpha = 0.7)
    axes[1].plot(energy_bins_centers, Dp_decay_analytic_lengths, color = '#A597B6', alpha = 0.7)

    axes[0].set_ylabel(r"$l_{\textrm{Free}}$")
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')

    axes[0].set_ylim(1e-3, 5000)
    axes[1].set_ylim(1e-3, 5000)

    legend0 = axes[0].legend(loc = 'upper left')
    for text in legend0.get_texts():
        text.set_color('#A597B6')

    legend1 = axes[1].legend(loc = 'upper left')
    for text in legend1.get_texts():
        text.set_color('#FEF3E8')

    fig.savefig("./plots/Free_Path_Distribution", bbox_inches = 'tight')
    return

plot_free_path_distribution()