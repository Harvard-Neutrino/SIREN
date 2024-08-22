import pandas as pd
import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

filename = "output/FullSim.parquet"

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

def mass(p):
    return np.sqrt(p[0] ** 2 - (p[1] ** 2 + p[2] ** 2 + p[3] ** 2))

def decay_length(v1, v2):
    return np.sqrt((v1[0] - v2[0]) ** 2  + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)

def extract_Q2(pe, pnu):
    return (pe[0] + pnu[0]) ** 2 - ((pe[1] + pnu[1]) ** 2 + (pe[2] + pnu[2]) ** 2 + (pe[3] + pnu[3]) ** 2)

def add_to_dict(list, dictionary):
    for item in list:
        if item in dictionary:
            dictionary[item] += 1
        else: dictionary[item] = 1

class event:
    def __init__(self, row) -> None:
        self.row = row
        self.num_interaction = row["num_interactions"]
    
    def DTYPE(self) -> int:
        # type of charmed meson is the second one
        return int(self.row["secondary_types"][1][1])
    
    def D_DECAY_LEN(self) -> float:
        # the first vertex is 0th, the decay vertex is the last one
        return decay_length(self.row["vertex"][0], self.row["vertex"][-1])

class analysis:
    def __init__(self, f) -> None:
        self.df = pd.read_parquet(f)
        self.num_events = len(self.df["event_weight"])
        # print("Initializing... There are {} events".format(f.attrs["num_events"]))
        self.num_interactions = {}
        self.secondary_types = {}
    
    def separation_analysis(self):
        D0_energies = []
        D0_separations = []
        Dp_energies = []
        Dp_separations = []
        D0_weights = []
        Dp_weights = []
        for i in range(self.num_events):
            cur_event = event(self.df.iloc[i])
            print("{}/{}".format(i, self.num_events), end = '\r')
            # check if current event is D0 or D+
            if cur_event.DTYPE() == 421: # This is D0
                # extract the vertex separations
                D0_separations.append(cur_event.D_DECAY_LEN())
                D0_energies.append(cur_event.row["primary_momentum"][2][0])
                D0_weights.append(cur_event.row["event_weight"])

            elif cur_event.DTYPE() == 411: # This is D+
                # extract the vertex separations
                Dp_separations.append(cur_event.D_DECAY_LEN())
                Dp_energies.append(cur_event.row["primary_momentum"][2][0])
                Dp_weights.append(cur_event.row["event_weight"])

        return D0_energies, D0_separations, Dp_energies, Dp_separations, D0_weights, Dp_weights

    def energy_loss_analysis_2d(self):
        E_D0 = []
        E_Dp = []
        n_D0 = []
        n_Dp = []
        w_D0 = []
        w_Dp = []
        for i in range(self.num_events):
            cur_event = event(self.df.iloc[i])
            print("{}/{}".format(i, self.num_events), end = '\r')
            if cur_event.DTYPE() == 421: # This is D0
                # extract the vertex separations
                n_D0.append(cur_event.row["num_interactions"] - 3)
                E_D0.append(cur_event.row["primary_momentum"][2][0])
                w_D0.append(cur_event.row["event_weight"])

            elif cur_event.DTYPE() == 411: # This is D+
                # extract the vertex separations
                n_Dp.append(cur_event.row["num_interactions"] - 3)
                E_Dp.append(cur_event.row["primary_momentum"][2][0])
                w_Dp.append(cur_event.row["event_weight"])

        return E_D0, E_Dp, n_D0, n_Dp, w_D0, w_Dp
        
