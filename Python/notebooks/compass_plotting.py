# -*- coding: utf-8 -*-

from __future__ import division, print_function
import numpy as np
import pandas as pd
from matplotlib import rc
# import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from fitlab.resman import RESMAN
from tools.tools import load_config

# %matplotlib inline
plt.ion()

rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})

conf = load_config("../fitlab/inputs/upol_compass.py")
conf["resman"] = RESMAN(conf)
conf["resman"].get_residuals(conf["parman"].par)

# Bins
x_bins = [0.003, 0.008, 0.013, 0.02, 0.032, 0.055, 0.1, 0.21, 0.4]
q2_bins = [1.0, 1.7, 3.0, 7.0, 16.0, 81.0]

raw = pd.read_excel('../database/sidis/expdata/5001.xlsx')

data = pd.concat(pd.DataFrame(d)
                 for d in conf["resman"].sidisres.tabs.values())

data = data[data["hadron"] == "pi+"]

z_func = lambda z: [0.2, 0.3, 0.4, 0.6].index(z)
k_max = 4

raw["qT"] = raw["pT"] / raw["z"]
data["qT"] = data["pT"] / data["z"]

# hor_lab = "pT"
hor_lab = "qT"

# function (raw, data, hor_lab, x_bins, q2_bins, z_func, k_max)
# possibly don't use z_func if "k(z)" already in raw and data
raw = raw.copy(deep=True)
data = data.copy(deep=True)

if "delta" not in raw.columns:
    raw["delta"] = np.sqrt(raw["stat_u"] ** 2 + raw["sys_u"] ** 2)

if "delta" not in data.columns:
    data["delta"] = np.sqrt(data["stat_u"] ** 2 + data["sys_u"] ** 2)

raw["k(z)"] = raw["z"].apply(z_func)
data["k(z)"] = data["z"].apply(z_func)

# Plotting
colors = plt.get_cmap("tab20").colors
colors = list(colors[::2]) + list(colors[1::2])

nrows = len(q2_bins) + 1
ncols = len(x_bins) + 1

subplot_kw = {"yscale": "log"}
gridspec_kw = {"wspace": 0.0, "hspace": 0.0}

fig, axs = plt.subplots(nrows, ncols, sharex="row", sharey="col",
                        subplot_kw=subplot_kw, gridspec_kw=gridspec_kw)

ax_bins = {
    (i + 1, j + 1): {
        "Q2": (q2_bins[-i - 2], q2_bins[-i - 1]),
        "x": (x_bins[j], x_bins[j + 1])
    }
    for i in range(len(q2_bins) - 1)
    for j in range(len(x_bins) - 1)
}
# ax_bins maps row, col to Q2 and x bins

raw_slices = {
    (i, j):
        raw[(raw["Q2"] <= ax_bins[i, j]["Q2"][1]) &
            (raw["x"] <= ax_bins[i, j]["x"][1]) &
            ((ax_bins[i, j]["Q2"][0] < raw["Q2"]) if (i != nrows - 2) else
             (ax_bins[i, j]["Q2"][0] <= raw["Q2"])) &
            ((ax_bins[i, j]["x"][0] < raw["x"]) if (j != 1) else
             (ax_bins[i, j]["x"][0] <= raw["x"]))]
    for i, j in ax_bins
}

for i, j in ax_bins:
    if len(raw_slices[i, j]) == 0:
        del raw_slices[i, j]
# raw_slices maps row, col to corresponding raw data slices
# raw_slices does not contain empty slices

data_slices = {
    (i, j):
        data[(data["Q2"] <= ax_bins[i, j]["Q2"][1]) &
             (data["x"] <= ax_bins[i, j]["x"][1]) &
             ((ax_bins[i, j]["Q2"][0] < data["Q2"]) if (i != nrows - 2) else
              (ax_bins[i, j]["Q2"][0] <= data["Q2"])) &
             ((ax_bins[i, j]["x"][0] < data["x"]) if (j != 1) else
              (ax_bins[i, j]["x"][0] <= data["x"]))]
    for i, j in ax_bins
}
# data_slices maps row, col to corresponding data slices

for i in range(nrows):
    for j in range(ncols):
        ax = axs[i, j]

        # If there's nothing to plot
        if (i, j) not in raw_slices:
            ax.set_axis_off()
            continue

        # Remove unnecessary axis labels
        if (i + 1, j) in raw_slices:
            ax.tick_params(labelbottom=False)
        if (i, j - 1) in raw_slices:
            ax.tick_params(labelleft=False)

        q2_min, q2_max = ax_bins[i, j]["Q2"]
        x_min, x_max = ax_bins[i, j]["x"]

        raw_sub = raw_slices[i, j]

        data_sub = data_slices[i, j]

        for k in range(k_max):
            raw_z = raw_sub[raw_sub["k(z)"] == k]
            data_z = data_sub[data_sub["k(z)"] == k]

            ax.errorbar(raw_z[hor_lab],
                        raw_z["value"],
                        raw_z["delta"],
                        color=colors[k],
                        marker="o",
                        linestyle="",
                        linewidth=1,
                        markersize=1,
                        alpha=0.2)

            ax.errorbar(data_z[hor_lab],
                        data_z["value"],
                        data_z["delta"],
                        color=colors[k],
                        marker="o",
                        linestyle="",
                        linewidth=1,
                        markersize=1)

        ax.relim()
        ax.autoscale_view()
