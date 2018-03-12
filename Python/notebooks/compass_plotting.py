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
z_cats = 4

raw["qT"] = raw["pT"] / raw["z"]
data["qT"] = data["pT"] / data["z"]

col_lab = "pT"
hor_lab = r"$p_T$ (GeV)"
# col_lab = "qT"
# hor_lab = r"$q_T$ (GeV)"

vert_lab = r"$M_{D}^{\pi^+}$"

# function (raw, data, col_lab, x_bins, q2_bins, z_func, z_cats,
#           hor_lab, vert_lab)
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
# tab20 colors are qualitative
# Every other color is a lighter version of the former

nrows = len(q2_bins) + 1
ncols = len(x_bins) + 1

subplot_kw = {"yscale": "log", "zorder": 1}
gridspec_kw = {"wspace": 0.0, "hspace": 0.0}

fig, axs = plt.subplots(nrows, ncols, sharex="col", sharey="row",
                        subplot_kw=subplot_kw, gridspec_kw=gridspec_kw)

bigax = fig.add_subplot(111, zorder=0)

bigax.spines["top"].set_visible(False)
bigax.spines["right"].set_visible(False)

bigax.set_xticks(np.linspace(0, 1, ncols + 1))
bigax.set_yticks(np.linspace(0, 1, nrows + 1))

bigax.set_xticklabels([""] + x_bins)
bigax.set_yticklabels([""] + q2_bins)

bigax.set_xlabel("$x$")
bigax.set_ylabel("$Q^2$", rotation="horizontal")

bigax.set_title("COMPASS")

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
# raw_slices does not contain empty slices (very useful)

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

        if (i, j) not in raw_slices:  # If there's nothing to plot
            ax.set_axis_off()
            continue

        ax.tick_params(direction="in",
                       labelsize="xx-small")

        # Remove unnecessary x-axis labels
        if (i + 1, j) in raw_slices:  # If not bottom
            ax.tick_params(labelbottom=False)
        else:  # If bottom
            ax.set_xlabel(hor_lab,
                          horizontalalignment="center",
                          fontsize="x-small")
            plt.setp(ax.get_xticklabels(), visible=True)

        # Remove unnecessary y-axis labels
        if (i, j - 1) in raw_slices:  # If not leftmost
            ax.tick_params(labelleft=False)
        else:  # If leftmost
            ax.set_ylabel(vert_lab,
                          rotation="horizontal",
                          verticalalignment="center",
                          fontsize="small",
                          labelpad=8)
            ax.tick_params(axis="y", pad=0.0)
            plt.setp(ax.get_yticklabels(), visible=True)

        q2_min, q2_max = ax_bins[i, j]["Q2"]
        x_min, x_max = ax_bins[i, j]["x"]

        raw_sub = raw_slices[i, j]

        data_sub = data_slices[i, j]

        # Plotting for each z bin
        for k in range(z_cats):
            raw_z = raw_sub[raw_sub["k(z)"] == k]
            data_z = data_sub[data_sub["k(z)"] == k]

            # Raw data points
            ax.errorbar(raw_z[col_lab],
                        raw_z["value"],
                        raw_z["delta"],
                        color=colors[2 * k + 1],
                        marker="o",
                        linestyle="",
                        linewidth=1,
                        markersize=1,
                        alpha=0.2,
                        zorder=2)

            # Selected data points
            ax.errorbar(data_z[col_lab],
                        data_z["value"],
                        data_z["delta"],
                        color=colors[2 * k + 1],
                        marker="o",
                        linestyle="",
                        linewidth=1,
                        markersize=1,
                        alpha=0.8,
                        zorder=3)

            # Theory
            ax.plot(data_z[col_lab],
                    data_z["thy"],
                    color=colors[2 * k],
                    zorder=4)

        ax.relim()
        ax.autoscale_view()
