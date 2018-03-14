# -*- coding: utf-8 -*-

from __future__ import division, print_function
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


class Plotter(object):
    def __init__(self,
                 raw=None,
                 data=None,
                 col_lab=None,
                 q2_bins=None,
                 x_bins=None,
                 z_func=None,
                 z_ids=None,
                 z_labs=None,
                 sub_xlabel=None,
                 sub_ylabel=None,
                 big_xlabel=None,
                 big_ylabel=None,
                 title=None,
                 raw_colors=None,
                 data_colors=None,
                 thy_colors=None,
                 subplot_kw=None,
                 gridspec_kw=None,
                 fig_kw=None,
                 bigax_kw=None,
                 raw_plot_kw=None,
                 data_plot_kw=None,
                 thy_plot_kw=None,
                 legend_kw=None
                 ):
        self._raw = None
        self._data = None
        self._col_lab = None
        self._q2_bins = None
        self._x_bins = None
        self._z_func = None
        self._z_ids = None
        self._subplot_kw = None
        self._gridspec_kw = None
        self._fig_kw = None
        self._bigax_kw = None
        self._raw_plot_kw = None
        self._data_plot_kw = None
        self._thy_plot_kw = None

        if raw is not None:
            self.raw = raw

        if data is not None:
            self.data = data

        if col_lab is not None:
            self.col_lab = col_lab

        self.nrows = None
        if q2_bins is not None:
            self.q2_bins = q2_bins

        self.ncols = None
        if x_bins is not None:
            self.x_bins = x_bins

        if z_func is not None:
            self.z_func = z_func

        if z_ids is not None:
            self.z_ids = z_ids

        if z_labs is not None:
            self.z_labs = z_labs

        self.sub_xlabel = sub_xlabel
        self.sub_ylabel = sub_ylabel
        self.big_xlabel = big_xlabel
        self.big_ylabel = big_ylabel
        self.title = title

        # Use tab20 colors as defaults
        tab20_colors = ("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
                        "#2ca02c", "#98df8a", "#d62728", "#ff9896",
                        "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
                        "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
                        "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")

        self.raw_colors = (tab20_colors[1::2] if (raw_colors is None)
                           else raw_colors)
        self.data_colors = (tab20_colors[1::2] if (data_colors is None)
                            else data_colors)
        self.thy_colors = (tab20_colors[::2] if (thy_colors is None)
                           else thy_colors)

        self.subplot_kw = subplot_kw
        self.gridspec_kw = gridspec_kw
        self.fig_kw = fig_kw
        self.bigax_kw = bigax_kw
        self.raw_plot_kw = raw_plot_kw
        self.data_plot_kw = data_plot_kw
        self.thy_plot_kw = thy_plot_kw

        self.legend_kw = legend_kw
        if self.legend_kw is None:
            self.legend_kw = {"loc": "upper left"}

        self._ax_bins = None  # Determined in self.plot()
        self._raw_slices = None  # Determined in self.plot()
        self._data_slices = None  # Determined in self.plot()

        self._fig = None  # Created in self.plot()
        self._axs = None  # Created in self.plot()
        self._bigax = None  # Created in self.plot()

    @property
    def raw(self):
        return self._raw

    @raw.setter
    def raw(self, value):
        self._raw = value.copy(deep=True)
        if "delta" not in self._raw.columns:
            if "sys_u" in self._raw.columns:
                self._raw["delta"] = np.sqrt(self._raw["stat_u"] ** 2 +
                                             self._raw["sys_u"] ** 2)
            else:
                self._raw["delta"] = self._raw["stat_u"]
        if "k(z)" not in self._raw.columns and self.z_func is not None:
            self._raw["k(z)"] = self._raw["z"].apply(self.z_func)
        if self.col_lab is not None:
            self._raw.sort_values(by=self.col_lab, inplace=True)

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value.copy(deep=True)
        if "delta" not in self._data.columns:
            if "sys_u" in self._data.columns:
                self._data["delta"] = np.sqrt(self._data["stat_u"] ** 2 +
                                              self._data["sys_u"] ** 2)
            else:
                self._data["delta"] = self._data["stat_u"]
        if "k(z)" not in self._data.columns and self.z_func is not None:
            self._data["k(z)"] = self._data["z"].apply(self.z_func)
        if self.col_lab is not None:
            self._data.sort_values(by=self.col_lab, inplace=True)

    @property
    def col_lab(self):
        return self._col_lab

    @col_lab.setter
    def col_lab(self, value):
        self._col_lab = value
        if self.raw is not None:
            self.raw.sort_values(by=self._col_lab, inplace=True)
        if self.data is not None:
            self.data.sort_values(by=self._col_lab, inplace=True)

    @property
    def q2_bins(self):
        return self._q2_bins

    @q2_bins.setter
    def q2_bins(self, value):
        self.nrows = len(value) + 1
        self._q2_bins = value

    @property
    def x_bins(self):
        return self._x_bins

    @x_bins.setter
    def x_bins(self, value):
        self.ncols = len(value) + 1
        self._x_bins = value

    @property
    def z_func(self):
        return self._z_func

    @z_func.setter
    def z_func(self, value):
        self._z_func = value
        if self.raw is not None and "k(z)" not in self.raw.columns:
            self.raw["k(z)"] = self.raw["z"].apply(self._z_func)
        if self.data is not None and "k(z)" not in self.data.columns:
            self.data["k(z)"] = self.data["z"].apply(self._z_func)

    @property
    def z_ids(self):
        return self._z_ids

    @z_ids.setter
    def z_ids(self, value):
        self._z_ids = value
        self.z_labs = {k: None for k in self._z_ids}

    @property
    def subplot_kw(self):
        return self._subplot_kw

    @subplot_kw.setter
    def subplot_kw(self, value):
        if value is None:
            self._subplot_kw = {"yscale": "log", "zorder": 1}  # Defaults
        else:
            self._subplot_kw = value
            if "zorder" not in self._subplot_kw:
                self._subplot_kw["zorder"] = 1

    @property
    def gridspec_kw(self):
        return self._gridspec_kw

    @gridspec_kw.setter
    def gridspec_kw(self, value):
        if value is None:
            self._gridspec_kw = {"wspace": 0.0, "hspace": 0.0}  # Defaults
        else:
            self._gridspec_kw = value

    @property
    def fig_kw(self):
        return self._fig_kw

    @fig_kw.setter
    def fig_kw(self, value):
        if value is None:
            self._fig_kw = {}  # Defaults
        else:
            self._fig_kw = value

    @property
    def bigax_kw(self):
        return self._bigax_kw

    @bigax_kw.setter
    def bigax_kw(self, value):
        if value is None:
            self._bigax_kw = {"zorder": 0}  # Defaults
        else:
            self._bigax_kw = value
            if "zorder" not in self._bigax_kw:
                self._bigax_kw["zorder"] = 0

    @property
    def raw_plot_kw(self):
        return self._raw_plot_kw

    @raw_plot_kw.setter
    def raw_plot_kw(self, value):
        if value is None:
            self._raw_plot_kw = {"marker": "o",
                                 "linestyle": "",
                                 "linewidth": 1,
                                 "markersize": 1,
                                 "alpha": 0.2,
                                 "zorder": 2
                                 }  # Defaults
        else:
            self._raw_plot_kw = value
            if "zorder" not in self._raw_plot_kw:
                self._raw_plot_kw["zorder"] = 2

    @property
    def data_plot_kw(self):
        return self._data_plot_kw

    @data_plot_kw.setter
    def data_plot_kw(self, value):
        if value is None:
            self._data_plot_kw = {"marker": "o",
                                  "linestyle": "",
                                  "linewidth": 1,
                                  "markersize": 1,
                                  "alpha": 0.8,
                                  "zorder": 3
                                  }  # Defaults
        else:
            self._data_plot_kw = value
            if "zorder" not in self._data_plot_kw:
                self._data_plot_kw["zorder"] = 3

    @property
    def thy_plot_kw(self):
        return self._thy_plot_kw

    @thy_plot_kw.setter
    def thy_plot_kw(self, value):
        if value is None:
            self._thy_plot_kw = {"zorder": 4}  # Defaults
        else:
            self._thy_plot_kw = value
            if "zorder" not in self._thy_plot_kw:
                self._thy_plot_kw["zorder"] = 4

    @property  # Read-only
    def ax_bins(self):
        # ax_bins maps (row, col) to Q2 and x bins
        return self._ax_bins

    @property  # Read-only
    def raw_slices(self):
        # raw_slices maps (row, col) to corresponding raw slices
        # raw_slices does not contain empty slices (very useful)
        return self._raw_slices

    @property  # Read-only
    def data_slices(self):
        # data_slices maps (row, col) to corresponding data slices
        return self._data_slices

    @property  # Read-only
    def fig(self):
        return self._fig

    @property  # Read-only
    def axs(self):
        return self._axs

    @property  # Read-only
    def bigax(self):
        return self._bigax

    def plot(self):
        # Preparation
        raw = self.raw
        data = self.data
        q2_bins = self.q2_bins
        x_bins = self.x_bins
        nrows = self.nrows
        ncols = self.ncols
        raw_colors = self.raw_colors
        nraw_colors = len(raw_colors)
        data_colors = self.data_colors
        ndata_colors = len(data_colors)
        thy_colors = self.thy_colors
        nthy_colors = len(thy_colors)

        self._ax_bins = {
            (i + 1, j + 1): {
                "Q2": (q2_bins[-i - 2], q2_bins[-i - 1]),
                "x": (x_bins[j], x_bins[j + 1])
            }
            for i in range(len(q2_bins) - 1)
            for j in range(len(x_bins) - 1)
        }
        ax_bins = self.ax_bins

        self._raw_slices = {
            (i, j):
                raw[(raw["Q2"] <= ax_bins[i, j]["Q2"][1]) &
                    (raw["x"] <= ax_bins[i, j]["x"][1]) &
                    ((ax_bins[i, j]["Q2"][0] < raw["Q2"]) if (i != nrows - 2)
                     else (ax_bins[i, j]["Q2"][0] <= raw["Q2"])) &
                    ((ax_bins[i, j]["x"][0] < raw["x"]) if (j != 1)
                     else (ax_bins[i, j]["x"][0] <= raw["x"]))]
            for i, j in ax_bins
        }

        for i, j in ax_bins:
            if len(self._raw_slices[i, j]) == 0:
                del self._raw_slices[i, j]

        raw_slices = self.raw_slices

        self._data_slices = {
            (i, j):
                data[(data["Q2"] <= ax_bins[i, j]["Q2"][1]) &
                     (data["x"] <= ax_bins[i, j]["x"][1]) &
                     ((ax_bins[i, j]["Q2"][0] < data["Q2"]) if (i != nrows - 2)
                      else (ax_bins[i, j]["Q2"][0] <= data["Q2"])) &
                     ((ax_bins[i, j]["x"][0] < data["x"]) if (j != 1)
                      else (ax_bins[i, j]["x"][0] <= data["x"]))]
            for i, j in ax_bins
        }
        data_slices = self.data_slices

        # Plotting
        self._fig, self._axs = plt.subplots(nrows, ncols,
                                            sharex="col", sharey="row",
                                            subplot_kw=self.subplot_kw,
                                            gridspec_kw=self.gridspec_kw,
                                            **self.fig_kw)

        fig = self.fig
        axs = self.axs

        self._bigax = fig.add_subplot(111, **self.bigax_kw)
        bigax = self.bigax

        bigax.spines["top"].set_visible(False)
        bigax.spines["right"].set_visible(False)

        bigax.set_xticks(np.linspace(0, 1, ncols + 1))
        bigax.set_yticks(np.linspace(0, 1, nrows + 1))

        bigax.set_xticklabels([""] + x_bins)
        bigax.set_yticklabels([""] + q2_bins)

        bigax.set_xlabel(self.big_xlabel)
        bigax.set_ylabel(self.big_ylabel, rotation="horizontal")

        bigax.set_title(self.title)

        for i in range(nrows):
            for j in range(ncols):
                ax = axs[i, j]

                if (i, j) not in raw_slices:  # If there's nothing to plot
                    ax.set_axis_off()
                    continue

                ax.tick_params(which="both",
                               direction="in",
                               labelsize="xx-small")

                # Remove unnecessary x-axis labels
                if (i + 1, j) in raw_slices:  # If not bottom
                    ax.tick_params(labelbottom=False)
                else:  # If bottom
                    ax.set_xlabel(self.sub_xlabel,
                                  horizontalalignment="center",
                                  fontsize="x-small")
                    plt.setp(ax.get_xticklabels(), visible=True)

                # Remove unnecessary y-axis labels
                if (i, j - 1) in raw_slices:  # If not leftmost
                    ax.tick_params(labelleft=False)
                else:  # If leftmost
                    ax.set_ylabel(self.sub_ylabel,
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
                for k in self.z_ids:
                    raw_z = raw_sub[raw_sub["k(z)"] == k]
                    data_z = data_sub[data_sub["k(z)"] == k]

                    # Raw data points
                    ax.errorbar(raw_z[self.col_lab],
                                raw_z["value"],
                                raw_z["delta"],
                                color=raw_colors[k % nraw_colors],
                                **self.raw_plot_kw)

                    # Selected data points
                    ax.errorbar(data_z[self.col_lab],
                                data_z["value"],
                                data_z["delta"],
                                color=data_colors[k % ndata_colors],
                                **self.data_plot_kw)

                    # Theory
                    ax.plot(data_z[self.col_lab],
                            data_z["thy"],
                            color=thy_colors[k % nthy_colors],
                            **self.thy_plot_kw)

                ax.relim()
                ax.autoscale_view()

        if set(self.z_labs.values()) != {None}:
            patches = [mpatches.Patch(color=thy_colors[k % nthy_colors],
                                      label=self.z_labs[k])
                       for k in self.z_ids]
            bigax.legend(handles=patches, **self.legend_kw)
