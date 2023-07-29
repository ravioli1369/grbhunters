#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
@File    :   final_script.py
@Time    :   2023/07/25 05:58:45
@Author  :   Ravi K.
@Desc    :   The final functions and script to run the whole project
"""

import os
import glob
import time
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.stats import sigma_clipped_stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from pipelinev3 import bindata

import warnings

warnings.simplefilter("ignore", np.RankWarning)

# Using LaTeX fonts
params = {
    "text.usetex": True,
    "font.family": "serif",
    "figure.dpi": 150,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.top": True,
    "ytick.left": True,
    "ytick.right": True,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.minor.size": 2.5,
    "xtick.major.size": 5,
    "ytick.minor.size": 2.5,
    "ytick.major.size": 5,
    "axes.axisbelow": True,
}
matplotlib.rcParams.update(params)
##############################################################################################
############ Functions for detrending, outlier detection, snr, energy binning ################
##############################################################################################


def saav2(a):
    """
    Return maximum runs of consecutive zeros in a 1D array.
    Equivalent to SAA start and end indices for lc created with v2 pipeline
    """
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    no_of_zeros = ranges[:, 1] - ranges[:, 0]
    if len(no_of_zeros) > 0:
        max_zeros = np.where(no_of_zeros == np.max(no_of_zeros))[0][0]
        return ranges[max_zeros]
    else:
        return [len(a)]


def saav3(a):
    """
    Returns the start and end indices of the South Atlantic Anomaly
    for lcs created with v3 pipeline
    """
    diff = np.diff(a)
    saa_start = np.argmax(diff)
    saa_end = saa_start + 1
    return saa_start, saa_end


def get_saa_indices(data):
    """
    Returns the final SAA start and end indices regardless of the pipeline used
    """
    if saav2(data["RATE"])[-1] == len(data["RATE"]):
        saa_start, saa_end = saav3(data["TIME"])
    else:
        length_saa_zeroes = (
            data["TIME"][saav2(data["RATE"])[-1]] - data["TIME"][saav2(data["RATE"])[0]]
        )
        length_saa_nans = (
            data["TIME"][saav3(data["TIME"])[-1]] - data["TIME"][saav3(data["TIME"])[0]]
        )
        if length_saa_zeroes > length_saa_nans:
            saa_start, saa_end = saav2(data["RATE"])
        else:
            saa_start, saa_end = saav3(data["TIME"])
    return saa_start, saa_end


def quadratic(x, a, b, c):
    return a * x**2 + b * x + c


def get_trigger_index(filename, trigger_time):
    """
    Returns the index of the trigger in the given file
    """
    data = fits.getdata(filename)
    diffs = np.abs(
        np.round(data["TIME"], 0).astype(int) - np.round(trigger_time, 0).astype(int)
    )
    trigger_index = np.where(diffs == np.min(diffs))[0][-1]
    return trigger_index


def quadratic_detrend_trigger(
    filename, trigger_index, polyorder=2, detrend_window=21, data=None
):
    """
    Detrends the given file using a quadratic fit to the data around the trigger
    """
    if data is None:
        data = fits.getdata(filename)
    saa_start, saa_end = get_saa_indices(data)
    timebin = data["TIME"][trigger_index + 1] - data["TIME"][trigger_index]
    detrend_window = (
        np.rint(detrend_window / timebin).astype(int) // 2 * 2 + 1
    )  # make it odd
    background_window = np.rint(500 / timebin).astype(int)
    if trigger_index < saa_start:
        if trigger_index > background_window:
            if trigger_index + background_window < saa_start:
                counts = data["RATE"][
                    trigger_index
                    - background_window : trigger_index
                    + background_window
                ]
                times = data["TIME"][
                    trigger_index
                    - background_window : trigger_index
                    + background_window
                ]
                new_trigger_index = background_window
            elif trigger_index + background_window >= saa_start:
                counts = data["RATE"][trigger_index - background_window : saa_start]
                times = data["TIME"][trigger_index - background_window : saa_start]
                new_trigger_index = background_window
        elif trigger_index <= background_window:
            if trigger_index + background_window < saa_start:
                counts = data["RATE"][: trigger_index + background_window]
                times = data["TIME"][: trigger_index + background_window]
                new_trigger_index = trigger_index
            elif trigger_index + background_window >= saa_start:
                counts = data["RATE"][:saa_start]
                times = data["TIME"][:saa_start]
                new_trigger_index = trigger_index
        else:
            raise ValueError("Please check manually, something is wrong")
    elif trigger_index > saa_end:
        if trigger_index - saa_end > background_window:
            if trigger_index + background_window < len(data):
                counts = data["RATE"][
                    trigger_index
                    - background_window : trigger_index
                    + background_window
                ]
                times = data["TIME"][
                    trigger_index
                    - background_window : trigger_index
                    + background_window
                ]
                new_trigger_index = background_window
            elif trigger_index + background_window >= len(data):
                counts = data["RATE"][trigger_index - background_window :]
                times = data["TIME"][trigger_index - background_window :]
                new_trigger_index = background_window
        elif trigger_index - saa_end <= background_window:
            if trigger_index + background_window < len(data):
                counts = data["RATE"][saa_end : trigger_index + background_window]
                times = data["TIME"][saa_end : trigger_index + background_window]
                new_trigger_index = trigger_index - saa_end
            elif trigger_index + background_window >= len(data):
                counts = data["RATE"][saa_end:]
                times = data["TIME"][saa_end:]
                new_trigger_index = trigger_index - saa_end
        else:
            raise ValueError("Please check manually, something is wrong")
    else:
        raise ValueError("Trigger index is in SAA")
    # clipping the outliers before fitting the quadratic
    mean, _, std = sigma_clipped_stats(counts)
    counts = np.copy(counts)
    counts[np.abs(counts - mean) > 3 * std] = np.nan
    # filtered = counts
    if detrend_window > 2:
        filtered = savgol_filter(counts, detrend_window, 2)
    else:
        filtered = np.copy(counts)
    idx = np.isfinite(filtered)
    x = times[idx]
    y = filtered[idx]
    popt = np.polyfit(x, y, polyorder)
    # popt, _ = curve_fit(quadratic, x, y)
    window_start, window_end = (
        np.where(data["TIME"] == times[0])[0][0],
        np.where(data["TIME"] == times[-1])[0][0],
    )
    # putting back the outliers before detrending
    counts = data["RATE"][window_start : window_end + 1]
    # detrending
    # trend = quadratic(times, *popt)
    trend = np.polyval(popt, times)
    detrended_counts = counts - trend
    detrended = QTable([times, detrended_counts], names=("TIME", "RATE"))
    raw = QTable([times, counts], names=("TIME", "RATE"))
    return detrended, raw, trend, filtered, new_trigger_index, popt


def create_master_lc(directory, timebin=1):
    """
    Creates the master light curve (20-200 keV) for the given directory
    """
    gti = glob.glob(f"{directory}/*bc.gti")[0]
    evt = glob.glob(f"{directory}/*bc.evt")[0]
    mkf = glob.glob(f"{directory}/*.mkf")[0]
    if not os.path.exists(f"{directory}/{timebin}s/master_lc"):
        print(
            "Creating master light curve at {}/{}s/master_lc".format(directory, timebin)
        )
        os.makedirs(f"{directory}/{timebin}s/master_lc")
        if not os.path.exists(gti):
            os.system("python3 pipelinev3.py -d {} -time {}".format(directory, timebin))
        else:
            bindata(evt, mkf, timebin, 20, 200)
        print(
            "\n\n Moving light curves to {}/{}s/master_lc \n\n".format(
                directory, timebin
            )
        )
        os.system("mv {}/*.lc {}/{}s/master_lc".format(directory, directory, timebin))
    else:
        print(
            "Master light curve already exists at {}/{}s/master_lc".format(
                directory, timebin
            )
        )


def gen_energy_bins(directory, timebin=1):
    """
    Generates energy bins for the given number of bins
    """
    n_bins = 3
    energy_ranges = [20, 60, 100, 200]
    lc_paths = []
    if not os.path.exists(f"{directory}/{timebin}s/{n_bins}_bins"):
        print(
            "{}s: Creating {} energy bins at {}/{}s/{}_bins".format(
                timebin, n_bins, directory, timebin, n_bins
            )
        )
        os.makedirs(f"{directory}/{timebin}s/{n_bins}_bins")
        for i in range(n_bins):
            emin = energy_ranges[i]
            emax = energy_ranges[i + 1]
            evt = glob.glob(f"{directory}/*bc.evt")[0]
            mkf = glob.glob(f"{directory}/*.mkf")[0]
            print(emin, emax)
            bindata(evt, mkf, timebin, emin, emax)
            if not os.path.exists(
                f"{directory}/{timebin}s/{n_bins}_bins/{int(emin)}-{int(emax)}"
            ):
                os.makedirs(
                    f"{directory}/{timebin}s/{n_bins}_bins/{int(emin)}-{int(emax)}"
                )
            print(
                f"\n\n Moving light curves to {directory}/{timebin}s/{n_bins}_bins/{int(emin)}-{int(emax)}/\n\n"
            )
            os.system(
                f"mv {directory}/*.lc {directory}/{timebin}s/{n_bins}_bins/{int(emin)}-{int(emax)}/"
            )
            lc_paths.append(
                f"{directory}/{timebin}s/{n_bins}_bins/{int(emin)}-{int(emax)}/"
            )
    else:
        print(
            "{}s: Requested energy bins already exist at {}/{}s/{}_bins".format(
                timebin, directory, timebin, n_bins
            )
        )
        for i in range(n_bins):
            emin = energy_ranges[i]
            emax = energy_ranges[i + 1]
            lc_paths.append(
                f"{directory}/{timebin}s/{n_bins}_bins/{int(emin)}-{int(emax)}/"
            )
    return lc_paths


def outlier(filename, trigger_index, detection_sigma=3):
    """
    Returns the indices of the outliers in the given file
    """
    t, *_ = quadratic_detrend_trigger(filename, trigger_index)
    mean, _, std = sigma_clipped_stats(t["RATE"])
    outliers = np.where(t["RATE"] > mean + detection_sigma * std)[0]
    return outliers


def snr_outlier(filename1, filename2, trigger_index, detection_sigma=3):
    """
    Returns the SNR of the outliers in the given file
    """
    outliers = outlier(filename1, trigger_index, detection_sigma)
    t, *_ = quadratic_detrend_trigger(filename2, trigger_index)
    mean, _, std = sigma_clipped_stats(t["RATE"])
    signal = t["RATE"][outliers] + mean
    noise = mean + std
    snr = signal / noise
    return snr, outliers


def find_outliers(directory, trigger_time, timebin=1, detection_sigma=3):
    """
    Finds the outliers and potential GRBs in each quadrant and returns their SNRs and indices
    """

    def each_quad(lc_paths, trigger_index, quadrant):
        """
        Returns the SNR and indices of the outliers in each quadrant
        """
        snr = []
        lcs = []
        for path in lc_paths:
            lcs.append(glob.glob(f"{path}/*{str(quadrant)}.lc"))
        for i in range(3):
            snr.append(
                snr_outlier(lcs[0][0], lcs[i][0], trigger_index, detection_sigma)[0]
            )
        return snr, snr_outlier(lcs[0][0], lcs[0][0], trigger_index, detection_sigma)[1]

    def snr_grb(master_lc, filtered_outliers, trigger_index):
        """
        Returns the SNR of the potential GRBs in the master light curve
        """
        t, *_ = quadratic_detrend_trigger(master_lc, trigger_index)
        mean, _, std = sigma_clipped_stats(t["RATE"])
        noise = mean + std
        signal = t["RATE"][filtered_outliers]
        snr = signal / noise
        return snr

    create_master_lc(directory, timebin)
    master_lcs = np.sort(glob.glob(f"{directory}/{timebin}s/master_lc/*.lc"))
    lc_paths = gen_energy_bins(directory, timebin)
    results = []
    for i in range(4):
        trigger_index = get_trigger_index(master_lcs[i], trigger_time)
        snr_outliers_each_quad, outliers_each_quad = each_quad(
            lc_paths, trigger_index, i
        )
        master_snr_each_quad = snr_grb(master_lcs[i], outliers_each_quad, trigger_index)
        filtered_outliers_mask_each_quad = np.logical_and(
            np.logical_or(snr_outliers_each_quad[1] > 3, snr_outliers_each_quad[2] > 3),
            master_snr_each_quad > 3,
        )
        filtered_outliers_each_quad = outliers_each_quad[
            filtered_outliers_mask_each_quad
        ]
        filtered_outliers_snr_each_quad = snr_grb(
            master_lcs[i], filtered_outliers_each_quad, trigger_index
        )
        results.append(
            [
                snr_outliers_each_quad,
                outliers_each_quad,
                filtered_outliers_mask_each_quad,
                filtered_outliers_snr_each_quad,
            ]
        )
    results.append(master_lcs)
    results.append(lc_paths)
    return results


def find_potential_grbs(
    master_lcs, lc_paths, trigger_time, results, timebin, grb_name, plot=False
):
    u = potential_grb_times(master_lcs, trigger_time, results, timebin)
    counter = 0
    potential_grb_snr = [0, 0, 0, 0]
    potential_grb_time = [0, 0, 0, 0]
    figs = []
    for i in range(4):
        trigger_index = get_trigger_index(master_lcs[i], trigger_time)
        detrended, *_ = quadratic_detrend_trigger(
            master_lcs[i], trigger_index, polyorder=2
        )
        _, outliers, filtered_outliers_mask, filtered_outliers_snr = results[i]
        filtered_outliers = np.zeros_like(detrended["RATE"])
        filtered_outliers[outliers[filtered_outliers_mask]] = filtered_outliers_snr
        potential_grbs = np.intersect1d(
            u, detrended["TIME"][outliers[filtered_outliers_mask]]
        )
        if len(potential_grbs) > 0:
            matched_times_mask = np.isin(detrended["TIME"], potential_grbs)
            u_mask = np.isin(u, potential_grbs)
            potential_grb_time[i] = u[u_mask][
                np.argmax(filtered_outliers[matched_times_mask])
            ]
            potential_grb_snr[i] = np.max(filtered_outliers[matched_times_mask])
            print(
                f"Potential GRB found in Quadrant {i} at {potential_grb_time[i]}s with SNR {np.round(potential_grb_snr[i], 2)}!!!!"
            )
            counter += 1
        if i == 3:
            if counter > 1:
                if plot:
                    figs = plot_a_bunch_of_stuff(
                        master_lcs,
                        lc_paths,
                        results,
                        u,
                        grb_name,
                        trigger_time,
                        timebin,
                    )
                print(f"Potential GRB found for trigger time {trigger_time}s.")
            else:
                print(f"No Potential GRB found for trigger time {trigger_time}s.")
    return potential_grb_snr, potential_grb_time, figs


def potential_grb_times(master_lcs, trigger_time, results, timebin=1):
    outlier_times = []
    for i in range(4):
        trigger_index = get_trigger_index(master_lcs[i], trigger_time)
        detrended, *_ = quadratic_detrend_trigger(
            master_lcs[i], trigger_index, polyorder=2
        )
        _, outliers, grb_mask, _ = results[i]
        quad_outlier_times = detrended["TIME"][outliers[grb_mask]]
        quad_outlier_times = np.concatenate(
            (
                quad_outlier_times,
                quad_outlier_times + timebin,
                quad_outlier_times - timebin,
            )
        )
        outlier_times = np.concatenate((outlier_times, np.unique(quad_outlier_times)))
    u, c = np.unique(outlier_times, return_counts=True)
    u = u[c > 1]
    return u


def plot_a_bunch_of_stuff(
    master_lcs, lc_paths, results, u, grb_name, trigger_time, timebin
):
    fig_raw, ax_raw = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=True)
    fig_raw.set_tight_layout(True)
    fig_detrended, ax_detrended = plt.subplots(
        2, 2, figsize=(15, 10), sharex=True, sharey=True
    )
    fig_detrended.set_tight_layout(True)
    fig_mark_outlier, ax_mark_outlier = plt.subplots(
        2, 2, figsize=(15, 10), sharex=True, sharey=True
    )
    fig_mark_outlier.set_tight_layout(True)
    fig_snr_outlier, ax_snr_outlier = plt.subplots(
        4, 2, figsize=(15, 10), sharex=True, sharey=True
    )
    fig_snr_outlier.set_tight_layout(True)
    fig_snrvsenergy, ax_snrvsenergy = plt.subplots(
        2, 2, figsize=(15, 10), sharex=True, sharey=True
    )
    fig_snrvsenergy.set_tight_layout(True)
    fig_marked_grb, ax_marked_grb = plt.subplots(
        4, 2, figsize=(15, 10), sharex=True, sharey=True
    )
    fig_marked_grb.set_tight_layout(True)
    figs = [
        fig_raw,
        fig_detrended,
        fig_mark_outlier,
        fig_snr_outlier,
        fig_snrvsenergy,
        fig_marked_grb,
    ]
    counter = 0
    for i in range(4):
        trigger_index = get_trigger_index(master_lcs[i], trigger_time)
        detrended, raw, trend, *_ = quadratic_detrend_trigger(
            master_lcs[i], trigger_index, polyorder=2
        )
        detrended20to60, *_ = quadratic_detrend_trigger(
            lc_paths[i], trigger_index, polyorder=2
        )
        snr, outliers, grb_mask, grb_snr = results[i]
        ax_raw[i // 2, i % 2].plot(
            raw["TIME"],
            raw["RATE"],
            color="slateblue",
            label="Raw Count Rate",
            alpha=0.85,
        )
        ax_raw[i // 2, i % 2].plot(
            raw["TIME"], trend, color="salmon", label="Trend", linewidth=2
        )
        ax_raw[i // 2, i % 2].fill_between(
            raw["TIME"], 0, raw["RATE"], color="slateblue", alpha=0.2
        )
        ax_raw[i // 2, i % 2].set_xlim(raw["TIME"][0], raw["TIME"][-1])
        ax_raw[i // 2, i % 2].set_xlabel("Time (s)")
        ax_raw[i // 2, i % 2].set_ylabel("Count Rate (counts/s)")
        ax_raw[i // 2, i % 2].set_title("Quadrant {}".format(i))
        if i == 1:
            ax_raw[i // 2, i % 2].legend()

        ax_detrended[i // 2, i % 2].plot(
            detrended["TIME"],
            detrended["RATE"],
            color="salmon",
            label="Detrended Count Rate",
        )
        ax_detrended[i // 2, i % 2].set_xlim(
            detrended["TIME"][0], detrended["TIME"][-1]
        )
        ax_detrended[i // 2, i % 2].set_xlabel("Time (s)")
        ax_detrended[i // 2, i % 2].set_ylabel("Count Rate (counts/s)")
        ax_detrended[i // 2, i % 2].set_title("Quadrant {}".format(i))
        if i == 1:
            ax_detrended[i // 2, i % 2].legend()

        ax_mark_outlier[i // 2, i % 2].plot(
            detrended20to60["TIME"],
            detrended20to60["RATE"],
            color="slateblue",
            label="20-60 keV",
        )
        ax_mark_outlier[i // 2, i % 2].scatter(
            detrended20to60["TIME"][outliers],
            detrended20to60["RATE"][outliers],
            color="red",
            alpha=0.6,
            s=10 * snr[0],
            label="Outliers",
        )
        ax_mark_outlier[i // 2, i % 2].set_xlim(
            detrended20to60["TIME"][0], detrended20to60["TIME"][-1]
        )
        ax_mark_outlier[i // 2, i % 2].set_xlabel("Time (s)")
        ax_mark_outlier[i // 2, i % 2].set_ylabel("Count Rate (counts/s)")
        ax_mark_outlier[i // 2, i % 2].set_title("Quadrant {}".format(i))
        if i == 1:
            ax_mark_outlier[i // 2, i % 2].legend()
        grb = np.zeros_like(detrended["RATE"])
        grb[outliers] = snr[0]
        ax_snr_outlier[i, 0].plot(
            detrended["TIME"], grb, alpha=0.6, color="red", label="20-60keV"
        )
        grb[outliers] = snr[1]
        ax_snr_outlier[i, 0].plot(
            detrended["TIME"], grb, alpha=0.6, color="blue", label="60-100keV"
        )
        grb[outliers] = snr[2]
        ax_snr_outlier[i, 0].plot(
            detrended["TIME"], grb, alpha=0.6, color="green", label="100-200keV"
        )
        ax_snr_outlier[i, 0].set_xlim(detrended["TIME"][0], detrended["TIME"][-1])
        ax_snr_outlier[i, 0].set_title("Quadrant {} Outliers".format(i))
        ax_snr_outlier[i, 0].set_xlabel("Outliers")
        ax_snr_outlier[i, 0].set_ylabel("SNR")
        grb[outliers] = 0
        grb[outliers[grb_mask]] = grb_snr
        ax_snr_outlier[i, 1].plot(
            detrended["TIME"], grb, alpha=0.6, color="slateblue", label="20-200keV"
        )
        ax_snr_outlier[i, 1].set_xlim(detrended["TIME"][0], detrended["TIME"][-1])
        ax_snr_outlier[i, 1].set_title("Quadrant {} Potential GRBs".format(i))
        ax_snr_outlier[i, 1].set_xlabel("Outliers")
        ax_snr_outlier[i, 1].set_ylabel("SNR")
        eranges = [40, 80, 150]
        for j in range(len(snr[0])):
            if j == 0:
                ax_snrvsenergy[i // 2, i % 2].plot(
                    eranges,
                    [snr[0][j], snr[1][j], snr[2][j]],
                    color="slateblue",
                    marker="o",
                    markersize=5,
                    alpha=0.05,
                    label="Outliers",
                    linestyle="--",
                )
            else:
                ax_snrvsenergy[i // 2, i % 2].plot(
                    eranges,
                    [snr[0][j], snr[1][j], snr[2][j]],
                    color="slateblue",
                    marker="o",
                    markersize=5,
                    alpha=0.05,
                    linestyle="--",
                )
        for j in range(len(snr[0][grb_mask])):
            if j == 0:
                ax_snrvsenergy[i // 2, i % 2].plot(
                    eranges,
                    [snr[0][grb_mask][j], snr[1][grb_mask][j], snr[2][grb_mask][j]],
                    color="tomato",
                    marker="o",
                    markersize=5,
                    alpha=0.3,
                    label="Filtered Outliers",
                    linestyle="--",
                )
            else:
                ax_snrvsenergy[i // 2, i % 2].plot(
                    eranges,
                    [snr[0][grb_mask][j], snr[1][grb_mask][j], snr[2][grb_mask][j]],
                    color="tomato",
                    marker="o",
                    markersize=5,
                    alpha=0.3,
                    linestyle="--",
                )
        potential_grbs = np.intersect1d(u, detrended["TIME"][outliers[grb_mask]])
        ax_marked_grb[i, 0].plot(
            raw["TIME"],
            raw["RATE"],
            color="slateblue",
            label="Raw Count Rate",
            alpha=0.85,
        )
        ax_marked_grb[i, 0].fill_between(
            raw["TIME"], 0, raw["RATE"], color="slateblue", alpha=0.2
        )
        ax_marked_grb[i, 1].plot(
            detrended["TIME"],
            detrended["RATE"],
            color="tomato",
            label="Detrended Count Rate",
            alpha=0.85,
        )
        ax_marked_grb[i, 0].set_xlabel("Time (s)")
        ax_marked_grb[i, 0].set_ylabel("Count Rate (counts/s)")
        ax_marked_grb[i, 1].set_xlabel("Time (s)")
        ax_marked_grb[i, 1].set_ylabel("Count Rate (counts/s)")
        ax_marked_grb[i, 0].set_title(
            "Quadrant {} Potential GRBs - Raw Counts".format(i)
        )
        ax_marked_grb[i, 1].set_title(
            "Quadrant {} Potential GRBs - Detrended Counts".format(i)
        )
        if len(potential_grbs) > 0:
            counter += 1
            matched_times_mask = np.isin(detrended["TIME"], potential_grbs)
            u_mask = np.isin(u, potential_grbs)
            snr_potential_grbs = []
            grb[outliers] = snr[0]
            snr_potential_grbs.append(grb[matched_times_mask])
            grb[outliers] = snr[1]
            snr_potential_grbs.append(grb[matched_times_mask])
            grb[outliers] = snr[2]
            snr_potential_grbs.append(grb[matched_times_mask])
            grb[outliers] = 0
            grb[outliers[grb_mask]] = grb_snr
            for j in range(len(snr_potential_grbs[0])):
                if j == 0:
                    ax_snrvsenergy[i // 2, i % 2].plot(
                        eranges,
                        [
                            snr_potential_grbs[0][j],
                            snr_potential_grbs[1][j],
                            snr_potential_grbs[2][j],
                        ],
                        color="forestgreen",
                        marker="o",
                        markersize=5,
                        alpha=1,
                        label="Potential GRBs",
                        linestyle="--",
                    )
                else:
                    ax_snrvsenergy[i // 2, i % 2].plot(
                        eranges,
                        [
                            snr_potential_grbs[0][j],
                            snr_potential_grbs[1][j],
                            snr_potential_grbs[2][j],
                        ],
                        color="forestgreen",
                        marker="o",
                        markersize=5,
                        alpha=1,
                        linestyle="--",
                    )
            ax_snr_outlier[i, 1].scatter(
                u[u_mask],
                grb[matched_times_mask],
                color="forestgreen",
                alpha=0.6,
                s=2 * grb[matched_times_mask],
                label="Potential GRBs",
            )
            ax_marked_grb[i, 0].scatter(
                raw["TIME"][matched_times_mask],
                raw["RATE"][matched_times_mask],
                color="forestgreen",
                alpha=0.6,
                s=2 * grb[matched_times_mask],
                label="Potential GRB",
            )
            ax_marked_grb[i, 1].scatter(
                detrended["TIME"][matched_times_mask],
                detrended["RATE"][matched_times_mask],
                color="forestgreen",
                alpha=0.6,
                s=2 * grb[matched_times_mask],
                label="Potential GRB",
            )
            if counter == 1:
                ax_snr_outlier[i, 0].legend()
                ax_snr_outlier[i, 1].legend()
                ax_marked_grb[i, 0].legend()
                ax_marked_grb[i, 1].legend()
            ax_snrvsenergy[i // 2, i % 2].legend(loc="best")
        ax_snrvsenergy[i // 2, i % 2].set_title("Quadrant {}".format(i))
        ax_snrvsenergy[i // 2, i % 2].set_xlabel("Energy (keV)")
        ax_snrvsenergy[i // 2, i % 2].set_ylabel("SNR")

    fig_detrended.suptitle(f"Detrended Count Rate for {grb_name} - {timebin}s Binsize")
    fig_raw.suptitle(f"Raw Count Rate and Trend for {grb_name} - {timebin}s Binsize")
    fig_mark_outlier.suptitle(
        f"Detrended Count Rate + Outliers for {grb_name} - {timebin}s Binsize"
    )
    fig_snr_outlier.suptitle(f"SNR vs Outliers for {grb_name} - {timebin}s Binsize")
    fig_snrvsenergy.suptitle(f"SNR vs Energy for {grb_name} - {timebin}s Binsize")
    fig_marked_grb.suptitle(
        f"Potential GRB Detections for {grb_name} - {timebin}s Binsize"
    )
    # pdf = PdfPages(f"output_for_{grb_name}.pdf")
    # for fig in figs:
    #     pdf.savefig(fig)
    #     plt.close()
    # pdf.close()
    return figs


def run_timebins(
    directory, trigger_time, grb_name, timebin, detection_sigma=1, plot=False
):
    results = find_outliers(directory, trigger_time, timebin, detection_sigma)
    master_lcs = results[4]
    lc_paths = np.sort(glob.glob(f"{results[5][0]}/*.lc"))

    potential_grb_snrs, potential_grb_times, figs = find_potential_grbs(
        master_lcs, lc_paths, trigger_time, results, timebin, grb_name, plot
    )
    return potential_grb_snrs, potential_grb_times, figs


def main(directory, trigger_time, grb_name, input_timebin=None, detection_sigma=1):
    if input_timebin is not None:
        potential_grb_snrs, potential_grb_times, figs = run_timebins(
            directory,
            trigger_time,
            grb_name,
            input_timebin,
            detection_sigma,
            plot=True,
        )
        if len(np.nonzero(potential_grb_snrs)[0]) > 0:
            pdf = PdfPages(f"output_for_{grb_name}.pdf")
            for fig in figs:
                pdf.savefig(fig)
                plt.close()
            pdf.close()
            quadrants = np.nonzero(potential_grb_snrs)[0]
            for i in quadrants:
                print(
                    f"Potential GRB found in Quadrant {i} at {potential_grb_times[i]}s with SNR {np.round(potential_grb_snrs[i], 2)}!!!!"
                )
    else:
        timebins = [
            0.2,
            0.4,
            0.8,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
        ]
        snrs = []
        for timebin in timebins:
            snr = np.array(
                run_timebins(
                    directory, trigger_time, grb_name, timebin, detection_sigma
                )[0]
            )
            snr = snr[snr > 0]
            if len(snr) > 0:
                snrs.append(np.mean(snr))
            else:
                snrs.append(0)
        if len(np.nonzero(snrs)[0]) > 0:
            optimal_timebin = timebins[np.argmax(snrs)]
            print("\nOptimal timebin found!!!")
            print(
                "Generating plots for the optimal timebin of {}s\n".format(
                    optimal_timebin
                )
            )
            potential_grb_snrs, potential_grb_times, figs = run_timebins(
                directory,
                trigger_time,
                grb_name,
                optimal_timebin,
                detection_sigma,
                plot=True,
            )

            fig_snrvstime, ax_snrvstime = plt.subplots(
                figsize=(15, 10), sharex=True, sharey=True
            )
            ax_snrvstime.plot(
                timebins,
                snrs,
                color="slateblue",
                marker="o",
                markersize=7,
                linewidth=2,
                linestyle="--",
            )
            ax_snrvstime.set_xlabel("Timebin (s)", fontsize=16, labelpad=10)
            ax_snrvstime.set_ylabel("SNR", fontsize=16, labelpad=10)
            ax_snrvstime.set_title(
                f"SNR vs Timebin for {grb_name}",
                fontsize=18,
                pad=10,
            )
            pdf = PdfPages(f"output_for_{grb_name}.pdf")
            for fig in figs:
                pdf.savefig(fig)
                plt.close()
            pdf.savefig(fig_snrvstime)
            plt.close()
            pdf.close()
            print(f"\n\nBEST TIMEBIN: {optimal_timebin}s")
            quadrants = np.nonzero(potential_grb_snrs)[0]
            for i in quadrants:
                print(
                    f"Potential GRB found in Quadrant {i} at {potential_grb_times[i]}s with SNR {np.round(potential_grb_snrs[i], 2)}!!!!"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to identify potential GRBs in the given directory"
    )
    parser.add_argument("-d", help="Input directory", type=str)
    parser.add_argument("-t", help="Trigger time", type=float)
    parser.add_argument("-s", help="Detection sigma", type=float, default=1)
    parser.add_argument("-n", help="GRB name", type=str)
    parser.add_argument(
        "--timebin", help="Enter Manual timebin", type=float, default=None
    )

    directory = parser.parse_args().d
    trigger_time = parser.parse_args().t
    detection_sigma = parser.parse_args().s
    grb_name = parser.parse_args().n
    input_timebin = parser.parse_args().timebin

    t = time.time()
    main(directory, trigger_time, grb_name, input_timebin, detection_sigma)
    print(f"Total Time taken: {time.time() - t}s")
