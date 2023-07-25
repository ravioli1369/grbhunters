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
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import QTable
from astropy.stats import sigma_clipped_stats
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from pipelinev3 import bindata

import warnings

warnings.simplefilter("ignore", np.RankWarning)

##############################################################################################
############ Functions for detrending, outlier detection, snr, energy binning ################
##############################################################################################


def zero_runs(a):
    """
    Return maximum runs of consecutive zeros in a 1D array.
    """
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    no_of_zeros = ranges[:, 1] - ranges[:, 0]
    max_zeros = np.where(no_of_zeros == np.max(no_of_zeros))[0][0]
    return ranges[max_zeros]


def saa(a):
    """
    Returns the start and end indices of the South Atlantic Anomaly
    """
    diff = np.diff(a)
    saa_start = np.argmax(diff)
    saa_end = saa_start + 1
    return saa_start, saa_end


def rebin(count_data, time_data, bin_size):
    n = int(bin_size * 10)
    rebin_count = [
        np.median(count_data[i : i + n]) for i in range(0, len(count_data), n)
    ]
    rebin_time = [np.mean(time_data[i : i + n]) for i in range(0, len(time_data), n)]
    return rebin_count, rebin_time


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
    filename, trigger_index, polyorder=2, detrend_window=101, data=None
):
    """
    Detrends the given file using a quadratic fit to the data around the trigger
    """
    if data is None:
        data = fits.getdata(filename)
    if zero_runs(data["RATE"])[-1] == len(data["RATE"]):
        saa_start, saa_end = saa(data["TIME"])
    else:
        length_saa_zeroes = data['TIME'][zero_runs(data["RATE"])[-1]] - data['TIME'][zero_runs(data["RATE"])[0]]
        length_saa_nans = data['TIME'][saa(data["TIME"])[-1]] - data['TIME'][saa(data["TIME"])[0]]
        if length_saa_zeroes > length_saa_nans:
            saa_start, saa_end = zero_runs(data["RATE"])
        else:
            saa_start, saa_end = saa(data["TIME"])
    # saa_start, saa_end = saa(data["TIME"])
    timebin = data["TIME"][trigger_index + 1] - data["TIME"][trigger_index]
    detrend_window = (
        np.rint(detrend_window / timebin).astype(int) // 2 * 2 + 1
    )  # make it odd
    background_window = np.rint(700 / timebin).astype(int)
    if trigger_index < saa_start:
        if (
            trigger_index > background_window
            and trigger_index + background_window < saa_start
        ):
            counts = data["RATE"][
                trigger_index - background_window : trigger_index + background_window
            ]
            times = data["TIME"][
                trigger_index - background_window : trigger_index + background_window
            ]
            new_trigger_index = background_window
        elif (
            trigger_index < background_window
            and trigger_index + background_window < saa_start
        ):
            counts = data["RATE"][: trigger_index + background_window]
            times = data["TIME"][: trigger_index + background_window]
            new_trigger_index = trigger_index
        elif (
            trigger_index > background_window
            and trigger_index + background_window > saa_start
        ):
            counts = data["RATE"][trigger_index - background_window :]
            times = data["TIME"][trigger_index - background_window :]
            new_trigger_index = background_window
        else:
            counts = data["RATE"]
            times = data["TIME"]
            new_trigger_index = trigger_index
    elif trigger_index > saa_end:
        if (
            trigger_index - saa_end > background_window
            and trigger_index + background_window < len(data)
        ):
            counts = data["RATE"][
                trigger_index - background_window : trigger_index + background_window
            ]
            times = data["TIME"][
                trigger_index - background_window : trigger_index + background_window
            ]
            new_trigger_index = background_window
        elif (
            trigger_index - saa_end < background_window
            and trigger_index + background_window < len(data)
        ):
            counts = data["RATE"][saa_end : trigger_index + background_window]
            times = data["TIME"][saa_end : trigger_index + background_window]
            new_trigger_index = trigger_index
        elif (
            trigger_index - saa_end > background_window
            and trigger_index + background_window > len(data)
        ):
            counts = data["RATE"][trigger_index - background_window :]
            times = data["TIME"][trigger_index - background_window :]
            new_trigger_index = background_window
        else:
            counts = data["RATE"]
            times = data["TIME"]
            new_trigger_index = trigger_index
    else:
        raise ValueError("Trigger index is in SAA")
    # clipping the outliers before fitting the quadratic
    mean, _, std = sigma_clipped_stats(counts)
    counts = np.copy(counts)
    counts[np.abs(counts - mean) > 3 * std] = np.nan
    filtered = counts
    # filtered = savgol_filter(counts, detrend_window, 3)
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
    return detrended, raw, trend, new_trigger_index, popt


def create_master_lc(directory):
    """
    Creates the master light curve (20-200 keV) for the given directory
    """
    emin, emax = 20, 200
    if not os.path.exists(f"{directory}/master_lc"):
        print("Creating master light curve at {}/master_lc".format(directory))
        os.mkdir(f"{directory}/master_lc")
        os.system(
            "python3 pipelinev3.py -d {} -emin {} -emax {}".format(
                directory, emin, emax
            )
        )
        os.system("mv {}/*.lc {}/master_lc".format(directory, directory))
    else:
        print("Master light curve already exists at {}/master_lc".format(directory))


def gen_energy_bins(directory, n_bins=3):
    """
    Generates energy bins for the given number of bins
    """
    if n_bins != 3:
        energy_ranges = np.linspace(20, 200, n_bins + 1)
    else:
        energy_ranges = [20, 60, 100, 200]
    lc_paths = []
    if not os.path.exists(f"{directory}/{n_bins}_bins"):
        print("Creating {} energy bins at {}/{}_bins".format(n_bins, directory, n_bins))
        os.mkdir(f"{directory}/{n_bins}_bins")
        for i in range(n_bins):
            emin = energy_ranges[i]
            emax = energy_ranges[i + 1]
            evt = glob.glob(f"{directory}/*bc.evt")[0]
            mkf = glob.glob(f"{directory}/*.mkf")[0]
            print(emin, emax)
            bindata(evt, mkf, 1, emin, emax)
            if not os.path.exists(f"{directory}/{n_bins}_bins/{int(emin)}-{int(emax)}"):
                os.mkdir(f"{directory}/{n_bins}_bins/{int(emin)}-{int(emax)}")
            os.system(
                f"mv {directory}/*.lc {directory}/{n_bins}_bins/{int(emin)}-{int(emax)}/"
            )
            lc_paths.append(f"{directory}/{n_bins}_bins/{int(emin)}-{int(emax)}/")
    else:
        print(
            "Requested energy bins already exist at {}/{}_bins".format(
                directory, n_bins
            )
        )
        for i in range(n_bins):
            emin = energy_ranges[i]
            emax = energy_ranges[i + 1]
            lc_paths.append(f"{directory}/{n_bins}_bins/{int(emin)}-{int(emax)}/")
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


def find_grb(directory, trigger_time, detection_sigma=3):
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

    def snr_grb(master_lc, possible_grb, trigger_index):
        """
        Returns the SNR of the potential GRBs in the master light curve
        """
        t, *_ = quadratic_detrend_trigger(master_lc, trigger_index)
        mean, _, std = sigma_clipped_stats(t["RATE"])
        noise = mean + std
        signal = t["RATE"][possible_grb]
        snr = signal / noise
        return snr

    create_master_lc(directory)
    master_lcs = np.sort(glob.glob(f"{directory}/master_lc/*.lc"))
    lc_paths = gen_energy_bins(directory)
    results = []
    for i in range(4):
        trigger_index = get_trigger_index(master_lcs[i], trigger_time)
        snr_each_quad, outliers_each_quad = each_quad(lc_paths, trigger_index, i)
        grb_mask_each_quad = np.logical_or(snr_each_quad[1] > 3, snr_each_quad[2] > 3)
        possible_grb_each_quad = outliers_each_quad[grb_mask_each_quad]
        grb_snr_each_quad = snr_grb(
            master_lcs[i], possible_grb_each_quad, trigger_index
        )
        results.append(
            [snr_each_quad, outliers_each_quad, grb_mask_each_quad, grb_snr_each_quad]
        )
    results.append(master_lcs)
    results.append(lc_paths)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to identify potential GRBs in the given directory"
    )
    parser.add_argument("-d", help="Input directory", type=str)
    parser.add_argument("-t", help="Trigger time", type=float)
    directory = parser.parse_args().d
    trigger_time = parser.parse_args().t
