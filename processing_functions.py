import os
import glob
import numpy as np
from astropy.io import fits
from astropy.table import QTable
from astropy.stats import sigma_clipped_stats, sigma_clip
from scipy.signal import savgol_filter
from scipy.stats import poisson
from scipy.stats import gamma
from scipy.stats import skewnorm
from scipy.optimize import curve_fit


def openlc(filename):
    hdul = fits.open(filename)
    data = hdul[1].data
    hdul.close()
    return data


def zero_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    no_of_zeros = ranges[:, 1] - ranges[:, 0]
    max_zeros = np.where(no_of_zeros == np.max(no_of_zeros))[0][0]
    return ranges[max_zeros]


def saa(a):
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


def filter_and_detrend(filename, start, end, polyorder=3, window=101, data=None):
    if data is None:
        data = openlc(filename)
    saa_start = zero_runs(data["RATE"])[0]
    saa_end = zero_runs(data["RATE"])[-1]
    if end < saa_start:
        data["RATE"][start:end] -= (
            np.mean([data["RATE"][:start]])
            + np.mean(data["RATE"][end:saa_start] + np.mean(data["RATE"][saa_end:]))
        ) / 3
        data["RATE"][:start] -= savgol_filter(data["RATE"][:start], window, polyorder)
        data["RATE"][end:saa_start] -= savgol_filter(
            data["RATE"][end:saa_start], window, polyorder
        )
        data["RATE"][saa_end:] -= savgol_filter(
            data["RATE"][saa_end:], window, polyorder
        )
    elif start > saa_end:
        data["RATE"][start:end] -= (
            np.mean([data["RATE"][:saa_start]])
            + np.mean(data["RATE"][saa_end:start] + np.mean(data["RATE"][end:]))
        ) / 3
        data["RATE"][:saa_start] -= savgol_filter(
            data["RATE"][:saa_start], window, polyorder
        )
        data["RATE"][saa_end:start] -= savgol_filter(
            data["RATE"][saa_end:start], window, polyorder
        )
        data["RATE"][end:] -= savgol_filter(data["RATE"][end:], window, polyorder)
    else:
        print("Inputted start and end times are not valid")
    return data, saa_start, saa_end


def quadratic(x, a, b, c):
    return a * x**2 + b * x + c


def quadratic_detrend(filename, start, end, polyorder=3, window=101, data=None):
    if data is None:
        data = openlc(filename)
    saa_start, saa_end = saa(data["TIME"])
    timebin = (data.TIME[end] - data.TIME[start]) / (end - start)
    background_window = np.rint(500 / timebin).astype(int)
    grb = np.ones_like(data["RATE"][start:end]) * np.nan
    if end < saa_start:
        if start > background_window and (saa_start - end) > background_window:
            counts = np.concatenate(
                (
                    data["RATE"][start - background_window : start],
                    grb,
                    data["RATE"][end : end + background_window],
                )
            )
            times = data["TIME"][start - background_window : end + background_window]
            new_start, new_end = background_window, background_window + end - start
        elif start < background_window and (saa_start - end) > background_window:
            counts = np.concatenate(
                (data["RATE"][:start], grb, data["RATE"][end : end + background_window])
            )
            times = data["TIME"][: end + background_window]
            new_start, new_end = start, end
        elif start > background_window and (saa_start - end) < background_window:
            counts = np.concatenate(
                (
                    data["RATE"][start - background_window : start],
                    grb,
                    data["RATE"][end:saa_start],
                )
            )
            times = data["TIME"][start - background_window : saa_start]
            new_start, new_end = background_window, background_window + end - start
        else:
            counts = np.concatenate((data["RATE"][:start], grb, data["RATE"][end:]))
            times = data["TIME"][:]
            new_start, new_end = start, end
    elif start > saa_end:
        if (start - saa_end) > background_window and (
            len(data.RATE) - end
        ) > background_window:
            counts = np.concatenate(
                (
                    data["RATE"][start - background_window : start],
                    grb,
                    data["RATE"][end : end + background_window],
                )
            )
            times = data["TIME"][start - background_window : end + background_window]
        elif (start - saa_end) < background_window and (
            len(data.RATE) - end
        ) > background_window:
            counts = np.concatenate(
                (
                    data["RATE"][saa_end:start],
                    grb,
                    data["RATE"][end : end + background_window],
                )
            )
            times = data["TIME"][saa_end : end + background_window]
        elif (start - saa_end) > background_window and (
            len(data.RATE) - end
        ) < background_window:
            counts = np.concatenate(
                (
                    data["RATE"][start - background_window : start],
                    grb,
                    data["RATE"][end:],
                )
            )
            times = data["TIME"][start - background_window :]
        else:
            counts = np.concatenate((data["RATE"][:start], grb, data["RATE"][end:]))
            times = data["TIME"][:]
    else:
        print("Inputted start and end times are not valid")
    filtered = savgol_filter(counts, window, polyorder)
    idx = np.isfinite(filtered)
    x = times[idx]
    y = filtered[idx]
    popt, _ = curve_fit(quadratic, x, y)
    nans = np.where(np.isnan(counts))[0]
    counts[nans[0] : nans[-1] + 1] = data["RATE"][start:end]
    detrended = counts - quadratic(times, *popt)
    t = QTable([times, detrended], names=("TIME", "RATE"))
    return t, saa_start, saa_end, new_start, new_end


def get_trigger_index(filename, trigger_time):
    data = openlc(filename)
    diffs = np.abs(
        np.round(data["TIME"], 0).astype(int) - np.round(trigger_time, 0).astype(int)
    )
    trigger_index = np.where(diffs == np.min(diffs))[0][-1]
    return trigger_index


def quadratic_detrend_trigger(
    filename, trigger_index, polyorder=3, window=101, data=None
):
    if data is None:
        data = openlc(filename)
    saa_start, saa_end = saa(data["TIME"])
    timebin = data["TIME"][trigger_index + 1] - data["TIME"][trigger_index]
    background_window = np.rint(700 / timebin).astype(int)
    mean, median, std = sigma_clipped_stats(data["RATE"], sigma=3)
    clipped_data = np.copy(data["RATE"])
    clipped_data[clipped_data > mean + 3 * std] = np.nan
    if trigger_index < saa_start:
        if (
            trigger_index > background_window
            and trigger_index + background_window < saa_start
        ):
            counts = clipped_data[
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
            counts = clipped_data[: trigger_index + background_window]
            times = data["TIME"][: trigger_index + background_window]
            new_trigger_index = trigger_index
        elif (
            trigger_index > background_window
            and trigger_index + background_window > saa_start
        ):
            counts = clipped_data[trigger_index - background_window :]
            times = data["TIME"][trigger_index - background_window :]
            new_trigger_index = background_window
        else:
            counts = clipped_data
            times = data["TIME"]
            new_trigger_index = trigger_index
    elif trigger_index > saa_end:
        if (
            trigger_index > background_window
            and trigger_index + background_window < len(data)
        ):
            counts = clipped_data[
                trigger_index - background_window : trigger_index + background_window
            ]
            times = data["TIME"][
                trigger_index - background_window : trigger_index + background_window
            ]
            new_trigger_index = background_window
        elif (
            trigger_index < background_window
            and trigger_index + background_window < len(data)
        ):
            counts = clipped_data[: trigger_index + background_window]
            times = data["TIME"][: trigger_index + background_window]
            new_trigger_index = trigger_index
        elif (
            trigger_index > background_window
            and trigger_index + background_window > len(data)
        ):
            counts = clipped_data[trigger_index - background_window :]
            times = data["TIME"][trigger_index - background_window :]
            new_trigger_index = background_window
        else:
            counts = clipped_data
            times = data["TIME"]
            new_trigger_index = trigger_index
    else:
        raise ValueError("Trigger index is in SAA")
    filtered = savgol_filter(counts, window, polyorder)
    idx = np.isfinite(filtered)
    x = times[idx]
    y = filtered[idx]
    popt, _ = curve_fit(quadratic, x, y)
    window_start, window_end = (
        np.where(data["TIME"] == times[0])[0][0],
        np.where(data["TIME"] == times[-1])[0][0],
    )
    counts = data["RATE"][window_start : window_end + 1]
    detrended = counts - quadratic(times, *popt)
    t = QTable([times, detrended], names=("TIME", "RATE"))
    return t, new_trigger_index


def snr_rms(filename, start, end, polyorder=3):
    data, saa_start, saa_end = filter_and_detrend(filename, start, end, polyorder)
    if end < saa_start:
        rms = np.sqrt(
            np.mean(data["RATE"][:start] ** 2)
        )  # ignoring the noise after the GRB
        signal = np.max(data["RATE"][start:end])
        snr = signal / rms
    elif start > saa_end:
        rms = np.sqrt(
            np.mean(data["RATE"][:saa_start] ** 2 + data["RATE"][saa_end:start] ** 2)
        )
        signal = np.max(data["RATE"][start:end])
        snr = signal / rms
    else:
        print("Inputted start and end times are not valid")
        snr = 0
    return snr


def snr_abs(filename, start, end, polyorder=3):
    data, saa_start, saa_end = filter_and_detrend(filename, start, end, polyorder)
    if end < saa_start:
        abs = np.mean(np.abs(data["RATE"][:start]))
        signal = np.max(data["RATE"][start:end])
        snr = signal / abs
    elif start > saa_end:
        abs = (
            np.mean(np.abs(data["RATE"][:saa_start]))
            + np.mean(np.abs(data["RATE"][saa_end:start]))
        ) / 2
        signal = np.max(data["RATE"][start:end])
        snr = signal / abs
    else:
        print("Inputted start and end times are not valid")
        snr = 0
    return snr


def gaussian(x, A, m, s, c):
    return A * np.exp(-((x - m) ** 2) / (2 * s**2)) + c


def snr_gauss(filename, start, end, polyorder=3, in_bins=100, window=101, d=None):
    data, saa_start, saa_end, start, end = quadratic_detrend(
        filename, start, end, polyorder, window, d
    )
    # if end<saa_start:
    #     total_noise = np.concatenate((data['RATE'][:start], data['RATE'][end:saa_start], data['RATE'][saa_end:]))
    # elif start>saa_end:
    #     total_noise = np.concatenate((data['RATE'][:saa_start], data['RATE'][saa_end:start], data['RATE'][end:]))
    # else:
    #     print('Inputted start and end times are not valid')

    total_noise = np.concatenate((data["RATE"][:start], data["RATE"][end:]))

    mu = np.mean(total_noise)
    bins = np.arange(int(mu - in_bins / 2), int(mu + in_bins / 2)) - 0.5
    n, bin_edges = np.histogram(total_noise, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(
        gaussian,
        bin_center,
        n,
        p0=[np.max(n), np.mean(total_noise), np.std(total_noise), 0],
    )
    signal = np.max(data["RATE"][start:end])
    noise = popt[1] + 3 * popt[2]
    snr = (signal + popt[1]) / noise
    return snr, n, bin_center, popt


def poisson_fit(k, lamb, c):
    return c * poisson.pmf(k, lamb)


def snr_poisson(filename, start, end, polyorder=3, in_bins=100, window=101, d=None):
    data, saa_start, saa_end, start, end = quadratic_detrend(
        filename, start, end, polyorder, window, d
    )
    # if end<saa_start:
    #     noise = np.concatenate((data['RATE'][:start], data['RATE'][end:saa_start], data['RATE'][saa_end:]))
    # elif start>saa_end:
    #     noise = np.concatenate((data['RATE'][:saa_start], data['RATE'][saa_end:start], data['RATE'][end:]))
    # else:
    #     print('Inputted start and end times are not valid')

    total_noise = np.concatenate((data["RATE"][:start], data["RATE"][end:]))

    lamb = np.std(noise) ** 2
    noise_for_fit = total_noise + lamb
    bins = np.arange(int(lamb - in_bins / 2), int(lamb + in_bins / 2)) - 0.5
    n, bin_edges = np.histogram(noise_for_fit, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(
        poisson_fit, bin_center, n, p0=[np.mean(noise_for_fit), 2000]
    )
    signal = np.max(data["RATE"][start:end]) + popt[0]
    noise = popt[0] + 3 * np.sqrt(popt[0])
    snr = signal / noise
    return snr, n, bin_center, popt


def snr_gamma(filename, start, end, polyorder=3, in_bins=100, window=101, d=None):
    def gamma_fit(x, c):
        return c * gamma.pdf(x, k, scale=theta)

    data, saa_start, saa_end, start, end = quadratic_detrend(
        filename, start, end, polyorder, window, d
    )
    # if end<saa_start:
    #     total_noise = np.concatenate((data['RATE'][:start], data['RATE'][end:saa_start], data['RATE'][saa_end:]))
    # elif start>saa_end:
    #     total_noise = np.concatenate((data['RATE'][:saa_start], data['RATE'][saa_end:start], data['RATE'][end:]))
    # else:
    #     print('Inputted start and end times are not valid')

    total_noise = np.concatenate((data["RATE"][:start], data["RATE"][end:]))

    params = gamma.fit(total_noise)
    k, loc, theta = params[0], params[1], params[2]
    noise_for_fit = total_noise - loc
    bins = np.arange(int(-loc - in_bins / 2), int(-loc + in_bins / 2)) - 0.5
    n, bin_edges = np.histogram(noise_for_fit, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(gamma_fit, bin_center, n)
    fit = gamma_fit(bin_center, *popt)
    signal = np.max(data["RATE"][start:end]) - loc
    noise = -loc + 3 * np.sqrt(k * theta**2)
    snr = signal / noise
    popt = np.append(popt, (k, loc, theta))
    return snr, n, bin_center, fit, popt


def snr_skewnorm(filename, start, end, polyorder=3, in_bins=100, window=101, d=None):
    def skewnorm_fit(x, k):
        return k * skewnorm.pdf(x, a, scale=scale)

    data, saa_start, saa_end, start, end = quadratic_detrend(
        filename, start, end, polyorder, window, d
    )
    # if end<saa_start:
    #     total_noise = np.concatenate((data['RATE'][:start], data['RATE'][end:saa_start], data['RATE'][saa_end:]))
    # elif start>saa_end:
    #     total_noise = np.concatenate((data['RATE'][:saa_start], data['RATE'][saa_end:start], data['RATE'][end:]))
    # else:
    #     print('Inputted start and end times are not valid')

    total_noise = np.concatenate((data["RATE"][:start], data["RATE"][end:]))

    params = skewnorm.fit(total_noise)
    a, loc, scale = params[0], params[1], params[2]
    noise_for_fit = total_noise - loc
    bins = np.arange(int(-loc - in_bins / 2), int(-loc + in_bins / 2)) - 0.5
    n, bin_edges = np.histogram(noise_for_fit, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(skewnorm_fit, bin_center, n)
    fit = skewnorm_fit(bin_center, *popt)
    signal = np.max(data["RATE"][start:end]) - loc
    noise = -loc + 3 * skewnorm.std(a, scale=scale)
    snr = signal / noise
    popt = np.append(popt, (a, loc, scale))
    return snr, n, bin_center, fit, popt


def skewnorm_trigger(filename, trigger, polyorder=3, in_bins=100, window=101, d=None):
    def skewnorm_fit(x, k):
        return k * skewnorm.pdf(x, a, scale=scale)

    data, trigger = quadratic_detrend_trigger(filename, trigger, polyorder, window, d)
    total_noise = sigma_clip(data["RATE"], sigma=3)

    params = skewnorm.fit(total_noise)
    a, loc, scale = params[0], params[1], params[2]
    noise_for_fit = total_noise - loc
    bins = np.arange(int(-loc - in_bins / 2), int(-loc + in_bins / 2)) - 0.5
    n, bin_edges = np.histogram(noise_for_fit, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(skewnorm_fit, bin_center, n)
    fit = skewnorm_fit(bin_center, *popt)
    noise = -loc + 3 * skewnorm.std(a, scale=scale)
    popt = np.append(popt, (a, loc, scale))
    return n, bin_center, fit, popt


def snr_counts(filename, start, end, polyorder=3, window=101):
    data, saa_start, saa_end = filter_and_detrend(
        filename, start, end, polyorder, window
    )
    duration_lc = data["TIME"][-1] - data["TIME"][0]
    duration_burst = data["TIME"][end] - data["TIME"][start]
    duration_sao = data["TIME"][saa_end] - data["TIME"][saa_start]
    duration_noise = duration_lc - duration_burst - duration_sao
    if end < saa_start:
        noise = (
            np.sum(np.abs(data["RATE"][:start]))
            + np.sum(np.abs(data["RATE"][end:saa_start]))
            + np.sum(np.abs(data["RATE"][saa_end:]))
        ) / duration_noise

    elif start > saa_end:
        noise = (
            np.sum(np.abs(data["RATE"][:saa_start]))
            + np.sum(np.abs(data["RATE"][saa_end:start]))
            + np.sum(np.abs(data["RATE"][end:]))
        ) / duration_noise
    else:
        print("Inputted start and end times are not valid")
        snr = 0

    signal = np.sum(np.abs(data["RATE"][start:end])) / duration_burst
    snr = signal / noise
    return snr


def outlier(filename, trigger_index, detection_sigma=3):
    # *_, popt = skewnorm_trigger(filename, trigger_index)
    t, *_ = quadratic_detrend_trigger(filename, trigger_index)
    mean, median, std = sigma_clipped_stats(t["RATE"])
    # _, std = -popt[2], skewnorm.std(popt[1], scale=popt[3])
    outliers = np.where(t["RATE"] > mean + detection_sigma * std)[0]
    return outliers


def snr_outlier(filename1, filename2, trigger_index, detection_sigma=3):
    outliers = outlier(filename1, trigger_index, detection_sigma)
    # *_, popt = skewnorm_trigger(filename2, trigger_index)
    # mean, std = -popt[2], skewnorm.std(popt[1], scale=popt[3])
    t, *_ = quadratic_detrend_trigger(filename2, trigger_index)
    mean, median, std = sigma_clipped_stats(t["RATE"])
    signal = t["RATE"][outliers] + mean
    noise = mean + std
    snr = signal / noise
    return snr, outliers


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
        os.mkdir(f"{directory}/{n_bins}_bins")
        for i in range(n_bins):
            emin = energy_ranges[i]
            emax = energy_ranges[i + 1]
            print(energy_ranges[i], energy_ranges[i + 1])
            os.system(
                "python3 pipelinev3.py -d {} -emin {} -emax {}".format(
                    directory, emin, emax
                )
            )
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


def create_master_lc(directory):
    emin, emax = 20, 200
    if not os.path.exists(f"{directory}/master_lc"):
        os.mkdir(f"{directory}/master_lc")
        os.system(
            "python3 pipelinev3.py -d {} -emin {} -emax {}".format(
                directory, emin, emax
            )
        )
        os.system("mv {}/*.lc {}/master_lc".format(directory, directory))
    else:
        print("Master light curve already exists at {}/master_lc".format(directory))


def find_grb(directory, trigger_time, detection_sigma=3):
    def each_quad(lc_paths, trigger_index, quadrant):
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
        t, *_ = quadratic_detrend_trigger(master_lc, trigger_index)
        mean, _, std = sigma_clipped_stats(t["RATE"])
        noise = mean + std
        signal = t["RATE"][possible_grb]
        snr = signal / noise
        return snr

    create_master_lc(directory)
    master_lcs = np.sort(glob.glob(f"{directory}/master_lc/*.lc"))
    trigger_index = get_trigger_index(master_lcs[0], trigger_time)
    lc_paths = gen_energy_bins(directory)
    results = []
    for i in range(4):
        snr_each_quad, outliers_each_quad = each_quad(lc_paths, trigger_index, i)
        grb_mask_each_quad = np.logical_or(snr_each_quad[1] > 3, snr_each_quad[2] > 3)
        possible_grb_each_quad = outliers_each_quad[grb_mask_each_quad]
        grb_snr_each_quad = snr_grb(
            master_lcs[0], possible_grb_each_quad, trigger_index
        )
        results.append(
            [snr_each_quad, outliers_each_quad, grb_mask_each_quad, grb_snr_each_quad]
        )
    return results
