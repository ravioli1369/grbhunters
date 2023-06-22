from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from scipy.special import factorial
from scipy.stats import poisson
from astropy.stats import sigma_clipped_stats
from scipy.optimize import curve_fit
import matplotlib
params={
    'text.usetex':True,
    'font.family':'serif',
    'xtick.minor.visible':True,
    'ytick.minor.visible':True,
    'xtick.top':True,
    'ytick.left':True,
    'ytick.right':True,
    'xtick.direction':'out',
    'ytick.direction':'out',
    'xtick.minor.size':2.5,
    'xtick.major.size':5,
    'ytick.minor.size':2.5,
    'ytick.major.size':5,
    'axes.axisbelow':True
}
matplotlib.rcParams.update(params)

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
    no_of_zeros = ranges[:,1] - ranges[:,0]
    max_zeros = np.where(no_of_zeros == np.max(no_of_zeros))[0][0]
    return ranges[max_zeros]

def filter_and_detrend(filename, start, end, polyorder=3, window=101):
    data = openlc(filename)
    south_atlantic_start = zero_runs(data['RATE'])[0]
    south_atlantic_end = zero_runs(data['RATE'])[-1]
    if end<south_atlantic_start:
        data['RATE'][start:end] -= (np.mean([data['RATE'][:start]]) + np.mean(data['RATE'][end:south_atlantic_start] + np.mean(data['RATE'][south_atlantic_end:])))/3
        data['RATE'][:start] -= savgol_filter(data['RATE'][:start], window, polyorder)
        data['RATE'][end:south_atlantic_start] -= savgol_filter(data['RATE'][end:south_atlantic_start], window, polyorder)
        data['RATE'][south_atlantic_end:] -= savgol_filter(data['RATE'][south_atlantic_end:], window, polyorder)
    elif start>south_atlantic_end:
        data['RATE'][start:end] -= (np.mean([data['RATE'][:south_atlantic_start]]) + np.mean(data['RATE'][south_atlantic_end:start] + np.mean(data['RATE'][end:])))/3
        data['RATE'][:south_atlantic_start] -= savgol_filter(data['RATE'][:south_atlantic_start], window, polyorder)
        data['RATE'][south_atlantic_end:start] -= savgol_filter(data['RATE'][south_atlantic_end:start], window, polyorder)
        data['RATE'][end:] -= savgol_filter(data['RATE'][end:], window, polyorder)
    else:
        print('Inputted start and end times are not valid')
    return data, south_atlantic_start, south_atlantic_end

def snr_rms(filename, start, end, polyorder=3):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder)
    if end<south_atlantic_start:
        rms = np.sqrt(np.mean(data['RATE'][:start]**2))   #ignoring the noise after the GRB
        signal = np.max(data['RATE'][start:end])
        snr = signal/rms
    elif start>south_atlantic_end:
        rms = np.sqrt(np.mean(data['RATE'][:south_atlantic_start]**2+data['RATE'][south_atlantic_end:start]**2))
        signal = np.max(data['RATE'][start:end])
        snr = signal/rms
    else:
        print('Inputted start and end times are not valid')
        snr = 0
    return snr

def snr_abs(filename, start, end, polyorder=3):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder)
    if end<south_atlantic_start:
        abs = np.mean(np.abs(data['RATE'][:start]))
        signal = np.max(data['RATE'][start:end])
        snr = signal/abs
    elif start>south_atlantic_end:
        abs = (np.mean(np.abs(data['RATE'][:south_atlantic_start]))+np.mean(np.abs(data['RATE'][south_atlantic_end:start])))/2
        signal = np.max(data['RATE'][start:end])
        snr = signal/abs
    else:
        print('Inputted start and end times are not valid')
        snr = 0
    return snr

def gaussian(x, A, m, s, c):
    return A*np.exp(-(x-m)**2/(2*s**2)) + c

def snr_gauss(filename, start, end, polyorder=3, in_bins=80):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder)
    if end<south_atlantic_start:
        total_noise = np.concatenate((data['RATE'][:start], data['RATE'][end:south_atlantic_start], data['RATE'][south_atlantic_end:]))
    elif start>south_atlantic_end:
        total_noise = np.concatenate((data['RATE'][:south_atlantic_start], data['RATE'][south_atlantic_end:start], data['RATE'][end:]))
    else:
        print('Inputted start and end times are not valid')
        snr = 0
        total_noise = 0
        popt = [0]
        bin_center = 0
    n, bins = np.histogram(total_noise, bins=in_bins)
    bin_center = np.array([0.5*(bins[i]+bins[i+1]) for i in range(len(bins)-1)])
    popt, pcov = curve_fit(gaussian, bin_center, n, p0=[np.max(n), np.mean(total_noise), np.std(total_noise), 0])
    signal = np.max(data['RATE'][start:end])
    noise = popt[1]+3*popt[2]
    snr = signal/noise
    return snr, total_noise, bin_center, popt

def poisson_fit(k, lamb, c):
    return c*poisson.pmf(k, lamb)

def snr_poisson(filename, start, end, in_bins=100):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, 3, 11)
    if end<south_atlantic_start:
        noise = np.concatenate((data['RATE'][:start], data['RATE'][end:south_atlantic_start], data['RATE'][south_atlantic_end:]))
    elif start>south_atlantic_end:
        noise = np.concatenate((data['RATE'][:south_atlantic_start], data['RATE'][south_atlantic_end:start], data['RATE'][end:]))
    else:
        print('Inputted start and end times are not valid')
    lamb = np.std(noise)**2
    noise_for_fit = noise+lamb
    bins = np.arange(int(lamb - in_bins/2), int(lamb + in_bins/2)) - 0.5
    n, bin_edges = np.histogram(noise_for_fit, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(poisson_fit, bin_center, n, p0=[np.mean(noise_for_fit), 2000])
    signal = np.max(data['RATE'][start:end])+lamb
    noise = popt[0]+3*np.sqrt(popt[0])
    snr = signal/noise
    return snr, n, bin_center, popt

def snr_counts(filename, start, end, polyorder=3):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder)
    duration_lc = data['TIME'][-1]-data['TIME'][0]
    duration_burst = data['TIME'][end]-data['TIME'][start]
    duration_sao = data['TIME'][south_atlantic_end]-data['TIME'][south_atlantic_start]
    duration_noise = duration_lc-duration_burst-duration_sao
    if end<south_atlantic_start:
        signal = np.sum(data['RATE'][start:end])/duration_burst
        noise = np.sum(np.abs(data['RATE'][:start]))/duration_noise+np.sum(np.abs(data['RATE'][end:south_atlantic_start]))/duration_noise+np.sum(np.abs(data['RATE'][south_atlantic_end:]))/duration_noise
        snr = signal/noise
    elif start>south_atlantic_end:
        signal = np.sum(data['RATE'][start:end])/duration_burst
        noise = np.sum(np.abs(data['RATE'][:south_atlantic_start]))/duration_noise+np.sum(np.abs(data['RATE'][south_atlantic_end:start]))/duration_noise+np.sum(np.abs(data['RATE'][end:]))/duration_noise
        snr = signal/noise
    else:
        print('Inputted start and end times are not valid')
        snr = 0
    return snr
