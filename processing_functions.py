from astropy.io import fits
import numpy as np
from scipy.signal import savgol_filter
from scipy.stats import poisson
from scipy.stats import gamma
from scipy.stats import exponweib
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

def snr_gauss(filename, start, end, polyorder=3, in_bins=100, window=101):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder, window)
    if end<south_atlantic_start:
        total_noise = np.concatenate((data['RATE'][:start], data['RATE'][end:south_atlantic_start], data['RATE'][south_atlantic_end:]))
    elif start>south_atlantic_end:
        total_noise = np.concatenate((data['RATE'][:south_atlantic_start], data['RATE'][south_atlantic_end:start], data['RATE'][end:]))
    else:
        print('Inputted start and end times are not valid')
    mu = np.mean(total_noise)
    bins = np.arange(int(mu-in_bins/2), int(mu+in_bins/2)) - 0.5
    n, bin_edges = np.histogram(total_noise, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(gaussian, bin_center, n, p0=[np.max(n), np.mean(total_noise), np.std(total_noise), 0])
    signal = np.max(data['RATE'][start:end])
    noise = popt[1]+3*popt[2]
    snr = (signal+popt[1])/noise
    return snr, n, bin_center, popt

def poisson_fit(k, lamb, c):
    return c*poisson.pmf(k, lamb)

def snr_poisson(filename, start, end, polyorder=3, in_bins=100, window=101):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder, window)
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
    signal = np.max(data['RATE'][start:end])+popt[0]
    noise = popt[0]+3*np.sqrt(popt[0])
    snr = signal/noise
    return snr, n, bin_center, popt

def snr_gamma(filename, start, end, polyorder=3, in_bins=100, window=101):
    
    def gamma_fit(x, c):
        return c*gamma.pdf(x, k, scale=theta)
    
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder, window)
    if end<south_atlantic_start:
        total_noise = np.concatenate((data['RATE'][:start], data['RATE'][end:south_atlantic_start], data['RATE'][south_atlantic_end:]))
    elif start>south_atlantic_end:
        total_noise = np.concatenate((data['RATE'][:south_atlantic_start], data['RATE'][south_atlantic_end:start], data['RATE'][end:]))
    else:
        print('Inputted start and end times are not valid')
    
    params = gamma.fit(total_noise)
    k, loc, theta = params[0], params[1], params[2]
    noise_for_fit = total_noise - loc
    bins = np.arange(int(-loc-in_bins/2), int(-loc+in_bins/2)) - 0.5
    n, bin_edges = np.histogram(noise_for_fit, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(gamma_fit, bin_center, n)
    fit = gamma_fit(bin_center, *popt)
    signal = np.max(data['RATE'][start:end])-loc
    noise = -loc + 3*np.sqrt(k*theta**2)
    snr = signal/noise
    return snr, n, bin_center, fit, popt

def snr_weibull(filename, start, end, polyorder=3, in_bins=100, window=101):
    
    def weibull_fit(x, k):
        return k*exponweib.pdf(x, a, c, scale=scale)

    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder, window)
    if end<south_atlantic_start:
        total_noise = np.concatenate((data['RATE'][:start], data['RATE'][end:south_atlantic_start], data['RATE'][south_atlantic_end:]))
    elif start>south_atlantic_end:
        total_noise = np.concatenate((data['RATE'][:south_atlantic_start], data['RATE'][south_atlantic_end:start], data['RATE'][end:]))
    else:
        print('Inputted start and end times are not valid')
    
    params = exponweib.fit(total_noise)
    a, c, loc, scale = params[0], params[1], params[2], params[3]
    noise_for_fit = total_noise - loc
    bins = np.arange(int(-loc-in_bins/2), int(-loc+in_bins/2)) - 0.5
    n, bin_edges = np.histogram(noise_for_fit, bins=bins)
    bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    popt, pcov = curve_fit(weibull_fit, bin_center, n)
    fit = weibull_fit(bin_center, *popt)
    signal = np.max(data['RATE'][start:end])-loc
    noise = -loc + 3*exponweib.std(a, c, scale=scale)
    snr = signal/noise
    return snr, n, bin_center, fit, popt

def snr_counts(filename, start, end, polyorder=3):
    data, south_atlantic_start, south_atlantic_end = filter_and_detrend(filename, start, end, polyorder)
    duration_lc = data['TIME'][-1]-data['TIME'][0]
    duration_burst = data['TIME'][end]-data['TIME'][start]
    duration_sao = data['TIME'][south_atlantic_end]-data['TIME'][south_atlantic_start]
    duration_noise = duration_lc-duration_burst-duration_sao
    if end<south_atlantic_start:
        signal = np.sum(data['RATE'][start:end])/duration_burst
        noise = (np.sum(np.abs(data['RATE'][:start]))+np.sum(np.abs(data['RATE'][end:south_atlantic_start]))+np.sum(np.abs(data['RATE'][south_atlantic_end:])))/duration_noise
        snr = signal/noise
    elif start>south_atlantic_end:
        signal = np.sum(data['RATE'][start:end])/duration_burst
        noise = (np.sum(np.abs(data['RATE'][:south_atlantic_start]))+np.sum(np.abs(data['RATE'][south_atlantic_end:start]))+np.sum(np.abs(data['RATE'][end:])))/duration_noise
        snr = signal/noise
    else:
        print('Inputted start and end times are not valid')
        snr = 0
    return snr
