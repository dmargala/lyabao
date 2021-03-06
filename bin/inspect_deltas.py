#!/usr/bin/env python
import argparse
from itertools import chain

import numpy as np
import numpy.ma as ma
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import h5py

from tqdm import tqdm

from uniform_grid import get_fiducial_pixel_index_offset
from uniform_grid import get_fiducial_wavelength


def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ## targets to fit
    parser.add_argument('--name', type=str, default=None,
        help='base filename to use for input and output')
    parser.add_argument('--max-rows', type=int, default=0,
        help='max rows')
    parser.add_argument('--max-ivar', type=float, default=200,
        help='maximum pipeline ivar for binning')
    parser.add_argument('--min-ivar', type=float, default=0,
        help='minimum pipeline ivar for binning')
    parser.add_argument('--num-ivar', type=int, default=50,
        help='number of ivar bins for lss measurement')
    parser.add_argument('--max-z', type=float, default=3.4,
        help='maximum redshift for binning')
    parser.add_argument('--min-z', type=float, default=2.0,
        help='minimum redshift for binning')
    parser.add_argument('--num-z', type=int, default=7,
        help='number of redshift bins for lss var measurement')
    parser.add_argument('--max-delta', type=float, default=1e3,
        help='filter pixels where delta is >= this value')
    args = parser.parse_args()

    # import data
    infile = h5py.File(args.name + '-delta.hdf5', 'r')
    lines_of_sight = infile['lines_of_sight']
    # read attributes with info on processing so far
    coeff0 = infile.attrs['coeff0']
    coeff1 = infile.attrs['coeff1']
    num_wave_obs = infile.attrs['max_fid_index']
    try:
        wave_lya = infile.attrs['wave_lya']
    except KeyError:
        wave_lya = 1216.0


    num_sightlines = len(lines_of_sight.keys())
    # if specified, only use max_rows number of targets
    if args.max_rows > 0:
        num_sightlines = args.max_rows
    print('num sightlines: {:d}'.format(num_sightlines))

    # loop over targets
    waves = []
    deltas = []
    ivars = []
    weights = []
    thing_ids = lines_of_sight.keys()[:num_sightlines]
    for i, thing_id in tqdm(enumerate(thing_ids), total=num_sightlines):
        los = lines_of_sight[thing_id]
        z = los.attrs['z']
        loglam = los['loglam'].value
        delta = los['delta'].value

        # if np.any(np.abs(delta) > 1e3):
        #     # print thing_id, los.attrs['plate'], los.attrs['mjd'], los.attrs['fiber']
        #     continue
        ivar = los['ivar'].value
        weight = los['weight'].value

        valid = (ivar > 0) & (np.abs(delta) < args.max_delta)
        waves.append(np.power(10, loglam[valid]))
        deltas.append(delta[valid])
        ivars.append(ivar[valid])
        weights.append(weight[valid])
    infile.close()

    # concatenate all pixels
    all_waves = np.fromiter(chain.from_iterable(waves), np.float)
    all_deltas = np.fromiter(chain.from_iterable(deltas), np.float)
    all_ivars = np.fromiter(chain.from_iterable(ivars), np.float)
    all_weights = np.fromiter(chain.from_iterable(weights), np.float)
    print('Total number of pixels: {:d}'.format(len(all_waves)))

    # set up binning for observed frame stats

    # subtract 0.5 so that these are bin edges instead of bin centers
    wave_bins = get_fiducial_wavelength(
        np.arange(num_wave_obs + 1) - 0.5, coeff1=coeff1)
    wave_bin_centers = get_fiducial_wavelength(
        np.arange(num_wave_obs), coeff1=coeff1)
    # determine observed wavelength bin indices
    wave_bin_indices = np.digitize(all_waves, wave_bins)
    # calculate stats as a function of observed wavelength
    counts_per_obs_pixel = ma.masked_all(num_wave_obs)
    mean_delta = ma.masked_all(num_wave_obs)
    wmean_delta = ma.masked_all(num_wave_obs)
    for i in np.unique(wave_bin_indices):
        i_indices = wave_bin_indices == i
        if i - 1 >= num_wave_obs:
            continue
        counts_per_obs_pixel[i - 1] = ma.count(all_waves[i_indices])
        mean_delta[i - 1] = ma.average(all_deltas[i_indices])
        wmean_delta[i - 1] = ma.average(
            all_deltas[i_indices], weights=all_weights[i_indices]
        )

    # set up binning for lss variance measurement
    ivar_bin_edges = np.linspace(
        args.min_ivar, args.max_ivar, args.num_ivar + 1)
    ivar_bin_centers = 0.5 * (ivar_bin_edges[:-1] + ivar_bin_edges[1:])
    # we want to perform this measurement in redshift slices
    z_bin_edges = np.linspace(args.min_z, args.max_z, args.num_z+1)
    z_bin_centers = 0.5 * (z_bin_edges[:-1] + z_bin_edges[1:])
    all_redshifts = all_waves/wave_lya - 1.0

    ivar_delta_per_ivar = ma.zeros((args.num_z, args.num_ivar))
    ivar_delta_per_ivar.mask = True
    counts_per_ivar = ma.zeros((args.num_z, args.num_ivar))
    counts_per_ivar.mask = True
    for j in range(args.num_z):
        # select pixels for the current redshift interval
        z_slice = (
            (all_redshifts >= z_bin_edges[j]) &
            (all_redshifts < z_bin_edges[j+1])
        )
        # determine pipeline ivar bins
        ivar_bin_indices = np.digitize(all_ivars[z_slice], ivar_bin_edges)
        # compute stats per ivar bin
        for i in np.unique(ivar_bin_indices):
            # skip ivar bins outside our grid
            if i - 1 >= args.num_ivar:
                continue
            ivar_slice = ivar_bin_indices == i
            mean = ma.average(
                all_deltas[z_slice][ivar_slice],
                weights=all_ivars[z_slice][ivar_slice])
            var = ma.average(
                (all_deltas[z_slice][ivar_slice]-mean)**2,
                weights=all_ivars[z_slice][ivar_slice])
            if var > 0:
                ivar_delta_per_ivar[j, i - 1] = 1.0/var
            counts_per_ivar[j, i - 1] = ma.count(all_deltas[z_slice][ivar_slice])

    print('Fitting lss ivar...')

    def lss_ivar_chisq(p):
        ivar_model = (
            p[:args.num_z, np.newaxis]/ivar_bin_centers[np.newaxis, :] +
            p[args.num_z:, np.newaxis]
        )
        dy = 1.0/ivar_delta_per_ivar - ivar_model
        return ma.average(dy*dy, weights=1.0/counts_per_ivar)

    p0 = np.concatenate([np.ones(args.num_z), np.zeros(args.num_z)])
    result = minimize(lss_ivar_chisq, p0)
    print(result.x)

    # create figures
    savefig_kwargs = {
        'bbox_inches': 'tight',
        'dpi': 100
    }

    plt.figure(figsize=(12, 9))
    for i in range(args.num_z):
        p, = plt.plot(
            ivar_bin_centers, ivar_delta_per_ivar[i],
            marker='o', label=z_bin_centers[i], lw=0
        )
        plt.plot(
            ivar_bin_centers,
            1.0/(result.x[i]/ivar_bin_centers + result.x[args.num_z + i]),
            color=p.get_color()
        )
    plt.xlim(ivar_bin_edges[0], ivar_bin_edges[-1])
    plt.xlabel('Pipeline ivar')
    plt.ylim(0, 30)
    plt.ylabel('Inverse Delta Variance')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(args.name + '-delta-ivar.png', **savefig_kwargs)

    plt.figure(figsize=(12, 9))
    plt.plot(wave_bin_centers, mean_delta,
        marker='+', lw=0, label='Mean')
    plt.plot(wave_bin_centers, wmean_delta,
        marker='+', lw=0, label='Weighted Mean')
    plt.xlim(wave_bins[0], wave_bins[-1])
    plt.ylim(-2, 2)
    plt.xlabel('Observed Wavelength')
    plt.ylabel('Mean Delta')
    plt.legend()
    plt.savefig(args.name + '-delta-mean-obs.png', **savefig_kwargs)

    plt.figure(figsize=(12, 9))
    plt.scatter(wave_bin_centers, counts_per_obs_pixel)
    plt.xlim(wave_bins[0], wave_bins[-1])
    plt.xlabel('Observed Wavelength')
    plt.ylabel('Counts per obs pixel')
    plt.savefig(args.name + '-wave-hist-obs.png', **savefig_kwargs)

    plt.figure(figsize=(12, 9))
    plt.hist(all_waves, bins=wave_bins, histtype='step')
    plt.xlim(wave_bins[0], wave_bins[-1])
    plt.xlabel('Observed Wavelength')
    plt.ylabel('Num pixels')
    plt.savefig(args.name + '-wave-hist.png', **savefig_kwargs)

    plt.figure(figsize=(12, 9))
    plt.scatter(all_waves, all_deltas, marker=',', lw=0, s=1)
    plt.xlim(wave_bins[0], wave_bins[-1])
    plt.xlabel('Observed Wavelength')
    plt.savefig(args.name + '-delta-scatter.png', **savefig_kwargs)

    plt.figure(figsize=(12, 9))
    plt.scatter(all_waves, all_ivars, marker=',', lw=0, s=1)
    plt.xlim(wave_bins[0], wave_bins[-1])
    plt.xlabel('Observed Wavelength')
    plt.savefig(args.name + '-ivar-scatter.png', **savefig_kwargs)

    plt.figure(figsize=(12, 9))
    plt.scatter(all_ivars, np.abs(all_deltas), marker=',', lw=0, s=1)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-8, 1e18)
    plt.xlim(1e-18, 1e4)
    plt.xlabel('ivar')
    plt.ylabel('delta')
    plt.savefig(args.name + '-ivar-delta-scatter.png', **savefig_kwargs)

    plt.figure(figsize=(12, 9))
    plt.scatter(all_weights, np.abs(all_deltas), marker=',', lw=0, s=1)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-8, 1e18)
    plt.xlim(1e-18, 1e4)
    plt.xlabel('weight')
    plt.ylabel('delta')
    plt.savefig(args.name + '-weight-delta-scatter.png', **savefig_kwargs)


if __name__ == '__main__':
    main()
