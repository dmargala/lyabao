#!/usr/bin/env python
import argparse
import glob

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from tqdm import tqdm


def smooth_cov(data, weight, rp, rt, drt=4, drp=4):
    """
    Based on https://github.com/igmhub/pyLyA/blob/master/py/pylya/utils.py
    """

    nsamples, nbins = data.shape

    w = weight.sum(axis=0)
    mean = (data * weight).sum(axis=0) / w
    resid = weight * (data - mean)

    cov = resid.T.dot(resid) / np.outer(w, w)

    var = np.diagonal(cov)
    sigma = np.sqrt(np.outer(var, var))
    cor = cov / sigma

    cor_smooth = np.zeros([nbins, nbins])
    dcor = {}
    dncor = {}

    for i in range(nbins):
        for j in range(i + 1, nbins):
            idrp = round(abs(rp[j] - rp[i]) / drp)
            idrt = round(abs(rt[i] - rt[j]) / drt)
            if not (idrp, idrt) in dcor:
                dcor[(idrp, idrt)] = 0.0
                dncor[(idrp, idrt)] = 0

            dcor[(idrp, idrt)] += cor[i, j]
            dncor[(idrp, idrt)] += 1

    for i in range(nbins):
        cor_smooth[i, i] = 1.0
        for j in range(i + 1, nbins):
            idrp = round(abs(rp[j] - rp[i]) / drp)
            idrt = round(abs(rt[i] - rt[j]) / drt)
            cor_smooth[i, j] = dcor[(idrp, idrt)] / dncor[(idrp, idrt)]
            cor_smooth[j, i] = cor_smooth[i, j]

    return cor_smooth * sigma


def weighted_cov(X, W):
    """
    x,w shape is (num_variables, num_observations)
    """
    xx0 = X - np.average(X, weights=W, axis=1)[:, np.newaxis]
    return np.cov(xx0*np.sqrt(W))


def explicit_cov(X, W):
    xw = (X * W)
    xi = xw.sum(axis=1) / W.sum(axis=1)
    first = xw[:, np.newaxis, :] * xw[np.newaxis, :, :]
    second = xi[:, np.newaxis, np.newaxis] * xi[:, np.newaxis, np.newaxis]
    numer = (first - second).sum(axis=2)
    denom = (W[:, np.newaxis, :] * W[np.newaxis, :, :]).sum(axis=2)
    cov = numer / denom
    return cov


def nested_for_cov(xi_s, w_s):
    assert xi_s.shape == w_s.shape
    N, M = xi_s.shape

    xi = np.zeros(N)
    w = np.zeros(N)
    for a in xrange(N):
        for s in xrange(M):
            xi[a] += w_s[a, s]*xi_s[a, s]
            w[a] += w_s[a, s]
        xi[a] /= w[a]

    cov = np.zeros((N, N))
    for a in xrange(N):
        for b in xrange(N):
            xi_sq = xi[a]*xi[b]
            for s in xrange(M):
                cov[a, b] += w_s[a, s]*w_s[b, s]*(xi_s[a, s]*xi_s[b, s] - xi_sq)
            cov[a, b] /= w[a]*w[b]

    return xi, cov


def save_xi(filename, xi, float_fmt='%.18e'):
    nbins = xi.shape[0]
    indices = np.arange(nbins)
    outdata = np.array([indices, xi])
    outfmt = ' '.join(['%d', float_fmt])
    np.savetxt(filename, outdata.transpose(), fmt=outfmt)


def save_cov(filename, cov, float_fmt='%.18e'):
    num_entries = cov.shape[0]
    indices = np.arange(num_entries)
    xi, yi = np.meshgrid(indices, indices)
    tri_indices = np.tril_indices(num_entries)
    outdata = np.array([xi[tri_indices], yi[tri_indices], cov[tri_indices]])
    outfmt = ' '.join(['%d', '%d', float_fmt])
    np.savetxt(filename, outdata.transpose(), fmt=outfmt)


def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', type=str, default=None,
        help='filename glob pattern')
    parser.add_argument('--output', type=str, default=None,
        help='output filename for combined data')
    args = parser.parse_args()


    filenames = glob.glob(args.input)
    nsamples = len(filenames)
    print('Reading {:d} files matching pattern: "{}"'.format(
        nsamples, args.input))

    first = np.loadtxt(filenames[0])
    nbins, ncols = first.shape

    # index, xi, weight, npairs
    assert ncols == 4

    subsamples = np.empty((nsamples, nbins, ncols))
    for i, filename in enumerate(tqdm(filenames, total=nsamples)):
        subsamples[i] = np.loadtxt(filename)

    xi_s = subsamples[:, :, 1]
    weight_s = subsamples[:, :, 2]
    npair_s = subsamples[:, :, 3]

    print('Computing xi...')

    xi = np.average(xi_s, weights=weight_s, axis=0)

    print('Saving xi...')
    save_xi(args.output + '.data', xi)

    # TODO: set via command line args, config, or something else...
    r = 2 + np.arange(0, 200, 4)
    xx, yy = np.meshgrid(r, r)
    rr = xx*xx + yy*yy

    print('Computing cov...')
    cov = smooth_cov(xi_s, weight_s, yy.ravel(), xx.ravel())

    # check pos def
    sign, det = np.linalg.slogdet(cov)
    print('sign, logdet: {}, {}'.format(sign, det))

    # check symmetric
    print('symmatric: {}'.format(np.allclose(cov.transpose(), cov)))

    print('Saving covariance matrix...')
    save_cov(args.output + '.cov', cov)



def plot_2d_xi():
    # plot 2d xi image
    fig, ax = plt.subplots(figsize=(8, 6))

    xx, yy = np.meshgrid(2 + np.arange(0, 200, 4), 2 + np.arange(0, 200, 4))
    rr = xx*xx + yy*yy

    rsq_xi = rr*avg_xi.reshape(50, 50)
    ma_rsq_xi = ma.masked_where((rr < 10**2) | (rr > 190**2), rsq_xi)

    mask_radius = 185
    circ = patches.Circle((0, 0), mask_radius, facecolor='none')
    ax.add_patch(circ)  # Plot the outline

    im = ax.imshow(ma_rsq_xi, origin='lower',
                   extent=(0, 200, 0, 200), interpolation='none',
                   cmap='bwr',
                   clip_path=circ, clip_on=True)

    ax.set_xlim(0, mask_radius)
    ax.set_ylim(0, mask_radius)

    ax.set_ylabel(r'$r_\parallel h^{-1} \mathrm{Mpc}$', fontsize=16)
    ax.set_xlabel(r'$r_\perp h^{-1} \mathrm{Mpc}$', fontsize=16)

    ax.tick_params(direction='in', labelsize=16)

    fig.colorbar(im, ax=ax)


def plot_npair_per_bin():
    # xi bin pair count
    fig, axes = plt.subplots(nrows=2, figsize=(12, 8))

    ax = axes[0]
    ax.plot(-np.log10(np.sum(npair_s, axis=1)))

    ax.set_xlabel('xi bin index', fontsize=16)
    ax.set_ylabel('-log(npair)', fontsize=16)

    ax.tick_params(direction='in', labelsize=14)

    ax = axes[1]
    ax.plot(-np.log10(np.sum(npair_s, axis=0)))

    ax.set_xlabel('healpix bin index', fontsize=16)
    ax.set_ylabel('-log(npair)', fontsize=16)

    ax.tick_params(direction='in', labelsize=14)


def plot_cov_diag():
    I, J, V = cov_dr11_raw.T
    cov = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr().toarray()

    fig, axes = plt.subplots(nrows=2, figsize=(12, 8), sharex=True)

    k = 0

    ax = axes[0]
    ax.plot(np.log10(np.diag(cov, k=k)), label='DR11')
    ax.plot(np.log10(np.diag(for_cov, k=k)), label='DM Mock-000')
    ax.set_ylabel('log(C)', fontsize=16)
    ax.legend()
    ax.tick_params(direction='in', labelsize=16)
    ax.grid(True)

    ax = axes[1]
    ax.plot(np.diag(for_cov, k=k) / np.diag(cov, k=k), color='C2', label='DM Mock-000 / DR11')
    ax.set_ylabel('Ratio', fontsize=16)
    ax.set_xlabel('Covariance diagonal index', fontsize=16)
    ax.legend()
    ax.tick_params(direction='in', labelsize=16)
    ax.grid(True)


def plot_cov():
    # see http://stackoverflow.com/questions/13784201/matplotlib-2-subplots-1-colorbar
    fig, axes = plt.subplots(
        nrows=2, ncols=2, figsize=(14, 14),
        sharex=True, sharey=True)

    imshow_kwargs = {
        'cmap': 'Purples',
        'origin': 'upper',
        'interpolation': 'none',
        'vmin': -13,
        'vmax': -7
    }

    s = slice(None, None)

    I, J, V = cov_dr11_raw.T
    cov = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr().toarray()
    print('min(cov): {}'.format(np.min(cov)))
    print('max(cov): {}'.format(np.max(cov)))

    ax = axes[0, 0]
    im = ax.imshow(ma.log10(cov[s, s]), **imshow_kwargs)
    ax.tick_params(direction='in')
    ax.set_ylabel('log(+C)', fontsize=18)
    # fig.colorbar(im, ax=ax)

    ax = axes[1, 0]
    im = ax.imshow(ma.log10(-cov[s, s]), **imshow_kwargs)
    ax.tick_params(direction='in')
    ax.set_ylabel('log(-C)', fontsize=18)
    ax.set_xlabel('DR11', fontsize=18)
    # fig.colorbar(im, ax=ax)

    # I, J, V = cov_mock000_raw.T
    # cov = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr().toarray()
    cov = np.tril(smo_cov)
    print('min(cov): {}'.format(np.min(cov)))
    print('max(cov): {}'.format(np.max(cov)))

    ax = axes[0, 1]
    im = ax.imshow(ma.log10(cov[s, s]), **imshow_kwargs)
    ax.tick_params(direction='in')
    # fig.colorbar(im, ax=ax)

    ax = axes[1, 1]
    im = ax.imshow(ma.log10(-cov[s, s]), **imshow_kwargs)
    ax.tick_params(direction='in')
    ax.set_xlabel('DM', fontsize=18)
    # fig.colorbar(im, ax=ax)

    plt.tight_layout()


if __name__ == '__main__':
    main()