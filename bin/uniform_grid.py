#!/usr/bin/env python
import argparse

import numpy as np
import numpy.ma as ma
import h5py

import bossdata.path
import bossdata.remote
import bossdata.spec
import bossdata.meta

import fitsio

from tqdm import tqdm


class Target(dict):
    """
    Represents a BOSS target.

    Args:
        args: Variable length argument list.
        kwargs: Arbitrary keyword arguments.

    Raises:
        AssertionError: if 'target' key is not specified

    """
    def __init__(self, *args, **kwargs):
        super(Target, self).__init__(*args, **kwargs)
        # assert 'target' in self, \
        #     'Target: must have plate-mjd-fiber identifier key'
        # plate, mjd, fiber = self['target'].split('-')
        self['target'] = '{}-{}-{}'.format(self['plate'], self['mjd'], self['fiber'])
        self['plate'] = int(self['plate'])
        self['mjd'] = int(self['mjd'])
        self['fiber'] = int(self['fiber'])

    def to_string(self):
        """
        Returns the standard plate-mjd-fiber string represntation of the target.

        Returns:
            plate-mjd-fiber string represntation of target
        """
        return self['target']

    @classmethod
    def from_string(cls, target_string):
        """
        Returns a Target object constructed from a target string identifier.

        Args:
            target_string (str): a target string identifier.

        Returns:
            :class:`Target` object
        """
        plate, mjd, fiber = target_string.split('-')
        return cls({'plate': plate, 'mjd': mjd, 'fiber': fiber})
        # return cls({'target':target_string})

    @classmethod
    def from_plate_mjd_fiber(cls, plate, mjd, fiber):
        """
        Returns a Target object constructed from plate, mjd, fiber.

        Args:
            plate (int): target's plate id
            mjd (int): mjd observation
            fiber (int): target's fiber id

        Returns:
            :class:`Target` object
        """
        return cls({'plate': plate, 'mjd': mjd, 'fiber': fiber})
        # target_string = '-'.join([str(field) for field in [plate, mjd, fiber]])
        # return cls.from_string(target_string)


def load_target_list(filename, fields=None, verbose=False):
    """
    Loads a target data from a text file.

    The first column must be plate-mjd-fiber target identifier.
    Use the fields argument to specify additional columns to
    read. Must specify a (name, type, column index) tuple for each field.

    Args:
        filename (str): The filename to load.
        fields (list, optional): A list of columns to read, see example.
            Defaults to None.
        verbose (bool, optional): Whether or not to print verbose output.
            Defaults to False.

    Returns:
        list of :class:`Target` objects.
    """
    if fields is None:
        fields = []
    fields = [('plate', 'S4', 0), ('mjd', 'S5', 1), ('fiber', 'S4', 2)] + fields
    names, formats, cols = zip(*fields)
    if verbose:
        print('Target list: {}'.format(filename))
        print('Reading fields: {}'.format(', '.join(names)))
    targets = np.genfromtxt(
        filename, dtype={'names':names, 'formats':formats}, usecols=cols, skip_header=1)

    return [Target(dict(zip(targets.dtype.names, t))) for t in targets]


def get_fiducial_wavelength(pixel_index,  coeff1=1e-4, log10lam0=np.log10(3500.26)):
    """
    Returns the wavelength at the center of the specified index
    of the BOSS co-add fiducial wavelength grid.

    Args:
        pixel_index (int): index of the BOSS co-add fiducial wavelength grid.

    Returns:
        wavelength (float): central wavelength of the specified index on the fiducial wavelength grid
    """
    return np.power(10.0, log10lam0 + coeff1*pixel_index)

def get_fiducial_pixel_index_offset(loglam, coeff1=1e-4, log10lam0=np.log10(3500.26)):
    """
    Returns the pixel index offset from the start of the
    BOSS co-add fiducial wavelength grid.

    Args:
        coeff0 (float): central wavelength (log10) of first pixel
        coeff1 (float, optional): log10 dispersion per pixel

    Returns:
        pixel index offset from the start of the fiducial wavelength grid.
    """
    return np.round((loglam - log10lam0)/coeff1).astype(int)


class UniformGridException(Exception):
    pass

class RedshiftError(UniformGridException):
    pass

class NormError(UniformGridException):
    pass

class ForestError(UniformGridException):
    pass


class FullSpecFile(object):
    '''
    Represents a full spec file object.

    Parameters
    ----------
    path: the filesystem path to the FITS file
    '''
    def __init__(self, path):
        self.hdulist = fitsio.FITS(path, mode=fitsio.READONLY)
        self.header = self.hdulist[0].read_header()

    def get_data(self):
        '''
        Returns the observed flux, ivar for this spectrum.
        Masks pixels indicated by `and_mask` and `ivar <= 0`

        '''
        # Look up the HDU for this spectrum and its pixel quality bitmap.
        hdu = self.hdulist[1]
        pixel_bits = hdu['and_mask'][:]
        num_pixels = len(pixel_bits)

        # Identify the pixels with valid data.
        ivar = hdu['ivar'][:]
        bad_pixels = (pixel_bits != 0) | (ivar <= 0.0)
        good_pixels = ~bad_pixels

        # Create and fill the unmasked structured array of data.
        dtype = [
            ('loglam', np.float32),
            ('flux', np.float32),
            ('ivar', np.float32)
        ]
        data = np.empty(num_pixels, dtype=dtype)
        data['loglam'][:] = hdu['loglam'][:]
        data['flux'][:] = hdu['flux'][:]
        data['ivar'][:] = ivar[:]

        return ma.MaskedArray(data, mask=bad_pixels)


class MetaDataManager(object):

    meta_dtype = [
        ('ra', np.float32), ('dec', np.float32), ('plate', np.int32),
        ('mjd', np.int32), ('fiber', np.int32), ('thing_id', np.int64)
    ]

    def __init__(self, ntargets):
        self.meta_data = np.empty(ntargets,
            dtype=MetaDataManager.meta_dtype)

    def save_meta(self, i, target, spec):
        '''
        Saves the meta data for (target, spec) at i
        '''
        self.meta_data[i] = MetaDataManager.get_meta_data(target, spec)

    @staticmethod
    def get_meta_data(target, spec):
        '''
        Extracts meta data from target and spec file
        '''
        ra = np.radians(np.asscalar(spec.hdulist[2]['RA'][:]))
        dec = np.radians(np.asscalar(spec.hdulist[2]['DEC'][:]))
        plate = target['plate']
        mjd = target['mjd']
        fiber = target['fiber']
        thing_id = np.asscalar(spec.hdulist[2]['THING_ID'][:])

        meta_data = (ra, dec, plate, mjd, fiber, thing_id)

        assert len(meta_data) == len(MetaDataManager.meta_dtype)

        return meta_data

    @staticmethod
    def get_redshift(target, spec, quasar_catalog=None):
        '''
        Returns the redshift of target from spec, unless a quasar catalog
        is provided, then that is used to look up `Z_VI`. If target is not
        in the quasar catalog, the value `None` is returned.
        '''
        if quasar_catalog is None:
            z = np.asscalar(spec.hdulist[2]['Z'][:])
        else:
            where = 'PLATE={} and MJD={} and FIBER={}'.format(
                target['plate'], target['mjd'], target['fiber'])
            result = quasar_catalog.select_all(
                where=where, what='Z_VI', max_rows=1)
            if result is None:
                raise RedshiftError(
                    '{}: No redshift available for target'.format(
                        target['target']))
            z = result['Z_VI'][0]
        return z

class UniformGrid(object):
    '''Represents a line of sight on a uniform grid.

    Args:
        forest_lim (float, float): A pair of floats indicating the wavelength
            limits of the lya forest
        norm_lim (float, float): A pair for floats indicating the wavelength
            limits of the wavelength range to use for normalization
    '''
    def __init__(self, forest_lim, norm_lim, max_fid_index):
        pass


def get_loglam_index(loglam):
    loglam_index = (
        get_fiducial_pixel_index_offset(loglam)
    )
    return np.round(loglam_index).astype(int)


def get_index_lim(rest_wave_lim, z, clip=None):
    log_lo = np.log10(rest_wave_lim[0] * (1.0 + z))
    log_hi = np.log10(rest_wave_lim[1] * (1.0 + z))
    if clip is not None:
        log_lo = max(log_lo, clip[0])
        log_hi = min(log_hi, clip[-1])
    lo_index = get_loglam_index(log_lo)
    hi_index = get_loglam_index(log_hi)
    return (lo_index, hi_index)

def get_slices(target, forest_lim, norm_lim, z, loglam, ivar):

    offset = get_fiducial_pixel_index_offset(
        loglam.data[0])

    # determine fiducial wavelength offsets of observed forest pixels
    forest_lo_index, forest_hi_index = get_index_lim(
        forest_lim, z, clip=loglam.data)

    # check to see if the observed wavelength range overlaps the forest
    if forest_lo_index >= forest_hi_index:
        raise ForestError(
            '{}: no forest pixels [{}:{}], z = {}'.format(
                target['target'],  forest_lim[0], forest_lim[1], z))

    # the uniform wavelength grid slice to use for this observation
    forest_slice = slice(forest_lo_index, forest_hi_index)

    spec_slice = slice(forest_lo_index - offset, forest_hi_index - offset)
    # check for unmasked pixels in forest window
    if ma.sum(ivar[spec_slice].mask) == len(ivar[spec_slice]):
        raise ForestError(
            '{}: no unmasked pixels in forest [{}:{}], z = {}'.format(
                target['target'], forest_lim[0], forest_lim[1], z))

    # find normalization window
    norm_lo_index, norm_hi_index = get_index_lim(
        norm_lim, z)
    norm_slice = slice(norm_lo_index - offset, norm_hi_index - offset)

    return forest_slice, spec_slice, norm_slice

def get_norm(target, flux, ivar, norm_slice):
    # check okay ivar sum in range
    if np.sum(ivar[norm_slice].data) <= 0:
        raise NormError(
            '{}: No good pixels in norm window [{}]'.format(
                target['target'],  norm_slice))

    # calculate normalization as ivar weighted mean flux
    norm = ma.average(flux[norm_slice].data, weights=ivar[norm_slice].data)

    # verify normalization is valid
    if norm <= 0:
        raise NormError(
            '{}: Norm must be positive, norm = {}'.format(
                target['target'], norm))

    return norm

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action='store_true',
        help='print verbose output')
    parser.add_argument('-i', '--input', type=str, default=None,
        help='target list')
    parser.add_argument('-n', '--ntargets', type=int, default=0,
        help='number of targets to use, 0 for all')
    parser.add_argument('--mock', action='store_true',
        help='mock input data')
    parser.add_argument('--name', type=str, default=None,
        help='output file name base')
    parser.add_argument('--dry-run', action='store_true',
        help='dont save anything')
    # analysis options
    parser.add_argument('--forest-lo', type=float, default=1040,
        help='min forest wavelength')
    parser.add_argument('--forest-hi', type=float, default=1200,
        help='max forest wavelength')
    parser.add_argument('--norm-lo', type=float, default=1275,
        help='min forest wavelength')
    parser.add_argument('--norm-hi', type=float, default=1285,
        help='max forest wavelength')
    parser.add_argument('--max-fid-index', type=int, default=1900,
        help='Maximum fiducial pixel index')
    parser.add_argument('--spall-redshift', action='store_true',
        help='Use redshift from hdu 2 of spec file, instead of quasar catalog')
    args = parser.parse_args()

    if not args.dry_run and args.name is None:
        raise RuntimeError('Must specify output file name base with --name or '
                           'explicity set --dry-run for a "dry run"')

    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()

    # read target data
    targets = load_target_list(args.input, verbose=args.verbose)

    # use the first n targets
    ntargets = args.ntargets if args.ntargets > 0 else len(targets)
    targets = targets[:ntargets]

    # determine maximum forest wavelength observed
    skim_redshift = np.empty(ntargets)
    quasar_catalog = None
    if not args.spall_redshift:
        quasar_catalog = bossdata.meta.Database(quasar_catalog=True, lite=False)

    # convert max observed wavelength to fiducial pixel index
    max_index = args.max_fid_index

    skim_loglam = np.log10(get_fiducial_wavelength(np.arange(max_index)))

    # arrays for skim data
    skim_flux = ma.zeros((ntargets, max_index))
    skim_flux.mask = True
    skim_ivar = ma.zeros((ntargets, max_index))
    skim_ivar.mask = True
    skim_norm = np.zeros(ntargets)

    # target meta data
    md_manager = MetaDataManager(ntargets)

    bad_targets = {}
    def log_bad_target(kind, index, msg):
        if args.verbose:
            print(kind, i, msg)
        bad_targets.get(kind, []).append((index, msg))


    for i, target in tqdm(enumerate(targets), total=ntargets):
        try:
            # get spec filename and load data
            get_spec_path_args = {
                'plate': target['plate'],
                'mjd': target['mjd'],
                'fiber': target['fiber'],
                'lite': False
            }
            if args.mock:
                get_spec_path_args['mock'] = True
                get_spec_path_args['compressed'] = True
            remote_path = finder.get_spec_path(**get_spec_path_args)
            local_path = mirror.get(remote_path, progress_min_size=0.1)
            spec = FullSpecFile(local_path)

            # save meta data
            md_manager.save_meta(i, target, spec)

            # look up redshift
            z = MetaDataManager.get_redshift(target, spec, quasar_catalog)

            skim_redshift[i] = z

            # process spectrum data
            data = spec.get_data()
            loglam = data['loglam'][:]
            flux = data['flux'][:]
            ivar = data['ivar'][:]

            # calculate array slices for various regions
            forest_slice, spec_slice, norm_slice = get_slices(
                target, (args.forest_lo, args.forest_hi),
                (args.norm_lo, args.norm_hi), z, loglam, ivar)

            # calculate normalization
            norm = get_norm(target, flux, ivar, norm_slice)

            # save normalization
            skim_norm[i] = norm

            # copy forest data to skim
            skim_flux[i, forest_slice] = flux[spec_slice]
            skim_ivar[i, forest_slice] = ivar[spec_slice]

        except RedshiftError as e:
            log_bad_target('redshift', i, str(e))
        except NormError as e:
            log_bad_target('norm', i, str(e))
        except ForestError as e:
            log_bad_target('forest', i, str(e))

    # print summary of "bad" targets
    for kind in bad_targets.keys():
        print('Number of spectra with bad {}: {:d}'.format(
            kind, len(bad_targets[kind])))

    # verify flux and ivar masks are equal
    assert np.all(skim_flux.mask == skim_ivar.mask)

    if not args.dry_run:

        # only save rows that have unmasked pixels
        masked_rows = (np.sum(skim_ivar.mask, axis=1) == max_index)
        save_rows = ~masked_rows

        print('Saving {:d} rows...'.format(np.sum(save_rows)))

        outfile = h5py.File(args.name + '-skim.hdf5', 'w')
        # save args
        outfile.attrs['forest_lo'] = args.forest_lo
        outfile.attrs['forest_hi'] = args.forest_hi
        outfile.attrs['norm_lo'] = args.norm_lo
        outfile.attrs['norm_hi'] = args.norm_hi
        # save the index of the maximum observed forest wavelength
        outfile.attrs['max_fid_index'] = max_index
        outfile.attrs['coeff0'] = skim_loglam[0]
        outfile.attrs['coeff1'] = 1e-4
        outfile.attrs['wave_min'] = 0

        eps = 1e-15
        skim_flux.data[np.abs(skim_flux.data) < eps] = 0
        skim_ivar.data[np.abs(skim_ivar.data) < eps] = 0

        compression = 'gzip'
        # save uniform wavelength grid
        outfile.create_dataset('loglam',
                               data=skim_loglam,
                               compression=compression)
        # save pixel flux, ivar, and mask
        outfile.create_dataset('flux',
                               data=skim_flux.data[save_rows],
                               compression=compression)
        outfile.create_dataset('ivar',
                               data=skim_ivar.data[save_rows],
                               compression=compression)
        outfile.create_dataset('mask',
                               data=skim_ivar.mask[save_rows],
                               compression=compression)
        # save redshifts from input target list
        outfile.create_dataset('z',
                               data=skim_redshift[save_rows],
                               compression=compression)
        # save additional quantities
        outfile.create_dataset('norm',
                               data=skim_norm[save_rows],
                               compression=compression)
        # save target meta data
        outfile.create_dataset('meta',
                               data=md_manager.meta_data[save_rows],
                               compression=compression)
        # save masked target meta data
        outfile.create_dataset('masked_meta',
                               data=md_manager.meta_data[masked_rows],
                               compression=compression)

        outfile.close()


if __name__ == '__main__':
    main()
