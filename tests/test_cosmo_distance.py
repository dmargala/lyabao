#!/usr/bin/env python

import numpy as np
from astropy.cosmology import FlatLambdaCDM

def main():

    z, cosmo_distance = np.loadtxt('test_cosmo_distance.txt', unpack=True)
    cosmology = FlatLambdaCDM(H0=70, Om0=0.27)
    astropy_distance = cosmology.comoving_distance(z) * cosmology.h
    np.savetxt('compare_distance.txt',
               np.array((z, cosmo_distance, astropy_distance)).T,
               fmt='%.8f', header='z cosmo astropy')

if __name__ == '__main__':
    main()
