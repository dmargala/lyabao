#!/usr/bin/env python

import math
import numpy as np

from hypothesis import given, assume
from hypothesis.strategies import floats

from astropy.cosmology import FlatLambdaCDM

cosmology = FlatLambdaCDM(H0=70, Om0=0.27)


class SkyObject(object):
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        self.sin_ra = math.sin(ra)
        self.cos_ra = math.cos(ra)
        self.sin_dec = math.sin(dec)
        self.cos_dec = math.cos(dec)

    def get_cos_separation(self, other):
        a = self.sin_dec * other.sin_dec
        b = self.cos_dec * other.cos_dec
        c = self.sin_ra * other.sin_ra
        d = self.cos_ra * other.cos_ra
        return a + b*(c + d)


class LineOfSightObject(SkyObject):
    def __init__(self, ra, dec, z):
        SkyObject.__init__(self, ra, dec)
        self.z = z
        self.r = cosmology.comoving_distance(z).value * cosmology.h

    def get_separation_components(self, other):

        cos_sep = self.get_cos_separation(other);

        r1_sq = self.r * self.r
        two_r1_cos_sep = 2 * self.r * cos_sep;
        r_sq = r1_sq + (other.r - two_r1_cos_sep) * other.r;

        dr_para = abs(self.r - other.r);
        dr_perp = math.sqrt(abs(r_sq - dr_para*dr_para));

        return dr_para, dr_perp


ra_min, ra_max = 0, 360
dec_min, dec_max = -90, 90
z_min, z_max = 1, 4

eps1_min, eps1_max = -1, 1
threshold1 = 1e-2
cases1 = []

@given(
    ra=floats(min_value=ra_min, max_value=ra_max),
    dec=floats(min_value=dec_min, max_value=dec_max),
    z=floats(min_value=z_min, max_value=z_max),
    eps=floats(min_value=eps1_min, max_value=eps1_max))
def test_trig1(ra, dec, z, eps):
    assume(abs(eps) > 1e-10)
    # ra2=ra1, dec2=dec1+eps, z1=z2: dr_par=0, dr_perp = dist * eps
    ra = math.radians(ra)
    dec = math.radians(dec)
    eps = math.radians(eps)

    p1 = LineOfSightObject(ra, dec, z)
    p2 = LineOfSightObject(ra, dec + eps, z)

    test_dr_para = 0
    test_dr_perp = abs(p1.r * eps)

    dr_para, dr_perp = p1.get_separation_components(p2)

    assert abs(dr_para - test_dr_para) < threshold1
    assert abs(dr_perp - test_dr_perp) < threshold1

    cases1.append((
        p1.ra, p1.dec, p1.z,
        p2.ra, p2.dec, p2.z,
        test_dr_para, test_dr_perp, dr_para, dr_perp))


eps2_min, eps2_max = -1, 1
threshold2 = 1e-2
cases2 = []

@given(
    ra=floats(min_value=ra_min, max_value=ra_max),
    dec=floats(min_value=dec_min, max_value=dec_max),
    z=floats(min_value=z_min, max_value=z_max),
    eps=floats(min_value=eps2_min, max_value=eps2_max))
def test_trig2(ra, dec, z, eps):
    assume(abs(eps) > 1e-10)
    # ra2=ra1+eps, dec2=dec1, z1=z2: dr_par=0, dr_perp = dist * eps * cos(dec)
    ra = math.radians(ra)
    dec = math.radians(dec)
    eps = math.radians(eps)

    p1 = LineOfSightObject(ra, dec, z)
    p2 = LineOfSightObject(ra + eps, dec, z)

    test_dr_para = 0
    test_dr_perp = abs(p1.r * eps * math.cos(dec))

    dr_para, dr_perp = p1.get_separation_components(p2)

    assert abs(dr_para - test_dr_para) < threshold2
    assert abs(dr_perp - test_dr_perp) < threshold2

    cases2.append((
        p1.ra, p1.dec, p1.z,
        p2.ra, p2.dec, p2.z,
        test_dr_para, test_dr_perp, dr_para, dr_perp))

eps3_min, eps3_max = -1, 1
threshold3 = 1e-2
cases3 = []

@given(
    ra=floats(min_value=ra_min, max_value=ra_max),
    dec=floats(min_value=dec_min, max_value=dec_max),
    z=floats(min_value=z_min, max_value=z_max),
    eps=floats(min_value=eps3_min, max_value=eps3_max))
def test_trig3(ra, dec, z, eps):
    assume(abs(eps) > 1e-10)
    # ra2=ra1, dec2=dec1, z2=z1+eps: dr_par=eps * d(dist)/dz, dr_perp = 0
    ra = math.radians(ra)
    dec = math.radians(dec)

    p1 = LineOfSightObject(ra, dec, z)
    p2 = LineOfSightObject(ra, dec, z + eps)

    test_dr_para = abs(eps * (p2.r - p1.r)/(p2.z - p1.z))
    test_dr_perp = 0

    dr_para, dr_perp = p1.get_separation_components(p2)

    assert abs(dr_para - test_dr_para) < threshold3
    assert abs(dr_perp - test_dr_perp) < threshold3

    cases3.append((
        p1.ra, p1.dec, p1.z,
        p2.ra, p2.dec, p2.z,
        test_dr_para, test_dr_perp, dr_para, dr_perp))


if __name__ == '__main__':
    test_trig1()
    test_trig2()
    test_trig3()

    header = (
        'ra1 dec1 z1 '
        'ra2 dec2 z2 '
        'test_dr_para test_dr_perp dr_para dr_perp'
        )

    np.savetxt('test_trig1.txt', cases1, fmt='%.8f', header=header)
    np.savetxt('test_trig2.txt', cases2, fmt='%.8f', header=header)
    np.savetxt('test_trig3.txt', cases3, fmt='%.8f', header=header)

