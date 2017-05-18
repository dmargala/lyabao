#!/bin/sh
set -ex

NAME="$1"

mkdir -p /data/boss/lya/${NAME}/xi-2d || true
mkdir -p /data/boss/lya/${NAME}/baofit || true

time ~/source/turbo-octo-spice/build/h5healxi \
    -i /data/boss/lya/${NAME}/${NAME}-delta.hdf5 \
    --verbose \
    --cart \
    --axis1 [0:200]*50 \
    --axis2 [0:200]*50 \
    --nthreads 10 \
    --order 5 \
    -o /data/boss/lya/${NAME}/xi-2d/xi-2d-${NAME}-healpix \
    --save-subsamples

cd /data/boss/lya/${NAME}/xi-2d

mv xi-2d-${NAME}-healpix.cov xi-2d-${NAME}-healpix.cov.orig

time ~/source/turbo-octo-spice/python/est_cov.py \
    --name "xi-2d-${NAME}-healpix-*" \
    --save xi-2d-${NAME}-healpix.cov

cd /data/boss/lya/${NAME}/baofit

time baofit \
    --data /data/boss/lya/${NAME}/xi-2d/xi-2d-${NAME}-healpix \
    -i /home/dmargala/source/turbo-octo-spice/comoving.ini \
    --modelroot /home/dmargala/source/baofit/models/ \
    --dilmax 4 --dilmin 0.1 \
    --cov-sample-size 2993

