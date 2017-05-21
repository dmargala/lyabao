#!/bin/sh
set -ex

NAME="$1"
ROOT_DIR=/data/boss/lya

mkdir -p ${ROOT_DIR}/${NAME}/xi-2d || true

time ~/source/turbo-octo-spice/build/h5healxi \
    -i "${ROOT_DIR}/${NAME}/${NAME}-delta.hdf5" \
    --verbose \
    --cart \
    --axis1 [0:200]*50 \
    --axis2 [0:200]*50 \
    --nthreads 10 \
    --order 5 \
    --save-subsamples \
    -o "${ROOT_DIR}/${NAME}/xi-2d/cart-healpix"

time ~/source/lyabao/bin/combine_subsamples.py \
    --input "${ROOT_DIR}/${NAME}/xi-2d/cart-healpix-*" \
    --output "${ROOT_DIR}/${NAME}/xi-2d/cart-healpix"
