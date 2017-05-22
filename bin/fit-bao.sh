#!/bin/sh
set -ex

NAME="$1"
ROOT_DIR=/data/boss/lya

mkdir -p ${ROOT_DIR}/${NAME}/baofit || true

time baofit -i ~/source/baofit/config/BOSSDR11LyaF_k.ini \
    --modelroot ~/source/baofit/models/ \
    --data "${ROOT_DIR}/${NAME}/xi-2d/cart-healpix" \
    --output-prefix "${ROOT_DIR}/${NAME}/baofit/BOSSDR11LyaF_k_"
