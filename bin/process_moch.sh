#!/bin/bash
#$ -N mock_000
#$ -q free64
#$ -m beas
#$ -o 000.log

pwd
date

source activate lyabao
make -f makefiles/mock.mk MOCK_INDEX=000

date
