# lyabao
lya bao analysis

## Environment setup

Create a new conda env with required dependencies

```
conda env create -f environment.yml
```

Activate conda environement

```
source activate lyabao
```

## Run pipeline

Make sure to set `bossdata` settings:

```
export BOSS_LOCAL_ROOT=/share/dm/all
export BOSS_SAS_PATH=/sas/dr12/boss
export BOSS_REDUX_VERSION=v5_7_0
export BOSS_DATA_URL=http://dr12.sdss3.org
```

Run test

```
make -f makefiles/test.mk
```

Run DR12 analysis

```
make -f makefiles/dr12.mk
```
