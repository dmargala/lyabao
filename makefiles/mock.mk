MOCK_INDEX=000
NAME=output/mock_$(MOCK_INDEX)/mock_$(MOCK_INDEX)
OUTDIR=output/mock_$(MOCK_INDEX)

all: delta_field inspect


TARGETLIST=data/mock_goodtargets_1_10.txt
NTARGETS=0
FORESTLO=1040
FORESTHI=1200
NORMLO=1275
NORMHI=1285
MAXFIDINDEX=1900
SPALLREDSHIFT="--spall-redshift"
MOCK="--mock"


$(OUTDIR):
	-mkdir -p $(OUTDIR)

$(NAME)-skim.hdf5: bin/uniform_grid.py $(TARGETLIST) $(OUTDIR)
	export BOSS_LOCAL_ROOT=/share/dm/all; \
	export BOSS_SAS_PATH=/sas/dr12/boss; \
	export BOSS_REDUX_VERSION=M3_0_0/$(MOCK_INDEX); \
	python bin/uniform_grid.py --name $(NAME) \
		-i $(TARGETLIST) \
		--verbose \
		--ntargets $(NTARGETS) \
		--forest-lo $(FORESTLO) \
		--forest-hi $(FORESTHI) \
		--norm-lo $(NORMLO) \
		--norm-hi $(NORMHI) \
		--max-fid-index $(MAXFIDINDEX) \
		$(SPALLREDSHIFT) \
		$(MOCK)

NUMCOMBINE=3
WAVEMIN=3600.0

$(NAME)-cskim.hdf5: bin/analysis_prep.py $(NAME)-skim.hdf5
	python bin/analysis_prep.py --name $(NAME) \
		--num-combine $(NUMCOMBINE) \
		--wave-min $(WAVEMIN)


SUBSAMPLESTEP=1000
WAVELYA=1216.0
FORESTMAXZ=3.5

$(NAME)-forest.hdf5: bin/restframe_work.py $(NAME)-cskim.hdf5
	python bin/restframe_work.py --name $(NAME) \
		--subsample-step $(SUBSAMPLESTEP) \
		--wave-lya $(WAVELYA) \
		--forest-max-z $(FORESTMAXZ)


ABSALPHA=0.0018
ABSBETA=3.92
FORESTWAVEREF=1185.0

$(NAME)-linear-continuum.hdf5: bin/linear_continuum.py $(NAME)-forest.hdf5
	python bin/linear_continuum.py --name $(NAME) \
		--abs-alpha $(ABSALPHA) \
		--abs-beta $(ABSBETA) \
		--forest-wave-ref $(FORESTWAVEREF)


$(NAME)-delta.hdf5: bin/save_deltas.py $(NAME)-linear-continuum.hdf5
	python bin/save_deltas.py --name $(NAME) \
		--subsample-step $(SUBSAMPLESTEP)


delta_field: $(NAME)-delta.hdf5

$(NAME)-delta-scatter.png: bin/inspect_deltas.py $(NAME)-delta.hdf5
	python bin/inspect_deltas.py --name $(NAME)

inspect: $(NAME)-delta-scatter.png

.PHONY: all delta_field inspect
