export SHELL:=/bin/bash
.ONESHELL:
UNAME:=$(shell uname)

define help
This is the Makefile to help with running Conpair on lists of tumor and normal pileups

# INSTALLATION

Install with:

```
make install
```

Run the test suite with

```
make test
```

# USAGE

Example usages of Conpair parallel concordance can be used with:

```
make run
```

Another example can be used with

```
make run2 TUMOR_FILE=tumor_pileups.txt NORMAL_FILE=normal_pileups.txt
```

Where `TUMOR_FILE` and `NORMAL_FILE` are simple text files with one filepath per line

Extra Makefile args should also be provided;

```
MARKERS=markers.txt
THREADS=<number of threads to use>
CONCORDANCE_FILE=concordance_output.tsv
```

NOTE: the Makefile includes some default values for these which should be reviewed if you are not supplying them explicitly

To make it run faster, you should run the pre-processing step to convert all the pileup files in your list to Python .pickle likelihoods files

```
make likelihoods NORMAL_FILE=list_of_bams.txt MARKERS=markers.txt OUTPUT_DIR=output
```

You can then use the paths to these .pickle files in the lists used for `TUMOR_FILE` and `NORMAL_FILE` with `make run`

## WORKFLOWS

To make it easier to process large numbers of files through all steps needed for Conpair, some workflows have been included

### PREPROCESSING

Run the pre-processing workflow to convert a directory of .bam files into Python .pickle likelihoods files to use with Conpair.

NOTE: this requires a gatk.jar file for GATK 3

```
make preprocessing-workflow
```

Extra args to supply;

```
BAM_DIR=bams/
GATK_JAR=gatk.jar
REF_FASTA=genome.fasta
MARKERS_BED=markers.bed
MARKERS_TXT=markers.txt
```

### CONCORDANCE

Run the concordance workflow to run concordance for each tumor against all supplied normals, each one submitted as a separate HPC job for higher throughput. Also includes downstream processing and aggregation of concordance for all comparisons.

```
make concordance-workflow
```

Extra args to supply;

```
TUMOR_FILE=tumor_pileups.txt # list of GATK pileup files or Python likelihoods pickle files
NORMAL_FILE=normal_pileups.txt # list of GATK pileup files or Python likelihoods pickle files
MARKERS=markers.txt
```

endef
export help
help:
	@printf "$$help"
.PHONY : help

# ~~~~~ Install Dependencies ~~~~~ #
export PATH:=$(CURDIR)/conda/bin:$(CURDIR)/scripts:$(CURDIR):$(PATH)
unexport PYTHONPATH
unexport PYTHONHOME

# need to use Python 2.7 because Python 3 gives different results
# TODO: figure out why we get different results with Python 3

# CONDASH:=Miniconda2-4.7.12.1-MacOSX-x86_64.sh
# NOTE: had to switch to Miniconda3 because Miniconda2 no longer works on newer macOS due to this;
# https://github.com/conda/conda/issues/10361
# https://stackoverflow.com/questions/65130080/attributeerror-running-django-site-on-mac-11-0-1
# we can still run Conpair under Py2.7 but we just need to use Py3 miniconda for setup

ifeq ($(UNAME), Darwin)
CONDASH:=Miniconda3-4.7.12.1-MacOSX-x86_64.sh
endif

ifeq ($(UNAME), Linux)
CONDASH:=Miniconda3-4.7.12.1-Linux-x86_64.sh
endif

CONDAURL:=https://repo.anaconda.com/miniconda/$(CONDASH)

conda:
	@printf ">>> Setting up conda...\n\n"
	wget "$(CONDAURL)"
	bash "$(CONDASH)" -b -p conda
	rm -f "$(CONDASH)"

# reference for installation;
# https://github.com/mskcc/roslin-variant/blob/2.6.x/build/containers/conpair/0.3.3/Dockerfile
install: conda
	. ./conda/bin/activate && \
	conda env create --name conpair --file environment.yml && \
	printf "\n>>> To activate, run:\nsource conda/bin/activate\nconda activate conpair\n"

# example for running with simple args
TUMOR:=data/example/pileup/*tumor*.pileup.txt
NORMAL:=data/example/pileup/*normal*.pileup.txt
run:
	python run.py concordance '$(TUMOR)' '$(NORMAL)'

# output;
# concordance     num_markers_used        num_total_markers       tumor   normal  tumor_pileup    normal_pileup
# 0.9993209289691701      7363    7387    NA12878_tumor80x        NA12878_normal40x       data/example/pileup/NA12878_tumor80x.gatk.pileup.txt    data/example/pileup/NA12878_normal40x.gatk.pileup.txt

# example for running with more complicated args
THREADS:=8
TUMOR_FILE:=tumors.txt
NORMAL_FILE:=normals.txt
MARKERS:=/juno/work/ci/kellys5/projects/conpair-dev/markers/IMPACT468/FP_tiling_genotypes_for_Conpair.txt
CONCORDANCE_FILE:=concordance.tsv
run2:
	python run.py concordance \
	--tumors-list "$(TUMOR_FILE)" \
	--normals-list "$(NORMAL_FILE)" \
	--markers "$(MARKERS)" \
	--output-file "$(CONCORDANCE_FILE)" \
	--threads "$(THREADS)" \
	--save-benchmarks


SUB_THREADS:=4 8 16 24 32
SUB_LOG_SUFFIX:=log
SUB_NUM_NORMALS:=95
SUB_NUM_TUMORS:=184
SUB_NUM_TUMORS_START:=1
SUB_STEP:=2
$(SUB_THREADS):
	bsub \
	-W 48:00 \
	-n $@ \
	-sla CMOPI \
	-oo lsf.$@.$(SUB_LOG_SUFFIX) \
	/bin/bash -c 'cd $(CURDIR); for i in $$(seq $(SUB_NUM_TUMORS_START) $(SUB_STEP) $(SUB_NUM_TUMORS)); do make run NUM_NORMALS=$(SUB_NUM_NORMALS) THREADS=$@ NUM_TUMORS=$$i; done'
.PHONY:=$(SUB_THREADS)
submit: $(SUB_THREADS)

bash:
	bash

test:
	python modules/test_concordance.py

OUTPUT_DIR:=output
$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)
likelihoods: $(OUTPUT_DIR)
	python scripts/make_genotype_likelihoods.py \
	--pileup-list "$(NORMAL_FILE)" \
	--output-dir "$(OUTPUT_DIR)" \
	--markers "$(MARKERS)"

# run a .bam -> .pileup -> .pickle workflow
# NOTE: put your path to a dir with .bam files here;
BAM_DIR:=$(CURDIR)/bams
# example;
# BAM_DIR=/juno/work/ci/helix_filters_01/fixtures/Fillout01/bam/
WORKFLOW_DIR:=$(CURDIR)/workflow
export NXF_WORK:=$(CURDIR)/work
export NXF_LOG:=$(CURDIR)/nextflow.log
REF_FASTA:=/juno/work/ci/resources/genomes/GRCh37/fasta/b37.fasta
MARKERS_BED:=$(CURDIR)/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed
MARKERS_TXT:=$(CURDIR)/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt
preprocessing-workflow:
	nextflow -log "$(NXF_LOG)" run \
	$(WORKFLOW_DIR)/preprocessing-workflow.nf \
	-profile preprocessing \
	--input_dir $(BAM_DIR) \
	--output_dir $(OUTPUT_DIR) \
	--ref_fasta $(REF_FASTA) \
	--markers_bed $(MARKERS_BED) \
	--markers_txt $(MARKERS_TXT)
.PHONY:workflow


# run the concordance workflow on all the tumors in paralle
concordance-workflow:
	nextflow -log "$(NXF_LOG)" run \
	-resume \
	$(WORKFLOW_DIR)/concordance-workflow.nf \
	-profile concordance \
	--tumors_list "$(TUMOR_FILE)" \
	--normals_list "$(NORMAL_FILE)" \
	--markers_txt "$(MARKERS)"

clean:
	rm -f $(WORKFLOW_DIR)/*.log.*
	rm -f *.log.*
	rm -f *.html.*
	rm -f *trace.txt.*

clean-all: clean
	rm -rf $(NXF_WORK)
	rm -rf .nextflow
	rm -f *.log*
	rm -f *.html*
	rm -f *trace.txt*



# ~~~~~ Container ~~~~~ #
# make the Docker container
GIT_NAME:=conpair
# GIT_TAG:=$(shell git describe --tags --abbrev=0)
GIT_TAG:=dev
DOCKER_TAG:=mskcc/$(GIT_NAME):$(GIT_TAG)
docker-build:
	docker build -t "$(DOCKER_TAG)" .

# shell into the container to check that it looks right
docker-bash:
	docker run --rm -ti "$(DOCKER_TAG)" bash

# push the container to Dockerhub
# $ docker login --username=<username>
docker-push:
	docker push "$(DOCKER_TAG)"

# pull the Dockerhub container and convert to Singularity container
# NOTE: you cannot use a filename with a ':' as a Makefile target
SINGULARITY_SIF:=mskcc_$(GIT_NAME):$(GIT_TAG).sif
singularity-pull:
	unset SINGULARITY_CACHEDIR && \
	module load singularity/3.3.0 && \
	singularity pull --force --name "$(SINGULARITY_SIF)" docker://$(DOCKER_TAG)
