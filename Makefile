export SHELL:=/bin/bash
.ONESHELL:
UNAME:=$(shell uname)

define help
endef
export help
help:
	@printf "$$help"
.PHONY : help

# ~~~~~ Install Dependencies ~~~~~ #
# need to use Python 2.7 because Python 3 gives different results
export PATH:=$(CURDIR)/conda/bin:$(CURDIR)/scripts:$(CURDIR):$(PATH)
unexport PYTHONPATH
unexport PYTHONHOME

ifeq ($(UNAME), Darwin)
CONDASH:=Miniconda2-4.7.12.1-MacOSX-x86_64.sh
endif

ifeq ($(UNAME), Linux)
CONDASH:=Miniconda2-4.7.12.1-Linux-x86_64.sh
endif

CONDAURL:=https://repo.anaconda.com/miniconda/$(CONDASH)

conda:
	echo ">>> Setting up conda..."
	wget "$(CONDAURL)"
	bash "$(CONDASH)" -b -p conda
	rm -f "$(CONDASH)"

export NXF_VER:=20.07.1
./nextflow:
	if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
	curl -fsSL get.nextflow.io | bash

# https://github.com/mskcc/roslin-variant/blob/2.6.x/build/containers/conpair/0.3.3/Dockerfile
install: conda ./nextflow
	pip install scipy==1.1.0 numpy==1.15.4

THREADS:=8
NUM_TUMORS:=1
NUM_NORMALS:=5
TUMOR_FILE:=tumors.txt
NORMAL_FILE:=normals.txt
MARKERS:=/juno/work/ci/kellys5/projects/conpair-dev/markers/IMPACT468/FP_tiling_genotypes_for_Conpair.txt
CONCORDANCE_FILE:=concordance.tsv
run:
	python run.py concordance \
	--tumors-list "$(TUMOR_FILE)" \
	--normals-list "$(NORMAL_FILE)" \
	--num-tumors "$(NUM_TUMORS)" \
	--num-normals "$(NUM_NORMALS)" \
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
	--pileup-file "$(NORMAL_FILE)" \
	--output-dir "$(OUTPUT_DIR)" \
	--markers "$(MARKERS)"

# run a .bam -> .pileup -> .pickle workflow
WORKFLOW_DIR:=$(CURDIR)/workflow
export NXF_WORK:=$(WORKFLOW_DIR)/work
export NXF_PID_FILE:=$(WORKFLOW_DIR)/.nextflow.pid
export NXF_LOG:=$(WORKFLOW_DIR)/nextflow.log
BAM_DIR:=$(CURDIR)/bams
GATK_JAR:=/juno/work/ci/kellys5/projects/conpair-dev/gatk.jar
REF_FASTA:=/juno/work/ci/resources/genomes/GRCh37/fasta/b37.fasta
MARKERS_BED:=$(CURDIR)/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed
MARKERS_TXT:=$(CURDIR)/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt
preprocessing-workflow:
	if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
	nextflow -log "$(NXF_LOG)" run \
	-resume \
	$(WORKFLOW_DIR)/preprocessing-workflow.nf \
	-profile preprocessing \
	--input_dir $(BAM_DIR) \
	--output_dir $(OUTPUT_DIR) \
	--gatk_jar $(GATK_JAR) \
	--ref_fasta $(REF_FASTA) \
	--markers_bed $(MARKERS_BED) \
	--markers_txt $(MARKERS_TXT)
.PHONY:workflow

concordance-workflow:
	if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
	nextflow -log "$(NXF_LOG)" run \
	-resume \
	$(WORKFLOW_DIR)/concordance-workflow.nf \
	-profile concordance \
	--tumors_list "$(TUMOR_FILE)" \
	--normals_list "$(NORMAL_FILE)" \
	--markers_txt "$(MARKERS)"

clean:
	rm -f $(WORKFLOW_DIR)/*.log
	rm -rf .nextflow
	rm -rf $(NXF_WORK)
	rm -f *.html.*
	rm -f *trace.txt.*

# ~~~~~ #
# python ../Conpair/scripts/verify_concordances.py -p pairing.txt -N ... -T ...
# python ../Conpair/scripts/estimate_tumor_normal_contaminations.py -p pairing.txt -N ... -T ...
# export GATK_JAR=/path/to/gatk.jar
# TUMOR_PILEUP:=$(shell head -1 $(TUMOR_FILE))
# NORMAL_PILEUP:=$(shell head -1 $(NORMAL_FILE))
# test-run1:
# 	export CONPAIR_DIR=$(CURDIR)
# 	python2.7 scripts/verify_concordance.py \
# 	--tumor_pileup $(TUMOR_PILEUP) \
# 	--normal_pileup $(NORMAL_PILEUP) && \
# 	python2.7 scripts/estimate_tumor_normal_contamination.py \
# 	--tumor_pileup $(TUMOR_PILEUP) \
# 	--normal_pileup $(NORMAL_PILEUP)
# # 0.782
# # Based on 188/7387 markers (coverage per marker threshold: 10 reads)
# # Minimum mappinq quality: 10
# # Minimum base quality: 20
# # Normal sample contamination level: 0.113%
# # Tumor sample contamination level: 0.0%
#
# test-run2:
# 	export CONPAIR_DIR=$(CURDIR)
# 	python2.7 scripts/verify_concordance2.py \
# 	--tumor_pileup $(TUMOR_PILEUP) \
# 	--normal_pileup $(NORMAL_PILEUP)
