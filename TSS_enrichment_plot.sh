#!/bin/bash
# This is the submission script for "TSS_enrichment_plot.py"

# submit for a specific sample
# bash TSS_enrichment_plot.sh $SAMPLE.bam

##### submit for all samples in CWD
# for x in `/bin/ls *.bam` ; do bash TSS_enrichment_plot.sh $x; done

##### set environment
source activate TSS_enrichment_py27

##### specify variables to pass to TSS_enrichment_plot.py
BAM_FILE=`echo $1`
## THRESHOLD is not an FDR, but score based?
PREFIX=`basename $BAM_FILE .bam`
## set annotation directory
ANNO_DIR="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/" #scg
# ANNO_DIR="/Users/jchap12/Downloads/" #local
TSS=`echo $ANNO_DIR\hg19_gencode_tss_unique.bed`
BP_EDGE=2000
CHROMSIZES=`echo $ANNO_DIR\hg19.chrom.sizes`
READ_LEN=102
BINS=400
PROCESSES=8
GREENLEAF_NORM="True"
#
NAME=`echo $PREFIX`
EXE_DIR="/srv/gsfs0/projects/snyder/chappell/scripts/ChIPseq/" #scg
# EXE_DIR="/Users/jchap12/Desktop/BIOINFORMATICS/ChIPseq/" #local

## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.TSSplot
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -pe shm 1
#$ -l h_rt=0:59:00
#$ -l s_rt=0:59:00

## add modules & source specific conda environment
source activate TSS_enrichment_py27

## run script
python $EXE_DIR\TSS_enrichment_plot.py \
--bam_file $BAM_FILE \
--prefix $PREFIX \
--tss $TSS \
--bp_edge $BP_EDGE \
--chromsizes $CHROMSIZES \
--read_len $READ_LEN \
--bins $BINS \
--processes $PROCESSES \
--greenleaf_norm $GREENLEAF_NORM
                    

## deactivate conda environment
source deactivate
EOF

## qsub then remove the tempscript
# bash $NAME.tempscript.sh #local
qsub $NAME.tempscript.sh #scg
sleep 1
rm $NAME.tempscript.sh
