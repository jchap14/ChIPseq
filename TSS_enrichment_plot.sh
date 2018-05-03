#!/bin/bash
# This is the submission script for "TSS_enrichment_plot.py"

# submit for a specific sample
# bash TSS_enrichment_plot.sh $SAMPLE.bam

##### submit for all samples in CWD
# for x in `/bin/ls *.bam` ; do bash TSS_enrichment_plot.sh $x; done

##### specify variables to pass to TSS_enrichment_plot.py
BAM_FILE=`echo $1`
## THRESHOLD is not an FDR, but score based?
PREFIX=`basename $BAM_FILE .bam`
NAME=`echo $PREFIX`
## set annotation directory
ANNO_DIR="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/" #scg
# ANNO_DIR="/Users/jchap12/Downloads/" #local
TSS=`echo $ANNO_DIR\hg19_gencode_tss_unique.bed`
BP_EDGE=2000
CHROMSIZES=`echo $ANNO_DIR\hg19.chrom.sizes`
READ_LEN=101
BINS=400
PROCESSES=8
GREENLEAF_NORM="True"
#
EXE_DIR="/srv/gsfs0/projects/snyder/chappell/scripts/ChIPseq/" #scg
# EXE_DIR="/Users/jchap12/Desktop/BIOINFORMATICS/ChIPseq/" #local

## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash -l
#SBATCH --job-name $NAME.TSSplot
#SBATCH --output=$NAME.TSSplot.out
#SBATCH --mail-user jchap14@stanford.edu
#SBATCH --mail-type=ALL
# Request run time & memory
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --account=mpsnyder
#SBATCH --nodes=1

## add modules & source specific conda environment
source activate TSS_enrichment_py27
echo "The current Conda environment is" $(conda env list | grep \* | cut -f1 -d ' ')

## index bam file
if [ -f $(echo $BAM_FILE.bai) ]
then
    echo "Index Exists"
else
    samtools index $BAM_FILE    
fi


## run python script
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
sbatch $NAME.tempscript.sh #scg
sleep 1
# rm $NAME.tempscript.sh
