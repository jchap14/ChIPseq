#!/bin/bash
# This is the submission script for "get_signal_to_noise.py"

# submit for a specific sample
# bash get_signal_to_noise.sh $SAMPLE.nmSort.bam

##### submit for all samples in CWD
# for x in `/bin/ls *.nmSort.bam` ; do bash get_signal_to_noise.sh $x; done

##### set environment
# source activate Spyder_Python3 #local only

##### specify variables to pass to get_signal_to_noise.py
BAMFILE=`echo $1`
OUTNAME=`basename $BAMFILE .nmSort.bam`
FINAL_BED=`echo $OUTNAME.bam.bed`
## set annotation & executables directory for either local or scg use
ANNO_DIR="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/" #scg
# ANNO_DIR="/Users/jchap12/Desktop/BIOINFORMATICS/Annotations/hg19/" #local
EXE_DIR="/srv/gsfs0/projects/snyder/chappell/scripts/ChIPseq/" #scg
# EXE_DIR="/Users/jchap12/Desktop/BIOINFORMATICS/ChIPseq/" #local
## Set variables to pass to python
DNASE_REGIONS=`echo $ANNO_DIR\reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz`
BLACKLIST_REGIONS=`echo $ANNO_DIR\hg19EncodeMapabilityBlacklist.bed`
PROM_REGIONS=`echo $ANNO_DIR\reg2map_honeybadger2_dnase_prom_p2.bed.gz`
ENH_REGIONS=`echo $ANNO_DIR\reg2map_honeybadger2_dnase_enh_p2.bed.gz`
PEAKS=`echo $OUTNAME.FDRe1.Peaks.bed`

## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $OUTNAME.SigToNoise
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -pe shm 1
#$ -l h_rt=0:59:00
#$ -l s_rt=0:59:00

## add modules & source specific conda environment
# source activate Spyder_Python3 #local only

## convert BAM file to BED (if test!)
if [ -f $(echo $FINAL_BED) ]
then
    echo "BED exists, continue ..."
else
    bedtools bamtobed -i $BAMFILE > $FINAL_BED
fi

## run python script
python $EXE_DIR\get_signal_to_noise.py \
--outname $OUTNAME \
--final_bed $FINAL_BED \
--dnase_regions $DNASE_REGIONS \
--blacklist_regions $BLACKLIST_REGIONS \
--prom_regions $PROM_REGIONS \
--enh_regions $ENH_REGIONS \
--peaks $PEAKS

## deactivate conda environment
# source deactivate
EOF

## qsub then remove the tempscript
#bash $NAME.tempscript.sh #local
qsub $NAME.tempscript.sh #scg
sleep 1
rm $NAME.tempscript.sh
