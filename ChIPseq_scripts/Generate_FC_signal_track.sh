#!/bin/bash

#$ -j y
#$ -cwd
#$ -pe shm 12
#$ -V
#$ -l h_vmem=4G

################ Generate FC signal track from MACS2 ChIP-seq output #####################
## qsub Generate_ChIP_FC_signal_track.sh $.treat_pileup.bdg $.control_lambda.bdg

module add MACS2
module add ucsc_tools
module add bedtools
treatment=$1
control=$2
chromsizes="/srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/hg19.chrom.sizes"
name=`basename $treatment .treat_pileup.bdg`


##### 
# echo "macs2 bdgcmp -t $treatment -c $control --o-prefix $name -m FE"
# macs2 bdgcmp -t $treatment -c $control --o-prefix $name -m FE

echo "slopBed -i $name_FE.bdg -g $chromsizes -b 0 | bedClip stdin $chromsizes $name.fc.signal.bedGraph"
slopBed -i ${name}_FE.bdg -g $chromsizes -b 0 | bedClip stdin $chromsizes $name.fc.signal.bedGraph

echo "sort -k1,1 -k2,2n $name.fc.signal.bedGraph  > $name.fc.signal.srt.bedGraph"
sort -k1,1 -k2,2n $name.fc.signal.bedGraph  > $name.fc.signal.srt.bedGraph

echo "bedGraphToBigWig $name.fc.signal.srt.bedGraph $chromsizes $name.fc.signal.bw"
bedGraphToBigWig $name.fc.signal.srt.bedGraph $chromsizes $name.fc.signal.bw

# rm -f $name.FE.bdg $name.fc.signal.bedGraph $name.fc.signal.srt.bedGraph