#!/bin/bash

##### Submission script for BDS controlled ChIPseq pipeline

##### for x in `/bin/ls *.trim.R1.fq.gz` ; do bash Kundaje_ChIPseq.sh $x; done

source activate aquas_chipseq_py3

##### specify variables to pass to chipseq.py
FASTQ1=$1
NAME=`basename $FASTQ1 .trim.R1.fq.gz`
FASTQ2==`echo $NAME.trim.R2.fq.gz`

## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash

## add modules & source specific conda environment
source activate aquas_chipseq_py3

##### run script
## -type can be histone or TF

python /srv/gsfs0/projects/snyder/chappell/TF_chipseq_pipeline/chipseq.py \
-type histone --screen $NAME -pe -species hg19 -nth 12 \
-fastq1 $FASTQ1 -fastq2 $FASTQ2

## deactivate conda environment
EOF

## qsub then remove the tempscript
bash $NAME.tempscript.sh
sleep 1
# rm $NAME.tempscript.sh
