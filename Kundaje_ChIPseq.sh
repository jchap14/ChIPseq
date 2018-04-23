#!/bin/bash

##### Submission script for BDS controlled ChIPseq pipeline
## for x in `/bin/ls *.trim.R1.fq.gz` ; do bash Kundaje_ChIPseq.sh $x; done

##### specify variables to pass to chipseq.py
FASTQR1=$1
# NAME=`basename $FASTQR1 .untrimmed.R1.fq.gz`
# FASTQR2=`echo $NAME.untrimmed.R2.fq.gz`
NAME=`basename $FASTQR1 .trim.R1.fq.gz`
FASTQR2=`echo $NAME.trim.R2.fq.gz`
CHIPTYPE="histone"
# CHIPTYPE="TF"


## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash

## conda environment
## don't source conda envs, they are auto-activated by the bds pipeline

##### run script
## -type can be histone or TF

python /srv/gsfs0/projects/snyder/chappell/TF_chipseq_pipeline/chipseq.py \
-type histone --screen $NAME -pe -species hg19 -nth 12 \
-fastq1_1 $FASTQR1 -fastq1_2 $FASTQR2 -out_dir $NAME

## deactivate conda environment
EOF

## qsub then remove the tempscript
bash $NAME.tempscript.sh
sleep 1
rm $NAME.tempscript.sh
