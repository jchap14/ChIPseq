#!/bin/bash
##### Submission script for BDS controlled ChIPseq pipeline

##### forloop batch submission command (don't use if input controls are different)
## for x in `/bin/ls *.repA.trim.R1.fq.gz` ; do bash Kundaje_ChIPseq_w_2reps_w_input.sh $x; done

##### INPUTs required
## the input control file name
FQ_INPUT_CNTL_PRE="GFP.Input.repAB" #caMEK5.Input.repAB
## type of factor for ChIP (-type in python script)
CHIPTYPE="histone"
# CHIPTYPE="TF"

##### input files to pass to chipseq.py (FQs should be trimmed)
## calc replicate 1 FQs (format == $NAME.repB)
FQ_REPA_R1=$1
NAME=`basename $FQ_REPA_R1 .repA.trim.R1.fq.gz`
FQ_REPA_R2=`echo $NAME.repA.trim.R2.fq.gz`
## calc replicate 2 FQs (format == $NAME.repB)
FQ_REPB_R1=`echo $NAME.repB.trim.R1.fq.gz`
FQ_REPB_R2=`echo $NAME.repB.trim.R2.fq.gz`
## calc input control FQs
FQ_INPUT_CNTL_R1=`echo $FQ_INPUT_CNTL_PRE.trim.R1.fq.gz`
FQ_INPUT_CNTL_R2=`echo $FQ_INPUT_CNTL_PRE.trim.R2.fq.gz`

##### create the tempscript for queue submission
cat > $NAME.tempscript.sh << EOF
#!/bin/bash

## conda environment
## don't source conda envs, they are auto-activated by the bds pipeline

##### run chipseq.py ('-peak_caller macs2' left out w inputs included)
python /srv/gsfs0/projects/snyder/chappell/TF_chipseq_pipeline/chipseq.py \
-type $CHIPTYPE --screen $NAME -pe -species hg19 -nth 12 \
-fastq1_1 $FQ_REPA_R1 -fastq1_2 $FQ_REPA_R2 \
-fastq2_1 $FQ_REPB_R1 -fastq2_2 $FQ_REPB_R2 \
-ctl_fastq1_1 $FQ_INPUT_CNTL_R1 -ctl_fastq1_2 $FQ_INPUT_CNTL_R2 \
-out_dir $NAME -mem_dedup 15G -q mpsnyder

## end tempscript
EOF

## bash the script (will submit BDS job), then remove the tempscript
bash $NAME.tempscript.sh
sleep 1
rm $NAME.tempscript.sh
