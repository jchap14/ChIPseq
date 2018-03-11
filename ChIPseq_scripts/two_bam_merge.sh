#!/bin/sh
#sub command - for x in `/bin/ls *_2.bam` ; do bash two_bam_merge.sh $x; done
jobname=`echo $1`
output=`basename $1 _2.bam`
cat > /tmp/tempscript.sh << EOF
#!/bin/bash
#$ -N $1
#$ -j y
cd .
samtools merge $output.merged.bam $output'_1.bam' $output'_2.bam'  
EOF
qsub -q rcc-30d /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh
rm *.o*
