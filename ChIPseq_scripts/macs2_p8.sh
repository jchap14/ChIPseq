#!/bin/sh
#sub command - for x in `/bin/ls *.merged.bam` ; do bash macs2_p2.sh $x; done
jobname=`echo $1`
output=`basename $1 .merged.bam`
cat > /tmp/tempscript.sh << EOF
#!/bin/bash
#$ -N $1
#$ -j y
cd .
time python2.7 /usr/local/macs2/2.0.10.09132012/bin/macs2 callpeak \
-t  $output.merged.bam -c WCE_wt.bam -f BAM -g mm --keep-dup 1 \
-n $output.pMinus2 -B --nomodel --shiftsize 200 -p 0.00000001 --broad --broad-cutoff 0.0000001
EOF
qsub -q rcc-30d /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh
