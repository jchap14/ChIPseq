#!/bin/sh
jobname=`echo $1`
output=`basename $1 .fastq.gz`
cat > /tmp/tempscript.sh << EOF
#!/bin/sh
#$ -N $1
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q rcc-30d
#$ -pe thread 4
cd `pwd`
gunzip $1
EOF
qsub -q rcc-30d /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh