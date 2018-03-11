#!/bin/sh
jobname=`echo $1`
output=`basename $1 .fastq`
cat > /tmp/tempscript.sh << EOF
#!/bin/sh
#$ -N $1
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q rcc-30d
#$ -pe thread 4
#$ -l mem_total=24g
cd `pwd`
bowtie2 -p 4 --local -N 1 --phred33 /db/bowtie2/11192013/hg19 $output.fastq -S $output.sam
samtools view -bS -h -F 4 $output.sam > $output.bam
samtools sort $output.bam $output.sorted
EOF
qsub -q rcc-30d /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh
