#!/bin/bash
#for x in `/bin/ls *.fastq` ; do bash bowtie2.sh $x; done
jobname=`echo $1`
output=`basename $1 .fastq`
cat > /tmp/tempscript.sh << EOF
#!/bin/bash
#$ -N $1
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q rcc-30d
#$ -pe thread 4
#$ -l mem_total=24g
cd .
bowtie2 -p 4 --local -N 1 --phred33 /db/bowtie2/11192013/mm9 $output.fastq -S $output.sam
rm $output.fastq
samtools view -bS -h -F 4 $output.sam > $output.bam
samtools sort $output.bam $output.sorted
rm $output.sam
rm $output.bam
EOF
qsub -q rcc-30d /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh