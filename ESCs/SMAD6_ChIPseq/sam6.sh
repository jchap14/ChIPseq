#!/bin/sh
#$ -N $1
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q rcc-30d
#$ -pe thread 4
#$ -l mem_total=24g
cd `pwd`
samtools merge IgG_DE.bam C1AHRACXX_s8_0_GSLv2-7_04_SL27363.sorted.bam C1AHRACXX_s7_0_GSLv2-7_04_SL27363.sorted.bam D1FT7ACXX_s6_0_GSLv2-7_04_SL27363.sorted.bam
cd `pwd`
