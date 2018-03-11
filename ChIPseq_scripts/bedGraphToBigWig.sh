#!/bin/bash
# for x in `/bin/ls *p.bdg` ; do bash bedGraphToBigWig.sh $x & done
bedGraphToBigWig $1 /home/jchap14/hi_quota_folder/Annotations/UCSC-hg19/hg19.genomeSizes_subset $1.bw