##### Script to annotate (ATAC/ChIP) peaks with nearest gene using HOMER ####
## instructions provided here: http://homer.ucsd.edu/homer/ngs/annotation.html

##### Submit this script locally
## for x in `find . -name "*.bed"`; do bash Peak_Annotation.sh $x; done

## define variables
BEDFILE=$1
NAME=`basename $BEDFILE .sorted.thr1.filt.Peaks.bed`

## peaks: ~/SAMPLE/peak/idr/conservative_set/SAMPLE_rep1-pr.IDR0.1.filt.narrowPeak.gz
## decompress peaks
#gunzip $BEDFILE

## add homer module
#module add homer

##### run annotatePeaks.pl with built-in hg19 RefSeq
annotatePeaks.pl $BEDFILE hg19 > $NAME.annoPeaks.txt