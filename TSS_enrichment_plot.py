##### load required modules
import os
import logging
import pybedtools
import metaseq
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import mlab
import timeit
from datetime import datetime
import argparse

##### set working directory
# os.chdir('/Users/jchap12/Downloads')

##### Define parseArguments to use positional variables in bash submission
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('--bam_file', type=str)
    parser.add_argument('--prefix', type=str)
    parser.add_argument('--tss', type=str)
    parser.add_argument('--bp_edge', type=int)
    parser.add_argument('--chromsizes', type=str)
    parser.add_argument('--read_len', type=int)
    parser.add_argument('--bins', type=int)
    parser.add_argument('--processes', type=int)
    parser.add_argument('--greenleaf_norm', type=bool)

    # Print version
    parser.add_argument('--version', action='version', version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()
    
    return args

##### set temporary function variables
#bam_file="SC3.chr10.bam"
#prefix='test'
#tss='/Users/jchap12/Downloads/hg19_gencode_tss_unique.bed'
#bp_edge=2000
#chromsizes="/Users/jchap12/Downloads/hg19.chrom.sizes"
#read_len=102
#bins=400
#processes=8
#greenleaf_norm=True

def make_tss_plot(bam_file, tss, prefix, chromsizes, read_len, bins=400, bp_edge=2000,
                  processes=8, greenleaf_norm=True):
    #### Set up the log file and timing
    start = timeit.default_timer()

    '''
    Take bootstraps, generate tss plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    logging.info('Generating tss plot...')
    tss_plot_file = '{0}_tss-enrich.png'.format(prefix)
    tss_plot_large_file = '{0}_large_tss-enrich.png'.format(prefix)

    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)

    # Load the bam file
    bam = metaseq.genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
    # Shift to center the read on the cut site
    bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, processes=processes, stranded=True)

    # Actually first build an "ends" file
    #get_ends = '''zcat {0} | awk -F '\t' 'BEGIN {{OFS="\t"}} {{if ($6 == "-") {{$2=$3-1; print}} else {{$3=$2+1; print}} }}' | gzip -c > {1}_ends.bed.gz'''.format(bed_file, prefix)
    #print(get_ends)
    #os.system(get_ends)

    #bed_reads = metaseq.genomic_signal('{0}_ends.bed.gz'.format(prefix), 'bed')
    #bam_array = bed_reads.array(tss_ext, bins=bins,
    #                      processes=processes, stranded=True)

    # Normalization (Greenleaf style): Find the avg height
    # at the end bins and take fold change over that
    if greenleaf_norm:
        # Use enough bins to cover 100 bp on either end
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise
    else:
        bam_array /= bam.mapped_read_count() / 1e6

    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    # Note the middle high point (TSS)
    tss_point_val = max(bam_array.mean(axis=0))

    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(tss_plot_file)

    # Print a more complicated plot with lots of info

    # Find a safe upper percentile - we can't use X if the Xth percentile is 0
    upper_prct = 99
    if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
        upper_prct = 100.0

    plt.rcParams['font.size'] = 8
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5, 10),
                                   vmin=5, vmax=upper_prct, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k', alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))

    # And save the file
    fig.savefig(tss_plot_large_file)

    return tss_plot_file, tss_plot_large_file, tss_point_val

    # finish
    stop = timeit.default_timer()
    print("Run time:", str(datetime.timedelta(seconds=int(stop - start))))
    
    return 'Finished'

if __name__ == '__main__':
    # Parse the arguments
    args = parseArguments()

    # Raw print arguments
    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))

    # Run function
    make_tss_plot(args.bam_file, args.tss, args.prefix, args.chromsizes,
                  args.read_len, args.bins, args.bp_edge, args.processes,
                  args.greenleaf_norm)







