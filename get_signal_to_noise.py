##### load required modules
import os
import logging
import subprocess
import signal
import pandas as pd
import argparse

## not sure if these necessary
# import matplotlib
# matplotlib.use('Agg')
# import pybedtools
# import metaseq
# from matplotlib import pyplot as plt
# import numpy as np
# from matplotlib import mlab
# import timeit
# from datetime import datetime

##### function to run shell commands in python (Taken from ENCODE DCC ATAC pipeline)
def run_shell_cmd(cmd):
    print(cmd)
    try:
        p = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid)
        pid = p.pid
        pgid = os.getpgid(pid)
        ret = ''
        while True:
            line = p.stdout.readline()
            if line=='' and p.poll() is not None:
                break
            # log.debug('PID={}: {}'.format(pid,line.strip('\n')))
            #print('PID={}: {}'.format(pid,line.strip('\n')))
            ret += line
        p.communicate() # wait here
        if p.returncode > 0:
            raise subprocess.CalledProcessError(
                p.returncode, cmd)
        return ret.strip('\n')
    except:
        # kill all child processes
        os.killpg(pgid, signal.SIGKILL)
        p.terminate()
        raise Exception('Unknown exception caught. PID={}'.format(pid))

##### function to get the number of lines from a file (Taken from ENCODE DCC ATAC pipeline)
def get_num_lines(f):
    cmd = 'cat {} | wc -l'.format(f)
    return int(run_shell_cmd(cmd))

##### Define parseArguments to use positional variables in bash submission
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('--outname', help='experiment name')
    parser.add_argument('--final_bed', help='Final filtered alignments in BED format')
    parser.add_argument('--dnase_regions', help='Open chromatin region file')
    parser.add_argument('--blacklist_regions', help='Blacklisted region file')
    parser.add_argument('--prom_regions', help='Promoter region file')
    parser.add_argument('--enh_regions', help='Enhancer region file')
    parser.add_argument('--peaks', help='Peak file')
    parser.add_argument('--jz_k27ac_regions', help='Judith Zaugg K27ac peaks')
    parser.add_argument('--jz_k4me1_regions', help='Judith Zaugg K4me1 peaks')
    parser.add_argument('--huvec_k27me3_regions', help='Merged HUVEC K27me3 peaks')

    # Print version
    parser.add_argument('--version', action='version', version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()
    
    return args

##### Function: take bed file of reads & bed file of regions & gets % of reads sitting in said regions
def get_fract_reads_in_regions(reads_bed, regions_bed):
    # uses new run_shell_cmd
    cmd = "bedtools sort -i {}  | "
    cmd += "bedtools merge -i stdin | "
    cmd += "bedtools intersect -u -nonamecheck -a {} -b stdin | "
    cmd += "wc -l"
    #cmd += "bedtools intersect -c -nonamecheck -a stdin -b {} | "
    #cmd += "awk '{{ sum+=$4 }} END {{ print sum }}'"
    cmd = cmd.format(regions_bed, reads_bed)
    intersect_read_count = int(run_shell_cmd(cmd))
    total_read_count = get_num_lines(reads_bed)
    fract_reads = float(intersect_read_count) / total_read_count
    return intersect_read_count, fract_reads

##### set temporary variables to test functions
# outname='GFP_K27me3.rep3'
# final_bed='GFP_K27me3.rep3.bam.bed'
# dnase_regions="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
# blacklist_regions="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/hg19EncodeMapabilityBlacklist.bed"
# prom_regions="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
# enh_regions="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
# peaks='GFP_K27me3.rep3.FDRe1.Peaks.bed'

##### Function: run get_fract_reads_in_regions on reginon
def get_signal_to_noise(final_bed, dnase_regions, blacklist_regions, 
                        prom_regions,enh_regions, peaks, jz_k27ac_regions,
                        jz_k4me1_regions, huvec_k27me3_regions):
    '''
    Given region sets, determine whether reads are
    falling in or outside these regions
    '''
    logging.info('signal to noise...')

    # Dnase regions
    reads_dnase, fract_dnase = get_fract_reads_in_regions(final_bed, dnase_regions)

    # Blacklist regions
    reads_blacklist, fract_blacklist = get_fract_reads_in_regions(final_bed, blacklist_regions)

    # Prom regions
    reads_prom, fract_prom = get_fract_reads_in_regions(final_bed, prom_regions)

    # Enh regions
    reads_enh, fract_enh = get_fract_reads_in_regions(final_bed, enh_regions)

    # Peak regions
    reads_peaks, fract_peaks = get_fract_reads_in_regions(final_bed, peaks)
    
    # JZ K27ac regions
    reads_jz_k27ac, fract_jz_k27ac = get_fract_reads_in_regions(final_bed, jz_k27ac_regions)

    # JZ K4me1 regions
    reads_jz_k4me1, fract_jz_k4me1 = get_fract_reads_in_regions(final_bed, jz_k4me1_regions)

    # HUVEC K27me3 regions
    reads_huvec_k27me3, fract_huvec_k27me3 = get_fract_reads_in_regions(final_bed, huvec_k27me3_regions)


    return reads_dnase, fract_dnase, reads_blacklist, fract_blacklist, \
        reads_prom, fract_prom, reads_enh, fract_enh, \
        reads_peaks, fract_peaks, reads_jz_k27ac, fract_jz_k27ac, \
        reads_jz_k4me1, fract_jz_k4me1, reads_huvec_k27me3, fract_huvec_k27me3


##### This is run in the main function
   
if __name__ == '__main__':
    ##### Parse the arguments
    args = parseArguments()

    ##### Raw print arguments
    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))

    ##### Run functions
    
    ##### determine number of reads in each of the following
#    ## DNase hypersensitivity regions
#    reads_dnase, fract_dnase = get_signal_to_noise(args.final_bed, args.dnase_regions)
#    ## blacklist regions    
#    reads_blacklist, fract_blacklist = get_signal_to_noise(args.final_bed, args.blacklist_regions)
#    ## promoter regions
#    reads_prom, fract_prom = get_signal_to_noise(args.final_bed, args.prom_regions)
#    ## enhancer regions
#    reads_enh, fract_enh = get_signal_to_noise(args.final_bed, args.enh_regions)
#    ## peaks from this experiment
#    reads_peaks, fract_peaks = get_signal_to_noise(args.final_bed, args.peaks)
    
    ## combines all 5 function calls above
    reads_dnase, fract_dnase, reads_blacklist, fract_blacklist, \
    reads_prom, fract_prom, reads_enh, fract_enh, \
    reads_peaks, fract_peaks, reads_jz_k27ac, fract_jz_k27ac, \
    reads_jz_k4me1, fract_jz_k4me1, reads_huvec_k27me3, fract_huvec_k27me3 \
    = get_signal_to_noise(args.final_bed, args.dnase_regions, args.blacklist_regions, \
    args.prom_regions, args.enh_regions, args.peaks, args.jz_k27ac_regions, \
    args.jz_k4me1_regions, args.huvec_k27me3_regions)
    
    ##### get total_read_count
    total_read_count = get_num_lines(args.final_bed)

    ##### print results to output file (or create specific output file)
    df = pd.DataFrame([[str(args.outname),total_read_count, reads_dnase, fract_dnase,
                        reads_blacklist, fract_blacklist, reads_prom, fract_prom,
                        reads_enh, fract_enh, reads_peaks,fract_peaks, reads_jz_k27ac,
                        fract_jz_k27ac, reads_jz_k4me1, fract_jz_k4me1,
                        reads_huvec_k27me3, fract_huvec_k27me3]],
        columns=['Condition','TotalReads','ReadsInDNAse','PercInDNAse',
                 'ReadsInBlacklist', 'PercInBlacklist',
                 'ReadsInPromoters', 'PercInPromoters',
                 'ReadsInEnhancers', 'PercInEnhancers',
                 'ReadsInPeaks', 'PercInPeaks',
                 'ReadsInJZk27ac', 'PercInJZk27ac',
                 'ReadsInJZk4me1', 'PercInJZk4me1',
                 'ReadsInHUVEC27me3', 'PercInHUVEC27me3'])
    
    df.to_csv(path_or_buf=str(args.outname + ".signalToNoise.csv"), index=False)
