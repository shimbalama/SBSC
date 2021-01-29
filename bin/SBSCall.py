#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from multiprocessing import Pool, TimeoutError
import argparse
import sys
import json
import pandas as pd
import os
from glob import glob

#import time
#from subprocess import run#call
import SBSC.parsers.parse_pile as pp
import SBSC.output.filt as ff
#import liams_simple_scaffold_pkg.prefil.prefil_files as fil

def parse_args(args):


    parser = argparse.ArgumentParser(description='somatic var caller')

    required = parser.add_argument_group(
        'Required',
        'ref and pileups')

    required.add_argument(
        '-r',
        '--ref',
        type=str,
        help='ref 37 or 38')

    required.add_argument(
        '-o',
        '--output',
        type=str,
        help='output name [output]',
        default = 'output')

    required.add_argument(
        '-c',
        '--cancer_pile',
        type=str,
        help='pileup of tumour/cancer')

    required.add_argument(
        '-n',
        '--normal_pile',
        type=str,
        help='pileup of normal/blood')

    optional = parser.add_argument_group(
        'Optional',
        'threads and chroms')

    optional.add_argument(
        '-x',
        '--chrom',
        type=str,
        help='Which chromosomes to query. comma,separated,list or [all]',
        default='all')

    optional.add_argument(
        '-t',
        '--threads',
        type=int,
        help='Threads [7]. Uses about 12GB RAM per thread with window of 1MB, assuming 30x norm and 60x tumour',
        default=5)

    optional.add_argument(
        '-b',
        '--base_qual',
        type=int,
        help='Minimum individual base qual to use for var calling. [7]',
        default=7)#10=10% error

    optional.add_argument(
        '-w',
        '--window_size',
        type=int,
        help='To use threading the genome is chunked into chunks of window size. the bigger the better but more RAM needed [100k]',
        default=10000)

    filters = parser.add_argument_group(
        'Filters',
        'P, depth etc')

    filters.add_argument(
        '-q',
        '--min_mean_base_quality',
        type=int,
        help='Minimum mean base quality per position [9]',
        default = 9)

    filters.add_argument(
        '-p',
        '--P_value',
        type=float,
        help='Exclude P values (from Fishers Exact test) larger than [0.001]',
        default = 0.001)

    filters.add_argument(
        '-m',
        '--min_depth_tumour',
        type=int,
        help='Minimun read depth for tumour [12]',
        default = 12)

    filters.add_argument(
        '-y',
        '--min_depth_normal',
        type=int,
        help='Minimun read depth for normal [10]',
        default = 10)

    filters.add_argument(
        '-j',
        '--json',
        action='store_false',
        help='Make a json of all raw results to allow re-run of filtering [True].',
        default = True)

    filters.add_argument(
        '-u',
        '--use_json',
        action='store_true',
        help='Use a json of all raw results to re-run filtering [False].',
        default = False)

    filters.add_argument(
        '-z',
        '--MNVs',
        action='store_false',
        help='Call MNVs (join contigous SNVs in output) [True].',
        default = True)

    # annotate = parser.add_argument_group(
    #     'Annotation',
    #     'annotate')

    # annotate.add_argument(
    #     '-g',
    #     '--gnomad',
    #     type=str,
    #     help='gnomad.genomes.r3.0.sites.vcf.bgz [38]',
    #     default = '/reference/data/gnomAD/gnomad-public/r3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz')

    # annotate.add_argument(
    #     '-k',
    #     '--onek',
    #     type=str,
    #     help='1k genomes')

    # annotate.add_argument(
    #     '-l',
    #     '--low_complexity',
    #     type=str,
    #     help='UCSC repeats [38]',
    #     default = '/working/lab_nicw/liamM/nano/COLO829/pile/repeats_GRCh38.bed')

    # annotate.add_argument(
    #     '-s',
    #     '--seg_dups',
    #     type=str,
    #     help='UCSC segmental duplications [38]',
    #     default = '/working/lab_nicw/liamM/nano/COLO829/pile/segmental_dups_38.bed')

    # annotate.add_argument(
    #     '-m',
    #     '--meres',
    #     type=str,
    #     help='telomeres_and_centromeres_38',
    #     default = '/working/lab_nicw/liamM/nano/COLO829/pile/telomeres_and_centromeres_38.txt')

    return parser.parse_args(args)

def main():

    '''
    Make a simple python package scaffold
    '''

    args = parse_args(sys.argv[1:])
    #pp.main(args)

    ####Several columns contain numeric quality values encoded as individual ASCII characters. Each character can range from “!” to “~” and is decoded by taking its ASCII value and subtracting 33; e.g., “A” encodes the numeric value 32.
    #start_total = time.time()
    if args.chrom == 'all':
            chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')
            
    if args.use_json:
        with open('result.json') as f:
            d = json.load(f)
    else:
        chunks = pp.chunk_ref(args, chroms)#[211:216] #testing
        print (f'len chunks {len(chunks)}')
        d={}
        with Pool(processes=args.threads) as pool:
            tmp = [(args, chunk) for chunk in chunks]
            res = pool.map(pp.doit, tmp)
            for variants in res:
                d = {**d, **variants}
        print ('d',len(d))
        if args.json:
            with open('result.json', 'w') as fout:
                json.dump(d, fout)
        
    #end = time.time()
    #p#rint('run time', end - start_total)
    #start_total = time.time()
    #ff.filt(d, 0.01, 9, 15, args, chroms, 12)
    ff.filt(d, args, chroms)
    #ff.filt(d, 0.0001, 9, 15, args, chroms, 12)

    df = pd.DataFrame.from_dict(d, orient='index')
    df.to_csv(args.output + '.tsv', sep='\t')
    #end = time.time()
    #print('filt time', end - start_total)



if __name__ == "__main__":
    main()


