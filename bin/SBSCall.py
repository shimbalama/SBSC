#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from functools import partial
from multiprocessing import Pool
import argparse
import sys
import json
import time
import pandas as pd
from Bio import SeqIO
import SBSC.parsers.parse_pile as pp
from SBSC.parsers.data_processor import process_genome_data
import SBSC.output.filt as ff
from SBSC import __version__


def parse_args(args):

    parser = argparse.ArgumentParser(
        description='somatic variant caller',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--version',
        action='version',
        version=__version__)

    parser.add_argument(
        '-o',
        '--output',
        type=str,
        help='output files base name',
        default='output')

    filters = parser.add_argument_group(
        'variant filtering options',
        'P, depth etc')

    filters.add_argument(
        '-q',
        '--min_mean_base_quality',
        type=int,
        help='Minimum mean base quality per position',
        default=9)

    filters.add_argument(
        '-p',
        '--P_value',
        type=float,
        help='Exclude P values (from Fishers Exact test) larger than this.',
        default=0.001)

    filters.add_argument(
        '-m',
        '--min_depth_tumour',
        type=int,
        help='Minimum read depth for tumour',
        default=12)

    filters.add_argument(
        '-x',
        '--chrom',
        type=str,
        help='Which chromosomes to query. comma,separated,list or "all"',
        default='all')

    filters.add_argument(
        '-y',
        '--min_depth_normal',
        type=int,
        help='Minimum read depth for normal',
        default=10)

    mnvgroup = filters.add_mutually_exclusive_group()

    mnvgroup.add_argument(
        '--MNVs',
        dest='MNVs',
        action='store_true',
        default=True,
        help='Call MNVs (join contigous SNVs in output).')
    
    mnvgroup.add_argument(
        '--no-MNVs',
        dest='MNVs',
        action='store_false')

    subparsers = parser.add_subparsers(
        required=True,
        dest='subparser_name',
        help=f'{parser.prog} [call|filt] --help for subcommand help',
        title='subcommands')

    call = subparsers.add_parser(
        'call',
        help='call and filter variants',
        description='Two output files are generated: a JSON format file '
        'containing all variants and a tsv format file containing only '
        'variants that pass filtering.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    call.add_argument(
        '-r',
        '--ref',
        required=True,
        type=str,
        help='reference fasta')

    call.add_argument(
        '-c',
        '--cancer_pile',
        required=True,
        type=str,
        help='pileup of tumour/cancer')

    call.add_argument(
        '-n',
        '--normal_pile',
        required=True,
        type=str,
        help='pileup of normal/blood')

    call.add_argument(
        '-b',
        '--min_base_qual',
        type=int,
        help='Minimum individual base qual to use for var calling.',
        default=10)  # 10=10% error

    call.add_argument(
        '-t',
        '--threads',
        type=int,
        help=('Number of threads. Uses about 6GB RAM per thread with window '
              'of 100k, assuming 30x norm and 60x tumour'),
        default=5)

    call.add_argument(
        '-w',
        '--window_size',
        type=int,
        help=('To use threading the genome is chunked into chunks of '
              'WINDOW_SIZE. The bigger the better but more RAM needed.'),
        default=10_000)

    filt = subparsers.add_parser(
        'filt',
        help='filter existing variants',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    filt.add_argument(
        'raw_results',
        help='JSON file containing unfiltered calls.')

    return parser.parse_args(args)

def chunk_ref(args, chroms):
    '''
    Split ref into chunks for parallel processing
    '''
    size = 100_000 #args.window_size - dont go lower than 100k, or super slow
    total_len = 0
    for record in SeqIO.parse(args.ref, 'fasta'):
        if record.id in chroms:
            seqlen = len(record.seq)
            for start in range(0, seqlen, size):
                if seqlen > start + size:
                    end = start + size
                else:  # end of chrom
                    end = seqlen
                seq_chunk = str(record.seq)[start:end]
                total_len += len(seq_chunk)
                chars = set(seq_chunk)
                if len(chars) == 1 and chars.pop().upper() == 'N':
                    continue
                else:
                    yield pp.GenomicRegion(
                        str(record.id).replace('chr', ''),
                        start,
                        end,
                        seq_chunk
                    )
            

def main():
    '''
    1. Parse pileup
    2. Call SNVs
    3. Filter SNVs
    '''

    args = parse_args(sys.argv[1:])
    start = time.time()
    if args.chrom == 'all':
        chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')

    if args.subparser_name == 'filt':
        with open(args.raw_results) as f:
            d = json.load(f)
    else:
        chunks = chunk_ref(args, chroms)
        dfs = []
        with Pool(processes=args.threads) as pool:
            process_genome_data_prefil = partial(process_genome_data, args)
            #tmp = [(args, chunk) for chunk in chunks]
            res = pool.imap_unordered(process_genome_data_prefil, chunks, chunksize=5)
            for df in res:
                dfs.append(df)
        df = pd.concat(dfs)
        # with open(f'{args.output}.json', 'w') as fout:
        #     json.dump(d, fout)
        print(df)
    df.to_csv('~/Downloads/regex.csv')
    #ff.filt(d, args)
    #df = pd.DataFrame.from_dict(d, orient='index')
    print(f'Total run time (seconds): {time.time() - start}')

if __name__ == "__main__":
    main()
