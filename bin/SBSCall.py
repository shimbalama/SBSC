#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from multiprocessing import Pool
import argparse
import sys
import json
import pandas as pd

import SBSC.parsers.parse_pile as pp
import SBSC.output.filt as ff
from SBSC import __version__


def parse_args(args):

    parser = argparse.ArgumentParser(
        description='somatic variant caller',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(
        'Required',
        'ref and pileups')

    required.add_argument(
        '-r',
        '--ref',
        required=True,
        type=str,
        help='reference fasta')

    required.add_argument(
        '-c',
        '--cancer_pile',
        required=True,
        type=str,
        help='pileup of tumour/cancer')

    required.add_argument(
        '-n',
        '--normal_pile',
        required=True,
        type=str,
        help='pileup of normal/blood')

    optional = parser.add_argument_group(
        'Optional',
        'output, threads and chroms')

    optional.add_argument(
        '-o',
        '--output',
        type=str,
        help='output name',
        default='output')

    optional.add_argument(
        '-x',
        '--chrom',
        type=str,
        help='Which chromosomes to query. comma,separated,list or "all"',
        default='all')

    optional.add_argument(
        '-t',
        '--threads',
        type=int,
        help=('Number of threads. Uses about 6GB RAM per thread with window '
              'of 100k, assuming 30x norm and 60x tumour'),
        default=5)

    optional.add_argument(
        '-b',
        '--base_qual',
        type=int,
        help='Minimum individual base qual to use for var calling.',
        default=7)  # 10=10% error

    optional.add_argument(
        '-w',
        '--window_size',
        type=int,
        help=('To use threading the genome is chunked into chunks of '
              'WINDOW_SIZE. The bigger the better but more RAM needed.'),
        default=100000)

    filters = parser.add_argument_group(
        'Filters',
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
        '-y',
        '--min_depth_normal',
        type=int,
        help='Minimum read depth for normal',
        default=10)

    filters.add_argument(
        '-j',
        '--json',
        action='store_false',
        help='Make a json of all raw results to allow re-run of filtering.',
        default=True)

    filters.add_argument(
        '-u',
        '--use_json',
        action='store_true',
        help='Use a json of all raw results to re-run filtering.',
        default=False)

    filters.add_argument(
        '-z',
        '--MNVs',
        action='store_false',
        help='Call MNVs (join contigous SNVs in output).',
        default=True)

    parser.add_argument(
        '--version',
        action='version',
        version=__version__)

    return parser.parse_args(args)


def main():
    '''
    Make a simple python package scaffold
    '''

    args = parse_args(sys.argv[1:])

    if args.chrom == 'all':
        chroms = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
    else:
        chroms = args.chrom.split(',')

    if args.use_json:
        with open('result.json') as f:
            d = json.load(f)
    else:
        chunks = pp.chunk_ref(args, chroms)
        d = {}
        with Pool(processes=args.threads) as pool:
            tmp = [(args, chunk) for chunk in chunks]
            res = pool.map(pp.doit, tmp)
            for variants in res:
                d.update(variants)
        if args.json:
            with open('result.json', 'w') as fout:
                json.dump(d, fout)

    ff.filt(d, args, chroms)
    df = pd.DataFrame.from_dict(d, orient='index')
    df.to_csv(args.output + '.tsv', sep='\t')


if __name__ == "__main__":
    main()
