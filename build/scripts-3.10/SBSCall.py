#!python
# -*- coding: utf-8 -*-

import argparse
import logging
import sys
import time
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from SBSC import __version__
from SBSC.output.filt import filt
from SBSC.parsers.parse_pileup import Pileups, process_genome_data
from SBSC.parsers.parse_reference import chunk_ref

logging.basicConfig(filename="log.txt", encoding="utf-8", level=logging.DEBUG)

CHROMOSOMES = ["chr" + str(i + 1) for i in range(22)] + ["chrX", "chrY"]


def parse_args(args):

    parser = argparse.ArgumentParser(
        description="somatic variant caller",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage='SBSCall --chrom chr17 --processes 4 call --cancer_pile HCC1937_t.pileup.gz --normal_pile HCC1937.pileup.gz --ref GRCh38.fa \n \
OR: SBSCall --chrom chr17 --processes 4 filt --raw_results results.pickle'
    )

    parser.add_argument("--version", action="version", version=__version__)

    parser.add_argument(
        "--output", type=lambda p: Path(p).absolute(), help="output VCF", default="output"
    )
    parser.add_argument(
        "--chrom",
        type=str,
        help="Which chromosome to query",
        default="standard",
        choices=CHROMOSOMES + ["all"],
    )  #        nargs='+',

    parser.add_argument(
        "--MNVs",
        dest="MNVs",
        action="store_false",
        default=True,
        help="Call MNVs (join contigous SNVs in output).",
    )

    parser.add_argument(
        "--processes", type=int, help=("Number of processes"), default=1
    )

    filters = parser.add_argument_group("variant filtering options", "P, depth etc")

    filters.add_argument(
        "--P_value",
        type=float,
        help="Exclude P values (from Fishers Exact test) larger than this.",
        default=0.001,
    )

    filters.add_argument(
        "--min_depth_tumour",
        type=int,
        help="Minimum read depth for tumour",
        default=12,
    )

    filters.add_argument(
        "--min_depth_normal",
        type=int,
        help="Minimum read depth for normal",
        default=10,
    )

    subparsers = parser.add_subparsers(
        required=True,
        dest="subparser_name",
        help=f"{parser.prog} [call|filt] --help for subcommand help",
        title="subcommands",
    )

    call = subparsers.add_parser(
        "call",
        help="call and filter variants",
        description="Two output files are generated: a JSON format file "
        "containing all variants and a tsv format file containing only "
        "variants that pass filtering.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=''
    )

    call.add_argument(
        "--ref",
        required=True,
        type=lambda p: Path(p).absolute(),
        help="reference fasta",
    )

    call.add_argument(
        "--cancer_pile",
        required=True,
        type=lambda p: Path(p).absolute(),
        help="pileup of tumour/cancer",
    )

    call.add_argument(
        "--normal_pile",
        required=True,
        type=lambda p: Path(p).absolute(),
        help="pileup of normal/blood",
    )

    call.add_argument(
        "--min_base_qual",
        type=int,
        help="Minimum individual base qual to use for var calling.",
        default=10,
    )  # 10=10% error

    call.add_argument(
        "--window_size",
        type=int,
        help=("To use multi processing the genome is chunked."),
        default=10_000,
    )

    filt = subparsers.add_parser(
        "filt",
        help="filter existing variants",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=''
    )

    filt.add_argument("--raw_results", help="JSON file containing unfiltered calls.")

    return parser.parse_args(args)


def main():
    """
    1. Parse pileup
    2. Call SNVs
    3. Filter SNVs
    """

    args = parse_args(sys.argv[1:])
    start = time.time()
    

    if args.chrom not in CHROMOSOMES + ['all']:
        raise ValueError(f'Please select chromosome from {CHROMOSOMES}')
    if args.chrom == 'all':
        chroms = CHROMOSOMES
    else:
        chroms = args.chrom

    if args.subparser_name == "filt":
        pd.read_pickle("results.pickle")
    else:
        chunks = chunk_ref(args.window_size, args.ref, chroms)
        dfs = []
        with Pool(processes=args.processes) as pool:
            process_genome_data_prefil = partial(
                process_genome_data,
                Pileups(
                    args.cancer_pile,
                    args.normal_pile,
                ),
                args.min_base_qual,
            )
            res = pool.imap_unordered(process_genome_data_prefil, chunks, chunksize=5)
            for df in res:
                dfs.append(df)
        df = pd.concat(dfs)
        df.to_pickle("results.pickle")

    df.to_csv("~/Downloads/regex.csv")
    # filt(df, args)
    # df = pd.DataFrame.from_dict(d, orient='index')
    print(f"Total run time (seconds): {time.time() - start}")


if __name__ == "__main__":
    main()
