#!/usr/bin/env python3
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
        usage="SBSFilter --chrom chr17 --processes 4 filt --raw_results results.pickle",
    )

    parser.add_argument("--version", action="version", version=__version__)

    parser.add_argument(
        "--output",
        type=lambda p: Path(p).absolute(),
        help="output VCF",
        default="output",
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
        help="Exclude P values (Fishers Exact test) larger than this.",
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

    filters.add_argument(
        "--raw_results",
        type=lambda p: Path(p).absolute(),
        help="Pickle file containing unfiltered calls.",
    )

    return parser.parse_args(args)


def main():
    """
    Filter SNVs
    """

    args = parse_args(sys.argv[1:])
    start = time.time()

    if args.chrom not in CHROMOSOMES + ["all"]:
        raise ValueError(f"Please select chromosome from {CHROMOSOMES}")
    if args.chrom == "all":
        chroms = CHROMOSOMES
    else:
        chroms = args.chrom

    
    # filt(df, args)
    # df = pd.DataFrame.from_dict(d, orient='index')
    print(f"Total run time (seconds): {time.time() - start}")


if __name__ == "__main__":
    main()
