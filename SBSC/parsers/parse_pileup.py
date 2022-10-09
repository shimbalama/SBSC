#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import re
from collections import Counter
from dataclasses import dataclass
from functools import partial, reduce
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import pysam

from .df_schema import Schema, remove_ref_from_tumour_res
from .parse_reference import GenomicRegion
from .variant_calling import CallAllVaraints


@dataclass
class Pileups:
    cancer_pileup: Path
    normal_pileup: Path


def create_df(region: GenomicRegion, pileup: Path) -> pd.DataFrame:
    """Converts a genomic region of a pileup to a df"""
    keys = [name for schema, name in vars(Schema).items() if "__" not in schema][:10]
    tabixfile = pysam.TabixFile(str(pileup))
    coords: list[Any] = ["chr" + str(region.chrom), region.start, region.end + 1]
    df = pd.DataFrame(
        [line.split("\t") for line in tabixfile.fetch(*coords)], columns=keys
    )
    index_names = Schema.CHROMOSOME + ":" + Schema.POSITION
    df[index_names] = df[Schema.CHROMOSOME] + ":" + df[Schema.POSITION]
    df.set_index(index_names, inplace=True)
    df = df.astype({"reads": "int64", "pos": "int64"})
    return df.copy(deep=True)


def fill_df(input_data: Pileups, region: GenomicRegion) -> pd.DataFrame:
    """"""
    logging.info(f"Creating a df for {str(region)}")
    df_tumour = create_df(region, input_data.cancer_pileup)
    df_normal = create_df(region, input_data.normal_pileup)
    df = df_normal.join(df_tumour, how="inner", lsuffix="_normal", rsuffix="_tumour")
    for col in [Schema.REFERENCE, Schema.CHROMOSOME, Schema.POSITION]:
        df = remove_redundant_column(df, col)
    check_the_ref_seq_matches_the_pileup_seq(df, region)
    df = remove_positions_with_little_support_for_somatic_var(df)
    df = add_homoploymers(df, region.homopolymer_lengths)

    return df.copy(deep=True)


def remove_positions_with_little_support_for_somatic_var(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Filters out positions with 4 or less putative somatic SNVs"""

    def putative_vars(res_t: str, res_n: str, ref: str) -> bool:
        res_t = remove_ref_from_tumour_res(res_n, res_t, ref)
        if res_t:
            return True if Counter(res_t).most_common()[0][1] > 4 else False
        return False

    return df.copy(deep=True)[
        df.apply(
            lambda x: putative_vars(
                x[f"{Schema.RESULTS}_tumour"],
                x[f"{Schema.RESULTS}_normal"],
                x[Schema.REFERENCE],
            ),
            axis=1,
        )
    ]


def remove_redundant_column(df: pd.DataFrame, col_name: str) -> pd.DataFrame:
    """Collapses columns that are identical in tumour and normal"""
    tumour_col = f"{col_name}_tumour"
    normal_col = f"{col_name}_normal"
    if not df[tumour_col].equals(df[normal_col]):
        raise ValueError(f"{col_name} differs!")
    df[col_name] = df[normal_col]
    del df[tumour_col]
    del df[normal_col]

    return df.copy(deep=True)


def check_the_ref_seq_matches_the_pileup_seq(
    df: pd.DataFrame, region: GenomicRegion
) -> pd.DataFrame:
    """Ensure internal consistency"""
    if len(df) == len(region.seq):
        df["seq_from_ref"] = list(region.seq)
        if not df["seq_from_ref"].equals(df["ref"]):
            raise ValueError(
                "The given reference sequence doesnt match the reference inteh pileup"
            )
        del df["seq_from_ref"]
    else:
        logging.info(
            "skipping ref check for %s due to length missmatch: \
            df=%s and len seq=%s",
            str(region),
            len(df),
            len(region.seq),
        )


def add_homoploymers(df: pd.DataFrame, region: dict[int, int]) -> pd.DataFrame:
    """Adds a homopolymer column to df"""
    df[Schema.HOMOPOLYMER] = df[Schema.POSITION].apply(lambda x: region.get(x, 0))

    return df.copy(deep=True)


def create_SV_column(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Adds SV column to df"""
    results_col = f"{Schema.RESULTS}_{sample_type}"
    SV_col = f"{Schema.STRUCTURAL_VARIANTS}_{sample_type}"
    df[SV_col] = df[results_col].str.count("^") + df[results_col].str.count("$")
    return df.copy(deep=True)


def remove_carrots_from_res(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """
    Removes caret from edn of results str to allow SNV calling
    ^ (caret) marks the start of a read segment and the ASCII
    of the character following `^' minus 33 gives the mapping quality
    """

    def slice_and_dice(res: str) -> str:
        if "^" in res:
            if res.startswith("^"):
                res = "".join([bit[1:] for bit in res.split("^")])
            else:
                res = "".join(
                    [bit[1:] if j != 0 else bit for j, bit in enumerate(res.split("^"))]
                )
        return res

    df[f"{Schema.RESULTS_NO_CARET}_{sample_type}"] = df[
        f"{Schema.RESULTS}_{sample_type}"
    ].apply(slice_and_dice)
    return df.copy(deep=True)



def correct_regex(regex_result: list[str]) -> list[str]:
    """Helper function for removing indels from str"""
    corected_regex_result: list[str] = []
    # can't find other way to stop regex adding unrelated nucs at end
    for indel in regex_result:
        indel_len = "".join(char for char in indel if char.isdigit())
        indel_nucs = "".join(char for char in indel if char.isalpha())[: int(indel_len)]
        corected_regex_result.append("".join([indel[0], indel_len, indel_nucs]))
    return corected_regex_result


def remove_indels_from_results_str(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Separates the indels out from the pre-processed results str"""

    def find_and_remove_indels(res: str) -> str:
        regex_result: list[str] = re.findall(
            "\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+", res
        )
        corected_regex_result = correct_regex(regex_result)
        indel_count = Counter(corected_regex_result)
        for indel in indel_count.keys():
            res = "".join(res.split(indel))
        return "".join(char for char in res if char in ".,ATCGatcg*")

    nucs = f"{Schema.RESULTS_NUCLEOTIDES}_{sample_type}"
    qual = f"{Schema.QUALITY}_{sample_type}"
    df[nucs] = df[f"{Schema.RESULTS_NO_CARET}_{sample_type}"].apply(
        find_and_remove_indels
    )
    assert df[nucs].str.len().equals(df[qual].str.len())
    f"{nucs} and {qual} should have same length"
    return df.copy(deep=True)


def find_indels_in_results_str(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Separates the indels out from the pre-processed results str"""

    def find_indels_in_res(res: str) -> tuple[str, int] | float:
        regex_result = re.findall("\+[0-9]+[ACGTN]+|\-[0-9]+[ACGTN]+", res.upper())
        corected_regex_result = correct_regex(regex_result)
        indel_count = Counter(corected_regex_result)
        return indel_count.most_common()[0] if indel_count else np.NaN

    df[f"{Schema.RESULTS_INDELS}_{sample_type}"] = df[
        f"{Schema.RESULTS_NO_CARET}_{sample_type}"
    ].apply(find_indels_in_res)
    return df.copy(deep=True)


def remove_low_quality_bases(
    sample_type: str, df: pd.DataFrame, phred: int
) -> pd.DataFrame:
    """Filters out any bases below given phred
    return filtered bases and coresponding phred"""

    def filter_on_base_qual(res: str, qual: int) -> pd.Series:
        new_res, new_qual = "", ""
        for base, quality in zip(res, qual):
            if (ord(quality) - 33) > phred:
                new_res += base
                new_qual += quality
        assert len(new_res) == len(new_qual)
        "resuts should match qualities"
        return pd.Series([new_res, new_qual])

    df[
        [
            f"{Schema.RESULTS_NUCLEOTIDES_FILTERED}_{sample_type}",
            f"{Schema.QUALITY_FILTERED}_{sample_type}",
        ]
    ] = df.apply(
        lambda x: filter_on_base_qual(
            x[f"{Schema.RESULTS_NUCLEOTIDES}_{sample_type}"],
            x[f"{Schema.QUALITY}_{sample_type}"],
        ),
        axis=1,
    )
    return df.copy(deep=True)


def get_read_depth_after_filtering(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Adds column for post filtering read depth"""
    df[f"{Schema.READ_DEPTH_POST_FILTER}_{sample_type}"] = df[
        f"{Schema.RESULTS_NUCLEOTIDES_FILTERED}_{sample_type}"
    ].str.len()
    return df.copy(deep=True)


Preprocessor = Callable[[str, pd.DataFrame], pd.DataFrame]


def compose(sample_type: str, *functions: Preprocessor) -> Preprocessor:
    """Helper function to call all df functions sequencially"""
    partially_filled_funcs = [partial(func, sample_type) for func in functions]
    return reduce(
        lambda func1, func2: lambda x: func2(func1(x)), partially_filled_funcs
    )


def process_genome_data(
    piles: Pileups, min_base_qual: int, region: GenomicRegion
) -> pd.DataFrame | None:
    # load the data from the pileup file
    df = fill_df(piles, region)
    if df.empty:
        logging.warn(f"Genomic region {str(region)} has no data in pileup")
        return None
    pipe = [
        create_SV_column,
        remove_carrots_from_res,
        remove_indels_from_results_str,
        find_indels_in_results_str,
        partial(remove_low_quality_bases, phred=min_base_qual),
        get_read_depth_after_filtering,
    ]
    for sample_type in ["tumour", "normal"]:
        preprocessor = compose(sample_type, *pipe)
        df = preprocessor(df)

    return CallAllVaraints.get_calls(df)
