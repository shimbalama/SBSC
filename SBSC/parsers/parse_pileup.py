import logging
import re
from abc import ABC, abstractmethod
from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
from functools import partial, reduce
from typing import Any, Callable

import numpy as np
import pandas as pd
import pysam
import scipy.stats as stats
from matplotlib.path import Path
from SBSC.parsers.parse_pile import GenomicRegion


#might slip this schema into tmp/final TODO
class Schema:  # TODO ???drop 'pos_in_read', 'read_names'????
    CHROMOSOME = "seq"
    POSITION = "pos"
    REFERENCE = "ref"
    READS = "reads"
    RESULTS = "res"
    QUALITY = "qual"
    POSITION_IN_READ = "pos_in_read"
    READ_NAMES = "read_names"
    FLAGS = "flags"
    MAPPING_QUALITY = "map"
    STRUCTURAL_VARIANTS = "SV"
    RESULTS_NO_CARET = "res_no_caret"
    RESULTS_NUCLEOTIDES = "res_nuc"
    RESULTS_NUCLEOTIDES_FILTERED = "res_nuc_filtered"
    QUALITY_FILTERED = "qual_filtred"
    RESULTS_INDELS = "res_indel"
    SINGLE_NUCLEOTIDE_CALLS = "SNV_calls"
    A_CALLS = "A_calls"
    T_CALLS = "T_calls"
    C_CALLS = "C_calls"
    G_CALLS = "G_calls"
    A_PVALUE = "A_pvalue"
    T_PVALUE = "T_pvalue"
    C_PVALUE = "C_pvalue"
    G_PVALUE = "G_pvalue"
    INDEL_CALLS = "INDEL_calls"
    STRUCTURAL_CALLS = "SV_calls"
    INDEL_PVALUE = "INDEL_pvalue"
    STRUCTURAL_PVALUE = "SV_pvalue"
    HOMOPOLYMER = "In_or_ajacent_to_homopolymer_of_length"
    READ_DEPTH_POST_FILTER = "read_depth_post_filtering"

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
    return df

def fill_df(input_data: Pileups, region: GenomicRegion) -> pd.DataFrame:
    logging.info(f"Creating a df for {str(region)}")
    df_tumour = create_df(region, input_data.cancer_pileup)
    df_normal = create_df(region, input_data.normal_pileup)
    df = df_normal.join(df_tumour, how="inner", lsuffix="_normal", rsuffix="_tumour")
    for col in [Schema.REFERENCE, Schema.CHROMOSOME, Schema.POSITION]:
        df = remove_redundant_column(df, col)
    check_the_ref_seq_matches_the_pileup_seq(df, region)
    df = remove_positions_with_little_support_for_somatic_var(df)
    df = add_homoploymers(df, region.homopolymer_lengths)

    return df

def remove_positions_with_little_support_for_somatic_var(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Filters out positions with 4 or less putative somatic SNVs"""

    def putative_vars(res_t: str, res_n: str, ref: str) -> bool:
        res_t = remove_ref_from_tumour_res(res_n, res_t, ref)
        if res_t:
            return True if Counter(res_t).most_common()[0][1] > 4 else False
        return False

    return df[
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

    return df


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
            df=%s and len seq=%s", str(region), len(df), len(region.seq)
        )

def add_homoploymers(df: pd.DataFrame, region: dict[int, int]) -> pd.DataFrame:
    """Adds a homopolymer column to df"""
    df[Schema.HOMOPOLYMER] = df[Schema.POSITION].apply(lambda x: region.get(x, 0))

    return df


def create_SV_column(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Adds SV column to df"""

    results_col = f"{Schema.RESULTS}_{sample_type}"
    SV_col = f"{Schema.STRUCTURAL_VARIANTS}_{sample_type}"
    df[SV_col] = df[results_col].str.count("^") + df[results_col].str.count("$")
    return df


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
    return df


def convert_to_upper(res: str, ref: str) -> str:
    """converts results str to upper case, replaces "," and "." to reference str
    and removes "*" to allow phred scores to align"""
    assert ref.isupper()
    return res.replace(",", ref).replace(".", ref).replace("*", "").upper()


def remove_ref_from_tumour_res(normal_results: str, tumour_results: str, ref: str) -> str:
    """Removes any reference seq from the tumour results str"""
    most_common_base_in_normal = Counter(convert_to_upper(normal_results, ref)).most_common()[0][0]
    return convert_to_upper(tumour_results, ref).replace(most_common_base_in_normal, "")

def correct_regex(regex_result: list[str]) -> list[str]:
    '''Helper function for removing indels from str'''
    corected_regex_result: list[str] = []  
    # can't find other way to stop regex adding unrelated nucs at end
    for indel in regex_result:
        indel_len = "".join(char for char in indel if char.isdigit())
        indel_nucs = "".join(char for char in indel if char.isalpha())[
            : int(indel_len)
        ]
        corected_regex_result.append("".join([indel[0], indel_len, indel_nucs]))
    return corected_regex_result

def remove_indels_from_results_str(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Separates the indels out from the pre-processed results str"""

    def find_and_remove_indels(res: str) -> str:
        regex_result: list[str] = re.findall("\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+", res)
        corected_regex_result = correct_regex(regex_result) 
        indel_count = Counter(corected_regex_result)
        for indel in indel_count.keys():
            res = "".join(res.split(indel))
        return "".join(char for char in res if char in ".,ATCGatcg*")

    nucs = f"{Schema.RESULTS_NUCLEOTIDES}_{sample_type}"
    qual = f"{Schema.QUALITY}_{sample_type}"
    df[nucs] = df[f"{Schema.RESULTS_NO_CARET}_{sample_type}"].apply(find_and_remove_indels)
    assert df[nucs].str.len().equals(df[qual].str.len()); f"{nucs} and {qual} should have same length"
    return df

def find_indels_in_results_str(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Separates the indels out from the pre-processed results str"""

    def find_indels_in_res(res: str) -> tuple[str, int] | float:
        regex_result = re.findall("\+[0-9]+[ACGTN]+|\-[0-9]+[ACGTN]+", res.upper())
        corected_regex_result = correct_regex(regex_result)
        indel_count = Counter(corected_regex_result)    
        return indel_count.most_common()[0] if indel_count else np.NaN

    df[f"{Schema.RESULTS_INDELS}_{sample_type}"] = df[f"{Schema.RESULTS_NO_CARET}_{sample_type}"].apply(find_indels_in_res)
    return df

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
        assert len(new_res) == len(new_qual); "resuts should match qualities"
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
    return df


def get_read_depth_after_filtering(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    """Adds column for post filtering read depth"""
    df[f"{Schema.READ_DEPTH_POST_FILTER}_{sample_type}"] = df[
        f"{Schema.RESULTS_NUCLEOTIDES_FILTERED}_{sample_type}"
    ].str.len()
    return df


Preprocessor = Callable[[str, pd.DataFrame], pd.DataFrame]


def compose(sample_type: str, *functions: Preprocessor) -> Preprocessor:
    """Helper function to call all df functions sequencially"""
    partially_filled_funcs = [partial(func, sample_type) for func in functions]
    return reduce(
        lambda func1, func2: lambda x: func2(func1(x)), partially_filled_funcs
    )


@dataclass
class Counts:
    tumour_var: int
    normal_var: int
    total_tumour_reads: int
    total_normal_reads: int

    @property
    def non_var_tumour_reads(self):
        return self.total_tumour_reads - self.tumour_var

    @property
    def non_var_normal_reads(self):
        return self.total_normal_reads - self.normal_var


class CallVariant(ABC):
    """Operations"""

    @abstractmethod
    def call() -> None:
        pass

    @staticmethod #this interface is now not applicable to all TODO
    def caller(var_counts: Counts) -> float:
        """Return p value from fisher exact test of given reads counts"""
        _, p_value = stats.fisher_exact(
            [
                [var_counts.tumour_var, var_counts.normal_var],
                [var_counts.non_var_tumour_reads, var_counts.non_var_normal_reads],
            ],
            alternative="greater",
        )

        return p_value


class CallSingleNucleotideVariant(CallVariant):
    def call(df: pd.DataFrame) -> pd.DataFrame:
        """Call SNVs"""

        def test_putative_somatic_SNVs(
            nuc: str, tumour_nucs: str, normal_nucs: str, ref: str
        ) -> tuple[str, float] | float:
            normal_nucs_upper = convert_to_upper(normal_nucs, ref)
            tumour_nucs_no_ref_upper = remove_ref_from_tumour_res(
                normal_nucs, tumour_nucs, ref
            )
            normal_nucs_most_common: dict[str, int] = dict(Counter(normal_nucs_upper).most_common())
            tumour_nucs_most_common: dict[str, int] = {
                t_var: count
                for t_var, count in Counter(tumour_nucs_no_ref_upper).items()
                if normal_nucs_most_common.get(t_var, 0) < 2
            }
            for putative_somatic_var, count in tumour_nucs_most_common.items():
                if putative_somatic_var == nuc and count > 4:
                    normal_count: int = normal_nucs_most_common.get(putative_somatic_var, 0)
                    p_value: float = CallVariant.caller(
                        Counts(count, normal_count, len(tumour_nucs), len(normal_nucs))
                    )
                    if p_value < 0.01:
                        return putative_somatic_var, p_value

            return np.NaN

        for col, nuc in [
            (Schema.A_CALLS, 'A'),
            (Schema.T_CALLS, 'T'),
            (Schema.C_CALLS, 'C'),
            (Schema.G_CALLS, 'G')
        ]:
            test = partial(test_putative_somatic_SNVs, nuc)
            df[col] = df.apply(
                lambda x: test(
                    x[f"{Schema.RESULTS_NUCLEOTIDES_FILTERED}_tumour"],
                    x[f"{Schema.RESULTS_NUCLEOTIDES_FILTERED}_normal"],
                    x[Schema.REFERENCE],
                ),
                axis=1,
            )
        return df


class CallInsertionOrDeletion(CallVariant):
    """Test if the most common INDEL is significantly more present in tumour"""

    def call(df: pd.DataFrame) -> pd.DataFrame:
        def test_top_putative_somatic_indel(
            tumour_indels: int, reads_t: int, reads_n: int, res_n: int
        ) -> tuple[str, float] | float:
            if not isinstance(tumour_indels, float):
                freq_of_most_common_tumour_indel_in_normal: int = res_n.upper().count(
                    tumour_indels[0]
                )
                if freq_of_most_common_tumour_indel_in_normal < 2:
                    normal_d: dict[str, int] = dict(
                        [(tumour_indels[0], freq_of_most_common_tumour_indel_in_normal)]
                    )
                    putative_somatic_INDEL, count = tumour_indels
                    if count > 8:
                        normal_count = normal_d.get(putative_somatic_INDEL, 0)
                        p_value: float = CallVariant.caller(
                            Counts(count, normal_count, reads_t, reads_n)
                        )
                        if p_value < 0.01:
                            return putative_somatic_INDEL, p_value
            return np.NaN

        df[Schema.INDEL_CALLS] = df.apply(
            lambda x: test_top_putative_somatic_indel(
                x[f"{Schema.RESULTS_INDELS}_tumour"],
                x[f"{Schema.READS}_tumour"],
                x[f"{Schema.READS}_normal"],
                x[f"{Schema.RESULTS}_normal"],
            ),
            axis=1,
        )
        return df


class CallStructuralVariant(CallVariant):
    """Test if there are significantly more SVs in tumour"""

    def call(df: pd.DataFrame) -> pd.DataFrame:
        def test_top_putative_somatic_SV(
            tumour_SVs: int, normal_SVs: int, reads_t: int, reads_n: int
        ) -> list[tuple[str, float]] | float:
            # no limit as some reads naturally start and stop....
            #TODO - can check flags and exclude those that aren't mapped chimerically
            if tumour_SVs > 8:
                p_value = CallVariant.caller(
                    Counts(tumour_SVs, normal_SVs, reads_t, reads_n)
                )
                if p_value < 0.01:
                    return "SV", p_value
            return np.NaN

        df[Schema.STRUCTURAL_CALLS] = df.apply(
            lambda x: test_top_putative_somatic_SV(
                x[f"{Schema.STRUCTURAL_VARIANTS}_tumour"],
                x[f"{Schema.STRUCTURAL_VARIANTS}_normal"],
                x[f"{Schema.READS}_tumour"],
                x[f"{Schema.READS}_normal"],
            ),
            axis=1,
        )
        return df


class removeNonCalls(CallVariant):
    """Remove all rows that are NaN for SNV, INDEL and SV calls"""

    def call(df: pd.DataFrame) -> pd.DataFrame:
        cols = [
            Schema.A_CALLS,
            Schema.T_CALLS,
            Schema.C_CALLS,
            Schema.G_CALLS,
            Schema.INDEL_CALLS,
            Schema.STRUCTURAL_CALLS,
        ]
        return df[df[cols].notna().any(1)]

class replaceNaN(CallVariant):
    def call(df: pd.DataFrame) -> pd.DataFrame:
        cols = [
            Schema.A_CALLS,
            Schema.T_CALLS,
            Schema.C_CALLS,
            Schema.G_CALLS,
            Schema.INDEL_CALLS,
            Schema.STRUCTURAL_CALLS,
        ]
        df.loc[:,cols] = df.loc[:,cols].applymap(lambda x: (np.NaN, np.NaN) if pd.isna(x) else x)
        return df

class reformatSNVCalls(CallVariant):
    def call(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        # https://stackoverflow.com/questions/35491274/split-a-pandas-column-of-lists-into-multiple-columns
        df[[Schema.A_CALLS, Schema.A_PVALUE]] = pd.DataFrame(df[Schema.A_CALLS].tolist(), index=df.index)
        df[[Schema.T_CALLS, Schema.T_PVALUE]] = pd.DataFrame(df[Schema.T_CALLS].tolist(), index=df.index)
        df[[Schema.C_CALLS, Schema.C_PVALUE]] = pd.DataFrame(df[Schema.C_CALLS].tolist(), index=df.index)
        df[[Schema.G_CALLS, Schema.G_PVALUE]] = pd.DataFrame(df[Schema.G_CALLS].tolist(), index=df.index)
        df[[Schema.INDEL_CALLS, Schema.INDEL_PVALUE]] = pd.DataFrame(df[Schema.INDEL_CALLS].tolist(), index=df.index)
        df[[Schema.STRUCTURAL_CALLS, Schema.STRUCTURAL_PVALUE]] = pd.DataFrame(df[Schema.STRUCTURAL_CALLS].tolist(), index=df.index)
        return df
        
class CallAllVaraints:
    """Call each variant type and update df"""

    def get_calls(df):
        for caller in CallVariant.__subclasses__():
            df = caller.call(df)
        return df


def process_genome_data(piles: Pileups, min_base_qual: int, region: GenomicRegion) -> pd.DataFrame | None:
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
