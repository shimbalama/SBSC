#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from functools import partial

import numpy as np
import pandas as pd
import scipy.stats as stats

from .df_schema import Schema, convert_to_upper, remove_ref_from_tumour_res


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

    @property
    def p_value(self):
        return self.caller()

    def caller(self) -> float:
        """Return p value from fisher exact test of given reads counts"""
        _, p_value = stats.fisher_exact(
            [
                [self.tumour_var, self.normal_var],
                [self.non_var_tumour_reads, self.non_var_normal_reads],
            ],
            alternative="greater",
        )
        return p_value


class CallVariant(ABC):
    """Operations"""

    @abstractmethod
    def call() -> None:
        pass


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
            normal_nucs_most_common: dict[str, int] = dict(
                Counter(normal_nucs_upper).most_common()
            )
            tumour_nucs_most_common: dict[str, int] = {
                t_var: count
                for t_var, count in Counter(tumour_nucs_no_ref_upper).items()
                if normal_nucs_most_common.get(t_var, 0) < 2
            }
            for putative_somatic_var, count in tumour_nucs_most_common.items():
                if putative_somatic_var == nuc and count > 4:
                    normal_count: int = normal_nucs_most_common.get(
                        putative_somatic_var, 0
                    )
                    counts = Counts(
                        count, normal_count, len(tumour_nucs), len(normal_nucs)
                    )
                    if counts.p_value < 0.01:
                        return putative_somatic_var, counts.p_value

            return np.NaN

        for col, nuc in [
            (Schema.A_CALLS, "A"),
            (Schema.T_CALLS, "T"),
            (Schema.C_CALLS, "C"),
            (Schema.G_CALLS, "G"),
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
        return df.copy(deep=True)


class CallInsertionOrDeletion(CallVariant):
    """Test if the most common INDEL is significantly more present in tumour"""

    def call(df: pd.DataFrame) -> pd.DataFrame:
        def test_top_putative_somatic_indel(
            tumour_indels: tuple[str, int], reads_t: int, reads_n: int, res_n: int
        ) -> tuple[str, float] | float:
            if isinstance(tumour_indels, float):
                return np.NaN
            putative_somatic_indel, count = tumour_indels
            normal_count: int = res_n.upper().count(putative_somatic_indel)
            if count > 8 and normal_count < 2:
                counts = Counts(count, normal_count, reads_t, reads_n)
                if counts.p_value < 0.01:
                    return putative_somatic_indel, counts.p_value
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
        return df.copy(deep=True)


class CallStructuralVariant(CallVariant):
    """Test if there are significantly more SVs in tumour"""

    def call(df: pd.DataFrame) -> pd.DataFrame:
        def test_top_putative_somatic_SV(
            tumour_SVs: int, normal_SVs: int, reads_t: int, reads_n: int
        ) -> list[tuple[str, float]] | float:
            # no limit as some reads naturally start and stop....
            # TODO - can check flags and exclude those that aren't mapped chimerically
            if tumour_SVs > 8:
                counts = Counts(tumour_SVs, normal_SVs, reads_t, reads_n)
                if counts.p_value < 0.01:
                    return "SV", counts.p_value
            return np.NaN

        df.loc[:, Schema.STRUCTURAL_CALLS] = df.apply(
            lambda x: test_top_putative_somatic_SV(
                x[f"{Schema.STRUCTURAL_VARIANTS}_tumour"],
                x[f"{Schema.STRUCTURAL_VARIANTS}_normal"],
                x[f"{Schema.READS}_tumour"],
                x[f"{Schema.READS}_normal"],
            ),
            axis=1,
        )
        return df.copy(deep=True)


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
        return df.copy(deep=True)[df[cols].notna().any(1)]


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
        df.loc[:, cols] = df.loc[:, cols].applymap(
            lambda x: (np.NaN, np.NaN) if pd.isna(x) else x
        )
        return df.copy(deep=True)


class reformatSNVCalls(CallVariant):
    def call(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        # https://stackoverflow.com/questions/35491274/split-a-pandas-column-of-lists-into-multiple-columns
        df[[Schema.A_CALLS, Schema.A_PVALUE]] = pd.DataFrame(
            df[Schema.A_CALLS].tolist(), index=df.index
        )
        df[[Schema.T_CALLS, Schema.T_PVALUE]] = pd.DataFrame(
            df[Schema.T_CALLS].tolist(), index=df.index
        )
        df[[Schema.C_CALLS, Schema.C_PVALUE]] = pd.DataFrame(
            df[Schema.C_CALLS].tolist(), index=df.index
        )
        df[[Schema.G_CALLS, Schema.G_PVALUE]] = pd.DataFrame(
            df[Schema.G_CALLS].tolist(), index=df.index
        )
        df[[Schema.INDEL_CALLS, Schema.INDEL_PVALUE]] = pd.DataFrame(
            df[Schema.INDEL_CALLS].tolist(), index=df.index
        )
        df[[Schema.STRUCTURAL_CALLS, Schema.STRUCTURAL_PVALUE]] = pd.DataFrame(
            df[Schema.STRUCTURAL_CALLS].tolist(), index=df.index
        )
        return df.copy(deep=True)


class CallAllVaraints:
    """Call each variant type and update df"""

    def get_calls(df: pd.DataFrame) -> pd.DataFrame:
        for caller in CallVariant.__subclasses__():
            df = caller.call(df)
        return df.copy(deep=True)
