from argparse import ArgumentParser
from functools import partial, reduce
from typing import Callable, Dict
import pysam
import numpy as np
import pandas as pd
from collections import Counter
import re
import scipy.stats as stats
from abc import ABC, abstractmethod
from SBSC.parsers.parse_pile import GenomicRegion


class Schema:#TODO ???drop 'pos_in_read', 'read_names'????
    CHROMOSOME = 'seq'
    POSITION = 'pos'
    REFERENCE = 'ref'
    READS = 'reads'
    RESULTS = 'res'
    QUALITY = 'qual'
    POSITION_IN_READ = 'pos_in_read'
    READ_NAMES = 'read_names'
    FLAGS = 'flags'
    MAPPING_QUALITY = 'map'
    STRUCTURAL_VARIANTS = 'SV'
    RESULTS_NO_CARET = 'res_no_caret'
    RESULTS_NUCLEOTIDES = 'res_nuc'
    RESULTS_NUCLEOTIDES_FILTERED = 'res_nuc_filtered'
    QUALITY_FILTERED = 'qual_filtred'
    RESULTS_INDELS = 'res_indel'
    SINGLE_NUCLEOTIDE_CALLS = 'SNV_calls'
    INDEL_CALLS = 'INDEL_calls'
    STRUCTURAL_CALLS = 'SV_calls'
    HOMOPOLYMER = 'In_or_ajacent_to_homopolymer_of_length'
    READ_DEPTH_POST_FILTER = 'read_depth_post_filtering'


def create_df(region: GenomicRegion, pileup: str) -> pd.DataFrame:
    '''Converts a genomic region of a pileup to a df'''
    keys = [name for schema, name in vars(Schema).items() if "__" not in schema][:10]
    tabixfile = pysam.TabixFile(pileup)
    coords = ['chr' + str(region.chrom), region.start, region.end]
    df = pd.DataFrame([line.split('\t') for line in tabixfile.fetch(*coords)], columns=keys)
    index_names = Schema.CHROMOSOME + ':' + Schema.POSITION
    df[index_names] = df[Schema.CHROMOSOME] + ':' + df[Schema.POSITION]
    df.set_index(index_names, inplace=True)
    df = df.astype({'reads': 'int64', 'pos': 'int64'})
    return df

def remove_positions_with_little_support_for_somatic_var(df: pd.DataFrame) -> pd.DataFrame:

    def putative_vars(res_t, res_n, ref):
        res_t = remove_ref_from_tumour_res(res_n, res_t, ref)
        if res_t:
            return True if Counter(res_t).most_common()[0][1] > 4 else False
        return False
    return df[df.apply(lambda x: putative_vars(
        x[f'{Schema.RESULTS}_tumour'],
        x[f'{Schema.RESULTS}_normal'],
        x[Schema.REFERENCE]),
        axis = 1
    )]

def remove_redundant_column(df: pd.DataFrame, col_name: str) -> pd.DataFrame:

    tumour_col = f'{col_name}_tumour'
    normal_col = f'{col_name}_normal'
    if not df[tumour_col].equals(df[normal_col]):
        raise ValueError(f'{col_name} differs!')
    df[col_name] = df[normal_col]
    del df[tumour_col]
    del df[normal_col]
   
    return df

def add_homoploymers(df: pd.DataFrame, region: Dict) -> pd.DataFrame:

    df[Schema.HOMOPOLYMER] = df[Schema.POSITION].apply(lambda x: region.get(x, 0))

    return df

def create_SV_column(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    col = f'{Schema.RESULTS}_{sample_type}'
    df[f'{Schema.STRUCTURAL_VARIANTS}_{sample_type}'] = df[col].str.count('^') + df[col].str.count('$')
    return df

def remove_carrots_from_res(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    # ^ (caret) marks the start of a read segment and the ASCII
    # of the character following `^' minus 33 gives the mapping quality
    def slice_and_dice(res):
        if '^' in res:
            if res.startswith('^'):
                res = ''.join([bit[1:] for bit in res.split('^')])
            else:
                res = ''.join([bit[1:] if j != 0 else bit for j, bit in enumerate(res.split('^'))])
        return res
    df[f'{Schema.RESULTS_NO_CARET}_{sample_type}'] = df[f'{Schema.RESULTS}_{sample_type}'].apply(slice_and_dice)
    return df

#def remove_low_quality(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
def convert_to_upper(res, ref):
    assert ref.isupper()
    return res.replace(',', ref).replace('.', ref).replace('*', '').upper()

def remove_ref_from_tumour_res(res_n, res_t, ref):

    most_common_base_in_normal = Counter(convert_to_upper(res_n, ref)).most_common()[0][0]
    return convert_to_upper(res_t, ref).replace(most_common_base_in_normal, '')

#tried complied regex, still slow
def separate_indels_from_res(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:
    def find_indels_in_res(res):
        regex_result = re.findall("\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+", res)
        corected_regex_result = []#can't find other way to stop regex adding unrelated nucs at end
        for indel in regex_result:
            indel_len = ''.join(char for char in indel if char.isdigit())
            indel_nucs = ''.join(char for char in indel if char.isalpha())[:int(indel_len)]
            corected_regex_result.append(''.join([indel[0], indel_len, indel_nucs]))
        indel_count = Counter(corected_regex_result)
        for indel in indel_count.keys():
            res = ''.join(res.split(indel))
        most_seen_indel = indel_count.most_common()[0] if indel_count else 'NA'
        most_seen_indel_upper = np.nan if most_seen_indel == 'NA' else (most_seen_indel[0].upper(), most_seen_indel[1])
        res = ''.join(char for char in res if char in '.,ATCGatcg*')
        return pd.Series([res, most_seen_indel_upper])
    nucs = f'{Schema.RESULTS_NUCLEOTIDES}_{sample_type}'
    qual = f'{Schema.QUALITY}_{sample_type}'
    df[[nucs, f'{Schema.RESULTS_INDELS}_{sample_type}']] = df[f'{Schema.RESULTS_NO_CARET}_{sample_type}'].apply(find_indels_in_res)
    assert df[nucs].str.len().equals(df[qual].str.len()); f'{nucs} and {qual} should have same length!! Exiting...'
    return df

def remove_low_quality_bases(sample_type: str, df: pd.DataFrame, phred: int) -> pd.DataFrame:
    def filter_on_base_qual(res, qual):
        new_res, new_qual = '', ''
        for base, quality in zip(res, qual):
            if (ord(quality) - 33) > phred:
                new_res += base
                new_qual += quality
        assert len(new_res) == len(new_qual); 'resuts should match qualities'
        return pd.Series([new_res, new_qual])
    df[[f'{Schema.RESULTS_NUCLEOTIDES_FILTERED}_{sample_type}',
        f'{Schema.QUALITY_FILTERED}_{sample_type}']] =  df.apply(
            lambda x: filter_on_base_qual(
                x[f'{Schema.RESULTS_NUCLEOTIDES}_{sample_type}'],
                x[f'{Schema.QUALITY}_{sample_type}']),
            axis = 1)
    return df

def get_read_depth_after_filtering(sample_type: str, df: pd.DataFrame) -> pd.DataFrame:

    df[f'{Schema.READ_DEPTH_POST_FILTER}_{sample_type}'] = df[f'{Schema.RESULTS_NUCLEOTIDES_FILTERED}_{sample_type}'].str.len()
    return df

Preprocessor = Callable[[str, pd.DataFrame], pd.DataFrame]

def compose(sample_type: str, *functions: Preprocessor) -> Preprocessor:
    partially_filled_funcs = [partial(func, sample_type) for func in functions]
    return reduce(lambda func1, func2: lambda x: func2(func1(x)), partially_filled_funcs)

class CallVariant(ABC):
    '''Operations'''
    @abstractmethod
    def call():
        pass

    def caller(count, normal_count, reads_t, reads_n, putative_somatic_var, vars):

        vals = [count, normal_count]
        nons = [int(reads_t) - count, int(reads_n) - normal_count]
        oddsratio, pvalue1 = stats.fisher_exact([vals, nons], alternative='greater')
        if pvalue1 < 0.01:
            vars.append((putative_somatic_var, pvalue1))
        return vars
    
class CallSingleNucleotideVariant(CallVariant):
    def call(df: pd.DataFrame, caller: CallVariant.caller) -> pd.DataFrame:
        def test_putative_somatic_SNVs(tumour_nucs, normal_nucs, ref):
            normal_nucs_upper = convert_to_upper(normal_nucs, ref)
            tumour_nucs_no_ref_upper = remove_ref_from_tumour_res(normal_nucs, tumour_nucs, ref)
            normal_nucs_most_common = dict(Counter(normal_nucs_upper).most_common())
            tumour_nucs_most_common = {
                t_var: count for t_var, count in Counter(tumour_nucs_no_ref_upper).items() if normal_nucs_most_common.get(t_var, 0) < 2
            }
            vars = []
            for putative_somatic_var, count in tumour_nucs_most_common.items():
                if count > 4:
                    normal_count = normal_nucs_most_common.get(putative_somatic_var, 0)
                    vars = caller(count, normal_count, len(tumour_nucs), len(normal_nucs), putative_somatic_var, vars)
            return vars if vars else np.nan
        df[Schema.SINGLE_NUCLEOTIDE_CALLS] = df.apply(lambda x: test_putative_somatic_SNVs(
            x[f'{Schema.RESULTS_NUCLEOTIDES_FILTERED}_tumour'],
            x[f'{Schema.RESULTS_NUCLEOTIDES_FILTERED}_normal'],
            x[Schema.REFERENCE]),
            axis = 1
        )
        return df
class CallInsertionOrDeletion(CallVariant):
    def call(df: pd.DataFrame, caller: CallVariant.caller) -> pd.DataFrame:
        def test_top_putative_somatic_indel(tumour_indels, reads_t, reads_n, res_n):
            vars = []
            if type(tumour_indels) != float:
                freq_of_most_common_tumour_indel_in_normal = res_n.upper().count(tumour_indels[0])
                if freq_of_most_common_tumour_indel_in_normal < 2:
                    normal_d = dict([(tumour_indels[0], freq_of_most_common_tumour_indel_in_normal)])
                    putative_somatic_var, count = tumour_indels
                    if count > 4:
                        normal_count = normal_d.get(putative_somatic_var, 0)
                        vars = caller(count, normal_count, reads_t, reads_n, putative_somatic_var, vars)
            return vars if vars else np.nan
        df[Schema.INDEL_CALLS] = df.apply(lambda x: test_top_putative_somatic_indel(
            x[f'{Schema.RESULTS_INDELS}_tumour'],
            x[f'{Schema.READS}_tumour'],
            x[f'{Schema.READS}_normal'],
            x[f'{Schema.RESULTS}_normal']),
            axis = 1
        )
        return df
class CallStructuralVariant(CallVariant):
    def call(df: pd.DataFrame, caller: CallVariant.caller) -> pd.DataFrame:
        def test_top_putative_somatic_SV(tumour_SVs, normal_SVs, reads_t, reads_n):
            #no limit as some reads naturally start and stop....
            vars = []
            if tumour_SVs > 4:
                vars = caller(tumour_SVs, normal_SVs, reads_t, reads_n, 'SV', vars)
            return vars if vars else np.nan
        df[Schema.STRUCTURAL_CALLS] = df.apply(lambda x: test_top_putative_somatic_SV(
            x[f'{Schema.STRUCTURAL_VARIANTS}_tumour'],
            x[f'{Schema.STRUCTURAL_VARIANTS}_normal'],
            x[f'{Schema.READS}_tumour'],
            x[f'{Schema.READS}_normal']),
            axis = 1
        )
        return df

class CallAllVaraints:
    '''Call each variant type and update df'''
    @abstractmethod
    def get_calls(df):
        for caller in CallVariant.__subclasses__():
            df = caller.call(df, CallVariant.caller)
        return df


def process_genome_data(args: ArgumentParser, region: GenomicRegion) -> pd.DataFrame:
    # load the data from the pileup file
    df_tumour = create_df(region, args.cancer_pile)
    df_normal = create_df(region, args.normal_pile)
    df = df_normal.join(df_tumour, how='inner', lsuffix='_normal', rsuffix='_tumour')
    if df.empty:
        return None
    # df.index = df.index.astype('int64')
    # df[Schema.HOMOPOLYMER] = df.index.isin(region.homopolymers)
    
    df = remove_redundant_column(df, Schema.REFERENCE)
    df = remove_redundant_column(df, Schema.CHROMOSOME)
    df = remove_redundant_column(df, Schema.POSITION)
    df = remove_positions_with_little_support_for_somatic_var(df)
    df = add_homoploymers(df, region.get_hom_pol_lengths())
    pipe = [
        create_SV_column,
        remove_carrots_from_res,
        separate_indels_from_res,
        partial(remove_low_quality_bases, phred=args.min_base_qual),
        get_read_depth_after_filtering
    ]
    for sample_type in ['tumour', 'normal']:
        preprocessor = compose(
            sample_type,
            *pipe
        )
        df = preprocessor(df)
    df = CallAllVaraints.get_calls(df)
    return df[df[['SNV_calls','INDEL_calls','SV_calls']].notna().any(1)]
