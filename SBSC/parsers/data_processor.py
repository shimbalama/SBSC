from dataclasses import dataclass
from functools import partial, reduce
from typing import Callable
import pysam
import numpy as np
import pandas as pd
from collections import Counter
import re
import scipy.stats as stats



Preprocessor = Callable[[pd.DataFrame], pd.DataFrame]


class Schema:#TODO ???drop 'pos_in_read', 'read_names'????
    CHROMOSOME = "seq"
    POSITION = "pos"
    REFERENCE = "ref"
    READS = "reads"
    RESULTS = "res"
    QUALITY = "qual"
    POSITION_IN_READ = "pos_in_read"
    FLAGS = "flags"
    MAPPING_QUALITY = "map"
    STRUCTURAL_VARIANTS = 'SV'
    RESULTS_NO_CARET = 'res_no_caret'
    RESULTS_NUCLEOTIDES = 'res_nuc'
    RESULTS_INDELS = 'res_indel'
    SINGLE_NUCLEOTIDE_CALLS = 'SNV_calls'
    INDEL_CALLS = 'INDEL_calls'
    STRUCTURAL_CALLS = 'SV_calls'


def create_df(pile):#: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    '''Converts a genomic region of a pileup to a df'''
    keys = [name for schema, name in vars(Schema).items() if "__" not in schema][:-1]
    tabixfile = pysam.TabixFile(pile)
    df = pd.DataFrame([line.split('\t') for line in tabixfile.fetch()], columns=keys)
    df.set_index('pos', inplace=True)
    return df

def remove_redundant_column(df: pd.DataFrame, col_name: str) -> pd.DataFrame:

    tumour_col = f'{col_name}_tumour'
    normal_col = f'{col_name}_normal'
    if not df[tumour_col].equals(df[normal_col]):
        raise ValueError(f'{col_name} differs!')
    df[col_name] = df[normal_col]
    del df[tumour_col]
    del df[normal_col]
   
    return df

def create_SV_column(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
    col = f'{Schema.RESULTS}_{sample_type}'
    df[f'{Schema.STRUCTURAL_VARIANTS}_{sample_type}'] = df[col].str.count('^') + df[col].str.count('$')
    return df

# def create_indel_column(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
#     df[f'SV_{sample_type}'] = df[f'res_{sample_type}'].str.count('^') + df[f'res_{sample_type}'].str.count('$')
#     return df

def remove_carrots_from_res(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
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

def remove_positions_with_little_support_for_somatic_var(df: pd.DataFrame) -> pd.DataFrame:
    def putative_vars(res_t, res_n, ref):
        res_t = remove_ref_from_tumour_res(res_n, res_t, ref)
        if res_t:
            return True if Counter(res_t).most_common()[0][1] > 4 else False
        return False
    return df[df.apply(lambda x: putative_vars(
        x[f'{Schema.RESULTS_NO_CARET}_tumour'],
        x[f'{Schema.RESULTS_NO_CARET}_normal'],
        x[Schema.REFERENCE]),
        axis = 1
    )]
    
#tried complied regex, still slow
def separate_indels_from_res(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
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

def caller(count, normal_count, reads_t, reads_n, putative_somatic_var, vars):

    vals = [count, normal_count]
    nons = [int(reads_t) - count, int(reads_n) - normal_count]
    oddsratio, pvalue1 = stats.fisher_exact([vals, nons], alternative='greater')
    if pvalue1 < 0.01:
        vars.append((putative_somatic_var, pvalue1))
    return vars
    
def call_var(df: pd.DataFrame) -> pd.DataFrame:
    def test_putative_somatic_SNVs(tumour_nucs, normal_nucs, reads_t, reads_n, ref, res_n):
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
                vars = caller(count, normal_count, reads_t, reads_n, putative_somatic_var, vars)
        return vars if vars else np.nan
    df[Schema.SNV_CALLS] = df.apply(lambda x: test_putative_somatic_SNVs(
        x[f'{Schema.RESULTS_NUCLEOTIDES}_tumour'],
        x[f'{Schema.RESULTS_NUCLEOTIDES}_normal'],
        x[f'{Schema.READS}_tumour'],
        x[f'{Schema.READS}_normal'],
        x[Schema.REFERENCE],
        x[f'{Schema.READS}_normal']),
        axis = 1
    )
    def test_top_putative_somatic_indel(tumour_indels, reads_t, reads_n, ref, res_n):
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
        x[Schema.REFERENCE],
        x[f'{Schema.READS}_normal']),
        axis = 1
    )
    def test_top_putative_somatic_SV(tumour_SVs, normal_SVs, reads_t, reads_n, ref, res_n):
        #no limit as some reads naturally start and stop....
        vars = []
        if tumour_SVs > 4:
            vars = caller(tumour_SVs, normal_SVs, reads_t, reads_n, 'SV', vars)
        return vars if vars else np.nan
    df[Schema.STRUCTURAL_CALLS] = df.apply(lambda x: test_top_putative_somatic_SV(
        x[f'{Schema.STRUCTURAL_VARIANTS}_tumour'],
        x[f'{Schema.STRUCTURAL_VARIANTS}_normal'],
        x[f'{Schema.READS}_tumour'],
        x[f'{Schema.READS}_normal'],
        x[Schema.REFERENCE],
        x[f'{Schema.READS}_normal']),
        axis = 1
    )
    return df


def compose(*functions: Preprocessor) -> Preprocessor:
    return reduce(lambda f, g: lambda x: g(f(x)), functions)


def process_genome_data(path: str) -> pd.DataFrame:
    # load the data from the pileup file
    df_tumour = create_df('/Users/liam/code/personal/SBSC/tests/resources/HCC1937_t_RNF157.pileup.gz')
    df_normal = create_df('/Users/liam/code/personal/SBSC/tests/resources/HCC1937_n_RNF157.pileup.gz')
    df = df_normal.join(df_tumour, how='inner', lsuffix='_normal', rsuffix='_tumour')
    preprocessor = compose(
        partial(remove_redundant_column, col_name=Schema.REFERENCE),
        partial(remove_redundant_column, col_name=Schema.SEQUENCE),
        partial(create_SV_column, sample_type='tumour'),
        partial(create_SV_column, sample_type='normal'),
        partial(remove_carrots_from_res, sample_type='tumour'),
        partial(remove_carrots_from_res, sample_type='normal'),
        remove_positions_with_little_support_for_somatic_var,
        partial(separate_indels_from_res, sample_type='tumour'),
        partial(separate_indels_from_res, sample_type='normal'),
        call_var
    )
    return preprocessor(df)



'''
df = create_SV_column(df, 't')
df = create_SV_column(df, 'n')
df = remove_carrots_from_res(df, 't')
df = remove_carrots_from_res(df, 'n')
df = remove_positions_with_little_support_for_somatic_var(df)
print(df.shape)
df = separate_indels_from_res(df, 't')
df = separate_indels_from_res(df, 'n')
df = call_var(df)
df = df[df[['SNV_calls', 'indel_call', 'SV_call']].notna().any(1)]
df.shape'''